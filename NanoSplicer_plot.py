import matplotlib.pyplot as plt
import concurrent.futures
import time
import sys
import h5py
import getopt
import timeit
import os
import numpy as np
import re
import fcntl
import pysam
from tqdm import tqdm
from intervaltree import Interval
from intervaltree import IntervalTree
from pathlib import Path
from ont_fast5_api.fast5_interface import get_fast5_file        
from scipy.stats.mstats import theilslopes


import helper
from junction_identification import find_candidate, canonical_site_finder, \
                                                        candidate_motif_generator
from dtw import dtw


# MedMad normalisation
def medmad(x):
    x = np.array(x)
    x_med = np.median(x)
    x_mad = np.median(np.abs(x - x_med))
    return x_med, x_mad

def adjust_by_event(path, cum_matrix, candidate_squiggle):
    match_event_array = candidate_squiggle[path[:,1] - 1]
    event_index = match_event_array[1:,] != match_event_array[:-1,]
    event_index = np.append(np.array([True]), event_index[: ,0])
    event_index = np.add.accumulate(event_index) - 1
    
    cum_path = cum_matrix[path[:, 0], path[:, 1]]
    score_path = np.append(cum_path[0], cum_path[1:] -  cum_path[:-1])
    
    adjusted_event_loglik = np.array([])
    for i in range(0, max(event_index) + 1):
        adjusted_event_loglik = np.append(adjusted_event_loglik,
                                        np.mean(score_path[event_index == i]))
    
    return adjusted_event_loglik
        



# parse command line arg
def parse_arg():
    def print_help():
        print("\n\nUsage: python {} [OPTIONS]".format(argv[0]))
        print("Options:\n\tIndexing:")
        print('\t\t-i \t.bam/.sam file (required)')
        print('\t\t-f \tpath to fast5s (directory path) (required)')
        print('\t\t-r \tGenome reference file (required)')
        print("\t\t-o \toutput path, default: 'NanoSplicer_out'")
        print('\t\t-T \tNumber of events trimmed from scrappie model')
        print('\t\t-t \tNumber of samples trimmed from raw signal')
        print('\t\t-F \tFlanking sequence size in each side')
        print('\t\t-w \twindow size for searching the candidate')
        print('\t\t-G \trun junctions in csv file')
        return None

    argv = sys.argv
    if len(argv) <= 2:     
        print_help()       # print help doc when no command line args provided
        sys.exit(0)
    
    try: 
        opts, args = getopt.getopt(argv[1:],"hi:f:r:o:T:t:ab:F:w:G:",
                    ["help=","input_alignment=","input_fast5_dir=",
                    "genome_ref=","output_path=", "trim_model=",
                    "trim_signal=","dtw_adj","bandwidth=",
                    "flank_size=", "window=","group="])
    
    except getopt.GetoptError:
        print_help()
        sys.exit(0)

    # DEFAULT VALUE
    alignment_file, fast5_dir, genome_ref = None, None, None
    output_path = 'NanoSplicer_out'
    trim_signal = 6
    trim_model = 2
    dtw_adj = False
    bandwidth = 0.4
    flank_size = 20
    window = 10


    for opt, arg in opts:
        if opt in ('-h', "--help"):
            print_help()
            sys.exit(0)
        elif opt in ("-i", "--input_alignment"):
            alignment_file = arg
        elif opt in ("-f", "--input_fast5_dir"):
            fast5_dir = arg
        elif opt in ("-r", "--genome_ref"):
            genome_ref = arg
        elif opt in ("-o", "--output_path"):
            output_path = arg
        elif opt in ("-T", "--trim_model"):
           trim_model = int(arg)
        elif opt in ("-t", "--trim_signal"):
           trim_signal = int(arg)
        elif opt in ("-a", "--dtw_adj"):
           dtw_adj = True
        elif opt in ("-b", "--bandwidth"):
           bandwidth = float(arg)
        elif opt in ("-F", "--flank_size"):
           flank_size = int(arg)
        elif opt in ("-w", "--window"):
           flank_size = int(arg)
        elif opt in ("-G", "--group"):
            group_filename = arg


    # check input
    if not alignment_file or not fast5_dir or not genome_ref:
        print("Error:Missing input files.")
        sys.exit(1)

    # choose the version of dtw (sum: minimize the sum of cost for the path)
    def dtw_local_alignment(candidate_squiggle, junction_squiggle, 
                            bandwidth = bandwidth, dist_type = None):
        
        return dtw(candidate_squiggle=candidate_squiggle,
                   junction_squiggle=junction_squiggle, 
                   band_prop = bandwidth,
                   dist_type = dist_type).dtw_local_alignment()
    
    return fast5_dir, output_path, alignment_file, genome_ref, \
            bandwidth, trim_model, trim_signal, flank_size, \
            window, group_filename
def count_in_tmp(filename):
    '''
        count once in a certain tmp file
    '''
    f = open(filename, "r+")
    fcntl.flock(f,fcntl.LOCK_EX)
    count = f.read()
    count = int(count) + 1 if count else 1
    f.seek(0)
    f.write(str(count))
    f.truncate()
    f.close()

def get_gaps_in_read(AlignedSegment):
    blocks = AlignedSegment.get_blocks()
    gap = set([(blocks[i-1][1], blocks[i][0]) for i in range(len(blocks))])
    return gap
    

def read_file_for_ploting(filename):
    ids = []
    minimap_site = []
    NanoSplicer_site = []
    true_site = []
    minimap_priors = []
    NanoSplicer_priors = []
    with open(filename, 'r') as f:
        for line in f:
            rid, mst, mnd, nst, nnd, tst, tnd = line.strip().split(',')
            mp ,np = 2,2
            mst, mnd, nst, nnd, mp ,np, tst, tnd \
                = int(mst), int(mnd), int(nst), int(nnd), int(mp), int(np), \
                    int(tst), int(tnd)
            ids.append(rid)
            minimap_site.append(Interval(mst, mnd))
            NanoSplicer_site.append(Interval(nst, nnd))
            true_site.append(Interval(tst, tnd))
            minimap_priors.append(mp)
            NanoSplicer_priors.append(np)
    return ids, minimap_site, NanoSplicer_site, true_site, \
                minimap_priors, NanoSplicer_priors

def run_multifast5(fast5_path, all_junctions, AlignmentFile, ref_FastaFile, 
                    window, chrID, flank_size, trim_model, trim_signal,
                    bandwidth, output_file, group_filename):
                    
    def dtw_local_alignment(candidate_squiggle, junction_squiggle, 
                            bandwidth = bandwidth, dist_type = None):
        
        return dtw(candidate_squiggle=candidate_squiggle,
                   junction_squiggle=junction_squiggle, 
                   band_prop = bandwidth,
                   dist_type = dist_type).dtw_local_alignment()
    AlignmentFile = pysam.AlignmentFile(AlignmentFile)
    ref_FastaFile = pysam.FastaFile(ref_FastaFile)
    multi_fast5 = get_fast5_file(fast5_path, 'r')
    reads_in_file = set(multi_fast5.get_read_ids())

    print("reading junction list file...")
    ids, minimap_site, NanoSplicer_site, true_site, \
        minimap_priors, NanoSplicer_priors \
            = read_file_for_ploting(group_filename)
    print("reading junction list file done")
    
    for i in range(len(minimap_site)):

        junction = minimap_site[i]
        
        if ids[i] not in  reads_in_file:
            continue
        
        print("searching for {}".format(ids[i]))

        overlap_read = AlignmentFile.fetch(chrID, junction.begin, junction.end)
        
        read = None
        for read1 in overlap_read:
            #print(read.qname, ids[i])
            if read1.qname == ids[i]:
                print("{} found".format(read1.qname))
                read = read1
        if not read:
            print("read not found")
            sys.exit()

        print("processing {}".format(read.qname))
        if read.qname != ids[i]:
            exit("weird")
            continue

        if (junction.begin, junction.end) not in get_gaps_in_read(read)\
            or read.qname not in reads_in_file:
            continue
        
        donor_lst, acceptor_lst = canonical_site_finder(junction, 
                                            ref_FastaFile, AlignmentFile, 
                                            window, chrID)
        
        candidates_pos, candidate_motif, motif_start, motif_end = \
                candidate_motif_generator(chrID, donor_lst, acceptor_lst, 
                                            flank_size, ref_FastaFile)

        print(minimap_site[i])
        index_m = candidates_pos.index((minimap_site[i].begin, minimap_site[i].end))
        index_n = candidates_pos.index((NanoSplicer_site[i].begin, NanoSplicer_site[i].end))
        index_t = candidates_pos.index((true_site[i].begin, true_site[i].end))
        
        print(index_m, index_n)
        if True:
            '''
            two candidates mode
            '''
            candidate_motif = [candidate_motif[index_t], candidate_motif[index_n]]
        if not candidate_motif:# or len(candidate_motif) == 1:
            continue

        candidate_motif_rev = [helper.reverse_complement(seq) 
                                for seq in candidate_motif]

        # trim signal (use fewer base)
        motif_start += trim_signal
        motif_end -= trim_signal

        # check reverse
        if read.is_reverse:
            candidates_for_read = candidate_motif_rev
        else:
            candidates_for_read = candidate_motif

        # tombo resquiggle
        try:
            tombo_results, tombo_start_clip, tombo_end_clip = \
                tombo_squiggle_to_basecalls(multi_fast5, read)
        except:
            print("tombo resquiggle failed!!")


        read_length = len(tombo_results.genome_seq) + tombo_start_clip + tombo_end_clip
        normalised_raw_signal = tombo_results.raw_signal/1.4826

        # genome pos to read pos mapping vector
        g_r_mapping = \
            genome_to_read_pos_conversion(read.cigarstring)

        # convert to read relative pos (forward direction)
        start_pos_rel_to_mapped_start = motif_start - read.reference_start
        end_pos_rel_to_mapped_start = motif_end - read.reference_start

        # check the start and end pos are valid
        if  start_pos_rel_to_mapped_start >= 0 and start_pos_rel_to_mapped_start < len(g_r_mapping):
            motif_start_read = g_r_mapping[start_pos_rel_to_mapped_start]
            # discard junction squiggle with the queried motif start/end mapped to gaps
            if motif_start_read == -1:
                print("Warning: Junction squiggle start index point to mapped intron, junction squiggle skipped.")
                #count_in_tmp(fbad_junction_mapping)
                continue
            elif g_r_mapping[start_pos_rel_to_mapped_start] == \
                g_r_mapping[start_pos_rel_to_mapped_start - 1]:
                motif_start_read += 1 
        else:
            print("candidate start pos out of bound.")
            continue
        
        if end_pos_rel_to_mapped_start < len(g_r_mapping):
            motif_end_read = g_r_mapping[end_pos_rel_to_mapped_start - 1] + 1
            # discard junction squiggle with the queried motif start/end mapped to gaps
            if motif_end_read == -1 + 1:
                print("Warning: Junction squiggle end index point to mapped intron, junction squiggle skipped.")
                continue
        else:
            print("candidate end pos out of bound.")
            continue

        # get signal
        if not read.is_reverse:
            seg_start = max(motif_start_read - tombo_start_clip, 0)
            seg_end = motif_end_read - tombo_start_clip
        else:
            seg_start = max(read_length - motif_end_read -1 - tombo_start_clip, 0)
            seg_end = read_length - motif_start_read - 1 - tombo_start_clip

        # take into account the end clip
        seg_end = min(seg_end, len(tombo_results.segs) - 1)

        # getting junction squiggle
        signal = normalised_raw_signal[
            tombo_results.segs[seg_start]:tombo_results.segs[seg_end]]
        if not len(signal):
            continue
        else:
            print("pass")
        junction_squiggle = np.array(signal, float)
        # outlier removal
        junction_squiggle = junction_squiggle[abs(junction_squiggle) < 3]

        model_dic = helper.expect_squiggle_dict(candidates_for_read, 
                                                trim = trim_model, uniform_dwell=4)


        cum_path = {}
        score_dict = {}
        squiggle_match = {}
        output_prefix = ''
        num_of_cand = len(candidates_for_read)
        SAVE_DATA = False
        aligned_base = {}


        # edited
        for j, candidate in enumerate(candidates_for_read):
            print(candidate)
            candidate_squiggle = np.array(model_dic[candidate],float)
            #dtw run
            path , score, cum_matrix = \
                dtw_local_alignment(candidate_squiggle=candidate_squiggle, 
                                    junction_squiggle = junction_squiggle)
            


            score_dict[j] = -1 * score
            
            # candidate squiggle in matched len of junction suqiggle
            squiggle_match[j] = candidate_squiggle[path[:,1] - 1, :]
            aligned_base[j] = candidate[trim_model + int((path[0,1] - 1)/4): trim_model+ int((path[-1,1] - 1)/4) + 1]
            cum_path[j] = cum_matrix[path[:, 0], path[:, 1]]
            
            def find_candidate_middle_pos(squiggle_match):
                out = [0]
                for i in range(1, len(squiggle_match)):
                    if (squiggle_match[i] == squiggle_match[i-1]).all():
                        out[-1] += 0.5
                    else:
                        out.append(i)
                return(out)

            # renormalisation (theilslopesrenor)
            if False:
                def likelihood_contribute(junction_squiggle,squiggle_match):    

                    diff = abs(junction_squiggle - squiggle_match[:,0])
                    z = diff/squiggle_match[:,1]
                    laplacc_b = squiggle_match[:,1]/np.sqrt(2)
                    #return 0.9189385 + z**2/2 #norm
                    #return 1.14473 + log(1+z**2) # t with df = 1
                    return np.log(2*laplacc_b) + diff/laplacc_b
                def abline(slope, intercept):
                    """Plot a line from slope and intercept"""
                    axes = plt.gca()
                    x_vals = np.array(axes.get_xlim())
                    y_vals = intercept + slope * x_vals
                    plt.plot(x_vals, y_vals, '--')

                #edited
                #if j == 0:
                if True:
                    medslope, medintercept = theilslopes(y = junction_squiggle,
                                            x = squiggle_match[j][:,0])[:2]
                
                plt.plot(squiggle_match[j][:,0], junction_squiggle,'o')
                
                abline(medslope, medintercept)
                plt.title("slope:{},intercept:{}".format(medslope, medintercept))
                plt.savefig("{}_regression{}.png".format(output_prefix, j))
                plt.close()
                squiggle_match[j][:,0] = medslope *  squiggle_match[j][:,0] + medintercept
                squiggle_match[j][:,1] = np.sqrt(medslope) *  squiggle_match[j][:,1]
                candidate_squiggle[:,0] = medslope *  candidate_squiggle[:,0] + medintercept
                candidate_squiggle[:,1] = np.sqrt(medslope) *  candidate_squiggle[:,1]
                cum_path[j] = np.add.accumulate(
                            likelihood_contribute(junction_squiggle,
                                                  squiggle_match[j]))

                # # redo dtw
                if False:
                    path , score, cum_matrix = \
                        dtw_local_alignment(candidate_squiggle=candidate_squiggle, 
                                        junction_squiggle = junction_squiggle)
                    # candidate squiggle in matched len of junction suqiggle
                    squiggle_match[j] = candidate_squiggle[path[:,1] - 1, :]
                    cum_path[j] = cum_matrix[path[:, 0], path[:, 1]]
            
            # LR plot
            if j == num_of_cand - 1:
                distinguish_points_idx = \
                    np.where(np.abs(squiggle_match[0][:,0] - squiggle_match[1][:,0]) > \
                     squiggle_match[1][:,1])[0]
                LR_cum = cum_path[1] - cum_path[0]
                LR_con = np.append(LR_cum[0], LR_cum[1:] - LR_cum[:-1])
                fig, axes = plt.subplots(nrows=2, figsize=(20,12))
                axes[0].plot(distinguish_points_idx,
                                np.repeat(max(junction_squiggle)+0.5,len(distinguish_points_idx)), 'o')
                axes[0].plot(squiggle_match[0][:,0], 'r',linewidth=0.4, 
                        label = "minimap2 candidate")
                axes[0].plot(squiggle_match[1][:,0], 'g',linewidth=0.4,
                        label = "NanoSplcier candidate")
                axes[0].plot(junction_squiggle, 'b',linewidth=0.4, 
                        label = "junction squiggle")
                axes[0].legend(frameon=False, loc='best')
                axes[0].set_title(
                    "LR_full = {}, # of distinguish points = {}, LR_middle = {}, adj_dist_full = {},{}, adj_dist_middle = {},{}\n {}\n{}".format(
                    score_dict[0] - score_dict[1], len(distinguish_points_idx),
                    np.sum(LR_con[distinguish_points_idx]),
                    -score_dict[0]/len(path),-score_dict[1]/len(path),
                    np.mean(np.append(cum_path[0][0], cum_path[0][1:] - cum_path[0][:-1])[distinguish_points_idx]),
                    np.mean(np.append(cum_path[1][0], cum_path[1][1:] - cum_path[1][:-1])[distinguish_points_idx]),
                    aligned_base[0], aligned_base[1]
                ))
                middle_pos = find_candidate_middle_pos(squiggle_match[0])
                for i in range(len(aligned_base[0])):
                    axes[0].annotate(aligned_base[0][i],
                        xy=(middle_pos[i], max(junction_squiggle) + 0.2), xycoords='data', ha='center')
                middle_pos = find_candidate_middle_pos(squiggle_match[1])
                for i in range(len(aligned_base[1])):
                    axes[0].annotate(aligned_base[1][i],
                        xy=(middle_pos[i], min(junction_squiggle) - 0.2), xycoords='data', ha='center')
                                
                ax1 = axes[1]
                ax2 = ax1.twinx()
                # for j_i in range(num_of_cand):
                #     ax1.plot(squiggle_match[j_i][:,0], 'r',linewidth=0.2)
                # plt.savefig("{}_cum_LR_{}.png".format(output_prefix, read.qname))
                # plt.close()

                ax1.plot(squiggle_match[0][:,0], 'r',linewidth=0.4, 
                        label = "minimap2 candidate")
                ax1.plot(squiggle_match[1][:,0], 'g',linewidth=0.4, 
                        label = "NanoSplicer candidate")

                # # plot contribution of each data point in LR
                if True:
                    LR_cum = cum_path[1] - cum_path[0]
                    ax2.plot(np.append(LR_cum[0], LR_cum[1:] - LR_cum[:-1]), 
                            color='#4b0082', linewidth=0.4, alpha=1, 
                            dash_capstyle='round', label = 'LR contribution')
                    ax1.set_xlabel('Index')
                    ax1.set_ylabel('Candidate squiggle', color='g')
                    ax2.set_ylabel('log likelihood ratio contribution', color='#4b0082')
                    ax1.legend(frameon=False, loc='best')
                    fig.savefig("{}.png".format(read.qname))
                    plt.close()
                
                #plot cumulative contribution of LR of each data point
                if False:    
                    ax2.plot(cum_path[1] - cum_path[0], 
                            color='#4b0082', linewidth=0.2, alpha=1, dash_capstyle='round')
                    #plot diff of candidate instead
                    #plt.plot(squiggle_match[1][:,1] - squiggle_match[0][:,1])
                    plt.savefig("{}_cum_LR_{}.png".format(output_prefix, read.qname))
                    plt.close()
                
                    # med_j, mad_j = medmad(junction_squiggle)
                    # junction_squiggle_norm = (junction_squiggle - med_j)/ mad_j
                    # candidate_squiggle[:,0] = (candidate_squiggle[:,0] - med_c) / mad_c
                    # candidate_squiggle[:,1] = candidate_squiggle[:,1] / mad_c
                
        # Plot DTW 
            if False:
                print(score)
                if SAVE_DATA:
                    np.savetxt("{}_candidate{}_{}.csv".format(output_prefix, read.qname,j), squiggle_match[j], delimiter=",")
                    np.savetxt("junction_squiggle{}.csv".format(read.qname), junction_squiggle, delimiter=",")
                fig, axes = plt.subplots(nrows=4, figsize=(30,40))
                axes[0].plot(junction_squiggle,linewidth = 10)
                axes[0].tick_params(labelsize=40)
                            #axes[1].figure(figsize=(10,5))
                axes[1].plot(candidate_squiggle[:,0], color = "orange",linewidth =7)
                axes[1].plot(candidate_squiggle[:,0] + candidate_squiggle[:,1], ":",color = "orange",linewidth = 4)
                axes[1].plot(candidate_squiggle[:,0] - candidate_squiggle[:,1], ":",color = "orange",linewidth = 4)
                axes[1].tick_params(labelsize=40)   
                axes[1].set_title("Scrappie model",fontsize=30, pad = 1)
                axes[2].plot(junction_squiggle,linewidth = 10)
                axes[2].tick_params(labelsize=40)
                #path = np.array(path[::-1])
                axes[2].plot(path[:,0]-1, candidate_squiggle[[path[:,1]-1]][:,0],'',color = "orange",linewidth = 7)
                axes[2].plot(path[:,0]-1, candidate_squiggle[[path[:,1]-1]][:,0]\
                                                + candidate_squiggle[[path[:,1]-1]][:,1],':',color = "orange",linewidth = 4)
                axes[2].plot(path[:,0]-1, candidate_squiggle[[path[:,1]-1]][:,0]\
                                                - candidate_squiggle[[path[:,1]-1]][:,1],':',color = "orange",linewidth = 4)
                # axes[2].set_title(
                #     str(i) + \
                #     "\nDist: {:.2f}, path length: {}, Adjusted dist:"\
                #     " {:.2f}, prior: {}".format(
                #             score,len(path),
                #             score/len(path), 
                #             [minimap_priors[i], NanoSplicer_priors[i]][j])
                #         ,fontsize=40, pad = 1)

                pos = axes[3].imshow(cum_matrix, cmap='hot', interpolation='nearest',aspect='auto')
                fig.colorbar(pos,ax = axes[3])
                axes[3].plot(path[:,1], path[:,0])
                axes[3].set_title("Alignment path",fontsize=30, pad = 20)
                if read.is_reverse:
                    strand = 1
                else:
                    strand = 0
                fig.savefig("{}_fig{}_{}.png".format(output_prefix,read.qname, j))
                plt.close()

'''
Functions to select junction squiggles without the requirement of transcript reference
Input:
    BED line
    Cigar string
    fast5 filename
output:
    Normalised junction squiggles
'''

def tombo_squiggle_to_basecalls(multi_fast5, AlignedSegment):
    '''
    script from tombo (https://gist.github.com/marcus1487/5cf083644143aafc770dceb213783346)
    Input:
        read_fast5: fast5 filename
    Returns:
        rsqgl_results <class>: Resquiggle result
           # useful rsqgl_results attributes:
                #   - genome_seq: clipped basecalls to model-able positions (can't get k-mers and thus expected levels for first and last few bases)
                #   - raw_signal: median normalized raw signal values (only for portion of raw signal mapped to model-able basecalls)
                #   - read_start_rel_to_raw: start of rsqgl_results.raw_signal and rsqgl_results.segs within original all_raw_signal
                #   - segs: start position for each model-able base (should be one longer than length of genome_seq to include last end point) 
    '''
    from tombo import tombo_helper, tombo_stats, resquiggle
    # extract read info
    #fast5s_f = get_fast5_file(read_fast5_fn, 'r')
    #fast5_data = h5py.File
    seq_samp_type = tombo_helper.seqSampleType('DNA', False)
    #seq_data = resquiggle.get_read_seq(fast5_data, 'Basecall_1D_000', 'BaseCalled_template', seq_samp_type, 0)
    
    if AlignedSegment.is_reverse:
        read_seq = helper.reverse_complement(AlignedSegment.seq)
    else:
        read_seq = AlignedSegment.seq

    
    seq_data = tombo_helper.sequenceData(seq=read_seq, 
                                         id=AlignedSegment.qname, 
                                         mean_q_score=np.mean(AlignedSegment.query_qualities))
    # prep tombo objects
    std_ref = tombo_stats.TomboModel(seq_samp_type=seq_samp_type)
    start_clip = std_ref.central_pos
    end_clip = std_ref.kmer_width - std_ref.central_pos - 1
    rsqgl_params = tombo_stats.load_resquiggle_parameters(seq_samp_type)

    # extract raw signal
    #all_raw_signal = tombo_helper.get_raw_read_slot(fast5_data)['Signal'][:]
    all_raw_signal = multi_fast5.get_read(AlignedSegment.qname).get_raw_data()
    # if this is a direct RNA read, flip raw signal to process from 5' to 3'
    if seq_samp_type.rev_sig:
        all_raw_signal = all_raw_signal[::-1]

    # spoof mapping read to its own basecalls (with pseudo-clipped bases in genome_location for model-able postions)
    map_results = tombo_helper.resquiggleResults(
        align_info=tombo_helper.alignInfo(ID=seq_data.id, Subgroup='BaseCalled_template', ClipStart=0, ClipEnd=0, Insertions=0,
                                        Deletions=0, Matches=len(seq_data.seq), Mismatches=0),
        genome_loc=tombo_helper.genomeLocation(Start=std_ref.central_pos, Strand='+', Chrom='NA'),
        genome_seq=seq_data.seq, raw_signal=all_raw_signal, mean_q_score=seq_data.mean_q_score)

    # align raw signal to basecalls
    try:
        rsqgl_results = resquiggle.resquiggle_read(
        map_results, std_ref, rsqgl_params, all_raw_signal=all_raw_signal)
    except:
        rsqgl_params = rsqgl_params ._replace(bandwidth=2000)
        rsqgl_results = resquiggle.resquiggle_read(
        map_results, std_ref, rsqgl_params, all_raw_signal=all_raw_signal)

    return rsqgl_results, start_clip, end_clip


def genome_to_read_pos_conversion(cigar):
    '''
    Input:
        g_pos: 0-base position within genome
        cigar: CIGAR string from BAM/SAM file
    
    Returns:
        r_pos: 0-base position within a read
    '''
    cigar_long = []
    for count, type in re.findall('(\d+)([A-Za-z])', cigar):
        cigar_long += int(count) * [type]
    r_index = -1
    g_r_mapping = []
    for i in cigar_long:
        if i in "MIS":
            r_index += 1
        if i in "MD":
            g_r_mapping.append(r_index)
        if i == "N":
            g_r_mapping.append(-1)
    return g_r_mapping


def main():

    #temp file output
    out_fn= "NanoSplicer_out/plot.tsv"

    # get command line input
    fast5_dir, output_path, alignment_file, genome_ref, \
        bandwidth, trim_model, trim_signal, \
        flank_size, window, group_filename = parse_arg()
    os.system("mkdir -p {}".format(output_path))

    # create tmp files
    os.system("mkdir -p {}/.tmp".format(output_path))
    fcount_pass = "{}/.tmp/pass".format(output_path)
    ftombo_fail = "{}/.tmp/tombo_fail".format(output_path)
    fabnormal_sam_flag = "{}/.tmp/abnormal_sam_flag".format(output_path)
    fbad_junction_mapping = "{}/.tmp/bad_junction_mapping".format(output_path)
    fnoncanonical_read = "{}/.tmp/noncanonical_read".format(output_path)
    
    os.system("touch {} {} {} {} {}".format(
        fcount_pass,
        ftombo_fail,
        fabnormal_sam_flag,
        fbad_junction_mapping,
        fnoncanonical_read
        ))
    
    # get fast5 filenames recursively
    fast5_paths = Path(fast5_dir).rglob('*.fast5') # iterator

    # get intron boundary position.
    fbam = pysam.AlignmentFile(alignment_file)
    fref = pysam.FastaFile(genome_ref)
    f_fetch = fbam.fetch() # reads iterator
    f_introns = fbam.find_introns(f_fetch) # dictionary
    chrID = os.path.basename(alignment_file)[:-4] # e.g 'NC_000001.11.bam' to chrID
    
    # junction identification (use IntervalTree for fast query)
    intron_tree = IntervalTree()
    for (begin, end), data in f_introns.items():
        intron_tree.addi(begin, end, data)

        ## built non overlapped range
    intron_tree_non_overlaps = intron_tree.copy()
    intron_tree_non_overlaps.merge_overlaps()

    # run specific junction
    all_candidates = []

    # start to processing the fast5 (multiread format)
    print("running process")

    for fast5_path in fast5_paths:
        run_multifast5(fast5_path, 
                      all_candidates,
                      alignment_file,
                      genome_ref,
                      window,
                      chrID,
                      flank_size,
                      trim_model,
                      trim_signal,
                      bandwidth,
                      out_fn, 
                      group_filename)

if __name__ == "__main__":
    main()
