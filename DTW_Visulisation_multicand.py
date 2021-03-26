'''
visulisation code design

Input:
    subset of pandas dataframe of JWR to be tested (in HDF5), using key = "data"

Output:
    Figure 1:
     Correct cases: comparison for the best two candidate (figure) which containing:
     Incorrect cases: camparison for the true one and the NanoSplicer one:
        1. two candidate squiggles in two difference colors
        2. Highlight distinguishing point
    Figure 2: log LR contribution and accumulative log LR
        1. Highlight distinguishing point
    Results:
        1. DTW score for both candidate
        2. Aligned length
        3. Number of distinguishing point
        4. log LR contribution from the distinguishing point         
    Option to output the candidate and junction squiggles
'''
import textwrap
import pandas as pd
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

# denosing
from skimage.restoration import denoise_wavelet

import helper
from junction_identification import find_candidate, canonical_site_finder, \
                                                        candidate_motif_generator
from dtw import dtw


def __log_likelihood(a, b_mean, b_std, 
                        truncate_quantile = 0.01):
    '''
     log likelihood by assuming normal distribution
    retrun:
        log density in standard normal distribution
    '''
    def laplace_log_density(x, mean, sd, max_diff = None):                
        diff = np.abs(mean - x)
        if len(max_diff):
            diff = np.minimum(diff, max_diff)
        b = sd/np.sqrt(2)
        return -np.log(2*b) - diff/b
    def laplace_quantile(mean, sd, q):
        b = sd/np.sqrt(2)
        if q > 0.5:
            return mean - b*np.log(2-2*q) 
        else:
            return mean + b*np.log(2*q) 

    if truncate_quantile:
        max_diff = np.abs(b_mean - 
            laplace_quantile(b_mean, b_std, q = truncate_quantile))
        return laplace_log_density(a, b_mean, b_std, max_diff)
    else:
        return laplace_log_density(a, b_mean, b_std)


# parse command line arg
def parse_arg():
    def print_help():
        help_message =\
        '''
        Usage: python {} [OPTIONS] <plot hdf5 file>
        Options:
            -i      .bam/.sam file (required)
            -f      path to fast5s (required)
            -r      Genome reference file (required)
            -o      output path <default: "NanoSplicer_out">
            -T      Number of events trimmed from scrappie model <default: 2>
            -t      Number of bases trimmed from raw signam, the raw samples are
                        match to basecalled bases by tombo <default :6>
            -F      Flanking sequence size in each side of candidate searching 
                        window <default: 20>
            -w      Candidate searching window size <default: 10>
            
        '''.format(argv[0])

        print(textwrap.dedent(help_message))

    argv = sys.argv
    if len(argv) <= 2:     
        print_help()       # print help doc when no command line args provided
        sys.exit(0)
    
    try: 
        opts, args = getopt.getopt(argv[1:],"hi:f:r:o:T:t:ab:F:w:",
                    ["help=","input_alignment=","input_fast5_dir=",
                    "genome_ref=","output_path=", "trim_model=",
                    "trim_signal=","dtw_adj","bandwidth=",
                    "flank_size=", "window="])
    except getopt.GetoptError:
        print("ERROR:Invalid input.")
        print_help()
        sys.exit(1)

    # DEFAULT VALUE
    alignment_file, fast5_dir, genome_ref = None, None, None
    output_path = 'NanoSplicer_out'
    trim_signal = 6
    trim_model = 2
    dtw_adj = False
    bandwidth = 0.4
    flank_size = 20
    window = 10
    print(opts)
    pd_file = args[0]

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
           window = int(arg)

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
            bandwidth, trim_model, trim_signal, flank_size, window, pd_file

def main():
    # temp para
    #chrID = "chrIS"
    chrID = "NC_000001.11"
    out_fn = ''

    # parse command line argument
    fast5_dir, output_path, alignment_file, genome_ref, \
        bandwidth, trim_model, trim_signal, \
        flank_size, window, pd_file =  parse_arg()
    
    # import data
    test_df = pd.read_hdf(pd_file, key = 'data')
    
    # get fast5 filenames recursively
    fast5_paths = Path(fast5_dir).rglob('*.fast5') # iterator

    # start to processing the fast5 (multiread format)
    print("running process")
    from itertools import repeat
    with concurrent.futures.ProcessPoolExecutor(32) as executor:
        results = list(tqdm(executor.map(run_multifast5, list(fast5_paths), 
                     repeat(test_df), # all_junction in NanoSplicer2.py
                     repeat(alignment_file), 
                     repeat(genome_ref),
                     repeat(window),
                     repeat(chrID),
                     repeat(flank_size),
                     repeat(trim_model),
                     repeat(trim_signal),
                     repeat(bandwidth),
                     repeat(out_fn)),total = len(list(fast5_paths))))

# running a single process 
def run_multifast5(fast5_path, plot_df, AlignmentFile, ref_FastaFile, 
                    window, chrID, flank_size, trim_model, trim_signal,
                    bandwidth, output_file):
    '''
    Description:
        Run a single process within a certain multi-read fast5
    Input:
        fast5_path:     path to the target fast5 (multi-read version only)
        plot_df:  subset of the NanoSplcier output to be plotted(set of splice junctions to be test, it is not the 
                            splice junction candidates but the splice junctions to
                            match the junction within read to be queried)
        AlignmentFile:  BAM filename
        ref_fastaFile:  reference genome fasta file
        window:         candidate searching window size
        chrID:          chromosome name, come along with the all_junctions to locate the 
                            junction within reads
        flank_size:     flank size attached to the junction searching window to get candidate
                            motifs
        trim_models:    # of bases to be trimmed on candidate squiggles
        trim_signals:   # of bases to be trimmed on juntion squiggles
        bandwidth:      bandwidth for dtw
        output_file:    outbupt file name (.csv)    
    Output:
        on each linein the output_file write NanoSplice output in the following columns

        <to be added>
    '''
    # dtw version (function alias)
    def dtw_local_alignment(candidate_squiggle, junction_squiggle, 
                            bandwidth = bandwidth, dist_type = None):
        return dtw(candidate_squiggle=candidate_squiggle,
                   junction_squiggle=junction_squiggle, 
                   band_prop = bandwidth,
                   dist_type = dist_type, 
                   truncate_quantile = 0.01).dtw_local_alignment()

    # prepare the inputs
    AlignmentFile = pysam.AlignmentFile(AlignmentFile) 
    ref_FastaFile = pysam.FastaFile(ref_FastaFile)
    multi_fast5 = get_fast5_file(fast5_path, 'r')
    reads_in_file = set(multi_fast5.get_read_ids())
   
    # loop though each row in output_df
    for row in plot_df.itertuples():
        if row.id not in reads_in_file:
            continue
    #for junc_id, junction in all_junctions:
        
        # add a function to parse the row get:
            # junction.begin
            # junction.end
            # read id
            # candidates_tuples
            # true site
            # minimap2 site
            # NanoSplicer site
        
        # search read in bam file
        overlap_read = AlignmentFile.fetch(chrID, row.site_minimap[0], row.site_minimap[1])
        read = None
        for i in overlap_read:
            if i.qname == row.id:
                print("{} found".format(i.qname))
                read = i
        if not read:
            print("read {} is not found!".format(row.id))
            continue
        
        # 
            # # find GT-AG pattern in a window nearby
            # donor_lst, acceptor_lst = canonical_site_finder(junction, 
            #                                   ref_FastaFile, AlignmentFile, 
            #                                   window, chrID)
            
            # if not donor_lst or not acceptor_lst:
            #     candidates_tuple = []
            # else:
            #     candidates_tuple = list(itertools.product(donor_lst, acceptor_lst))
            
            # #  add best supported junction if it is not GT-AG
            # if (junction.begin, junction.end) not in candidates_tuple:
            #     candidates_tuple.append((junction.begin, junction.end))

        # get candidate motifs
        candidates_pos, candidate_motif, motif_start, motif_end = \
                candidate_motif_generator(chrID, row.candidates, 
                                          flank_size, ref_FastaFile)

        if not candidate_motif:# or len(candidate_motif) == 1:
            continue
        
        # Subset candidates  (two candidate) incomplete
        index_t = candidates_pos.index(row.junc_supported[0])
        index_n = candidates_pos.index(row.site_NanoSplicer)
        index_m = candidates_pos.index(row.site_minimap)
        if False:
            if len(row.junc_supported) > 1:
                print("read {} skipped because more than 1 true"
                "candidate were found".format(row.id))
                continue
            # minimap2 correct but NanoSplicer wrong
            if index_m == index_t and index_n != index_t:
                candidate_motif = [candidate_motif[index_m], candidate_motif[index_n]]
            # minimap2 wrong but NanoSplicer correct
            if index_m != index_t and index_n == index_t:
                candidate_motif = [candidate_motif[index_m], candidate_motif[index_n]]

        # check reverse
        candidate_motif_rev = [helper.reverse_complement(seq) 
                                for seq in candidate_motif]
        if read.is_reverse:
            candidates_for_read = candidate_motif_rev
        else:
            candidates_for_read = candidate_motif     
 
        # trim signal (use fewer base)
        motif_start += trim_signal
        motif_end -= trim_signal
            
        # tombo resquiggle
        try:
            tombo_results, tombo_start_clip, tombo_end_clip = \
                tombo_squiggle_to_basecalls(multi_fast5, read)
        except:
            print("tombo resquiggle failed!!")

        read_length = len(tombo_results.genome_seq) \
                                + tombo_start_clip + tombo_end_clip
        normalised_raw_signal = tombo_results.raw_signal/1.4826

        # genome pos to read pos mapping vector
        g_r_mapping = \
            genome_to_read_pos_conversion(read.cigarstring)

        # convert to read relative pos (forward direction)
        start_pos_rel_to_mapped_start = motif_start - read.reference_start
        end_pos_rel_to_mapped_start = motif_end - read.reference_start

        if  start_pos_rel_to_mapped_start >= 0 \
                    and start_pos_rel_to_mapped_start < len(g_r_mapping):
            motif_start_read = g_r_mapping[start_pos_rel_to_mapped_start]

            # discard junction squiggle with the queried motif start/end 
            # mapped to gaps
            if motif_start_read == -1:
                continue
            elif g_r_mapping[start_pos_rel_to_mapped_start] == \
                g_r_mapping[start_pos_rel_to_mapped_start - 1]:
                motif_start_read += 1 
        else:
            continue

        if end_pos_rel_to_mapped_start < len(g_r_mapping):
            motif_end_read = g_r_mapping[end_pos_rel_to_mapped_start - 1] + 1
            if motif_end_read == -1 + 1:
                continue
        else:
            continue

        # get signal
        if not read.is_reverse:
            seg_start = max(motif_start_read - tombo_start_clip, 0)
            seg_end = motif_end_read - tombo_start_clip
        else:
            seg_start = \
                max(read_length - motif_end_read -1 - tombo_start_clip, 0)
            seg_end = read_length - motif_start_read - 1 - tombo_start_clip
        # take into account the end clip
        seg_end = min(seg_end, len(tombo_results.segs) - 1)
        
        # getting junction squiggle
        signal = normalised_raw_signal[
            tombo_results.segs[seg_start]:tombo_results.segs[seg_end]]

        if not len(signal):
            continue
        else:
            junction_squiggle = np.array(signal, float)

            # outlier removal
        junction_squiggle = junction_squiggle[abs(junction_squiggle) < 3]        
        
        # candidate squiggle
        uniform_dwell = 8
        minimum_point_for_dist_seg = 4
        model_dic = helper.expect_squiggle_dict(seqs=candidates_for_read, 
                                                trim=trim_model,
                                                model='squiggle_r94',
                                                uniform_dwell=uniform_dwell)
        
        cum_path = {}
        score_dict = {}
        squiggle_match = {}
        output_prefix = ''
        num_of_cand = len(candidates_for_read)
        SAVE_DATA = False
        WAVELET_DENOISING = True
        aligned_base = {}
        prior_ratio = 9
        candidate_preference = np.array(row.candidate_preference)
        # edited

        for j, candidate in enumerate(candidates_for_read):
            candidate_squiggle = np.array(model_dic[candidate],float)
            
            if WAVELET_DENOISING:
                junction_squiggle_wav = \
                    denoise_wavelet(junction_squiggle, method='BayesShrink', 
                              mode='soft', wavelet_levels=1, wavelet='sym8', 
                              rescale_sigma='True')
                path , score, cum_matrix = \
                    dtw_local_alignment(candidate_squiggle = candidate_squiggle, 
                                        junction_squiggle = junction_squiggle_wav)   
            #dtw run
            else:
                path , score, cum_matrix = \
                    dtw_local_alignment(candidate_squiggle = candidate_squiggle, 
                                    junction_squiggle = junction_squiggle)

            # log likelihood
            score_dict[j] = -1 * score
            
            # candidate squiggle in matched len of junction suqiggle
            squiggle_match[j] = candidate_squiggle[path[:,1] - 1, :]
            aligned_base[j] = candidate[trim_model + int((path[0,1] - 1)/uniform_dwell): trim_model+ int((path[-1,1] - 1)/uniform_dwell) + 1]
            #cum_path[j] = cum_matrix[path[:, 0], path[:, 1]]
            
            def find_candidate_middle_pos(squiggle_match):
                out = [0]
                for i in range(1, len(squiggle_match)):
                    if (squiggle_match[i] == squiggle_match[i-1]).all():
                        out[-1] += 0.5
                    else:
                        out.append(i)
                return(out)

            def est_sd_from_flanking_seq(junction_squiggle, squiggle_match, mid_buffer = 8):
                # get sd from frank region
                n_base = 0
                pre_mean = 0
                diff = []
                for ix in range(len(junction_squiggle)):
                    if squiggle_match[ix,0] != pre_mean:
                        n_base += 1
                    if n_base < flank_size - mid_buffer:
                        diff.append(squiggle_match[ix,0] - junction_squiggle[ix])
                    else:
                        break
                for ix in range(len(junction_squiggle)):
                    if squiggle_match[::-1][ix,0] != pre_mean:
                        n_base += 1
                    if n_base < flank_size - mid_buffer:
                        diff.append(squiggle_match[::-1][ix,0] - junction_squiggle[::-1][ix])
                    else:
                        break
                    
                return np.std(np.array(diff))
                            
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

            # # redo dtw calculate sd
            if True:
                # get sd from frank region
                if WAVELET_DENOISING:
                    new_sd = est_sd_from_flanking_seq(junction_squiggle_wav, 
                                                            squiggle_match[j])
                    candidate_squiggle[:,1] = np.array([new_sd] * len(candidate_squiggle[:,1]))
                
                    path , score, cum_matrix = \
                        dtw_local_alignment(candidate_squiggle=candidate_squiggle, 
                                    junction_squiggle = junction_squiggle_wav)
                else:
                    new_sd = est_sd_from_flanking_seq(junction_squiggle, 
                                                            squiggle_match[j])
                    candidate_squiggle[:,1] = np.array([new_sd] * len(candidate_squiggle[:,1]))
                
                    path , score, cum_matrix = \
                        dtw_local_alignment(candidate_squiggle=candidate_squiggle, 
                                    junction_squiggle = junction_squiggle)

                # candidate squiggle in matched len of junction suqiggle
                squiggle_match[j] = candidate_squiggle[path[:,1] - 1, :]
                cum_path[j] = cum_matrix[path[:, 0], path[:, 1]]

            
            def get_distinguishing_segment(matched_candidate_ref, matched_candidate_ls):
                '''
                Input:
                    matched_candidate_ref: For the candidate that used as reference,  mean and sd from 
                                                                            that matched to junction squiggle
                    matched_candidate_ls: a list of queried candidate squiggle
                Output:
                    list of distinguishing segment: [(seg1_start, seg1_end), ...]
                '''
                matched_candidate_ls.append(matched_candidate_ref)
                all_means = np.array([x[:,0] for x in matched_candidate_ls])
                is_end_boundary = np.any(all_means[:,1:] != all_means[:,:-1], axis = 0)
                end_index = np.arange(1,len(is_end_boundary) + 1)[is_end_boundary]
                segment = list(zip(np.append(0, end_index[:-1]), end_index))
                seg_start = np.array(segment)[:,0]
                is_dist_seg = \
                    np.any(
                        np.abs(
                            all_means[:,seg_start] - matched_candidate_ref[seg_start,0]) \
                                    > matched_candidate_ref[seg_start,1], axis = 0)
                
                return np.array(segment), is_dist_seg

###################### LR plot
            if False and j == num_of_cand - 1:
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
                axes[0].plot(squiggle_match[0][:,0] + squiggle_match[0][:,1], 'y--',linewidth=0.4, 
                        label = "sd")
                axes[0].plot(squiggle_match[0][:,0] - squiggle_match[0][:,1], 'y--',linewidth=0.4)
                axes[0].plot(squiggle_match[1][:,0], 'g',linewidth=0.4,
                        label = "NanoSplicer candidate")
                axes[0].plot(junction_squiggle, 'b',linewidth=0.4, 
                        label = "junction squiggle")
                axes[0].legend(frameon=False, loc='best')
                axes[0].set_title(
                    "LR_full = {:.2f}, # of distinguish points = {}, LR_middle = {:.2f}, adj_dist_full = {:.2f},{:.2f}, adj_dist_middle = {:.2f},{:.2f}\n {}\n{}".format(
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
                    fig.savefig("{}_dtime8.png".format(read.qname))
                    plt.close()
                
################################## wavelet denoising

            # LR plot
            if False and j == num_of_cand - 1:
                distinguish_points_idx = \
                    np.where(np.abs(squiggle_match[0][:,0] - squiggle_match[1][:,0]) > \
                     squiggle_match[1][:,1])[0]
                junction_squiggle_wav = denoise_wavelet(junction_squiggle, method='BayesShrink', mode='soft', wavelet_levels=1, wavelet='sym8', rescale_sigma='True')
                LR_con = laplac_log_likelihood(junction_squiggle_wav, squiggle_match[0]) \
                         - laplac_log_likelihood(junction_squiggle_wav, squiggle_match[1])
                fig, axes = plt.subplots(nrows=2, figsize=(20,12))

                axes[0].plot(distinguish_points_idx,
                                np.repeat(max(junction_squiggle_wav)+0.5,len(distinguish_points_idx)), 'o')
                axes[0].plot(squiggle_match[0][:,0], 'r',linewidth=0.4, 
                        label = "minimap2 candidate")
                axes[0].plot(squiggle_match[0][:,0] + squiggle_match[0][:,1], 'y--',linewidth=0.4, 
                        label = "sd")
                axes[0].plot(squiggle_match[0][:,0] - squiggle_match[0][:,1], 'y--',linewidth=0.4)
                axes[0].plot(squiggle_match[1][:,0], 'g',linewidth=0.4,
                        label = "NanoSplicer candidate")
                axes[0].plot(junction_squiggle_wav, 'b',linewidth=0.4, 
                        label = "junction squiggle")
                axes[0].legend(frameon=False, loc='best')
                axes[0].set_title(
                    "LR_full = {:.2f}, # of distinguish points = {}, LR_middle = {:.2f}, adj_dist_full = {:.2f},{:.2f}, adj_dist_middle = {:.2f},{:.2f}\n {}\n{}".format(
                    np.sum(LR_con), len(distinguish_points_idx),
                    np.sum(LR_con[distinguish_points_idx]),
                    -score_dict[0]/len(path),-score_dict[1]/len(path),
                    np.mean(np.append(cum_path[0][0], cum_path[0][1:] - cum_path[0][:-1])[distinguish_points_idx]),
                    np.mean(np.append(cum_path[1][0], cum_path[1][1:] - cum_path[1][:-1])[distinguish_points_idx]),
                    aligned_base[0], aligned_base[1]
                ))
                middle_pos = find_candidate_middle_pos(squiggle_match[0])
                for i in range(len(aligned_base[0])):
                    axes[0].annotate(aligned_base[0][i],
                        xy=(middle_pos[i], max(junction_squiggle_wav) + 0.2), xycoords='data', ha='center')
                middle_pos = find_candidate_middle_pos(squiggle_match[1])
                for i in range(len(aligned_base[1])):
                    axes[0].annotate(aligned_base[1][i],
                        xy=(middle_pos[i], min(junction_squiggle_wav) - 0.2), xycoords='data', ha='center')
                                
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
                    ax2.plot(laplac_log_likelihood(junction_squiggle_wav, squiggle_match[0]) - laplac_log_likelihood(junction_squiggle_wav, squiggle_match[1]), 
                            color='#4b0082', linewidth=0.4, alpha=1, 
                            dash_capstyle='round', label = 'LR contribution')
                    ax1.set_xlabel('Index')
                    ax1.set_ylabel('Candidate squiggle', color='g')
                    ax2.set_ylabel('log likelihood ratio contribution', color='#4b0082')
                    ax1.legend(frameon=False, loc='best')
                    fig.savefig("{}_wavelet.png".format(read.qname))
                    plt.close()

################################## mean denoising (all candidate seg)
            # LR plot
            if True and j == num_of_cand - 1:
                print(1)
                def mean_denoise(junction_squiggle, squiggle_match_list):
                    mean_junction_squiggle = np.array([])
                    event = junction_squiggle[0:1]
                    
                    for i in range(1, len(junction_squiggle)):
                        if all([(x[i] == x[i-1]).all() for x in squiggle_match_list]):
                            event = np.append(event, junction_squiggle[i])
                        else:
                            mean_junction_squiggle = np.append(mean_junction_squiggle,[np.mean(event)] * len(event))
                            event = np.array([junction_squiggle[i]])
                    mean_junction_squiggle = np.append(mean_junction_squiggle,[np.mean(event)] * len(event))
                    return(mean_junction_squiggle)

                # minimap2 candidate as reference
                matched_candidate_ref = squiggle_match[index_m]
                squiggle_match_list = [squiggle_match[x] for x in range(num_of_cand)]
                junction_squiggle_mean = mean_denoise(junction_squiggle, squiggle_match_list)
                segment, is_dist_seg = \
                    get_distinguishing_segment(matched_candidate_ref, squiggle_match_list)
                
                is_dist_seg = is_dist_seg & (segment[:,1] - segment[:,0] >= minimum_point_for_dist_seg)
                print(2)

                # LR calculation (distinguishing segment only)
                dist_seg_LR = []
                p_wise_LR_list = []
                for cand in range(num_of_cand):
                    p_wise_LR = __log_likelihood(junction_squiggle_mean, 
                                    squiggle_match[cand][:,0],squiggle_match[cand][:,1])
                    dist_seg_LR.append(sum(p_wise_LR[segment[is_dist_seg][:,0]]))
                    p_wise_LR_list.append(p_wise_LR)
                post_prob = np.exp(dist_seg_LR)/sum(np.exp(dist_seg_LR))
                post_prob_prior = post_prob * (prior_ratio **candidate_preference)
                post_prob_prior = post_prob_prior/sum(post_prob_prior)

                
                print(3)
                fig = plt.figure(figsize=(20,5 * num_of_cand))
                subplot_index = 0
                for cand in range(num_of_cand):
                    if cand == index_m:
                        continue
                    
                    # frame subplot
                    subplot_index += 1
                    ax = fig.add_subplot(num_of_cand -1,1,subplot_index)
                    for i, seg in enumerate(segment):
                        ax.axvspan(seg[0], seg[1]-1, -0.2,1.2, 
                            clip_on=True, alpha = 0.2 if is_dist_seg[i] else 0.1 , edgecolor = 'black', 
                            lw = 0.2,  facecolor = 'red' if is_dist_seg[i] else '#E2D3CA')

                    # # data-point wise LR

                    # LR_con = __log_likelihood(junction_squiggle_mean, 
                    #             matched_candidate_ref[:,0],matched_candidate_ref[:,1]) - \
                    #         __log_likelihood(junction_squiggle_mean, 
                    #             squiggle_match[cand][:,0],squiggle_match[cand][:,1])

                    #ax.plot(distinguish_points_idx,
                    #                np.repeat(max(junction_squiggle_mean)+0.5,len(distinguish_points_idx)), 'o')
                    ax.plot(matched_candidate_ref[:,0], 'r',linewidth=0.6, 
                            label = "minimap2 candidate")
                    ax.plot(matched_candidate_ref[:,0] + matched_candidate_ref[:,1], 'y--',linewidth=0.4, 
                            label = "sd")
                    ax.plot(matched_candidate_ref[:,0] - matched_candidate_ref[:,1], 'y--',linewidth=0.4)
                    ax.plot(squiggle_match[cand][:,0], 'g',linewidth=0.6,
                            label = "Candidate")
                    ax.plot(junction_squiggle_mean, 'b',linewidth=0.4, 
                            label = "junction squiggle")
                    ax.plot(junction_squiggle, 'b',linewidth=0.4, 
                            label = "junction squiggle")

                    middle_pos = find_candidate_middle_pos(squiggle_match[index_m])
                    for i in range(len(aligned_base[index_m])):
                        ax.annotate(aligned_base[index_m][i],
                            xy=(middle_pos[i], max(junction_squiggle) + 0.05), xycoords='data', ha='center')
                    middle_pos = find_candidate_middle_pos(squiggle_match[cand])
                    for i in range(len(aligned_base[cand])):
                        ax.annotate(aligned_base[cand][i],
                            xy=(middle_pos[i], min(junction_squiggle) - 0.2), xycoords='data', ha='center')
                
                    ax.legend(frameon=False, loc='best')
                    
                    ax.set_title('Log LR: {:.2f}, post probability: {:.3f}, candidate_preference: {}, post probability(seq prior ratio = 9): {:.3f} '.format(
                        dist_seg_LR[cand] - dist_seg_LR[index_m], post_prob[cand], candidate_preference[cand], post_prob_prior[cand]), y=1.0, pad=-14)
                
                fig.suptitle('minimap2 candidate likelihood: {:.2f}, post probability: {:.3f}, candidate_preference: {}, post probability(seq prior ratio = 9): {:.3f}'.format(
                        dist_seg_LR[index_m] - dist_seg_LR[index_m], post_prob[index_m],candidate_preference[index_m], post_prob_prior[index_m]))

                fig.savefig("{}_mean_all_cand.png".format(read.qname))

################################## mean denoising (all candidate seg) LR contri
            # LR plot
            if True and j == num_of_cand - 1:
                print(1)
                def mean_denoise(junction_squiggle, squiggle_match_list):
                    mean_junction_squiggle = np.array([])
                    event = junction_squiggle[0:1]
                    
                    for i in range(1, len(junction_squiggle)):
                        if all([(x[i] == x[i-1]).all() for x in squiggle_match_list]):
                            event = np.append(event, junction_squiggle[i])
                        else:
                            mean_junction_squiggle = np.append(mean_junction_squiggle,[np.mean(event)] * len(event))
                            event = np.array([junction_squiggle[i]])
                    mean_junction_squiggle = np.append(mean_junction_squiggle,[np.mean(event)] * len(event))
                    return(mean_junction_squiggle)

                # minimap2 candidate as reference
                matched_candidate_ref = squiggle_match[index_m]
                squiggle_match_list = [squiggle_match[x] for x in range(num_of_cand)]
                junction_squiggle_mean = mean_denoise(junction_squiggle, squiggle_match_list)
                segment, is_dist_seg = \
                    get_distinguishing_segment(matched_candidate_ref, squiggle_match_list)
                
                is_dist_seg = is_dist_seg & (segment[:,1] - segment[:,0] >= minimum_point_for_dist_seg)
                print(2)

                # LR calculation (distinguishing segment only)
                dist_seg_LR = []
                p_wise_LR_list = []
                for cand in range(num_of_cand):
                    p_wise_LR = __log_likelihood(junction_squiggle_mean, 
                                    squiggle_match[cand][:,0],squiggle_match[cand][:,1])
                    dist_seg_LR.append(sum(p_wise_LR[segment[is_dist_seg][:,0]]))
                    p_wise_LR_list.append(p_wise_LR)
                post_prob = np.exp(dist_seg_LR)/sum(np.exp(dist_seg_LR))
                post_prob_prior = post_prob * (prior_ratio**candidate_preference)
                post_prob_prior = post_prob_prior/sum(post_prob_prior)
                
                print(3)
                fig = plt.figure(figsize=(20,5 * num_of_cand))
                subplot_index = 0
                for cand in range(num_of_cand):
                    if cand == index_m:
                        continue
                    
                    # frame subplot
                    subplot_index += 1
                    ax = fig.add_subplot(num_of_cand -1,1,subplot_index)
                    for i, seg in enumerate(segment):
                        ax.axvspan(seg[0], seg[1]-1, -0.2,1.2, 
                            clip_on=True, alpha = 0.2 if is_dist_seg[i] else 0.1 , edgecolor = 'black', 
                            lw = 0.2,  facecolor = 'red' if is_dist_seg[i] else '#E2D3CA')

                    # # data-point wise LR

                    # LR_con = __log_likelihood(junction_squiggle_mean, 
                    #             matched_candidate_ref[:,0],matched_candidate_ref[:,1]) - \
                    #         __log_likelihood(junction_squiggle_mean, 
                    #             squiggle_match[cand][:,0],squiggle_match[cand][:,1])

                    #ax.plot(distinguish_points_idx,
                    #                np.repeat(max(junction_squiggle_mean)+0.5,len(distinguish_points_idx)), 'o')
                    ax.plot(matched_candidate_ref[:,0], 'r',linewidth=0.6, 
                            label = "minimap2 candidate")
                    ax.plot(matched_candidate_ref[:,0] + matched_candidate_ref[:,1], 'y--',linewidth=0.4, 
                            label = "sd")
                    ax.plot(matched_candidate_ref[:,0] - matched_candidate_ref[:,1], 'y--',linewidth=0.4)
                    ax.plot(squiggle_match[cand][:,0], 'g',linewidth=0.6,
                            label = "Candidate")

                    if True:
                        ax2 = ax.twinx()
                        ax2.plot( p_wise_LR_list[cand] - p_wise_LR_list[index_m], 
                                color='#4b0082', linewidth=0.4, alpha=1, 
                                dash_capstyle='round', label = 'LR contribution')
                        ax2.set_ylabel('log likelihood ratio contribution', color='#4b0082')
                        ax2.axhline(0, linestyle = '--', lw = 0.6)

                    middle_pos = find_candidate_middle_pos(squiggle_match[index_m])
                    for i in range(len(aligned_base[index_m])):
                        ax.annotate(aligned_base[index_m][i],
                            xy=(middle_pos[i], max(junction_squiggle) + 0.05), xycoords='data', ha='center')
                    middle_pos = find_candidate_middle_pos(squiggle_match[cand])
                    for i in range(len(aligned_base[cand])):
                        ax.annotate(aligned_base[cand][i],
                            xy=(middle_pos[i], min(junction_squiggle) - 0.2), xycoords='data', ha='center')
                
                    ax.legend(frameon=False, loc='best')
                    
                    ax.set_title('Log LR: {:.2f}, post probability: {:.3f}, candidate_preference: {}, post probability(seq prior ratio = 9): {:.3f} '.format(
                        dist_seg_LR[cand] - dist_seg_LR[index_m], post_prob[cand], candidate_preference[cand], post_prob_prior[cand]), y=1.0, pad=-14)
                
                fig.suptitle('minimap2 candidate likelihood: {:.2f}, post probability: {:.3f}, candidate_preference: {}, post probability(seq prior ratio = 9): {:.3f}'.format(
                        dist_seg_LR[index_m] - dist_seg_LR[index_m], post_prob[index_m],candidate_preference[index_m], post_prob_prior[index_m]))
                fig.savefig("{}_LR_mean_all_cand.png".format(read.qname))

################################# segment mean (pairwise version)
            # LR plot
            if True and j == num_of_cand - 1:
                print(1)
                def mean_denoise(junction_squiggle, squiggle_match_list):
                    mean_junction_squiggle = np.array([])
                    event = junction_squiggle[0:1]
                    
                    for i in range(1, len(junction_squiggle)):
                        if all([(x[i] == x[i-1]).all() for x in squiggle_match_list]):
                            event = np.append(event, junction_squiggle[i])
                        else:
                            mean_junction_squiggle = np.append(mean_junction_squiggle,[np.mean(event)] * len(event))
                            event = np.array([junction_squiggle[i]])
                    mean_junction_squiggle = np.append(mean_junction_squiggle,[np.mean(event)] * len(event))
                    return(mean_junction_squiggle)

                # minimap2 candidate as reference
                matched_candidate_ref = squiggle_match[index_m]
                squiggle_match_list = [squiggle_match[x] for x in range(num_of_cand)]
                junction_squiggle_mean_list = \
                    [mean_denoise(junction_squiggle, [matched_candidate_ref, x]) for x in squiggle_match_list]
                segment_list = \
                    [get_distinguishing_segment(matched_candidate_ref, [x]) for x in squiggle_match_list]
                
                
                print(2)

                # LR calculation (distinguishing segment only)
                dist_seg_LR = []
                
                for cand in range(num_of_cand):
                    junction_squiggle_mean = junction_squiggle_mean_list[cand]
                    segment, is_dist_seg = segment_list[cand]
                    is_dist_seg = is_dist_seg & (segment[:,1] - segment[:,0] >= minimum_point_for_dist_seg)
                    p_wise_LR = __log_likelihood(junction_squiggle_mean, 
                                    squiggle_match[cand][:,0],squiggle_match[cand][:,1]) -\
                                 __log_likelihood(junction_squiggle_mean, 
                                    matched_candidate_ref[:,0],matched_candidate_ref[:,1])

                    dist_seg_LR.append(sum(p_wise_LR[segment[is_dist_seg][:,0]]))
                
                rel_to_ref_LR = np.exp(dist_seg_LR)/np.exp(dist_seg_LR[index_m])
                post_prob = rel_to_ref_LR/sum(rel_to_ref_LR)
                post_prob_prior = post_prob * (prior_ratio **candidate_preference)
                post_prob_prior = post_prob_prior/sum(post_prob_prior)
                
                print(3)
                fig = plt.figure(figsize=(20,5 * num_of_cand))
                subplot_index = 0
                for cand in range(num_of_cand):
                    junction_squiggle_mean = junction_squiggle_mean_list[cand]
                    segment, is_dist_seg = segment_list[cand]
                    is_dist_seg = is_dist_seg & (segment[:,1] - segment[:,0] >= minimum_point_for_dist_seg)

                    if cand == index_m:
                        continue
                    
                    # frame subplot
                    subplot_index += 1
                    ax = fig.add_subplot(num_of_cand -1,1,subplot_index)
                    for i, seg in enumerate(segment):
                        ax.axvspan(seg[0], seg[1]-1, -0.2,1.2, 
                            clip_on=True, alpha = 0.2 if is_dist_seg[i] else 0.1 , edgecolor = 'black', 
                            lw = 0.2,  facecolor = 'red' if is_dist_seg[i] else '#E2D3CA')

                    # # data-point wise LR

                    # LR_con = __log_likelihood(junction_squiggle_mean, 
                    #             matched_candidate_ref[:,0],matched_candidate_ref[:,1]) - \
                    #         __log_likelihood(junction_squiggle_mean, 
                    #             squiggle_match[cand][:,0],squiggle_match[cand][:,1])

                    #ax.plot(distinguish_points_idx,
                    #                np.repeat(max(junction_squiggle_mean)+0.5,len(distinguish_points_idx)), 'o')
                    ax.plot(matched_candidate_ref[:,0], 'r',linewidth=0.6, 
                            label = "minimap2 candidate")
                    ax.plot(matched_candidate_ref[:,0] + matched_candidate_ref[:,1], 'y--',linewidth=0.4, 
                            label = "sd")
                    ax.plot(matched_candidate_ref[:,0] - matched_candidate_ref[:,1], 'y--',linewidth=0.4)
                    ax.plot(squiggle_match[cand][:,0], 'g',linewidth=0.6,
                            label = "Candidate")
                    ax.plot(junction_squiggle_mean, 'b',linewidth=0.4, 
                            label = "junction squiggle")
                    ax.plot(junction_squiggle, 'b',linewidth=0.4, 
                            label = "junction squiggle")
                    
                    middle_pos = find_candidate_middle_pos(squiggle_match[index_m])
                    for i in range(len(aligned_base[index_m])):
                        ax.annotate(aligned_base[index_m][i],
                            xy=(middle_pos[i], max(junction_squiggle) + 0.05), xycoords='data', ha='center')
                    middle_pos = find_candidate_middle_pos(squiggle_match[cand])
                    for i in range(len(aligned_base[cand])):
                        ax.annotate(aligned_base[cand][i],
                            xy=(middle_pos[i], min(junction_squiggle) - 0.2), xycoords='data', ha='center')
                    ax.legend(frameon=False, loc='best')
                    
                    ax.set_title('Log LR: {:.2f}, post probability: {:.3f}, candidate_preference: {}, post probability(seq prior ratio = 9): {:.3f} '.format(
                        dist_seg_LR[cand], post_prob[cand], candidate_preference[cand], post_prob_prior[cand]), y=1.0, pad=-14)
                fig.suptitle('minimap2 candidate likelihood: {:.2f}, post probability: {:.3f}, candidate_preference: {}, post probability(seq prior ratio = 9): {:.3f}'.format(
                        dist_seg_LR[index_m], post_prob[index_m],candidate_preference[index_m], post_prob_prior[index_m]))

                fig.savefig("{}_mean_pairwise.png".format(read.qname))

                    # middle_pos = find_candidate_middle_pos(squiggle_match[0])
                    # for i in range(len(aligned_base[0])):
                    #     axes[0].annotate(aligned_base[0][i],
                    #         xy=(middle_pos[i], max(junction_squiggle_mean) + 0.2), xycoords='data', ha='center')
                    # middle_pos = find_candidate_middle_pos(squiggle_match[1])
                    # for i in range(len(aligned_base[1])):
                    #     axes[0].annotate(aligned_base[1][i],
                    #         xy=(middle_pos[i], min(junction_squiggle_mean) - 0.2), xycoords='data', ha='center')
                                    
                    # ax1 = axes[1]
                    # ax2 = ax1.twinx()
                    # # for j_i in range(num_of_cand):
                    # #     ax1.plot(squiggle_match[j_i][:,0], 'r',linewidth=0.2)
                    # # plt.savefig("{}_cum_LR_{}.png".format(output_prefix, read.qname))
                    # # plt.close()

                    # ax1.plot(squiggle_match[0][:,0], 'r',linewidth=0.4, 
                    #         label = "minimap2 candidate")
                    # ax1.plot(squiggle_match[1][:,0], 'g',linewidth=0.4, 
                    #         label = "NanoSplicer candidate")

                    # # # plot contribution of each data point in LR
                    # if True:
                    #     ax2.plot(LR_con, 
                    #             color='#4b0082', linewidth=0.4, alpha=1, 
                    #             dash_capstyle='round', label = 'LR contribution')
                    #     ax1.set_xlabel('Index')
                    #     ax1.set_ylabel('Candidate squiggle', color='g')
                    #     ax2.set_ylabel('log likelihood ratio contribution', color='#4b0082')
                    #     ax1.legend(frameon=False, loc='best')
                    #     fig.savefig("{}_mean_dtime8.png".format(read.qname))
                    #     plt.close()

################################# segment mean LR (pairwise version) LR contribution
            # LR plot
            if True and j == num_of_cand - 1:
                print(1)
                def mean_denoise(junction_squiggle, squiggle_match_list):
                    mean_junction_squiggle = np.array([])
                    event = junction_squiggle[0:1]
                    
                    for i in range(1, len(junction_squiggle)):
                        if all([(x[i] == x[i-1]).all() for x in squiggle_match_list]):
                            event = np.append(event, junction_squiggle[i])
                        else:
                            mean_junction_squiggle = np.append(mean_junction_squiggle,[np.mean(event)] * len(event))
                            event = np.array([junction_squiggle[i]])
                    mean_junction_squiggle = np.append(mean_junction_squiggle,[np.mean(event)] * len(event))
                    return(mean_junction_squiggle)

                # minimap2 candidate as reference
                matched_candidate_ref = squiggle_match[index_m]
                squiggle_match_list = [squiggle_match[x] for x in range(num_of_cand)]
                junction_squiggle_mean_list = \
                    [mean_denoise(junction_squiggle, [matched_candidate_ref, x]) for x in squiggle_match_list]
                segment_list = \
                    [get_distinguishing_segment(matched_candidate_ref, [x]) for x in squiggle_match_list]

                print(2)

                # LR calculation (distinguishing segment only)
                dist_seg_LR = []
                p_wise_LR_list = []
                
                for cand in range(num_of_cand):
                    junction_squiggle_mean = junction_squiggle_mean_list[cand]
                    segment, is_dist_seg = segment_list[cand]
                    is_dist_seg = is_dist_seg & (segment[:,1] - segment[:,0] >= minimum_point_for_dist_seg)
                    p_wise_LR = __log_likelihood(junction_squiggle_mean, 
                                    squiggle_match[cand][:,0],squiggle_match[cand][:,1]) -\
                                 __log_likelihood(junction_squiggle_mean, 
                                    matched_candidate_ref[:,0],matched_candidate_ref[:,1])
                    p_wise_LR_list.append(p_wise_LR)

                    dist_seg_LR.append(sum(p_wise_LR[segment[is_dist_seg][:,0]]))
                
                rel_to_ref_LR = np.exp(dist_seg_LR)/np.exp(dist_seg_LR[index_m])
                post_prob = rel_to_ref_LR/sum(rel_to_ref_LR)
                post_prob_prior = post_prob * (prior_ratio **candidate_preference)
                post_prob_prior = post_prob_prior/sum(post_prob_prior)
                
                print(3)
                fig = plt.figure(figsize=(20,5 * num_of_cand))
                subplot_index = 0
                for cand in range(num_of_cand):
                    junction_squiggle_mean = junction_squiggle_mean_list[cand]
                    segment, is_dist_seg = segment_list[cand]
                    is_dist_seg = is_dist_seg & (segment[:,1] - segment[:,0] >= minimum_point_for_dist_seg)

                    if cand == index_m:
                        continue
                    
                    # frame subplot
                    subplot_index += 1
                    ax = fig.add_subplot(num_of_cand -1,1,subplot_index)
                    for i, seg in enumerate(segment):
                        ax.axvspan(seg[0], seg[1]-1, -0.2,1.2, 
                            clip_on=True, alpha = 0.2 if is_dist_seg[i] else 0.1 , edgecolor = 'black', 
                            lw = 0.2,  facecolor = 'red' if is_dist_seg[i] else '#E2D3CA')

                    # # data-point wise LR

                    # LR_con = __log_likelihood(junction_squiggle_mean, 
                    #             matched_candidate_ref[:,0],matched_candidate_ref[:,1]) - \
                    #         __log_likelihood(junction_squiggle_mean, 
                    #             squiggle_match[cand][:,0],squiggle_match[cand][:,1])

                    #ax.plot(distinguish_points_idx,
                    #                np.repeat(max(junction_squiggle_mean)+0.5,len(distinguish_points_idx)), 'o')
                    ax.plot(matched_candidate_ref[:,0], 'r',linewidth=0.6, 
                            label = "minimap2 candidate")
                    ax.plot(matched_candidate_ref[:,0] + matched_candidate_ref[:,1], 'y--',linewidth=0.4, 
                            label = "sd")
                    ax.plot(matched_candidate_ref[:,0] - matched_candidate_ref[:,1], 'y--',linewidth=0.4)
                    ax.plot(squiggle_match[cand][:,0], 'g',linewidth=0.6,
                            label = "Candidate")

                    if True:
                        ax2 = ax.twinx()
                        ax2.plot(p_wise_LR_list[cand], 
                                color='#4b0082', linewidth=0.4, alpha=1, 
                                dash_capstyle='round', label = 'LR contribution')
                        ax2.set_ylabel('log likelihood ratio contribution', color='#4b0082')
                        ax2.axhline(0, linestyle = '--', lw = 0.6)
                    
                    middle_pos = find_candidate_middle_pos(squiggle_match[index_m])
                    for i in range(len(aligned_base[index_m])):
                        ax.annotate(aligned_base[index_m][i],
                            xy=(middle_pos[i], max(junction_squiggle) + 0.05), xycoords='data', ha='center')
                    middle_pos = find_candidate_middle_pos(squiggle_match[cand])
                    for i in range(len(aligned_base[cand])):
                        ax.annotate(aligned_base[cand][i],
                            xy=(middle_pos[i], min(junction_squiggle) - 0.2), xycoords='data', ha='center')
                    ax.legend(frameon=False, loc='best')
                    
                    ax.set_title('Log LR: {:.2f}, post probability: {:.3f}, candidate_preference: {}, post probability(seq prior ratio = 9): {:.3f} '.format(
                        dist_seg_LR[cand], post_prob[cand], candidate_preference[cand], post_prob_prior[cand]), y=1.0, pad=-14)
                fig.suptitle('minimap2 candidate likelihood: {:.2f}, post probability: {:.3f}, candidate_preference: {}, post probability(seq prior ratio = 9): {:.3f}'.format(
                        dist_seg_LR[index_m], post_prob[index_m],candidate_preference[index_m], post_prob_prior[index_m]))

                fig.savefig("{}_LR_mean_pairwise.png".format(read.qname))

################################## median denoising (all candidate seg)
            # LR plot
            if True and j == num_of_cand - 1:
                print(1)
                def median_denoise(junction_squiggle, squiggle_match_list):
                    median_junction_squiggle = np.array([])
                    event = junction_squiggle[0:1]
                    
                    for i in range(1, len(junction_squiggle)):
                        if all([(x[i] == x[i-1]).all() for x in squiggle_match_list]):
                            event = np.append(event, junction_squiggle[i])
                        else:
                            median_junction_squiggle = np.append(median_junction_squiggle,[np.median(event)] * len(event))
                            event = np.array([junction_squiggle[i]])
                    median_junction_squiggle = np.append(median_junction_squiggle,[np.median(event)] * len(event))
                    return(median_junction_squiggle)

                # minimap2 candidate as reference
                matched_candidate_ref = squiggle_match[index_m]
                squiggle_match_list = [squiggle_match[x] for x in range(num_of_cand)]
                junction_squiggle_median = median_denoise(junction_squiggle, squiggle_match_list)
                segment, is_dist_seg = \
                    get_distinguishing_segment(matched_candidate_ref, squiggle_match_list)
                
                is_dist_seg = is_dist_seg & (segment[:,1] - segment[:,0] >= minimum_point_for_dist_seg)
                print(2)

                # LR calculation (distinguishing segment only)
                dist_seg_LR = []
                p_wise_LR_list = []
                for cand in range(num_of_cand):
                    p_wise_LR = __log_likelihood(junction_squiggle_median, 
                                    squiggle_match[cand][:,0],squiggle_match[cand][:,1])
                    dist_seg_LR.append(sum(p_wise_LR[segment[is_dist_seg][:,0]]))
                    p_wise_LR_list.append(p_wise_LR)
                post_prob = np.exp(dist_seg_LR)/sum(np.exp(dist_seg_LR))
                post_prob_prior = post_prob * (prior_ratio **candidate_preference)
                post_prob_prior = post_prob_prior/sum(post_prob_prior)
                print(3)
                fig = plt.figure(figsize=(20,5 * num_of_cand))
                subplot_index = 0
                for cand in range(num_of_cand):
                    if cand == index_m:
                        continue
                    
                    # frame subplot
                    subplot_index += 1
                    ax = fig.add_subplot(num_of_cand -1,1,subplot_index)
                    for i, seg in enumerate(segment):
                        ax.axvspan(seg[0], seg[1]-1, -0.2,1.2, 
                            clip_on=True, alpha = 0.2 if is_dist_seg[i] else 0.1 , edgecolor = 'black', 
                            lw = 0.2,  facecolor = 'red' if is_dist_seg[i] else '#E2D3CA')

                    # # data-point wise LR

                    # LR_con = __log_likelihood(junction_squiggle_median, 
                    #             matched_candidate_ref[:,0],matched_candidate_ref[:,1]) - \
                    #         __log_likelihood(junction_squiggle_median, 
                    #             squiggle_match[cand][:,0],squiggle_match[cand][:,1])

                    #ax.plot(distinguish_points_idx,
                    #                np.repeat(max(junction_squiggle_median)+0.5,len(distinguish_points_idx)), 'o')
                    ax.plot(matched_candidate_ref[:,0], 'r',linewidth=0.6, 
                            label = "minimap2 candidate")
                    ax.plot(matched_candidate_ref[:,0] + matched_candidate_ref[:,1], 'y--',linewidth=0.4, 
                            label = "sd")
                    ax.plot(matched_candidate_ref[:,0] - matched_candidate_ref[:,1], 'y--',linewidth=0.4)
                    ax.plot(squiggle_match[cand][:,0], 'g',linewidth=0.6,
                            label = "Candidate")
                    ax.plot(junction_squiggle_median, 'b',linewidth=0.4, 
                            label = "junction squiggle")
                    ax.plot(junction_squiggle, 'b',linewidth=0.4, 
                            label = "junction squiggle")

                    middle_pos = find_candidate_middle_pos(squiggle_match[index_m])
                    for i in range(len(aligned_base[index_m])):
                        ax.annotate(aligned_base[index_m][i],
                            xy=(middle_pos[i], max(junction_squiggle) + 0.05), xycoords='data', ha='center')
                    middle_pos = find_candidate_middle_pos(squiggle_match[cand])
                    for i in range(len(aligned_base[cand])):
                        ax.annotate(aligned_base[cand][i],
                            xy=(middle_pos[i], min(junction_squiggle) - 0.2), xycoords='data', ha='center')
                
                    ax.legend(frameon=False, loc='best')
                    
                    ax.set_title('Log LR: {:.2f}, post probability: {:.3f}, candidate_preference: {}, post probability(seq prior ratio = 9): {:.3f} '.format(
                        dist_seg_LR[cand] - dist_seg_LR[index_m], post_prob[cand], candidate_preference[cand], post_prob_prior[cand]), y=1.0, pad=-14)
                
                fig.suptitle('minimap2 candidate likelihood: {:.2f}, post probability: {:.3f}, candidate_preference: {}, post probability(seq prior ratio = 9): {:.3f}'.format(
                        dist_seg_LR[index_m] - dist_seg_LR[index_m], post_prob[index_m],candidate_preference[index_m], post_prob_prior[index_m]))
                fig.savefig("{}_median_all_cand.png".format(read.qname))

################################## median denoising (all candidate seg) LR contri
            # LR plot
            if True and j == num_of_cand - 1:
                print(1)
                def median_denoise(junction_squiggle, squiggle_match_list):
                    median_junction_squiggle = np.array([])
                    event = junction_squiggle[0:1]
                    
                    for i in range(1, len(junction_squiggle)):
                        if all([(x[i] == x[i-1]).all() for x in squiggle_match_list]):
                            event = np.append(event, junction_squiggle[i])
                        else:
                            median_junction_squiggle = np.append(median_junction_squiggle,[np.median(event)] * len(event))
                            event = np.array([junction_squiggle[i]])
                    median_junction_squiggle = np.append(median_junction_squiggle,[np.median(event)] * len(event))
                    return(median_junction_squiggle)

                # minimap2 candidate as reference
                matched_candidate_ref = squiggle_match[index_m]
                squiggle_match_list = [squiggle_match[x] for x in range(num_of_cand)]
                junction_squiggle_median = median_denoise(junction_squiggle, squiggle_match_list)
                segment, is_dist_seg = \
                    get_distinguishing_segment(matched_candidate_ref, squiggle_match_list)
                
                is_dist_seg = is_dist_seg & (segment[:,1] - segment[:,0] >= minimum_point_for_dist_seg)
                print(2)

                # LR calculation (distinguishing segment only)
                dist_seg_LR = []
                p_wise_LR_list = []
                for cand in range(num_of_cand):
                    p_wise_LR = __log_likelihood(junction_squiggle_median, 
                                    squiggle_match[cand][:,0],squiggle_match[cand][:,1])
                    dist_seg_LR.append(sum(p_wise_LR[segment[is_dist_seg][:,0]]))
                    p_wise_LR_list.append(p_wise_LR)
                post_prob = np.exp(dist_seg_LR)/sum(np.exp(dist_seg_LR))
                post_prob_prior = post_prob * (prior_ratio **candidate_preference)
                post_prob_prior = post_prob_prior/sum(post_prob_prior)
                
                print(3)
                fig = plt.figure(figsize=(20,5 * num_of_cand))
                subplot_index = 0
                for cand in range(num_of_cand):
                    if cand == index_m:
                        continue
                    
                    # frame subplot
                    subplot_index += 1
                    ax = fig.add_subplot(num_of_cand -1,1,subplot_index)
                    for i, seg in enumerate(segment):
                        ax.axvspan(seg[0], seg[1]-1, -0.2,1.2, 
                            clip_on=True, alpha = 0.2 if is_dist_seg[i] else 0.1 , edgecolor = 'black', 
                            lw = 0.2,  facecolor = 'red' if is_dist_seg[i] else '#E2D3CA')

                    # # data-point wise LR

                    # LR_con = __log_likelihood(junction_squiggle_median, 
                    #             matched_candidate_ref[:,0],matched_candidate_ref[:,1]) - \
                    #         __log_likelihood(junction_squiggle_median, 
                    #             squiggle_match[cand][:,0],squiggle_match[cand][:,1])

                    #ax.plot(distinguish_points_idx,
                    #                np.repeat(max(junction_squiggle_median)+0.5,len(distinguish_points_idx)), 'o')
                    ax.plot(matched_candidate_ref[:,0], 'r',linewidth=0.6, 
                            label = "minimap2 candidate")
                    ax.plot(matched_candidate_ref[:,0] + matched_candidate_ref[:,1], 'y--',linewidth=0.4, 
                            label = "sd")
                    ax.plot(matched_candidate_ref[:,0] - matched_candidate_ref[:,1], 'y--',linewidth=0.4)
                    ax.plot(squiggle_match[cand][:,0], 'g',linewidth=0.6,
                            label = "Candidate")

                    if True:
                        ax2 = ax.twinx()
                        ax2.plot(p_wise_LR_list[cand] - p_wise_LR_list[index_m], 
                                color='#4b0082', linewidth=0.4, alpha=1, 
                                dash_capstyle='round', label = 'LR contribution')
                        ax2.set_ylabel('log likelihood ratio contribution', color='#4b0082')
                        ax2.axhline(0, linestyle = '--', lw = 0.6)

                    middle_pos = find_candidate_middle_pos(squiggle_match[index_m])
                    for i in range(len(aligned_base[index_m])):
                        ax.annotate(aligned_base[index_m][i],
                            xy=(middle_pos[i], max(junction_squiggle) + 0.05), xycoords='data', ha='center')
                    middle_pos = find_candidate_middle_pos(squiggle_match[cand])
                    for i in range(len(aligned_base[cand])):
                        ax.annotate(aligned_base[cand][i],
                            xy=(middle_pos[i], min(junction_squiggle) - 0.2), xycoords='data', ha='center')
                
                    ax.legend(frameon=False, loc='best')
                    
                    ax.set_title('Log LR: {:.2f}, post probability: {:.3f}, candidate_preference: {}, post probability(seq prior ratio = 9): {:.3f} '.format(
                        dist_seg_LR[cand] - dist_seg_LR[index_m], post_prob[cand], candidate_preference[cand], post_prob_prior[cand]), y=1.0, pad=-14)
                
                fig.suptitle('minimap2 candidate likelihood: {:.2f}, post probability: {:.3f}, candidate_preference: {}, post probability(seq prior ratio = 9): {:.3f}'.format(
                        dist_seg_LR[index_m] - dist_seg_LR[index_m], post_prob[index_m],candidate_preference[index_m], post_prob_prior[index_m]))
                fig.savefig("{}_LR_median_all_cand.png".format(read.qname))

################################# segment median (pairwise version)
            # LR plot
            if True and j == num_of_cand - 1:
                print(1)
                def median_denoise(junction_squiggle, squiggle_match_list):
                    median_junction_squiggle = np.array([])
                    event = junction_squiggle[0:1]
                    
                    for i in range(1, len(junction_squiggle)):
                        if all([(x[i] == x[i-1]).all() for x in squiggle_match_list]):
                            event = np.append(event, junction_squiggle[i])
                        else:
                            median_junction_squiggle = np.append(median_junction_squiggle,[np.median(event)] * len(event))
                            event = np.array([junction_squiggle[i]])
                    median_junction_squiggle = np.append(median_junction_squiggle,[np.median(event)] * len(event))
                    return(median_junction_squiggle)

                # minimap2 candidate as reference
                matched_candidate_ref = squiggle_match[index_m]
                squiggle_match_list = [squiggle_match[x] for x in range(num_of_cand)]
                junction_squiggle_median_list = \
                    [median_denoise(junction_squiggle, [matched_candidate_ref, x]) for x in squiggle_match_list]
                segment_list = \
                    [get_distinguishing_segment(matched_candidate_ref, [x]) for x in squiggle_match_list]
                
                
                print(2)

                # LR calculation (distinguishing segment only)
                dist_seg_LR = []
                
                for cand in range(num_of_cand):
                    junction_squiggle_median = junction_squiggle_median_list[cand]
                    segment, is_dist_seg = segment_list[cand]
                    is_dist_seg = is_dist_seg & (segment[:,1] - segment[:,0] >= minimum_point_for_dist_seg)
                    p_wise_LR = __log_likelihood(junction_squiggle_median, 
                                    squiggle_match[cand][:,0],squiggle_match[cand][:,1]) -\
                                 __log_likelihood(junction_squiggle_median, 
                                    matched_candidate_ref[:,0],matched_candidate_ref[:,1])

                    dist_seg_LR.append(sum(p_wise_LR[segment[is_dist_seg][:,0]]))
                
                rel_to_ref_LR = np.exp(dist_seg_LR)/np.exp(dist_seg_LR[index_m])
                post_prob = rel_to_ref_LR/sum(rel_to_ref_LR)
                post_prob_prior = post_prob * (prior_ratio **candidate_preference)
                post_prob_prior = post_prob_prior/sum(post_prob_prior)
                print(3)
                fig = plt.figure(figsize=(20,5 * num_of_cand))
                subplot_index = 0
                for cand in range(num_of_cand):
                    junction_squiggle_median = junction_squiggle_median_list[cand]
                    segment, is_dist_seg = segment_list[cand]
                    is_dist_seg = is_dist_seg & (segment[:,1] - segment[:,0] >= minimum_point_for_dist_seg)

                    if cand == index_m:
                        continue
                    
                    # frame subplot
                    subplot_index += 1
                    ax = fig.add_subplot(num_of_cand -1,1,subplot_index)
                    for i, seg in enumerate(segment):
                        ax.axvspan(seg[0], seg[1]-1, -0.2,1.2, 
                            clip_on=True, alpha = 0.2 if is_dist_seg[i] else 0.1 , edgecolor = 'black', 
                            lw = 0.2,  facecolor = 'red' if is_dist_seg[i] else '#E2D3CA')

                    # # data-point wise LR

                    # LR_con = __log_likelihood(junction_squiggle_median, 
                    #             matched_candidate_ref[:,0],matched_candidate_ref[:,1]) - \
                    #         __log_likelihood(junction_squiggle_median, 
                    #             squiggle_match[cand][:,0],squiggle_match[cand][:,1])

                    #ax.plot(distinguish_points_idx,
                    #                np.repeat(max(junction_squiggle_median)+0.5,len(distinguish_points_idx)), 'o')
                    ax.plot(matched_candidate_ref[:,0], 'r',linewidth=0.6, 
                            label = "minimap2 candidate")
                    ax.plot(matched_candidate_ref[:,0] + matched_candidate_ref[:,1], 'y--',linewidth=0.4, 
                            label = "sd")
                    ax.plot(matched_candidate_ref[:,0] - matched_candidate_ref[:,1], 'y--',linewidth=0.4)
                    ax.plot(squiggle_match[cand][:,0], 'g',linewidth=0.6,
                            label = "Candidate")
                    ax.plot(junction_squiggle_median, 'b',linewidth=0.4, 
                            label = "junction squiggle")
                    ax.plot(junction_squiggle, 'b',linewidth=0.4, 
                            label = "junction squiggle")
                    
                    middle_pos = find_candidate_middle_pos(squiggle_match[index_m])
                    for i in range(len(aligned_base[index_m])):
                        ax.annotate(aligned_base[index_m][i],
                            xy=(middle_pos[i], max(junction_squiggle) + 0.05), xycoords='data', ha='center')
                    middle_pos = find_candidate_middle_pos(squiggle_match[cand])
                    for i in range(len(aligned_base[cand])):
                        ax.annotate(aligned_base[cand][i],
                            xy=(middle_pos[i], min(junction_squiggle) - 0.2), xycoords='data', ha='center')
                    ax.legend(frameon=False, loc='best')
                    
                    ax.set_title('Log LR: {:.2f}, post probability: {:.3f}, candidate_preference: {}, post probability(seq prior ratio = 9): {:.3f} '.format(
                        dist_seg_LR[cand], post_prob[cand], candidate_preference[cand], post_prob_prior[cand]), y=1.0, pad=-14)
                fig.suptitle('minimap2 candidate likelihood: {:.2f}, post probability: {:.3f}, candidate_preference: {}, post probability(seq prior ratio = 9): {:.3f}'.format(
                        dist_seg_LR[index_m], post_prob[index_m],candidate_preference[index_m], post_prob_prior[index_m]))

                fig.savefig("{}_median_pairwise.png".format(read.qname))

                    # middle_pos = find_candidate_middle_pos(squiggle_match[0])
                    # for i in range(len(aligned_base[0])):
                    #     axes[0].annotate(aligned_base[0][i],
                    #         xy=(middle_pos[i], max(junction_squiggle_median) + 0.2), xycoords='data', ha='center')
                    # middle_pos = find_candidate_middle_pos(squiggle_match[1])
                    # for i in range(len(aligned_base[1])):
                    #     axes[0].annotate(aligned_base[1][i],
                    #         xy=(middle_pos[i], min(junction_squiggle_median) - 0.2), xycoords='data', ha='center')
                                    
                    # ax1 = axes[1]
                    # ax2 = ax1.twinx()
                    # # for j_i in range(num_of_cand):
                    # #     ax1.plot(squiggle_match[j_i][:,0], 'r',linewidth=0.2)
                    # # plt.savefig("{}_cum_LR_{}.png".format(output_prefix, read.qname))
                    # # plt.close()

                    # ax1.plot(squiggle_match[0][:,0], 'r',linewidth=0.4, 
                    #         label = "minimap2 candidate")
                    # ax1.plot(squiggle_match[1][:,0], 'g',linewidth=0.4, 
                    #         label = "NanoSplicer candidate")

                    # # # plot contribution of each data point in LR
                    # if True:
                    #     ax2.plot(LR_con, 
                    #             color='#4b0082', linewidth=0.4, alpha=1, 
                    #             dash_capstyle='round', label = 'LR contribution')
                    #     ax1.set_xlabel('Index')
                    #     ax1.set_ylabel('Candidate squiggle', color='g')
                    #     ax2.set_ylabel('log likelihood ratio contribution', color='#4b0082')
                    #     ax1.legend(frameon=False, loc='best')
                    #     fig.savefig("{}_median_dtime8.png".format(read.qname))
                    #     plt.close()

################################# segment median LR (pairwise version) LR contribution
            # LR plot
            if True and j == num_of_cand - 1:
                print(1)
                def median_denoise(junction_squiggle, squiggle_match_list):
                    median_junction_squiggle = np.array([])
                    event = junction_squiggle[0:1]
                    
                    for i in range(1, len(junction_squiggle)):
                        if all([(x[i] == x[i-1]).all() for x in squiggle_match_list]):
                            event = np.append(event, junction_squiggle[i])
                        else:
                            median_junction_squiggle = np.append(median_junction_squiggle,[np.median(event)] * len(event))
                            event = np.array([junction_squiggle[i]])
                    median_junction_squiggle = np.append(median_junction_squiggle,[np.median(event)] * len(event))
                    return(median_junction_squiggle)

                # minimap2 candidate as reference
                matched_candidate_ref = squiggle_match[index_m]
                squiggle_match_list = [squiggle_match[x] for x in range(num_of_cand)]
                junction_squiggle_median_list = \
                    [median_denoise(junction_squiggle, [matched_candidate_ref, x]) for x in squiggle_match_list]
                segment_list = \
                    [get_distinguishing_segment(matched_candidate_ref, [x]) for x in squiggle_match_list]

                print(2)

                # LR calculation (distinguishing segment only)
                dist_seg_LR = []
                p_wise_LR_list = []
                
                for cand in range(num_of_cand):
                    junction_squiggle_median = junction_squiggle_median_list[cand]
                    segment, is_dist_seg = segment_list[cand]
                    is_dist_seg = is_dist_seg & (segment[:,1] - segment[:,0] >= minimum_point_for_dist_seg)
                    p_wise_LR = __log_likelihood(junction_squiggle_median, 
                                    squiggle_match[cand][:,0],squiggle_match[cand][:,1]) -\
                                 __log_likelihood(junction_squiggle_median, 
                                    matched_candidate_ref[:,0],matched_candidate_ref[:,1])
                    p_wise_LR_list.append(p_wise_LR)

                    dist_seg_LR.append(sum(p_wise_LR[segment[is_dist_seg][:,0]]))
                
                rel_to_ref_LR = np.exp(dist_seg_LR)/np.exp(dist_seg_LR[index_m])
                post_prob = rel_to_ref_LR/sum(rel_to_ref_LR)
                post_prob_prior = post_prob * (prior_ratio **candidate_preference)
                post_prob_prior = post_prob_prior/sum(post_prob_prior)
                print(3)
                fig = plt.figure(figsize=(20,5 * num_of_cand))
                subplot_index = 0
                for cand in range(num_of_cand):
                    junction_squiggle_median = junction_squiggle_median_list[cand]
                    segment, is_dist_seg = segment_list[cand]
                    is_dist_seg = is_dist_seg & (segment[:,1] - segment[:,0] >= minimum_point_for_dist_seg)

                    if cand == index_m:
                        continue
                    
                    # frame subplot
                    subplot_index += 1
                    ax = fig.add_subplot(num_of_cand -1,1,subplot_index)
                    for i, seg in enumerate(segment):
                        ax.axvspan(seg[0], seg[1]-1, -0.2,1.2, 
                            clip_on=True, alpha = 0.2 if is_dist_seg[i] else 0.1 , edgecolor = 'black', 
                            lw = 0.2,  facecolor = 'red' if is_dist_seg[i] else '#E2D3CA')

                    # # data-point wise LR

                    # LR_con = __log_likelihood(junction_squiggle_median, 
                    #             matched_candidate_ref[:,0],matched_candidate_ref[:,1]) - \
                    #         __log_likelihood(junction_squiggle_median, 
                    #             squiggle_match[cand][:,0],squiggle_match[cand][:,1])

                    #ax.plot(distinguish_points_idx,
                    #                np.repeat(max(junction_squiggle_median)+0.5,len(distinguish_points_idx)), 'o')
                    ax.plot(matched_candidate_ref[:,0], 'r',linewidth=0.6, 
                            label = "minimap2 candidate")
                    ax.plot(matched_candidate_ref[:,0] + matched_candidate_ref[:,1], 'y--',linewidth=0.4, 
                            label = "sd")
                    ax.plot(matched_candidate_ref[:,0] - matched_candidate_ref[:,1], 'y--',linewidth=0.4)
                    ax.plot(squiggle_match[cand][:,0], 'g',linewidth=0.6,
                            label = "Candidate")

                    if True:
                        ax2 = ax.twinx()
                        ax2.plot(p_wise_LR_list[cand], 
                                color='#4b0082', linewidth=0.4, alpha=1, 
                                dash_capstyle='round', label = 'LR contribution')
                        ax2.set_ylabel('log likelihood ratio contribution', color='#4b0082')
                        ax2.axhline(0, linestyle = '--', lw = 0.6)
                    
                    middle_pos = find_candidate_middle_pos(squiggle_match[index_m])
                    for i in range(len(aligned_base[index_m])):
                        ax.annotate(aligned_base[index_m][i],
                            xy=(middle_pos[i], max(junction_squiggle) + 0.05), xycoords='data', ha='center')
                    middle_pos = find_candidate_middle_pos(squiggle_match[cand])
                    for i in range(len(aligned_base[cand])):
                        ax.annotate(aligned_base[cand][i],
                            xy=(middle_pos[i], min(junction_squiggle) - 0.2), xycoords='data', ha='center')
                    ax.legend(frameon=False, loc='best')
                    
                    ax.set_title('Log LR: {:.2f}, post probability: {:.3f}, candidate_preference: {}, post probability(seq prior ratio = 9): {:.3f} '.format(
                        dist_seg_LR[cand], post_prob[cand], candidate_preference[cand], post_prob_prior[cand]), y=1.0, pad=-14)
                fig.suptitle('minimap2 candidate likelihood: {:.2f}, post probability: {:.3f}, candidate_preference: {}, post probability(seq prior ratio = 9): {:.3f}'.format(
                        dist_seg_LR[index_m], post_prob[index_m],candidate_preference[index_m], post_prob_prior[index_m]))

                fig.savefig("{}_LR_median_pairwise.png".format(read.qname))

##################################
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
                
        
            if SAVE_DATA:
                np.savetxt("{}_candidate{}_{}.csv".format(output_prefix, read.qname,j), squiggle_match[j], delimiter=",")
                np.savetxt("junction_squiggle{}.csv".format(read.qname), junction_squiggle, delimiter=",")
            # Plot DTW 
            if False:
                print(score)
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
        
        
        
        
        ###########################################
        
        # score_output = []
        # score_trimmed = []
        # score_trimmed2 = []
        # score_trimmed_ref = []
        # cum_path = {}
        # squiggle_match = {}
        # element_wise_logL = {}
        # output_suffix = '_c4'
        # num_of_cand = len(candidates_for_read)
        
        # # loop over candidate for each junction within read
        # for j, candidiate in enumerate(candidates_for_read):
        #     print(candidate)
        #     candidate_squiggle = np.array(model_dic[candidiate],float)
            
        #     # dtw
        #     path , score, cum_matrix = \
        #         dtw_local_alignment(candidate_squiggle=candidate_squiggle, 
        #                             junction_squiggle = junction_squiggle)
            
        #     # candidate squiggle in matched len of junction suqiggle
        #     squiggle_match[j] = candidate_squiggle[path[:,1] - 1, :]
        #     cum_path[j] = cum_matrix[path[:, 0], path[:, 1]]            
        #     element_wise_logL[j] = np.append(cum_path[j][0],
        #                                     cum_path[j][1:] - cum_path[j][:-1])
            
        #     # plotting 
        #     if False:
        #         fig, axes = plt.subplots(nrows=4, figsize=(30,40))
        #         axes[0].plot(junction_squiggle,linewidth = 10)
        #         axes[0].tick_params(labelsize=40)
        #                     #axes[1].figure(figsize=(10,5))
        #         axes[1].plot(candidate_squiggle[:,0], color = "orange",linewidth =7)
        #         axes[1].plot(candidate_squiggle[:,0] + candidate_squiggle[:,1], ":",color = "orange",linewidth = 4)
        #         axes[1].plot(candidate_squiggle[:,0] - candidate_squiggle[:,1], ":",color = "orange",linewidth = 4)
        #         axes[1].tick_params(labelsize=40)
        #         axes[1].set_title("Scrappie model",fontsize=30, pad = 1)
        #         axes[2].plot(junction_squiggle,linewidth = 10)
        #         axes[2].tick_params(labelsize=40)
        #         #path = np.array(path[::-1])
        #         axes[2].plot(path[:,0]-1, candidate_squiggle[[path[:,1]-1]][:,0],'',color = "orange",linewidth = 7)
        #         axes[2].plot(path[:,0]-1, candidate_squiggle[[path[:,1]-1]][:,0]\
        #                                         + candidate_squiggle[[path[:,1]-1]][:,1],':',color = "orange",linewidth = 4)
        #         axes[2].plot(path[:,0]-1, candidate_squiggle[[path[:,1]-1]][:,0]\
        #                                         - candidate_squiggle[[path[:,1]-1]][:,1],':',color = "orange",linewidth = 4)
        #         axes[2].set_title(str(i)+"\nDist: {:.2f}, path length: {}, Adjusted dist: {:.2f}".format(score,len(path),score/len(path)),fontsize=40, pad = 1)
        #         pos = axes[3].imshow(cum_matrix, cmap='hot', interpolation='nearest',aspect='auto')
        #         fig.colorbar(pos,ax = axes[3])
        #         axes[3].plot(path[:,1], path[:,0])
        #         axes[3].set_title("Alignment path",fontsize=30, pad = 20)
        #         if read.is_reverse:
        #             strand = 1
        #         else:
        #             strand = 0
        #         fig.savefig("fig{}_{}.png".format('555f59a2-784b-4fbc-8392-6ff32fd5dcb4', i))

        # '''
        # number of points contribute to LR
        # '''   
        # #use best matched squiggle as reference
        # ref_cand = squiggle_match[np.argmin(score_output)]
        # ref_element_wise_logL = element_wise_logL[np.argmin(score_output)]
        # out_of_sd = np.zeros(len(junction_squiggle) ,dtype = int)
        # out_of_sd_ind = {}
        # for j in range(num_of_cand):
        #     out_of_sd_ind[j] = np.zeros(len(junction_squiggle) ,dtype = int)

        
        # for x in range(len(junction_squiggle)):
        #     for j in range(num_of_cand):
        #         if squiggle_match[j][x,0] - ref_cand[x,0] > ref_cand[x,1]:
        #             out_of_sd[x] = 1
        #             out_of_sd_ind[j][x] = 1
        
        # for j in range(num_of_cand):
        #     score_trimmed.append(sum(out_of_sd * element_wise_logL[j]))
        #     score_trimmed2.append(sum(out_of_sd_ind[j] * element_wise_logL[j]))
        #     score_trimmed_ref.append(sum(out_of_sd_ind[j] * ref_element_wise_logL))

def get_gaps_in_read(AlignedSegment):
    blocks = AlignedSegment.get_blocks()
    gap = set([(blocks[i-1][1], blocks[i][0]) for i in range(len(blocks))])
    return gap

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

if __name__ == "__main__":
    main()