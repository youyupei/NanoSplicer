import sys
import h5py
import getopt
import timeit
import os
import numpy as np
import re


from math import sqrt

#sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import helper
from helper import expect_squiggle_dict
from helper import parse_candidate_file
import junction_squiggle_selection

from dtw import dtw_local_alignment_max_mean as dtw_mean
#from dtw_cy.dtw import dtw_local_alignment_max_sum as dtw_sum
from dtw import dtw_local_alignment_max_sum as dtw_sum
from dtw import dtw_global_alignment_max_sum as dtw_global
from dtw import dtw_local_alignment_max_sum_band_flipped as dtw_band
from dtw import dist_to_likelihood, dist_to_likelihood_flipped_time_serie,dist_to_likelihood_flipped_new_path,dist_to_likelihood_flipped

def parse_arg():
    def print_help():
        print("\n\nUsage: python {} [OPTIONS] <fast5 filename>".format(argv[0]))
        print("Options:\n\tIndexing:")
        print('\t\t-c INT\tcandidate file name (required)')
        print("\t\t-o INT\toutput csv name, default: 'Untitled'")
        print('\t\t-T \tNumber of events trimmed from scrappie model')
        print('\t\t-t \tNumber of samples trimmed from raw signal')
     
        return None

    argv = sys.argv
    if len(argv) <= 2:     
        print_help()       # print help doc when no command line args provided
        sys.exit(0)
    
    try: 
        opts, args = getopt.getopt(argv[1:],"ho:c:T:t:ab:s:",
                    ["help=","output_csv=", "candidate_file=",\
                    "trim_model=","trim_signal=","dtw_adj","bandwidth=","SAM="])
    except getopt.GetoptError:
        print_help()
        sys.exit(0)

    output_file = "Untiled"
    try:
        fast5_filename = args[0]
    except:
        print("InputError: missing fast5 file name!")
        sys.exit(0)

    # DEFAULT VALUE
    trim_signal = 0
    trim_model = 4
    dtw_adj = False
    bandwidth = False

    for opt, arg in opts:
        if opt in ('-h', "--help"):
            print_help()
            sys.exit(0)
        elif opt in ("-o", "--output_csv"):
            output_file = arg
        elif opt in ("-c", "--candidate_file"):
           candidate_file = arg
        elif opt in ("-s", "--SAM"):
           sam_line = arg
        elif opt in ("-T", "--trim_model"):
           trim_model = int(arg)
        elif opt in ("-t", "--trim_signal"):
           trim_signal = int(arg)
        elif opt in ("-a", "--dtw_adj"):
           dtw_adj = True
        elif opt in ("-b", "--bandwidth"):
           bandwidth = float(arg)


    # choose the version of dtw (sum: minimize the sum of cost for the path)
    dtw_local_alignment = dtw_mean if dtw_adj else dtw_sum
    
    if bandwidth:
        def dtw_local_alignment(long, short, dist_type = None, upper = np.inf):
            return dtw_band(long=long, short=short, band_prop = bandwidth,\
            dist_type = dist_type, upper = upper)
    
    return fast5_filename, output_file, candidate_file,\
     dtw_local_alignment, trim_model, trim_signal, sam_line

def main():
    fast5_filename, output_file, candidate_file, \
        dtw_local_alignment, trim_model, trim_signal, sam_line = parse_arg()
    
    # get read id
    #########################################################
    #read_id  = helper.Fast5Class(fast5_filename).get_read_id()
    #########################################################
   
    outf = open(output_file,'w')
    '''
    outf.write("read id")
    for i in range(len(sequence)):
        outf.write(",z_score{},log_likelihood{},manhattan{}".format(i+1,i+1,i+1))
    
    outf.write('\n'+ fast5_filename)
    '''
    candidates = parse_candidate_file(candidate_file)

    # read line in sam file
#
#   print(sam_line)
    sam_flag, mapped_pos0, cigar = [sam_line.strip().split('\t')[index] for index in [1,3,5]]
    mapped_pos0 = int(mapped_pos0)
    # 1-based to 0-based
    mapped_pos0 -= 1

    # determine strand
    if sam_flag == "0":
        strand = "+"
    elif sam_flag == "16":
        strand = "-"
    else:
        print("Error: Abnormal mapped reads.")
        sys.exit(0)


    # tombo resquiggle
    try:
        tombo_results, tombo_start_clip, tombo_end_clip = \
            junction_squiggle_selection.tombo_squiggle_to_basecalls(fast5_filename)
    except:
        print("tombo resquiggle failed!!")
        sys.exit(0)
    
    read_length = len(tombo_results.genome_seq) + tombo_start_clip + tombo_end_clip
    normalised_raw_signal = tombo_results.raw_signal/1.4826
    
    # genome pos to read pos mapping vector
    g_r_mapping = junction_squiggle_selection.genome_to_read_pos_conversion(cigar)

    for candidate in candidates:
        outf.write(fast5_filename + ',' + strand)

        # take reverse compliment seq if nesseccary
        if strand == "-":
            candidate.sequences = [helper.reverse_complement(s)
             for s in candidate.sequences]

        # trim signal
        candidate.start += trim_signal
#   
#        
#        print(candidate.start, candidate.end, mapped_pos0, len(g_r_mapping),cigar)
        candidate.end -= trim_signal
        # convert to read relative pos (forward direction)
        start_pos_rel_to_mapped_start = candidate.start - mapped_pos0
        end_pos_rel_to_mapped_start = candidate.end- mapped_pos0
        
        if  start_pos_rel_to_mapped_start >= 0 and start_pos_rel_to_mapped_start < len(g_r_mapping):
            candidate.start = g_r_mapping[start_pos_rel_to_mapped_start]
            
            # discard junction squiggle with the queried motif start/end mapped to gaps
            if candidate.start == -1:
                print("Warning: Abnormal mapping, junction squiggle skipped.")
                outf.write("\n")
                continue

            elif g_r_mapping[start_pos_rel_to_mapped_start] == \
                g_r_mapping[start_pos_rel_to_mapped_start - 1]:
                candidate.start += 1 
        else:
            print("candidate start pos out of bound.")
            outf.write("\n")
            continue

        
        if end_pos_rel_to_mapped_start < len(g_r_mapping):
            candidate.end = g_r_mapping[end_pos_rel_to_mapped_start - 1] + 1
            
            # discard junction squiggle with the queried motif start/end mapped to gaps
            if candidate.end == -1 + 1:
                print("Warning: Abnormal mapping, junction squiggle skipped.")
                outf.write("\n")
                continue
        else:
            print("candidate end pos out of bound.")
#            print(candidate.start, candidate.end, mapped_pos0, g_r_mapping[-1])
            outf.write("\n")
            continue


        # get signal
        if strand == "+":
            seg_start = max(candidate.start - tombo_start_clip, 0)
            seg_end = candidate.end - tombo_start_clip
            
        elif strand == "-":
            seg_start = max(read_length - candidate.end -1 - tombo_start_clip, 0)
            seg_end = read_length - candidate.start - 1 - tombo_start_clip
        
        # take into account the end clip
        seg_end = min(seg_end, len(tombo_results.segs) - 1)

        signal = normalised_raw_signal[tombo_results.segs[seg_start]:tombo_results.segs[seg_end]]

        if not len(signal):
            for i in range(len(candidate.sequences)):
                outf.write(",NA,NA,NA")
            outf.write("\n")
            print("read discarded")
            continue
        
        else:
            print("pass")


########################################delete later######################################
 #       with open("squiggle.csv", 'a') as squiggle_f:
   #         squiggle_f.write(','.join([str(x) for x in signal]) + '\n')
   #     with open("candidate.csv", 'a') as candidate_f:
   #         candidate_f.write(read_id + ',' + ','.join(candidate.sequences) + '\n')
   #     continue
##########################################################################################

        # Normalisation
        #signal = helper.normalization(signal,"z_score") # "median_shift" or "z_score"

        model_dic = expect_squiggle_dict(candidate.sequences, trim = trim_model)
        
        dtw_long = np.array(signal, float)

        for key in candidate.sequences:

            dtw_short = np.array(model_dic[key],float)
            #dtw_short[:,0] = helper.normalization(dtw_short[:,0], "z_score")
            #dtw_short[:,1] = dtw_short[:,1]/sqrt(np.std(dtw_short[:,0]))

            #print("dtw_short")
            dtw_short = np.array(dtw_short)
            #np.random.shuffle(dtw_short)
            #print("Input queried signal: " + '\t'.join([str(i) for i in dtw_long]))

            #print("\n\n\nInput model: " + '\t'.join([str(i) for i in dtw_short]))

            #print("\n\n\nRunning DTW...")

            timer_start = timeit.default_timer()
            #dtw_long = np.repeat(dtw_long,3)
            #dtw_long = dtw_long[abs(dtw_long)-3 < 0]
            path1 , score1 = 'NA','NA'
            #path1 , score1 = dtw_local_alignment(dtw_long, dtw_short, dist_type = "z_score")[0:2]
            #path2 , score2 = 'NA','NA'
            path2 , score2 = dtw_local_alignment(dtw_long, dtw_short, dist_type = "log_likelihood")[0:2]
            #likelihood, unmatched = dist_to_likelihood(dtw_long, dtw_short, path2, dist_type = "log_likelihood")
            #likelihood, unmatched = dist_to_likelihood_flipped(dtw_short, dtw_long, path2, dist_type = "log_likelihood")    
            #new_path2, likelihood=dist_to_likelihood_flipped_new_path(dtw_short, dtw_long, path2, dist_type = "log_likelihood")
            path3 , score3 = 'NA','NA'
            #path3 , score3 = dtw_local_alignment(dtw_long, dtw_short, dist_type = 'manhattan')[0:2]
            outf.write(',{},{},{}'.format(score1,score2,score3))
            #outf.write(',{},{},{}'.format(likelihood,score2,unmatched))
            timer_stop = timeit.default_timer()
            runtime = timer_stop - timer_start

            # ploting 
            if False:
                for score, path in (score1, path1), (score2, path2), \
                (score3, path3):
                    helper.plot_dtw_alignment(figure_name = "Untitled", \
                    figure_title = "Untitled" , long_seq = dtw_long, \
                    short_seq = dtw_short, dtw_path = path, dtw_score = score, \
                    show_sd = True, figsize=(10,7))

        outf.write('\n')

if __name__ == "__main__":
    main()
