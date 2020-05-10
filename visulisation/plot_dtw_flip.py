import sys, h5py, getopt, timeit, os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from collections import defaultdict
from math import sqrt
import scrappy

sys.path.append(os.path.dirname(os.path.realpath(__file__))+"/..")
import helper
import junction_squiggle_selection
from dtw import dtw_local_alignment_max_sum_band_flipped as dtw_local_alignment
from dtw import dist_to_likelihood_flipped

TRIM_SIGNAL = 6

def parse_arg():
    '''
    parse the command line input
    '''
    def print_help():
        print("\n\nUsage: python {} [OPTIONS] <fast5 filename>".format(argv[0]))
        print("Options:\n\tIndexing:")
        print('\t\t-q STR\tqueried candidate sequence (required)')
        print('\t\t-l STR\talignment for the read in SAM format (required)')
        
        print('\t\t-s INT\tstart base of the signal')
        print('\t\t-e INT\tend base of the signal')
        print("\t\t-o INT\toutput figure name, default 'Untitled'\n"
        "\t\t\t if 'mean' is included in the figure name, this script will optimise"
        " step-num-adjucted versionof dtw")
        print('\t\t-p INT\tbase position to be fetched')
        print('\t\t-w INT\twindow size around the -p, default 20')
        print('\t\t-v \tverbose mode')
        print("\t\t-m \tscrappie models: 'squiggle_r94' (D),'squiggle_r94_rna','squiggle_r10'")
        print('\t\t-T \tNumber of events trimmed from scrappie model')
        print("\tNote, either '-p','-w' or '-s', '-e' are required!")
        return None

    argv = sys.argv
    if len(argv) < 3:
        print_help()
        sys.exit()
    
    opts, args = getopt.getopt(argv[1:],"hs:e:o:q:p:w:vT:m:l:",
                    ["help","start=","end=","output_figure="
                    ,"sequence=", "position=", "window=","verbose","trim=","scrappie_model=","SAM="])

    # default setting
    start_pos, end_pos,sequence,t_position\
        = None, None, None, None
    verbose = False
    window = 20
    trim_model = 2
    trim_signal = TRIM_SIGNAL
    figure_name = "Untiled"
    scrappie_model = 'squiggle_r94'
    
    for opt, arg in opts:
        if opt in ('-h', "--help"):
            print_help()
            sys.exit()
        elif opt in ("-s", "--start"):
            start_pos = int(arg)
        elif opt in ("-l", "--SAM"):
            sam_line = arg
        elif opt in ("-e", "--end"):
            end_pos= int(arg)
        elif opt in ("-o", "--output_figure"):
            figure_name = arg
        elif opt in ("-q", "--sequence"):
            sequences = arg
        elif opt in ("-p", "--position"):
            t_position = int(arg)
        elif opt in ("-w", "--window"):
            window = int(arg)
        elif opt in ("-v", "--verbose"):
            verbose = True
        elif opt in ("-T", "--trim"):
            trim_model = int(arg)
        elif opt in ("-m", "--scrappie_model"):
            scrappie_model = arg

    fast5_filename = args[0]

    return start_pos, end_pos,sequence,t_position, verbose, trim_model, fast5_filename, window, figure_name, scrappie_model,sequences, sam_line,trim_signal



def main():

    # parse command line args:
    start_pos, end_pos,sequence,t_position, verbose, trim_model,\
    fast5_filename, window, figure_name,\
        scrappie_model,sequences, sam_line, trim_signal = parse_arg()
    

    # get alignment info
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

    # Get candidate motifs
    sequences = sequences.strip().split(",")
    if strand == "-":
        sequences = [helper.reverse_complement(s) for s in sequences]


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

    # get signal
    start_pos += trim_signal
    end_pos -= trim_signal
    start_pos_rel_to_mapped_start = start_pos - mapped_pos0
    end_pos_rel_to_mapped_start = end_pos- mapped_pos0

    if  start_pos_rel_to_mapped_start >= 0 and start_pos_rel_to_mapped_start < len(g_r_mapping):
        candidate_start = g_r_mapping[start_pos_rel_to_mapped_start]
        
        # discard junction squiggle with the queried motif start/end mapped to gaps
        if candidate_start == -1:
            print("Warning: Abnormal mapping, junction squiggle skipped.")

        elif g_r_mapping[start_pos_rel_to_mapped_start] == \
            g_r_mapping[start_pos_rel_to_mapped_start - 1]:
            candidate_start += 1 
    else:
        print("candidate start pos out of bound.")


    if end_pos_rel_to_mapped_start < len(g_r_mapping):
        candidate_end = g_r_mapping[end_pos_rel_to_mapped_start - 1] + 1
        
        # discard junction squiggle with the queried motif start/end mapped to gaps
        if candidate_end == -1 + 1:
            print("Warning: Abnormal mapping, junction squiggle skipped.")

    else:
        print("candidate end pos out of bound.")



    # get signal
    if strand == "+":
        seg_start = max(candidate_start - tombo_start_clip, 0)
        seg_end = candidate_end - tombo_start_clip
        
    elif strand == "-":
        seg_start = max(read_length - candidate_end -1 - tombo_start_clip, 0)
        seg_end = read_length - candidate_start - 1 - tombo_start_clip
    
    # take into account the end clip
    seg_end = min(seg_end, len(tombo_results.segs) - 1)

    signal = normalised_raw_signal[tombo_results.segs[seg_start]:tombo_results.segs[seg_end]]
    #spike removal
    signal = signal[abs(signal) < 3]

    # running DTW
    model_dic = helper.expect_squiggle_dict(sequences,model = scrappie_model,trim = trim_model)
    dtw_long = np.array(signal, float)
    
    candidate_ID = 0
    for key in model_dic.keys():
        candidate_ID += 1
        dtw_short = np.array(model_dic[key],float)

        if verbose:
            print("Input queried signal: " + '\t'.join([str(i) for i in dtw_long]))
            print("\n\n\nInput model: " + '\t'.join([str(i) for i in dtw_short]))

        #running DTW
        timer_start = timeit.default_timer()
        print("\n\n\nRunning DTW...")
        #path , score, matrix = dtw_local_alignment(dtw_long, dtw_short, dist_type = "z_score")
        path , score, matrix = dtw_local_alignment(dtw_long, dtw_short, dist_type = "log_likelihood")
        #path , score, matrix = dtw_local_alignment(dtw_long, dtw_short, dist_type = 'manhattan')
        timer_stop = timeit.default_timer()
        print("\n\nDTW finished, runtime: {} sec".format(timer_stop - timer_start))
        #print("\n\nAlignment distance: {}".format(score))

        likelihood, unmatched = dist_to_likelihood_flipped(dtw_short, dtw_long, path, dist_type = "log_likelihood")   
        print("#############################################################")
        print(likelihood)
        print("#############################################################")
        plt.figure(figsize=(300,30))
        matplotlib.rc('xtick', labelsize=20)     
        matplotlib.rc('ytick', labelsize=20)
        fig, axes = plt.subplots(nrows=3, figsize=(30,20))


###################################################################
        axes[0].plot(dtw_long,linewidth = 10)
        #axes[0].plot(dtw_short[:,0],color = "orange",linewidth = 6)
        #axes[0].plot(dtw_short[:,0] + dtw_short[:,1],':',color = "orange",linewidth = 4)
        #axes[0].plot(dtw_short[:,0] - dtw_short[:,1],':',color = "orange",linewidth = 4)

        
        axes[0].tick_params(labelsize=40)
        path = np.array(path[::-1])
        #print("\n\n\nBest path start and end:\n{} {}".format(np.array(path)[0,1],np.array(path)[-1,1]))
        #print("\n\n\nBest path length:\n{}".format(len(np.array(path))))
        
        axes[0].plot(np.array(path)[:,0]-1, dtw_short[[np.array(path)[:,1]-1]][:,0],'',color = "orange",linewidth = 7)
        
        #axes[0].plot(np.array(path)[:,1]-1, dtw_long[[np.array(path)[:,0]-1]],'',linewidth = 4)
        
        axes[0].plot(np.array(path)[:,0]-1, dtw_short[[np.array(path)[:,1]-1]][:,0]\
                                        + dtw_short[[np.array(path)[:,1]-1]][:,1],':',color = "orange",linewidth = 4)

        axes[0].plot(np.array(path)[:,0]-1, dtw_short[[np.array(path)[:,1]-1]][:,0]\
                                        - dtw_short[[np.array(path)[:,1]-1]][:,1],':',color = "orange",linewidth = 4)

        #axes[0].set_title(figure_name+"_"+key+"({})\nDist: {:.2f}, path length: {}, Adjusted dist: {:.2f}".format(str(candidate_ID == 1),score*len(np.array(path)),len(np.array(path)),score),fontsize=40, pad = 30)
        axes[0].set_title(figure_name+"_"+key+"({})\nDist: {:.2f}, path length: {}, Adjusted dist: {:.2f}".format(str(candidate_ID == 1),score,len(np.array(path)),score/len(np.array(path))),fontsize=40, pad = 30)
        print("candidate " + key + "finished!")
        

        #plt.savefig(figure_name+"_"+str(candidate_ID)+".png")
        
        #plot simulated squiggle?

        if True:
            #axes[1].figure(figsize=(10,5))
            axes[1].plot(dtw_short[:,0], color = "orange",linewidth =7)
            axes[1].plot(dtw_short[:,0] + dtw_short[:,1], ":",color = "orange",linewidth = 4)
            axes[1].plot(dtw_short[:,0] - dtw_short[:,1], ":",color = "orange",linewidth = 4)
            axes[1].tick_params(labelsize=40)
            axes[1].set_title("Scrappie model",fontsize=30, pad = 20)
            #plt.savefig(figure_name+"_"+str(candidate_ID)+"simulated_squiggle.png")

        # plot path heatmap
        if True:
            #axes[2].figure(figsize=(10,5))
            #fig.figsize=(10,5)
            pos = axes[2].imshow(matrix, cmap='hot', interpolation='nearest',aspect='auto')
            fig.colorbar(pos,ax = axes[2])
            axes[2].plot(path[:,1], path[:,0])
            axes[2].set_title("Alignment path",fontsize=30, pad = 20)

        
        # plot squiggle
        if False:
            axes[2].plot(dtw_long,linewidth = 4)
        
        
        plt.subplots_adjust(hspace=0.5)
        plt.savefig(figure_name+"_" + "{}({}).png".format(str(candidate_ID), str(candidate_ID == 1)))


if __name__ == "__main__":
    main()
