import sys, h5py, getopt, timeit, os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from collections import defaultdict
from math import sqrt
import scrappy

sys.path.append(os.path.dirname(os.path.realpath(__file__))+"/..")
import helper

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
    try:
        argv = sys.argv
        if len(argv) < 3:
            print_help()
            sys.exit()
        
        opts, args = getopt.getopt(argv[1:],"hs:e:o:q:p:w:vT:m:",
                        ["help","start=","end=","output_figure="
                        ,"sequence=", "position=", "window=","verbose","trim=","scrappie_model="])

        # default setting
        start_pos, end_pos,sequence,t_position\
            = None, None, None, None
        verbose = False
        window = 20
        trim_model = 2
        figure_name = "Untiled"
        scrappie_model = 'squiggle_r94'
        
        for opt, arg in opts:
            if opt in ('-h', "--help"):
                print_help()
                sys.exit()
            elif opt in ("-s", "--start"):
                start_pos = int(arg)
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

        return start_pos, end_pos,sequence,t_position, verbose, trim_model, fast5_filename, window, figure_name, scrappie_model,sequences

    except:
        print("Error,failed to parse command line arguments!")
        print_help()
        sys.exit(1)


def main():


    # parse command line args:
    start_pos, end_pos,sequence,t_position, verbose, trim_model,\
    fast5_filename, window, figure_name,\
        scrappie_model,sequences = parse_arg()
    

    # get signal
    if start_pos and end_pos:
        signal = helper.get_junction_signal_by_pos(fast5_filename, \
         start_pos = start_pos + TRIM_SIGNAL, end_pos = end_pos - TRIM_SIGNAL)


    elif t_position:
        signal = helper.get_junction_signal_by_pos(fast5_filename, 
                            junction_pos = t_position, window = window)
    else:
        print("Invalid command line input:"
        "\t\tEither '-p','-w' or '-s', '-e' are required!")
        sys.exit(0)

    if not sequences:
        print("-q is required!!!")
        sys.exit(0)
    else:
        sequences = sequences.strip().split(',')
    
    if verbose:
        print("Candidate sequence:")
        print(sequences)

    # Get strand(in fast5)
    try:
        strand = helper.Fast5Class(fast5_filename).get_alignment(
            "mapping_info")["mapped_strand"]
        if strand == "-":
            sequences = [helper.reverse_complement(s) for s in sequences]
    except:
        sys.exit(0)

    model_dic = helper.expect_squiggle_dict(sequences,model = scrappie_model,trim = trim_model)
    dtw_long = np.array(signal, float)
    
    candidate_ID = 0
    for key in model_dic.keys():
        candidate_ID += 1
        dtw_short = np.array(model_dic[key],float)
        
        #normalise the simulated squiggle
        #dtw_short[:,0] = helper.normalization(dtw_short[:,0], "median_shift")
        #dtw_short[:,1] = dtw_short[:,1]/sqrt(np.std(dtw_short[:,0]))

        #print("dtw_short")
        dtw_short = np.array(dtw_short)
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

        axes[0].set_title(figure_name+"_"+key+"({})\nDist: {:.2f}, path length: {}, Adjusted dist: {:.2f}".format(str(candidate_ID == 1),score*len(np.array(path)),len(np.array(path)),score),fontsize=40, pad = 30)
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
        if False:
            #axes[2].figure(figsize=(10,5))
            #fig.figsize=(10,5)
            pos = axes[2].imshow(matrix, cmap='hot', interpolation='nearest',aspect='auto')
            fig.colorbar(pos,ax = axes[2])
            axes[2].plot(path[:,1], path[:,0])
            axes[2].set_title("Alignment path",fontsize=30, pad = 20)

        
        # plot squiggle
        if True:
            axes[2].plot(dtw_long,linewidth = 4)
        
        
        plt.subplots_adjust(hspace=0.5)
        plt.savefig(figure_name+"_" + "{}({}).png".format(str(candidate_ID), str(candidate_ID == 1)))


if __name__ == "__main__":
    main()
