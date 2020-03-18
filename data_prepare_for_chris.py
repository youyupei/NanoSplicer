import sys, h5py, getopt, timeit, os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from collections import defaultdict
from math import sqrt
import scrappy

sys.path.append(os.path.dirname(os.path.realpath(__file__))+"/..")
import helper

from dtw import dtw_local_alignment_max_mean as dtw_mean
from dtw import dtw_local_alignment_max_sum as dtw_sum
from dtw import dtw_local_alignment_max_sum_band as dtw_band


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
        trim_model = 4
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

        # choose the version of dtw (sum: minimize the sum of cost for the path)
        dtw_local_alignment = dtw_sum
        if "mean" in figure_name:
            dtw_local_alignment = dtw_mean
        elif "band" in figure_name:
            dtw_local_alignment = dtw_band
        fast5_filename = args[0]
        suffix = args[1]

        return start_pos, end_pos,sequence,t_position, verbose, trim_model,\
            dtw_local_alignment,fast5_filename, window, figure_name, scrappie_model,sequences,suffix

    except:
        print("Error,failed to parse command line arguments!")
        print_help()
        sys.exit(1)


def main():

    # parse command line args:
    start_pos, end_pos,sequence,t_position, verbose, trim_model,\
        dtw_local_alignment,fast5_filename, window, figure_name,\
        scrappie_model,sequences, suffix = parse_arg()

    # get signal
    if start_pos and end_pos:
        signal = helper.get_junction_signal_by_pos(fast5_filename, \
         start_pos = start_pos, end_pos = end_pos)

    elif t_position:
        signal = helper.get_junction_signal_by_pos(fast5_filename, 
                            junction_pos = t_position, window = window)
    else:
        print("Invalid command line input:"
        "\t\tEither '-p','-w' or '-s', '-e' are required!")
        sys.exit(0)

    with open("squiggle_"+suffix+".csv", 'a') as squiggle_f:
        squiggle_f.write(','.join([str(x) for x in signal]) + '\n')
    if not sequences:
        print("-q is required!!!")
        sys.exit(0)
    else:
        sequences = sequences.strip().split(',')

    # Get strand(in fast5)
    try:
        strand = helper.Fast5Class(fast5_filename).get_alignment(
            "mapping_info")["mapped_strand"]
        if strand == "-":
            sequences = [helper.reverse_complement(s) for s in sequences]
    except:
        sys.exit(0)

    with open("candidate_"+suffix+".csv", 'a') as candidate_f:
        candidate_f.write(','.join(sequences) + '\n')



if __name__ == "__main__":
    main()
