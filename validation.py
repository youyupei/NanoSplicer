import sys
import h5py
import getopt
import timeit
import os
import numpy as np


from math import sqrt

#sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import helper
from helper import expect_squiggle_dict
from helper import parse_candidate_file
from dtw import dtw_local_alignment as dtw_sum
from dtw_2 import dtw_local_alignment as dtw_mean

# choose the version of dtw (sum: minimize the sum of cost for the path)
dtw_local_alignment = dtw_sum



def main():
    def print_help():
        '''
        print command line instruction
        '''
        print("\n\nUsage: python {} [OPTIONS] <fast5 filename>".format(argv[0]))
        print("Options:\n\tIndexing:")
        print('\t\t-c INT\tcandidate file name')
        print('\t\t-o INT\toutput csv name')
        return None
      
    argv = sys.argv
    
    
    if len(argv) <= 2:     
        print_help()       # print help doc when no command line args provided
        exit()
    
    try: 
        opts, args = getopt.getopt(argv[1:],"hs:e:o:m:f:q:p:w:c:",
                    ["help=","start=","end=","output_figure=",
                    "model=","fasta=","sequence=", "position=", "window=",
                    "candidate_file="])
    except getopt.GetoptError:
        print_help()
        exit()

    start_pos, end_pos,read_from_model,read_from_fasta,sequence,candidate_file \
        = None, None, None, None, None, None

    t_position = None

    #window = 20
    output_file = "Untiled"
  
    try:
        fast5_filename = args[0]
    except:
        print("InputError: missing fast5 file name!")
        exit()

    for opt, arg in opts:
        if opt in ('-h', "--help"):
            print_help()
            exit(0)
        elif opt in ("-s", "--start"):
            start_pos = int(arg)
        elif opt in ("-e", "--end"):
            end_pos= int(arg)
        elif opt in ("-o", "--output_csv"):
            output_file = arg
        elif opt in ("-m", "--model"):
            read_from_model = arg
        elif opt in ("-f", "--fasta"):
            read_from_fasta = arg
        elif opt in ("-q", "--sequence"):
            sequence = arg.strip().split(',')
        elif opt in ("-p", "--position"):
            t_position = int(arg)
        elif opt in ("-w", "--window"): 
            window = int(arg)
        elif opt in ("-c", "--candidate_file"):
           candidate_file = arg


    outf = open(output_file,'w')
    '''
    outf.write("read id")
    for i in range(len(sequence)):
        outf.write(",z_score{},log_likelihood{},manhattan{}".format(i+1,i+1,i+1))
    
    outf.write('\n'+ fast5_filename)
    '''
    if "candidate_file" in locals():
        candidates = parse_candidate_file(candidate_file)


    strand = helper.Fast5Class(fast5_filename).get_alignment(
        "mapping_info")["mapped_strand"]
    for candidate in candidates:
        outf.write(fast5_filename + ',' + strand + ',')
        
        # get signal
        signal = helper.get_junction_signal_by_pos(
            fast5 = fast5_filename, start_pos = candidate.start ,
             end_pos = candidate.end)
        
        if not signal:
            for i in range(len(candidate.sequences)):
                outf.write(",NA,NA,NA")
            outf.write("\n")
            print("read discarded")
            continue

        # take reverse compliment seq if nesseccary
        
        if strand == "-":
            candidate.sequences = [helper.reverse_complement(s)
             for s in candidate.sequences]

        # Normalisation
        signal = helper.normalization(signal,"z_score") # "median_shift" or "z_score"

        model_dic = expect_squiggle_dict(candidate.sequences)
        
        dtw_long = signal

        for key in candidate.sequences:

            dtw_short = np.array(model_dic[key])
            dtw_short[:,0] = helper.normalization(dtw_short[:,0], "z_score")
            dtw_short[:,1] = dtw_short[:,1]/sqrt(np.std(dtw_short[:,0]))

            #print("dtw_short")
            dtw_short = np.array(dtw_short)
            #print("Input queried signal: " + '\t'.join([str(i) for i in dtw_long]))

            #print("\n\n\nInput model: " + '\t'.join([str(i) for i in dtw_short]))

            #print("\n\n\nRunning DTW...")

            timer_start = timeit.default_timer()
            #dtw_long = np.repeat(dtw_long,3)
            #dtw_long = dtw_long[abs(dtw_long)-3 < 0]
            
            path1 , score1 = dtw_local_alignment(
                dtw_long, dtw_short, dist_type = "z_score")
            path2 , score2 = dtw_local_alignment(
                dtw_long, dtw_short, dist_type = "log_likelihood")   
            path3 , score3 = dtw_local_alignment(
                dtw_long, dtw_short, dist_type = 'manhattan')
            
            outf.write(',{},{},{}'.format(
                score1/len(path1),score2/len(path2),score3/len(path3)))
            
            timer_stop = timeit.default_timer()
            runtime = timer_stop - timer_start

            # ploting 
            if True:
                for score, path in (score1, path1), (score2, path2), \
                (score3, path3):
                    helper.plot_dtw_alignment(figure_name = "Untitled", \
                    figure_title = "Untitled" , long_seq = dtw_long, \
                    short_seq = dtw_short, dtw_path = path, dtw_score = score, \
                    show_sd = True, figsize=(10,7))

        outf.write('\n')

if __name__ == "__main__":
    main()
