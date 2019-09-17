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
#from dtw_cy.dtw import dtw_local_alignment_max_sum as dtw_sum





#def sequence_to_squiggle(seq, trim = True):
#    '''input:
#        sequence <str>
#        output:
#        numpy array:
#            [[mean, std, dwell]...]
#    '''
#    simulated_seq = scrappy.sequence_to_squiggle(seq,rescale =True).data(
#                            as_numpy = True, sloika = False)
#    if trim:
#
#        return simulated_seq[2:-2]
#    else:
#        return simulated_seq
#
#def expect_squiggle_dict(seq=None, read_from_fasta = None, read_from_model = None):
#    '''
#    read squiggle data from scrappie, ready for dtw
#    '''
#    expect_squiggle_dic = defaultdict(list)
#    
#    if not seq and not read_from_fasta and not read_from_model:
#        sys.exit("No valid input detected when generating expect squiggle")
#    
#    if seq:
#        assert(isinstance(seq, str))
#        seq = seq.strip()
#        squiggle = sequence_to_squiggle(seq)
#        for mean, std, dwell_time in squiggle:
#            expect_squiggle_dic['seq'] += [[mean, std]] *int(round(dwell_time))
#            
#
#    if read_from_fasta:
#        with open(read_from_fasta,'r') as f:
#            for l in f:
#                l = l.strip('\n')
#                if l[0] == '>':
#                    model = l[1:]
#                else:
#                    squiggle = sequence_to_squiggle(l)
#                    for mean, std, dwell_time in squiggle:
#                        expect_squiggle_dic['seq'] += [[mean, std]] * dwell_time
#   
#    if read_from_model:
#        with open(read_from_model, 'r') as r:
#            for l in r:
#                l = l.strip('\n')
#                if l[0] == '#':
#                    model = l[1:]
#                elif l[:3] == "pos":
#                    continue
#                else:
#                    l = l.split()
#                    mean = float(l[2])
#                    std = float(l[3])
#                    dwell_time = int(round(float(l[4])))
#                    expect_squiggle_dic[model] += [[mean, std]] * dwell_time
#    
#    return expect_squiggle_dic


def main():
    # choose the version of dtw (sum: minimize the sum of cost for the path)
    dtw_local_alignment = dtw_sum

    def print_help():
        print("\n\nUsage: python {} [OPTIONS] <fast5 filename>".format(argv[0]))
        print("Options:\n\tIndexing:")
        print('\t\t-s INT\tstart base of the signal')
        print('\t\t-e INT\tend base of the signal')
        print('\t\t-o INT\toutput figure name')
        print('\t\t-p INT\tbase position to be fetched')
        print('\t\t-w INT\twindow size around the -p, default 20')
        print('\t\t-m STR\tscrappie squiggle model filename')
        print('\t\t-f STR\tqueried sequence FASTA filename')
        print('\t\t-q STR\tqueried candidate sequence')
        

        return None
      

    argv = sys.argv
    if len(argv) < 3:
        print_help()
        sys.exit()
    
    try: 
        opts, args = getopt.getopt(argv[1:],"hs:e:o:m:f:q:p:w:",
                    ["help=","start=","end=","output_figure=",
                    "model=","fasta=","sequence=", "position=", "window="])
    
    except getopt.GetoptError:
        print_help()
        sys.exit()
    
    
    start_pos, end_pos,read_from_model,read_from_fasta,sequence \
        = None, None, None, None, None

    t_position = None
    window = 20
    figure_name = "Untiled"
    
    
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
        elif opt in ("-m", "--model"):
            read_from_model = arg
        elif opt in ("-f", "--fasta"):
            read_from_fasta = arg
        elif opt in ("-q", "--sequence"):
            sequences = arg
        elif opt in ("-p", "--position"):
            t_position = int(arg)
        elif opt in ("-w", "--window"):
            window = int(arg)
    
    if "mean" in figure_name:
        dtw_local_alignment = dtw_mean


    fast5_filename = args[0]
    # get signal
    if start_pos and end_pos:
        signal = helper.get_junction_signal_by_pos(fast5_filename, \
         start_pos = start_pos, end_pos = end_pos)

    elif t_position:
        signal = helper.get_junction_signal_by_pos(fast5_filename, 
                            t_position, window)

    if not sequences:
        print("-q is required!!!")
        sys.exit(0)
    else:
        sequences = sequences.strip().split(',')
    print(sequences)


    
    # Get strand(in fast5)
    try:
        strand = helper.Fast5Class(fast5_filename).get_alignment(
            "mapping_info")["mapped_strand"]
        if strand == "-":
            sequences = [helper.reverse_complement(s) for s in sequences]
    except:
        sys.exit(0)

    # Normalisation
    #signal = helper.normalization(signal,"median_shift") # "median_shift" or "z_score"

    model_dic = helper.expect_squiggle_dict(sequences)
    dtw_long = np.array(signal, float)
    
    candidate_ID = 0
    for key in model_dic.keys():
        candidate_ID += 1
        dtw_short = np.array(model_dic[key],float)
        
        
        #normalise the simulated squiggle
        #dtw_short[:,0] = helper.normalization(dtw_short[:,0], "median_shift")
        
        
        dtw_short[:,1] = dtw_short[:,1]/sqrt(np.std(dtw_short[:,0]))

        #print("dtw_short")
        dtw_short = np.array(dtw_short)
        print("Input queried signal: " + '\t'.join([str(i) for i in dtw_long]))

        #print("\n\n\nInput model: " + '\t'.join([str(i) for i in dtw_short]))

        print("\n\n\nRunning DTW...")

        timer_start = timeit.default_timer()
        
        
        #dtw_long = np.repeat(dtw_long,3)
        #dtw_long = dtw_long[abs(dtw_long)-3 < 0]
        #path , score, matrix = dtw_local_alignment(dtw_long, dtw_short, dist_type = "z_score")
        path , score, matrix = dtw_local_alignment(dtw_long, dtw_short, dist_type = "log_likelihood")
        #path , score, matrix = dtw_local_alignment(dtw_long, dtw_short, dist_type = 'manhattan')


        timer_stop = timeit.default_timer()
        print("\n\nDTW finished, runtime: {} sec".format(timer_stop - timer_start))
        print("\n\nAlignment distance: {}".format(score))


        plt.figure(figsize=(300,30))
        matplotlib.rc('xtick', labelsize=20)     
        matplotlib.rc('ytick', labelsize=20)
        fig, axes = plt.subplots(nrows=3, figsize=(30,20))

        axes[0].plot(dtw_long)
        path = np.array(path[::-1])
        print("\n\n\nBest path start and end:\n{} {}".format(np.array(path)[0,1],np.array(path)[-1,1]))
        print("\n\n\nBest path length:\n{}".format(len(np.array(path))))
        
        axes[0].plot(np.array(path)[:,1]-1, dtw_short[[np.array(path)[:,0]-1]][:,0],'g')
        axes[0].plot(np.array(path)[:,1]-1, dtw_short[[np.array(path)[:,0]-1]][:,0]\
                                        + dtw_short[[np.array(path)[:,0]-1]][:,1],'g--')

        axes[0].plot(np.array(path)[:,1]-1, dtw_short[[np.array(path)[:,0]-1]][:,0]\
                                        - dtw_short[[np.array(path)[:,0]-1]][:,1],'g--')

        axes[0].set_title(figure_name+"_"+key+"({})\nDist: {:.2f}, path length: {}, Adjusted dist: {:.2f}".format(str(candidate_ID == 1),score*len(np.array(path)),len(np.array(path)),score),fontsize=40)
        
        
        #plt.savefig(figure_name+"_"+str(candidate_ID)+".png")
        
        #plot simulated squiggle only?
        if True:
            #axes[1].figure(figsize=(10,5))
            axes[1].plot(dtw_short[:,0], 'g')
            axes[1].set_title("Scrappie model",fontsize=30)
            #plt.savefig(figure_name+"_"+str(candidate_ID)+"simulated_squiggle.png")

        # plot path heatmap
        if True:
            #axes[2].figure(figsize=(10,5))
            #fig.figsize=(10,5)
            pos = axes[2].imshow(matrix, cmap='hot', interpolation='nearest',aspect='auto')
            fig.colorbar(pos,ax = axes[2])
            axes[2].plot(path[:,1], path[:,0])
            axes[2].set_title("Alignment path",fontsize=30)
            plt.savefig(figure_name+"_" + "{}({}).png".format(str(candidate_ID), str(candidate_ID == 1)))

if __name__ == "__main__":
    main()
