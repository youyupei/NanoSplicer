'''
 Plot the squiggle from fast5 file
'''

import sys, h5py, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 30})

sys.path.append(os.path.dirname(os.path.realpath(__file__))+"/..")
import helper


def main():
    argv = sys.argv
    if len(argv) < 3:
        exit("Usage: python {} <fast5 filename> <save fig as>\
        <optianl: base_postision_in_transcript>\
        <optional: window_size>".format(argv[0]))
    filename = argv[1]
    
    if len(argv) == 5:
        signal = helper.get_junction_signal_by_pos(filename, int(argv[3]), 
                                                                int(argv[4]))
    else:
        signal = helper.read_raw_signal(filename)
        med_sig = np.median(signal)
        mad_sig = np.median(abs(med_sig - signal)) * 1.4826
        signal = (signal - med_sig)/mad_sig

    # Normalisation
    #signal = helper.normalization(signal, "z_score")

    fig, ax = plt.subplots()
    fig.set_size_inches(20, 9.5)
    ax.plot(np.array(signal),linewidth = 1)
    ax.plot([3]*len(signal), ":",color = "orange",linewidth = 4)
    ax.plot([-3]*len(signal), ":",color = "orange",linewidth = 4)

    #ax.set_ylabel("Current levle", fontsize = 25)
    #ax.set_xlabel("index", fontsize = 25)
    #ax.grid()
    
    fig.savefig(argv[2])
    print("\nFigure saved as {}. \n".format(argv[2]))


if __name__ == "__main__":
    main()
    
