'''
 Plot the squiggle from fast5 file
'''
import sys, h5py, os
import numpy as np
import matplotlib.pyplot as plt

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

    # Normalisation
    signal = helper.normalization(signal, "z_score")

    fig, ax = plt.subplots()
    fig.set_size_inches(20, 9.5)
    ax.plot(np.array(signal))

    ax.set_ylabel("Current levle", fontsize = 25)
    ax.set_xlabel("index", fontsize = 25)
    ax.grid()
    
    fig.savefig(argv[2])
    print("\nFigure saved as {}. \n".format(argv[2]))


if __name__ == "__main__":
    main()
    
