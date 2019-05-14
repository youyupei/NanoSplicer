
import matplotlib.pyplot as plt
import numpy as np
import sys



def mean_normalization(signal):
    '''
    input: signal dataset from fast5
    output: nomalised signal <np.array>
    '''
    signal = np.array(signal)
    mad_scale = np.median(abs(signal - np.median(signal)))

    norm_signal = (signal - np.median(signal)) / mad_scale

    return(norm_signal)


def read_synth_model(filename):
    '''
    read squiggle data from scrappie, old and new version, ready for dtw
    '''
    dic = {}
    with open(filename, 'r') as r:
        for l in r:
            l = l.strip('\n')
            if l[0] == '#':
                name = l[1:]
                dic[name] = []
            elif l[:3] == "pos":
                continue
            else:
                l = l.split()
                dic[name] = dic[name] + [float(l[2])] * int(round(float(l[4])))
    # print >> sys.stderr, len(dic['adapter'])
    return dic

def main():
    argv = sys.argv
    if len(argv) < 3:
        exit("Usage: python {} <.model file> <save fig as>".format(argv[0]))
    filename = argv[1]

    signal_dic = read_synth_model(filename)

    # Normalisation
    for key in signal_dic.keys():
        signal = mean_normalization(signal_dic[key])
        fig, ax = plt.subplots()
        fig.set_size_inches(20, 9.5)
        ax.plot(np.array(signal))

        ax.set_ylabel("Current levle", fontsize = 25)
        ax.set_xlabel("index", fontsize = 25)
        ax.grid()

        fig.savefig(argv[2]+key+'.png')

if __name__ == "__main__":
    main()