'''
plot the scrappie model
'''
import matplotlib.pyplot as plt
import numpy as np
import sys
from collections import defaultdict
import scrappy, os
sys.path.append(os.path.dirname(os.path.realpath(__file__))+"/..")
import helper


def read_synth_model(filename):
    '''
    read squiggle data from scrappie, old and new version, ready for dtw
    '''
    expect_squiggle = defaultdict(list)
    mean = defaultdict(list)
    std = defaultdict(list)
    with open(filename, 'r') as r:
        for l in r:
            l = l.strip('\n')
            if l[0] == '#':
                name = l[1:]
            elif l[:3] == "pos":
                continue
            else:
                l = l.split()
                mean[name].append(float(l[2]))
                std[name].append(float(l[3]))
                expect_squiggle[name] += [float(l[2])] * int(round(float(l[4])))
    # print >> sys.stderr, len(dic['adapter'])
    return mean, std, expect_squiggle

def main():
    argv = sys.argv
    if len(argv) < 3:
        exit("Usage: python {} <.model file/or sequence> <save fig as>".format(argv[0]))
    
    sequence, filename = None, None 

    if '.model' not in argv[1]:
        sequence = argv[1]
    if True:
        sequence = helper.reverse_complement(sequence)
    else:
        filename = argv[1]


    if filename:
        synth_mean,synth_std, signal_dic = read_synth_model(filename)
    if sequence:
        sim_seq = helper.sequence_to_squiggle(sequence)
        signal_dic = defaultdict(list)
        for mean, std, dwell_time in sim_seq:
            signal_dic['seq'] += [mean] *int(round(dwell_time))

        


    # Normalisation
    for key in signal_dic.keys():
        signal = helper.normalization(signal_dic[key], 'z_score')
        fig, ax = plt.subplots()
        fig.set_size_inches(20, 9.5)
        ax.plot(np.array(signal))


        ax.set_ylabel("Current levle", fontsize = 25)
        ax.set_xlabel("index", fontsize = 25)
        ax.grid()

        fig.savefig(argv[2]+key+'.png')
        print("Figure saved!")

if __name__ == "__main__":
    main()