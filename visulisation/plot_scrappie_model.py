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
from dtw import dtw_local_alignment


def sequence_to_squiggle(seq):
    '''input:
        sequence <str>
        output:
        numpy array:
            [[mean, std, dwell]...]
    '''
    simulated_seq = scrappy.sequence_to_squiggle(seq,rescale =True).data(
                            as_numpy = True, sloika = False)
    return simulated_seq

def expect_squiggle_dict(seq=None, read_from_fasta = None, read_from_model = None):
    '''
    read squiggle data from scrappie, ready for dtw
    '''
    expect_squiggle_dic = defaultdict(list)
    
    if not seq and not read_from_fasta and not read_from_model:
        exit("No valid input detected when generating expect squiggle")
    
    if seq:
        assert(isinstance(seq, str))
        seq = seq.strip()
        squiggle = sequence_to_squiggle(seq)
        for mean, std, dwell_time in squiggle:
            expect_squiggle_dic['seq'] += [[mean, std]] *int(round(dwell_time))
            

    if read_from_fasta:
        with open(read_from_fasta,'r') as f:
            for l in f:
                l = l.strip('\n')
                if l[0] == '>':
                    model = l[1:]
                else:
                    squiggle = sequence_to_squiggle(l)
                    for mean, std, dwell_time in squiggle:
                        dwell_time = int(round(dwell_time))
                        expect_squiggle_dic['seq'] += [[mean, std]] * dwell_time
   
    if read_from_model:
        with open(read_from_model, 'r') as r:
            for l in r:
                l = l.strip('\n')
                if l[0] == '#':
                    model = l[1:]
                elif l[:3] == "pos":
                    continue
                else:
                    l = l.split()
                    mean = float(l[2])
                    std = float(l[3])
                    dwell_time = int(round(float(l[4])))
                    expect_squiggle_dic[model] += [[mean, std]] * dwell_time
    
    return expect_squiggle_dic

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
    else:
        filename = argv[1]


    if filename:
        synth_mean,synth_std, signal_dic = read_synth_model(filename)
    if sequence:
        sim_seq = sequence_to_squiggle(sequence)
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