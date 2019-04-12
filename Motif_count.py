# count randomly generated motif in a fasta/ fastq file

import sys
import numpy.random as random
import regex as re

# argument parsing
'''
arg = sys.argv
if len(arg) == 1:
    print("Usage: python {} <filename> ["DNA"/"RNA"]".format(arg[0]))
    exit()

'''


def count_in_file(motif, filename):
    occurence_overlap = 0
    occurence_no_overlap = 0
    with open(filename,'r') as f:
        for line in f:
            occurence_overlap += len(re.findall(motif, line, overlapped = True))
            occurence_no_overlap += len(re.findall(motif, line, overlapped = False))
    return occurence_overlap, occurence_no_overlap


def main(): 
    #FILENAME = arg[1] 
    FILENAME = "/home/youyupei/PhD_proj/cDNA_dataset_Analyses/nanopolish/intermedia_data/exp_reads.fasta"
    necleotide = 'DNA'
    n_choice = ["A", "G", "C", "T"] if necleotide == 'DNA' else ["A", "G", "C", "U"]
    # sample from n_choice with replacement
    with open("motif_count.csv",'w') as f:
        f.write("motif_len, count allow overlapping, count without overlapping\n")
        for motif_len in range(5, 21):
            motif = [''.join(random.choice(n_choice, size=motif_len)) for i in range(30)]
            print(motif)
            occurence = [count_in_file(i, FILENAME) for i in motif]
            for overlap,  no_overlap in occurence:
                f.write(str(motif_len)+','+str(overlap) +',' + str(no_overlap) + '\n')
                print("..debug")





if __name__ == "__main__":
    main()
