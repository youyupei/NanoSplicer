'''
extract signal from fast5 and write it to
 .tsv file(SquiggleKit Format)
'''
import sys
from raw_data_processing import *


def main():
    arg = sys.argv
    if len(arg) < 3:
        exit('Usage: python {}\
            <fast5_file> <tsv filename> <start> <end>'.format(arg[0]))

    if len(arg) == 5:
        signal = read_raw_segnal(arg[1], 
                start = arg[3], end = arg[4])
    else:
        signal = read_raw_segnal(arg[1])

    print(signal[:10])

    with open(arg[2], 'w') as tsv_f:
        tsv_f.write(arg[1]+'\t')
        tsv_f.write(str(len(signal)) + '\tall\t')
        tsv_f.write('\t'.join([str(i) for i in signal]))

    return None

if __name__ == '__main__':
    main()
        

    
