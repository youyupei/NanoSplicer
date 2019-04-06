
"""
This script is to select the mapped reads from a bam file and then write
it to a new bam file.
"""

import pysam, sys, getopt
from collections import defaultdict as dd

def parse_arg(argv):
    """ 
    parsing flagged command line argument
    return inputfile then outputfile
    """
    try: 
        opts, args = getopt.getopt(argv,"hi:o:q:",
                    ["help=","input=","output=","quality-thresh="])
    except getopt.GetoptError:
        print('Usage: {} -i <inputfile> -o <outputfile>'.format(sys.argv[0]))
        sys.exit
    
    # init mapq threshold
    mapq_thres = -1
    
    for opt, arg in opts:
        if opt in ('-h', "--help"):
            print('Usage: {} -i <inputfile> -o <outputfile>'.format(sys.argv[0]))
            sys.exit()
        elif opt in ("-i", "--input"):
            inputfile = arg
        elif opt in ("-o", "--output"):
            outputfile = arg
        elif opt in ("-q", "--quality-thresh"):
            mapq_thres = int(arg)   
    
    return inputfile, outputfile, mapq_thres

def main():

    """
    Read Argument from the command line:

        BAMFILE:
            <string> filename of input .bam file
        
        OUTFILENAME:
            <string> filename of output .bam file
        
        MAPQ_THRES(optional):
            <int> Threshold of mapping qualitis, read with mapping qulity 
            lower than this value will not be written out, default -1
    """
    BAMFILE, OUTFILENAME, MAPQ_THRES = parse_arg(sys.argv[1:])
    samfile = pysam.AlignmentFile(BAMFILE,"rb")


    # write out mapped reads to new bam file
    mappedreads = pysam.AlignmentFile(OUTFILENAME,"wb",
                                        template = samfile)
    for read in samfile:
        if not read.is_unmapped and read.mapq >= MAPQ_THRES:
            mappedreads.write(read)

    mappedreads.close()
    samfile.close()

if __name__ == "__main__":
    main()