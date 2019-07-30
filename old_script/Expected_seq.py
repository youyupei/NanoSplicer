'''
construct the expected_sequence given the .bam file
Expected Sequence = left soft clipping + mapped reference+ right soft clipping 
'''
import pysam, getopt, pickle, sys
from collections import defaultdict

def parse_arg(argv):
    """ 
    parsing flagged command line argument
    return inputfile then outputfile
    """
    try: 
        opts, args = getopt.getopt(argv,"hi:r:f:p:",
                    ["help=","input=","reference=","fasta=","pickle="])
    except getopt.GetoptError:
        print('Usage: {} -i <inputfile> -r <ref.fa> -f <fastafile> -p <picklefile>'.format(sys.argv[0]))
        sys.exit()
    
    fastafile = None 
    picklefile = None
    for opt, arg in opts:
        if opt in ('-h', "--help"):
            print('Usage: {} -i <inputfile> -r <ref.fa> -o <outputfile>'.format(sys.argv[0]))
            sys.exit()
        elif opt in ("-i", "--input"):
            inputfile = arg

        elif opt in ("-r", "--reference"):
            reference = arg

        elif opt in ("-f", "--fasta"):
            fastafile = arg  
        elif opt in ("-p", "--pickle"):
            picklefile = arg  
    return inputfile,reference, fastafile, picklefile
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

def pickle_dump(py_object, output_filename):
    # Pickle the object to <output.pickle> 
    # using the highest protocol available.
    with open(output_filename, 'wb') as f:
        pickle.dump(py_object, f, pickle.HIGHEST_PROTOCOL)
    print("Python Object Saved.")


if __name__ == "__main__":
    # Read BAM File
    BAMFILE,REFERENCE, FASTA_OUT, PICKLE_OUT = parse_arg(sys.argv[1:])
    print(BAMFILE,REFERENCE, FASTA_OUT, PICKLE_OUT)
    samfile = pysam.AlignmentFile(BAMFILE,"rb")
    ref = pysam.FastaFile(REFERENCE) 
    expected_seq_dict = defaultdict(list)
    if FASTA_OUT:
        fasta = open(FASTA_OUT, "w")
    for read in samfile:
        query_seq = read.query_sequence
        ref_seq = ref.fetch(read.reference_name)
        expected_seq = query_seq[:read.qstart] +ref_seq[read.reference_start: read.reference_end] + \
                       query_seq[read.qend:]
 
        expected_seq_dict[read.query_name] = expected_seq

        if FASTA_OUT:
            fasta.write(">" + read.query_name + "\n" + expected_seq + "\n")
    if FASTA_OUT:
        fasta.close()
        print(".fasta file saved!")
    if PICKLE_OUT:
        pickle_dump(expected_seq_dict, PICKLE_OUT)

        