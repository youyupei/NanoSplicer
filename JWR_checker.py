'''
    Finding junctions within reads (JWRs) from a bam file
    Input:
        bam_filename <str>:
            filename of the BAM
        output_name <str> <default: "JWR_full_list.hd5">:
            the filename of the output table.
    Output:
        A HDF5 with given output_name will be generated with the following 
        columns:
            1. read_id: <str> the id of original read of the JWR
            2. JWR_intron_start: <int> location of intron start in the 
                                    reference genome supported by the JWR.
            3. JWR_intron_end: <int> location of intron end in the 
                                    reference genome supported by the JWR.
            4. JAQ: <float> junction alignment quality

    Output:
        BED file
        1.chrom 
        2.chromStart: mapped 5' side splice site
        3.chromEnd: mapped 5' side splice site
        4.name: read_id
        5.score: JAQ
        6.

'''
import h5py
import pandas as pd
import pysam
import os
import re
import numpy as np
import sys
import getopt
# modify the output HDF5 file
class JWR_to_hdf(h5py.File):    
    def __init__(self, Filename):
        super().__init__(Filename, 'r+')
    def add_setting_info(self, name, value):
        name = 'setting/' + name
        self[name] = value

# Get information from bam
class JWR_from_reads:
    def __init__(self, read_iter):
        '''
        read_iter: read iterator from pysam. e.g. from .fetch function
        '''
        self.read_iter = read_iter
    def get_JWR(self):
        """
        Input:
            read: AlignedSegment object from pysam
        Return:
            list of introns [(start, stop),...]
            Listing the intronic sites in the reads (identified by 
            'N' in the cigar strings).
        """
        for read in self.read_iter:
            match_or_deletion = {0, 2, 7, 8} # only M/=/X (0/7/8) and D (2) are related to genome position
            ref_skip = 3
            base_position = read.pos
            for op, nt in read.cigartuples:
                if op in match_or_deletion:
                    base_position += nt
                elif op == ref_skip:
                    junc_start = base_position
                    base_position += nt
                    yield JWR_class(read, (junc_start, base_position))
        #             introns.append((junc_start, base_position))
        # return introns

class JWR_class:
    def __init__(self, read, loc):
        '''
        Define a junction within read (JWR)
        Input:
            read: pysam AlignedSegment
            loc: tuple of intron location corresponding to the JWR
        '''
        self.read = read
        self.loc = loc
    def get_JAQ(self, half_motif_size=25):
        '''
        Report the junction alignment quality
        '''
        junction_cigar = \
            self.get_junction_cigar(half_motif_size=25)
        junction_alignment_quality =\
            self.get_junc_map_quality(junction_cigar)
        return junction_alignment_quality
    def get_junction_cigar(self, half_motif_size):
        '''
        Return the cigar letter of around the JWR
        '''
        cigar_long = []
        for count, op in re.findall('(\d+)([A-Za-z=])', self.read.cigarstring):
            cigar_long += int(count) * [op]
        junc1_rel_read_start = self.loc[0] - self.read.reference_start
        junc2_rel_read_start = self.loc[1] - self.read.reference_start
        junction_cigar1 = ''
        junction_cigar2 = ''
        ref_index = -1
        for i in cigar_long:
            if i in 'MND=X':
                ref_index += 1
            if ref_index >= max(0, junc1_rel_read_start - half_motif_size) \
                    and ref_index < junc1_rel_read_start \
                    and i != 'N':
                junction_cigar1 +=  i
            if  ref_index >= junc1_rel_read_start \
                    and min(junc2_rel_read_start + half_motif_size,
                        self.read.reference_end - self.read.reference_start) \
                    and i != 'N':
                junction_cigar1 +=  i
            if ref_index >= min(junc2_rel_read_start + half_motif_size,
                                    self.read.reference_end - self.read.reference_start):
                break
        return junction_cigar1[-25:] + junction_cigar1[:25]
    def get_junc_map_quality(self,cigar):
        '''
        The junc map quality is simply as the proportion of 'M's within the cigar string
        '''
        if not cigar:
            return np.nan
        elif 'M' in cigar:
            return cigar.count('M')/len(cigar)
        elif '=' in cigar:
            return cigar.count('=')/len(cigar)
        else:
            return 0


# import os 
# PWD = '/data/cephfs/punim0614/yupei/deep_sequins/analysis/pipeline/NanoSplicer_11.3/NanoSplicer_run/seg_version'
# os.chdir(PWD)

# a = pysam.AlignmentFile('genome_map_sequins_barcode01.sorted.bam')
# a  = JWR_from_bam('genome_map_sequins_barcode01.sorted.bam')

# if 'N' in r.cigarstring:
#     last_read_pos = False
#     for read_loc, genome_loc in r.get_aligned_pairs():
#         if read_loc is None and last_read_pos:
#             start = genome_loc
#         elif read_loc and last_read_pos is None:
#             stop = genome_loc  # we are right exclusive ,so this is correct
#             res[(start, stop)] += 1
#             del start
#             del stop
#         last_read_pos = read_loc

# def introns_in_read(read):
#     """
#     Input:
#         read: AlignedSegment object from pysam
#     Return:
#         list of introns [(start, stop),...]
#         Listing the intronic sites in the reads (identified by 
#         'N' in the cigar strings).
#     """
#     match_or_deletion = {0, 2, 7, 8} # only M/=/X (0/7/8) and D (2) are related to genome position
#     ref_skip = 3
#     base_position = read.pos
 
#     introns = []
#     for op, nt in read.cigartuples:
#         if op in match_or_deletion:
#             base_position += nt
#         elif op == ref_skip:
#             junc_start = base_position
#             base_position += nt
#             introns.append((junc_start, base_position))
#     return introns

class JWR_checker_param:
    def __init__(self, arg = sys.argv):
        self.arg = arg
            argv = sys.argv
        
        try: 
            opts, args = getopt.getopt(argv[1:],"hw:",
                        ["help", 
                         "window=",
                         "chrID=",
                         "genome-loc="])
        except getopt.GetoptError:
            print("ERROR:Invalid input.")
            print_help()
            sys.exit(1)
        

        # DEFAULT VALUE
        self.junc_cigar_win = 25
        self.bamfile, self_outfile = args
        self.chrID, self.g_loc = None, None


        for opt, arg in opts:
            if opt in ("-h", "--help"):                                                                                                            ('-h', "--help"):
                print_help()
                sys.exit(0)
            elif opt in ("-w", "--window"):
                self.junc_cigar_win = int(arg)
            elif opt == "chrID=":
                self.chrID = arg
            elif opt == "genome-loc=":
                self.g_loc =\
                    tuple([int(i) for i in arg.split('-')])

    def print_help(self):
        help_message =\
        '''
        Usage: python {} [OPTIONS] <BAM file> <output file>
        Options:
            -h/--help        Print this help text
            -w/--window      Candidate searching window size <default: 25>
            --chrID          Target on specific chromosome, chrID should match
                                the chromosome name in the BAM
            --genome-loc     Target on specific genome region, chrID should be 
                                specified. e.g. --genome-loc=0-10000
        '''.format(argv[0])

        print(textwrap.dedent(help_message))

def main():
    param = JWR_checker_param()
    algn_file = pysam.AlignmentFile(param.bamfile)
    reads_fetch = algn_file.fetch(param.chrID,
                                    param.g_loc[0], 
                                    param.g_loc[1])
    JWR_fetch = JWR_from_reads(reads_fetch)
    JWRs = JWR_fetch.get_JWR()

    for JWR in JWRs:
        JWR.get_JAQ(param.junc_cigar_win)

if __name__ == "__main__":
    main()