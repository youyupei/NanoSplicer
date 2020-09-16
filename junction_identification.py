import pysam
import intervaltree
from intervaltree import IntervalTree
import re
import itertools

# read .bam file
def get_intron_tree(pysamAlignment, chrID):
    '''
        Read bam file and generate IntervalTrees with format: 
        (intron start, intron end, num of support)
        Parameter:
            bam_filename:
                <string> filename for the bam file.
    '''
    #pysamAlignment = pysam.AlignmentFile(bam_filename)
    f_fetch = pysamAlignmentClass.fetch(chrID) # reads iterator
    f_introns = pysamAlignmentClass.find_introns(f_fetch) # dictionary

    # converty dictionary to IntervalTree
    intron_tree = IntervalTree()
    for (begin, end), data in f_introns.items():
        intron_tree.addi(begin, end, data)

    #built non overlapped range
    intron_tree_non_overlaps = intron_tree.copy()
    intron_tree_non_overlaps.merge_overlaps()

    count = 0
    single_support = 0
    total_candidate = 0
    total_junc_with_in_read = 0
    for interval in intron_tree_non_overlaps:
        candidates = find_candidate(interval.begin, interval.end, intron_tree, 10)
        
        # some statistics
        total_candidate += len(candidates)
        total_junc_with_in_read += sum([x.data for x in candidates])
        single_support += sum([x.data for x in candidates if x.data < 3])
        count += 1 
        if count < 0:
            break

# get alignment given genome pos
def find_candidate(begin, end, tree, window=10, min_primary = 0, 
                   min_support=0, secondary_thres=0.0, primary_thres=1.0):
    '''
    Find candidate exon boundary (i.e. intron boundary) within a given range.
    Parameter:
        begin:
            start (left-most) position of the range to be searched (0-based)
        end:
            end (right-most) possition of the range to be searched (0-based)
        tree:
            IntervalTree containing all boundary pairs 
        window: 
            window size for group surrounding boundaries (difference 
            of boundary in each size of the intron will be grouped together if 
            the absolute difference < window size)
        min_support:
            The best supported boundary need will be included only when the num
            of support reaches the minimum
        secondary_thres:
            only the junctions with multiple well supported boundary will
            be included. Well supported junction is defined as 
            secondary_thres * support num of the most supported boundary.
    '''
    # get boundaries with in searching window, sorted by the number of support
    
    intervals_tree = IntervalTree(tree.envelop(begin, end))
    candidate_boundaries = []
    while intervals_tree:
        interval = max(intervals_tree, key = lambda x: x.data)
        best_support = interval.data
        if interval.data < min_primary: # lower bound of the support
            return candidate_boundaries
            
        #candidate_boundaries.append(interval)
        intervals_tree.remove(interval)
        
        # include surrounding boundaries
        enveloped_interval = intervals_tree.envelop(interval.begin - window, 
                                          interval.end + window)
        neighbour_found = []
        for i in enveloped_interval:
            if i.begin <= interval.begin + window and \
                    i.end >= interval.end - window:
                if i.data > secondary_thres * best_support:        
                    neighbour_found.append(i)
                intervals_tree.remove(i)
        if neighbour_found:
            neighbour_found.append(interval)
            count = sum([x.data for x in neighbour_found])
            if count >= min_support and best_support/count <= primary_thres:
                candidate_boundaries += neighbour_found
    return candidate_boundaries

def canonical_site_finder(candidate_Interval, ref_FastaFile,
                              AlignmentFile, window, chrID):
    '''
        Generate candidate motif given the position of intron boundary.
        Input:
            candidate_Interval:
                tuple of intron boundary position
            ref_FastaFile:
                pysam.Fastafile class read from reference genom
            AlignmentFile:
                pysam.AlignmentFile class read from bam/sam file
            window:
                size of window in which NanoSplicer search for the candidate
            chrID: 
                Chromosome ID name.
    '''
    start, end, data = candidate_Interval
    # test only
    # return [start], [end]
    
    # check determined transcript direction (by minimap2)
    donor_pattern = ref_FastaFile.fetch(chrID, start - window, \
                                        start + window).upper()
    acceptor_pattern = ref_FastaFile.fetch(chrID, end - window, \
                                        end + window).upper()
    
    if donor_pattern[window:window + 2] == "GT" and \
        acceptor_pattern[window -2 : window] == "AG":
        intron_start_candidate = [start - window + m.start() 
                                    for m in re.finditer("GT",donor_pattern)]
        
        intron_end_candidate = [end - window + m.start() + 2
                    for m in re.finditer("AG",acceptor_pattern)]
    
    elif donor_pattern[window:window + 2] == "CT" and \
        acceptor_pattern[window -2 : window] == "AC":
        intron_start_candidate = [start - window + m.start() 
                                    for m in re.finditer("CT",donor_pattern)]

        intron_end_candidate = [end - window + m.start() + 2
                    for m in re.finditer("AC",acceptor_pattern)]
        
    else:
        return None, None
    #print(intron_start_candidate, intron_end_candidate) 
    return intron_start_candidate, intron_end_candidate



def candidate_motif_generator(chrID, intron_start_candidate, 
                              intron_end_candidate, 
                              flank_size, ref_FastaFile):
    '''
        generate candidate motif (genome forward strand) with flanking seq
    '''
    if not intron_start_candidate or not intron_end_candidate:
        return None, None, None, None
    motif_start_pos = min(intron_start_candidate) - flank_size
    motif_end_pos = max(intron_end_candidate) + flank_size

    candidates_tuple = list(itertools.product(intron_start_candidate, 
                                              intron_end_candidate))

    motif_list = []
    for donor, acceptor in candidates_tuple:
        motif_list.append(ref_FastaFile.fetch( 
                                    chrID, motif_start_pos, donor).upper() + 
                          ref_FastaFile.fetch( 
                                    chrID, acceptor, motif_end_pos).upper())

    return candidates_tuple, motif_list, motif_start_pos, motif_end_pos


def main():
    return None

## test part
#import pysam
#
#BAM_FN = 'BAM/Chr_ID/NC_000001.11.bam'
#REF_FN = '/data/cephfs/punim0614/shared/shared_data/external_public/' \
#            'RNASeqMixology/splice_site_analysis/GRCh38_latest_genomic.fna'
#REF_FN_IDX = '/data/cephfs/punim0614/shared/shared_data/external_public/' \
#            'RNASeqMixology/splice_site_analysis/GRCh38_latest_genomic.fna.fai'
#
#fbam = pysam.AlignmentFile(BAM_FN)
#fref = pysam.FastaFile(REF_FN)
#
#chrID = 'NC_000001.11'
#test = fref.fetch(chrID, 100000,100500)
#
#fbam.check_index()
#fbam.count()
#f_fetch = fbam.fetch() # reads iterator
#f_introns = fbam.find_introns(f_fetch) # dictionary
#
#intron_tree = IntervalTree()
#for (begin, end), data in f_introns.items():
#    intron_tree.addi(begin, end, data)
#
## test transcript direction
#
#
#
##built non overlapped range
#intron_tree_non_overlaps = intron_tree.copy()
#intron_tree_non_overlaps.merge_overlaps()
#
#count = 0
#single_support = 0
#total_candidate = 0
#total_junc_within_read = 0
#for interval in intron_tree_non_overlaps:
#    candidates = find_candidate(interval.begin, interval.end, intron_tree, 10)
#    if candidates:
#        print(candidates)
#    # statistics
#    total_candidate += len(candidates)
#    total_junc_within_read += sum([x.data for x in candidates])
#    single_support += len([x.data for x in candidates if x.data < 30])
#    count += 1 
#    if count < 0:
#        break
#
#
#f_fetch = fbam.fetch('NC_000001.11',start = 29200000, stop = 30000000)
#f_pileup = fbam.pileup('NC_000001.11',start = 14000, stop = 15000, truncate=True)
#f_pileup = fbam.pileup(contig='NC_000001.11',start = 14275, stop = 15000,truncate=True, stepper = "nofilter")
#
#read_a = next(f_fetch,False)
#pileup_a = next(f_pileup,False)
#pileup_a.reference_pos
#pileup_a.get_num_aligned()


if __name__ == "__main__":
    main()