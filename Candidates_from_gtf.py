'''
Extract a list of exon candidates
input:
    1 ref genome
    2 annotation bed file
    3 splic bed file from mapping
'''
import numpy as np
import itertools, sys, os, re
from collections import defaultdict
import helper

def ReadBedLine(bedline):
    """
    Take a bed12 format line as input and return the corresponding splicing
    junction position (0-based and genome related)

    return:
        name: transcript name
        strand
        splice_sites list of tuple (genome 0-based position)
        junction_pos list (transript 0-based position)
    """
    bedline = bedline.strip().split()
    chrom, chromStart, chromEnd, name, score, strand, thickStart,\
     thickEnd, itemRgb, blockCount, blockSizes,\
      blockStarts = bedline
    
    blockSizes = [int(x) for x in blockSizes.strip(',').split(',')]
    blockStarts =[int(chromStart) + int(x) for x in blockStarts.strip(',').split(',')]
    blockEnds = [x + y for x,y in zip(blockSizes, blockStarts)]

    splice_sites = \
        [(blockEnds[i], blockStarts[i + 1]) for i in range(int(blockCount) - 1)]

    junction_pos = []
    for start_pos in blockStarts:
        exon_transcript_pos = start_pos - int(chromStart)
        for intron_pos in splice_sites:
            if start_pos >=  intron_pos[1]:
                exon_transcript_pos -= intron_pos[1] - intron_pos[0]
        junction_pos.append(exon_transcript_pos)
    
    junction_pos = junction_pos[1:]

    return name, strand, splice_sites, junction_pos
        


def splice_site_to_seq(splice_site, genome_ref, window = 30 , strand = "+", start = None, end = None): 
    '''
    input:
        splice_site: a splice site <tuple>
        genome_ref: <str>
    '''

    junction_bases_start = splice_site[0] - int(np.floor(window/2))
    junction_bases_end = splice_site[1] + int(np.ceil(window/2))

    if start:
        junction_bases_start = start
    if end: 
        junction_bases_end = end

    if junction_bases_start < 0 or junction_bases_end > len(genome_ref) + 1:
        print("Error: Specified Exon Junction goes outside of the Genome!\n")
        return 1
    return genome_ref[junction_bases_start:splice_site[0]] + \
                genome_ref[splice_site[1] : junction_bases_end]
    '''
    if strand == "+":
        return genome_ref[junction_bases_start:splice_site[0]] + \
                genome_ref[splice_site[1] : junction_bases_end]
        
    elif strand == "-":
        return helper.reverse_complement(
            genome_ref[junction_bases_start:splice_site[0]] + \
            genome_ref[splice_site[1] : junction_bases_end])
    '''
def check_exon(pos1, bedline):
    bedline = bedline.strip().split()
    chrom, chromStart, chromEnd, name, score, strand, thickStart,\
     thickEnd, itemRgb, blockCount, blockSizes,\
      blockStarts = bedline
    
    blockSizes = [int(x) for x in blockSizes.strip(',').split(',')]
    blockStarts =[int(chromStart) + int(x) for x in blockStarts.strip(',').split(',')]
    blockEnds = [x + y for x,y in zip(blockSizes, blockStarts)]

   
    junction_pos = []
    for exon_start, exon_end in zip(blockStarts, blockEnds):
        if pos1 <= exon_end and pos1 >= exon_start:
            return exon_start, exon_end

    print("Warning!!")

def search_potential_canonical_sites(splice_site, window, genome_ref, strand):
    search_start_0 = splice_site[0] - int(np.floor(window/2))
    search_end_0 = splice_site[0] + int(np.ceil(window/2))
    search_start_1 = splice_site[1] - int(np.floor(window/2))
    search_end_1 = splice_site[1] + int(np.ceil(window/2))
    motif0 = genome_ref[search_start_0: search_end_0]
    motif1 = genome_ref[search_start_1: search_end_1]
    if strand == "+":
        list0 = [search_start_0 + m.start() for m in re.finditer("GT",motif0)]
        list1 = [search_start_1 + m.start()+2 for m in re.finditer("AG",motif1)]

    if strand == "-":
        list0 = [search_start_0 + m.start() for m in re.finditer("CT",motif0)]
        list1 = [search_start_1 + m.start()+2 for m in re.finditer("AC",motif1)]
    
    return list(itertools.product(list0,list1))


def count_splice_site_supports(bedfile, ref_trans_strand):
    '''
    input:
        bedfile converted from bam/sam
    return:
        mapped splicing sites and number of supports
    '''
    
    count_intron_start = defaultdict(int)
    count_intron_end = defaultdict(int)
    count_intron_mapping = defaultdict(int)
    with open(bedfile, "r") as bf:
        for line in bf:
            name, mapped_strand, splice_sites, junction_pos = ReadBedLine(line)
            #if ref_trans_strand != mapped_strand:
            #    continue
            for site in splice_sites:
                count_intron_start[site[0]] += 1
                count_intron_end[site[1]] += 1
                count_intron_mapping[site] += 1

    return count_intron_start, count_intron_end, count_intron_mapping


def genome_pos_to_transcript_pos(genome_pos, bedline):
    bedline = bedline.strip().split()
    chrom, chromStart, chromEnd, name, score, strand, thickStart,\
     thickEnd, itemRgb, blockCount, blockSizes,\
      blockStarts = bedline
    
    blockSizes = [int(x) for x in blockSizes.strip(',').split(',')]
    blockStarts =[int(chromStart) + int(x) for x in blockStarts.strip(',').split(',')]
    blockEnds = [x + y for x,y in zip(blockSizes, blockStarts)]

    splice_sites = \
        [(blockEnds[i], blockStarts[i + 1]) for i in range(int(blockCount) - 1)]

    transcript_pos = genome_pos - int(chromStart) 
    for intron_pos in splice_sites:
        if genome_pos >=  intron_pos[1]:
            transcript_pos -= intron_pos[1] - intron_pos[0]
 
    
    return transcript_pos

class candidate_class:
        def __init__(self, sequences=None, start=-1, end = -1, num_of_correct_supports = 0):
            self.sequences = sequences
            self.start = int(start)
            self.end = int(end)
            self.num_of_correct_supports = num_of_correct_supports
 



def main():
    args = sys.argv
    if len(args) < 2:
        print("\n\nUsage: python3 {} <ref filename> <bed filename from bam>".format(argv[0]))
        exit()
    
    #gtfbedline = sys.stdin
    annotation_bedline = next(sys.stdin)
    
    genome_ref_filename, bambedfile = args[1:3]
    with open(genome_ref_filename, 'r') as gf:
        next(gf)
        genome_ref = next(gf)

      
    name, strand, annotated_sites, junction_pos = ReadBedLine(annotation_bedline)

    count_intron_start, count_intron_end, count_intron_mapping=\
        count_splice_site_supports(bambedfile, strand)

    candidates = []
    for site in annotated_sites:
        candidate_splice_site_list = [site]
        # best supported sites:
        search_win = 20
        flank_size = 10 #each side
        accept_thres = 3

        best_supported_count = accept_thres
        best_supported_site = []
        
        
        for counted_site in itertools.product(range(site[0] - int(np.floor(search_win/2)),site[0] + 1),\
              range(site[1], site[1] + int(np.ceil(search_win/2)) + 1)):
        
            if counted_site == site:
                continue
            
            if count_intron_mapping[counted_site] >= best_supported_count:
                best_supported_count = count_intron_mapping[counted_site]
                best_supported_site.append(counted_site)


        candidate_splice_site_list += best_supported_site +\
         search_potential_canonical_sites(site, window = search_win, genome_ref = genome_ref, strand = strand)


        candidate_splice_site_list = list(set(candidate_splice_site_list))
        true_index = candidate_splice_site_list.index(site)
        candidate_splice_site_list[0], candidate_splice_site_list[true_index] =\
        candidate_splice_site_list[true_index], candidate_splice_site_list[0]



        # generate candidate
        
        candidate = candidate_class(num_of_correct_supports = count_intron_mapping[site])
        candidate.start = min([x for x,y in candidate_splice_site_list]) - flank_size
        candidate.end = max([y for x,y in candidate_splice_site_list]) + flank_size

        if check_exon(candidate.start, annotation_bedline) != check_exon(site[0], annotation_bedline):
            candidate.start = check_exon(site[0], annotation_bedline)[0]
        
        if check_exon(candidate.end, annotation_bedline) != check_exon(site[1], annotation_bedline):
            candidate.end = check_exon(site[1], annotation_bedline)[1]

        candidate.sequences = [splice_site_to_seq(x, genome_ref, \
        start = candidate.start, end = candidate.end, strand = strand)\
         for x in candidate_splice_site_list]
        
        # convert to transcript related pos
        candidate.start = genome_pos_to_transcript_pos(candidate.start, annotation_bedline)
        candidate.end = genome_pos_to_transcript_pos(candidate.end, annotation_bedline)
        
        candidates.append(candidate)

    if strand  == '+':

        for candidate in candidates:
            print(','.join(candidate.sequences) + ",{},{},{}".format( candidate.start,\
                candidate.end,candidate.num_of_correct_supports))
    
    elif strand == '-':
        for candidate in candidates:
            print(','.join([helper.reverse_complement(s) for s in candidate.sequences]) + ",{},{},{}".format( -candidate.start,\
                -candidate.end,candidate.num_of_correct_supports))

    
    return None


if __name__ == '__main__':
    main()