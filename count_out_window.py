'''
count number of reads with mapped junction outside the sepcified window
'''

import sys
import os
import re
#bash command
#for i in $(ls); do transID=$(echo $i | sed 's/\.bed//g'); grep $transID ../annotation.bed | python3 ~/PhD_proj/Python/count_out_window.py $transID\.bed; done

annotation = next(sys.stdin)

#R1_102_1
#annotation = "chrIS   6790170 6829759 R1_102_1||R1_102        1000    +       6790170 6829759 196,196,196     7       126,332,103,172,124,192,441,    0,4197,8685,9514,12143,31935,39148,"
annotation = annotation.split()
chrom, chromStart, chromEnd, name, score, strand, thickStart,\
     thickEnd, itemRgb, blockCount, blockSizes,\
      blockStarts = annotation
blockSizes = [int(x) for x in blockSizes.strip(',').split(',')]
blockStarts =[int(chromStart) + int(x) for x in blockStarts.strip(',').split(',')]
blockEnds = [x + y for x,y in zip(blockSizes, blockStarts)]

exon_pos = list(zip(blockStarts,blockEnds))
annotated_sites = \
        [(blockEnds[i], blockStarts[i + 1]) for i in range(int(blockCount) - 1)]

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
    transcript_length = sum(blockSizes)
    return name, strand, splice_sites, junction_pos, transcript_length

def boundary_assignment(splice_site, exon_pos):
    site_start, site_end = splice_site
    num_of_boundaries = len(exon_pos) - 1
    if num_of_boundaries == 1:
        return 0
    if num_of_boundaries == 0:
        return -1
    for i in range(num_of_boundaries):
        if i == 0:
            if site_start <= exon_pos[1][0] and site_end >= exon_pos[0][1] and site_end <= exon_pos[1][1]:
                return i
        elif i == num_of_boundaries - 1:
            if site_start <= exon_pos[i+1][0] and \
                site_start >= exon_pos[i][0] and\
                site_end >= exon_pos[i][1]:
                return i
            else:
                return -1
        else:
            if site_start <= exon_pos[i+1][0] and \
                site_start >= exon_pos[i][0] and\
                site_end >= exon_pos[i][1] and \
                site_end >= exon_pos[i+1][1]:
                return i




#bedfile = 'R1_102_1.bed'
bedfile = sys.argv[1]
with open(bedfile,"r") as bedf:
    count_outside = 0
    count_all = 0
    for line in bedf:
        name, mapped_strand, splice_sites, junction_pos, transcript_length \
             = ReadBedLine(line)
        for splice_site in splice_sites:
            site_start, site_end = splice_site
            count_all += 1
            assigned_boundary = boundary_assignment(splice_site,exon_pos)
            if assigned_boundary == -1:
                continue
            try:
                annotated_start, annotated_end = annotated_sites[assigned_boundary]
            except:
                print(annotated_sites, exon_pos)
            if abs(site_start - annotated_start) > 10 or abs(site_end - annotated_end) > 10:
                #print(splice_site,assigned_boundary)
                count_outside += 1

if count_all:
    print(count_outside/count_all)
