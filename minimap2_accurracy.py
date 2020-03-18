# Script for calculating the splicing site mapping accuracy

import sys 
import os 
import numpy as np 
import regex
import csv



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


def add_col_to_csv(col, input, output):
    """
    adding a column at the end of a csv file
    aug:
        col: list with length = nrow(input)
        input: str, filename
        output: str filename
    return:
        None
    """
    try:
        with open(input,'r') as csvinput:
            with open(output, 'w') as csvoutput:
                writer = csv.writer(csvoutput, lineterminator='\n')
                reader = csv.reader(csvinput)
                try:
                    all = [ [b] + a for a,b in zip(reader,col)]
                    writer.writerows(all)
                except:
                    print("ERROR in adding column to csv!")
                    return(1)
    except:
        print("File missing")
        return(1)

def main():
    '''
    red bedline from transID.bed
    campare transID.bed with annotated.bed

    return:
        0/1 vector for each read
    '''
    bed_files = os.listdir("transcript_ID")
    annotated_bed = "annotation.bed"
    csv_dir = "/home/ubuntu/PhD_proj/pipeline_1.21/Validation/score_f20_ts6_t2_spikeT3_band_0.4_stand_likelihood_more_info_2_19"
    out_dir = "/home/ubuntu/PhD_proj/pipeline_1.21/Validation/score_f20_ts6_t2_spikeT3_band_0.4_stand_likelihood_more_info_2_19_minimap"
    
    for filename in bed_files:
        transID= filename[:-4]
        filename = "transcript_ID/" + filename

        annotated_junctions = []
        with open(annotated_bed, 'r') as anno_f:
            for line in anno_f:
                if transID in line:
                    annotated_junctions, strand = ReadBedLine(line)[2], ReadBedLine(line)[1]



        with open(filename, "r") as f:
            for line in f:
                mapping_correctness = []

                name, read_junctions =ReadBedLine(line)[0],  ReadBedLine(line)[2]
                for site in annotated_junctions:
                    mapping_correctness.append(int(site in read_junctions))
                    if site not in annotated_junctions:
                        pass

                #if strand == "-":
                #    mapping_correctness = mapping_correctness[::-1]

                add_col_to_csv(col = mapping_correctness,\
                 input = "{}/{}_{}.csv".format(csv_dir,transID,name),\
                 output = "{}/{}_{}.csv".format(out_dir,transID,name))
    
    return None

if __name__ == "__main__":
    main()
