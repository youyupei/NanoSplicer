from collections import defaultdict


def transcript_dict_generator(filename):
    with open(filename, 'r') as f:
        gene_dict = defaultdict(dict)
        transcript_id = defaultdict(list)
        for line in f:
            if 'exon' not in line:
                continue
            line = line.split()
            gene_dict[line[9]][line[11]] = transcript_id[line[11]]
            transcript_id[line[11]].append((int(line[3]), int(line[4])))
    return(gene_dict)
    '''
    Arg:
        filename: gtf file
    returns:
        dictionary: {gene_id:{transcript_id:[(exon1_start, exon1_end),...]}}
    '''

def main():
    gtf_filename = "/home/ubuntu/data/Sequin_resources/rnasequin_annotation_2.2.gtf"
    gene_dict = transcript_dict_generator(gtf_filename)
    for gene in gene_dict.keys():
        exon_boundary = []
        for trans in gene_dict[gene]:
            exon_boundary += gene_dict[gene][trans]
        uniq_exon_boundary = sorted(set(exon_boundary))
        for i in range(len(uniq_exon_boundary)-1):
            if uniq_exon_boundary[i][1] > uniq_exon_boundary[i+1][0]:
                print(gene,uniq_exon_boundary[i][1]) 



if __name__ == '__main__':
    main()
    
