
import regex as re

def motif_finder(motif, ref_file):
    '''
    Input: <spike-in reference> and <queried motif>
    Output: 0-based index of occurances of motif in reference
    '''
    with open(ref_file, "r") as ref:
        index_dic = {}
        '''
        Key: Chromosome 
        Value: list 0-based start indexes of the motif in ref
        ''' 
        ref = ref.readline()
        for line in ref:
            if line[0] == ">":
                chrom_name = line[1:].strip()
                ref_seq = next(ref).strip()
                index_dic[chrom_name]=[match.start() for match 
                            in re.finditer(motif, ref_seq,
                             overlapped=False)]         
            else:
                print("Warning! Chromosome skipped")
                next(ref)
    return index_dic

def main():
    import sys
    argv = sys.argv
    if len(argv < 3):
        exit("Usage: python {} <motif(string)> <ref_file>".format(argv[0]))
    else:
        motif = argv[1]
        ref_file = argv[2]
        motif_finder(motif, ref_file)


if __name__ == "__main__":
    main()

    
