# some states

import pysam
from intervaltree import IntervalTree
from collections import defaultdict

result_fn = "NanoSplicer_out/output.tsv"

report_f = open("NanoSplicer_out/Analysis_report.txt", 'w')
result_f = open("NanoSplicer_out/output.tsv", 'r')


def parse_result_line(line):
    minimap, candidates, NanoSplicer, min_score = line.split()
    minimap = tuple([int(x) for x in minimap.split(',')])
    #minimap = (minimap[0], minimap[1]-2)
    candidates = [int(i) for i in candidates.split(',')]
    candidates = [(candidates[i], candidates[i+1]) for i in range(0,len(candidates),2)]
    NanoSplicer = int(NanoSplicer)
    min_score = float(min_score)
    return minimap, candidates, NanoSplicer, min_score

def minimap_nanosplicer_campare(result_fn):
    with open(result_fn, 'r') as f:
        total = 0
        match = 0
        match_score = []
        mismatch_score =[]
        for line in f:
            total += 1
            minimap, candidates, NanoSplicer, min_score\
                                         = parse_result_line(line)
            if minimap == candidates[NanoSplicer]:
                match += 1
                if min_score < 100:
                    match_score.append(min_score)
            elif min_score < 100:
                mismatch_score.append(min_score)
    
    prop = match/total
    return total, prop, match_score,mismatch_score


# plot min score
import matplotlib.pyplot as plt
total, match, match_score,mismatch_score = \
                    minimap_nanosplicer_campare(result_fn)

plt.hist(match_score + mismatch_score, bins = 10000)
plt.hist(match_score , bins = 10000)
plt.xlim((0,10))
plt.legend(['All min distance (Si)','Si for results that match minimap2'])
plt.xlabel('Si')
plt.ylabel('count')
plt.title('Distribution of Si')
plt.savefig("histgram_of_S_new.png")
plt.close()

def to_short_read_compare(result_fn, find_intron_result, min_score_t = 4, short_thres = 0):
    with open(result_fn, 'r') as f:
        total = 0
        minimap_match = 0
        NanoSplicer_match = 0
        number_of_short_support = []
        for line in f:
            minimap, candidates, NanoSplicer,min_score\
                                 = parse_result_line(line)
            if min_score > min_score_t:
                continue
            total += 1
            if find_intron_result.get(minimap, 0) > short_thres:
                minimap_match += 1
            if find_intron_result.get(candidates[NanoSplicer], 0) > short_thres:
                NanoSplicer_match += 1
            number_of_short_support.append(len(set(candidates).intersection(set(find_intron_result.keys()))))
    return total, minimap_match, minimap_match/total, NanoSplicer_match #,number_of_short_support






BAM_FN = 'BAM/Chr_ID/NC_000001.11.bam'
SHORT_READ_no_anno = '/data/cephfs/punim0614/shared/shared_data/'\
                    'external_public/RNASeqMixology/splice_site_analysis/'\
                    'GRCh38_bam_no_annotation/SRR1706863_1Aligned.sortedByCoord.out.bam'

SHORT_READ_anno = '/data/cephfs/punim0614/shared/shared_data/external_public/'\
                    'RNASeqMixology/splice_site_analysis/GRCh38_bam_with'\
                    '_ncbi_annotation/SRR1706783_1Aligned.sortedByCoord.out.bam'
REF_FN = '/data/cephfs/punim0614/shared/shared_data/external_public/' \
            'RNASeqMixology/splice_site_analysis/GRCh38_latest_genomic.fna'
REF_FN_IDX = '/data/cephfs/punim0614/shared/shared_data/external_public/' \
            'RNASeqMixology/splice_site_analysis/GRCh38_latest_genomic.fna.fai'

fbam = pysam.AlignmentFile(BAM_FN)
f_short_anno = pysam.AlignmentFile(SHORT_READ_anno)
f_short_no_anno = pysam.AlignmentFile(SHORT_READ_no_anno)
fref = pysam.FastaFile(REF_FN)

chrID = 'NC_000001.11'
f_fetch = fbam.fetch() # reads iterator
f_short_anno_fetch = f_short_anno.fetch(chrID)
f_short_no_anno_fetch = f_short_no_anno.fetch(chrID)

f_introns = fbam.find_introns(f_fetch) # dictionary
f_short_anno_fetch_intron = f_short_anno.find_introns(f_short_anno_fetch)
f_short_no_anno_fetch_intron = f_short_no_anno.find_introns(f_short_no_anno_fetch)


to_short_read_compare(result_fn, f_short_anno_fetch_intron, 4,15)


# acc vs S_i threshold
import numpy as np
total_arr = []
minimap_arr = []
NanoSplicer_arr = []
sum_of_squiggle =  len(list(open(result_fn, 'r')))
for i in np.arange(1.5, 10, 0.5):
    total, minimap_match, minimap_prop, NanoSplicer_match = \
            to_short_read_compare(result_fn, f_short_anno_fetch_intron,i)
    total_arr.append(total/sum_of_squiggle)
    minimap_arr.append(minimap_match/total)
    NanoSplicer_arr.append(NanoSplicer_match/total)

plt.plot(np.arange(1.5, 10, 0.5), total_arr)
plt.plot(np.arange(1.5, 10, 0.5), minimap_arr)
plt.plot(np.arange(1.5, 10, 0.5), NanoSplicer_arr)
plt.legend(['total', 'minimap', 'NanoSplicer'])
plt.xlabel("S_i threshold")
plt.ylabel("proportion of identified junctions have short read support")
plt.savefig("acc_plot.png")
plt.close()

# acc vs num short read support threshold

import numpy as np
from tqdm import tqdm
total_arr = []
minimap_arr = []
NanoSplicer_arr = []
sum_of_squiggle =  len(list(open(result_fn, 'r')))
for i in np.arange(0, 20, 2 ):
    total, minimap_match, minimap_prop, NanoSplicer_match = \
            to_short_read_compare(result_fn, f_short_anno_fetch_intron,4, i)
    total_arr.append(total/sum_of_squiggle)
    minimap_arr.append(minimap_match/total)
    NanoSplicer_arr.append(NanoSplicer_match/total)

plt.plot(np.arange(0, 20), total_arr)
plt.plot(np.arange(0, 20), minimap_arr)
plt.plot(np.arange(0, 20), NanoSplicer_arr)
plt.legend(['total', 'minimap', 'NanoSplicer'])
plt.xlabel("munimum short read support")
plt.ylabel("proportion of identified junctions have short read support")
plt.savefig("acc_plot_short_support.png")
plt.close()


# proportion correlation

def to_short_read_abundence(result_fn, find_intron_result, min_score_t = 4, short_thres = 0):
    with open(result_fn, 'r') as f:
        total = 0
        minimap_match = defaultdict(0)
        NanoSplicer_match = defaultdict(0)
        total_number_of_short_support = 
        
        for line in f:
            minimap, candidates, NanoSplicer,min_score\
                                 = parse_result_line(line)
            if min_score > min_score_t:
                continue
            total += 1

            minimap_match[minimap] += find_intron_result.get(minimap, 0)
            NanoSplicer_match[candidates[NanoSplicer]] += find_intron_result.get(candidates[NanoSplicer], 0)

            number_of_short_support.append(len(set(candidates).intersection(set(find_intron_result.keys()))))
    return total, minimap_match, minimap_match/total, NanoSplicer_match #,number_of_short_support








f_short_n_intron_filter = dict()
for key, value in f_short_no_anno_fetch_intron.items():
    if value > 10:
        f_short_n_intron_filter[key] = value

f_short_intron_filter = dict()
for key, value in f_short_anno_fetch_intron.items():
    if value > 10:
        f_short_intron_filter[key] = value


count_short_anno_match = 0
count_short_no_anno_match = 0
for i in f_introns.keys():
    if f_short_intron_filter.get(i, None):
        count_short_anno_match += f_introns[i]
    if f_short_intron_filter.get(i, None):
        count_short_no_anno_match += f_introns[i]
        

intron_tree = IntervalTree()
for (begin, end), data in f_introns.items():
    intron_tree.addi(begin, end, data)

# write total number
total_jwr = sum([interval.data for interval in intron_tree])
report_f.write("Total number of junction within read:{}\n".format(
               total_jwr ))

report_f.write("Total number of jwr match short read (annotation guided):{}({}%)\n".format(
                count_short_anno_match, 100* count_short_anno_match/total_jwr))

report_f.write("Total number of jwr match short read (no guided):{}({}%)\n".format(
                count_short_no_anno_match, 100* count_short_no_anno_match/total_jwr))

report_f.close()
test = fref.fetch(chrID, 100000,100500)

fbam.check_index()
fbam.count()   




# test transcript direction



#built non overlapped range
intron_tree_non_overlaps = intron_tree.copy()
intron_tree_non_overlaps.merge_overlaps()

count = 0
single_support = 0
total_candidate = 0
total_junc_within_read = 0
for interval in intron_tree_non_overlaps:
    candidates = find_candidate(interval.begin, interval.end, intron_tree, 10)
    if candidates:
        print(candidates)
    # statistics
    total_candidate += len(candidates)
    total_junc_within_read += sum([x.data for x in candidates])
    single_support += len([x.data for x in candidates if x.data < 30])
    count += 1 
    if count < 0:
        break


f_fetch = fbam.fetch('NC_000001.11',start = 29200000, stop = 30000000)
f_pileup = fbam.pileup('NC_000001.11',start = 14000, stop = 15000, truncate=True)
f_pileup = fbam.pileup(contig='NC_000001.11',start = 14275, stop = 15000,truncate=True, stepper = "nofilter")

read_a = next(f_fetch,False)
pileup_a = next(f_pileup,False)
pileup_a.reference_pos
pileup_a.get_num_aligned()



report_f.close()