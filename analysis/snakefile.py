# this is a Snakefile
import os
from config import *

rule all:
    input:
        "BAM/genome_map{}.sorted.bam".format(SUFFIX),
        "BAM/transcript_map.sorted_Q60.bam"


# STEP 1 Index reference genome
rule index_genome_ref:
    input:
        "{fa_file}"
    output:
        "{fa_file}.fai"
    shell:
        "samtools faidx {input}"


# STEP 2 map to reference genome (minimap2)
rule map_genome:
    input:
        ref = REF_GENOME,
        fastq =FASTQ_PATH,
        ref_idx = REF_GENOME + ".fai"
    params:
        splice_param = MAP_CANONICAL
    output:
        "BAM/genome_map{}.sorted.bam".format(SUFFIX)
    threads: THREAD
    shell:
        '''
        mkdir -p BAM
        minimap2 -ax splice -t {threads} -k 14 {params.splice_param} --splice-flank=no --MD --eqx --secondary=no --sam-hit-only {input.ref} {input.fastq}/* | 
        samtools sort -o {output} &&
        samtools index {output}
        '''

# Identify the truth
    # 1. map to transcript ref
    # 2. filter based on mapQ
rule map_transcript:
   input:
       ref = REF_TRANSCRIPT,
       fastq = expand(FASTQ_PATH + '/{fq}', fq=os.listdir(FASTQ_PATH)),
       ref_idx = REF_TRANSCRIPT + ".fai"
   output:
       primary_align = "BAM/transcript_map.sorted.bam",
       Q60_filtered = "BAM/transcript_map.sorted_Q60.bam"
   threads: THREAD
   shell:
       "mkdir -p BAM;"
       "minimap2 -ax map-ont -t {threads} --secondary=no  {input.ref} {input.fastq} | "
       "samtools sort -o {output.primary_align} -T reads.tmp &&"
       "samtools index {output.primary_align};"
       "samtools view -b -q60 {output.primary_align} > {output.Q60_filtered};"
       "samtools index {output.Q60_filtered};"

rule bam_mapped_filter:
    input:
        "BAM/genome_map{}.sorted.bam".format(SUFFIX)
    output:
        mapped = "BAM/{}.sorted_mapped_subset.bam".format(SUFFIX)
    params:
        map_qual_thres = MIN_MAPQ
    shell:
        '''
        samtools view -b -F 4 -q 60 {input} > {output.mapped}
        samtools index {output.mapped}
        #rm {input}
        #rm {input}.bai
        '''

# # convert annotation into bed format
# rule gtf_to_bed:
#     input:
#         gtf_fn = TRANSCRIPT_ANNOTATION
#     output:
#         bed_fn = "BED/annotation.bed"
#     shell:
#         "mkdir -p BED;" +
#         PAFTOOLS + "gff2bed {input.gtf_fn} > {output.bed_fn}"

# # convert bam into bed format
# rule bam_to_bed:
#     input:
#         genome_bam = "BAM/genome_map.sorted_mapped_subset.bam"
#     output:
#         bed = "BED/BAM2Chr_ID_BED/{Chr_ID}.bed",
#         bam = "BAM/Chr_ID/{Chr_ID}.bam",
#         bai = "BAM/Chr_ID/{Chr_ID}.bam.bai"
#     shell:
#         '''
#         mkdir -p BED/BAM2Chr_ID_BED
#         mkdir -p BAM/Chr_ID
#         mkdir -p tmp

#         # splite bam file by chr id
#         samtools view -b {input.genome_bam} {wildcards.Chr_ID} > {output.bam}
#         samtools index {output.bam}

#         # convert echo .bam to .bed using paftools.js
#         samtools view {output.bam} > tmp/{wildcards.Chr_ID}.sam
#         '''
#         + PAFTOOLS + "splice2bed tmp/{wildcards.Chr_ID}.sam > {output.bed};"
#         "rm tmp/{wildcards.Chr_ID}.sam"

# BAM_FN='/data/cephfs/punim0614/yupei/deep_sequins/analysis/pipeline/BAM/genome_map.sorted_mapped_subset.bam'
# REF_FN='/data/cephfs/punim0614/yupei/deep_sequins/Sequin_resources/rnasequin_decoychr_2.2.fa'
# FAST_DIR='/data/cephfs/punim0614/yupei/deep_sequins/fast5_pass'
# python3 ../script/NanoSplicer2.py -i $BAM_FN -f $FAST_DIR -r $REF_FN