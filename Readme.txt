NanoSplicer

Usage example:

BAM_FN='/data/cephfs/punim0614/yupei/bulk_nanopore_Ritchie/analysis/PCR2_pool1/BAM/Chr_ID/NC_000001.11.bam'
REF_FN='/data/cephfs/punim0614/shared/shared_data/external_public/RNASeqMixology/splice_site_analysis/GRCh38_latest_genomic.fna'
FAST_DIR='/data/cephfs/punim0614/yupei/bulk_nanopore_Ritchie/PCR2_pool1/fast5_pass/barcode01'
SCRIPT='/data/cephfs/punim0614/yupei/bulk_nanopore_Ritchie/analysis/script'

# runing 
python3 $SCRIPT/NanoSplicer2.py -i $BAM_FN -f $FAST_DIR -r $REF_FN 
# ploting
python3 $SCRIPT/NanoSplicer_plot.py -i $BAM_FN -f $FAST_DIR -r $REF_FN -G

# subset run
python3 $SCRIPT/NanoSplicer_group.py -i $BAM_FN -f $FAST_DIR -r $REF_FN -G

Notes:
  - NanoSplicer2 is the multiprocess version and is more actively mainteined.  

  