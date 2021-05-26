
SUFFIX = '_sequins_barcode01'

#Data path
FASTQ_PATH = "/data/cephfs/punim0614/yupei/deep_sequins/fastq3.6.1/pass/barcode01"
FAST5_DIR = "/data/cephfs/punim0614/yupei/deep_sequins/fast5_pass/barcode01"

#Reference
#TRANSCRIPT_ANNOTATION = "/data/cephfs/punim0614/yupei/deep_sequins/Sequin_resources/rnasequin_annotation_2.2.gtf"
REF_GENOME = "/data/cephfs/punim0614/yupei/deep_sequins/Sequin_resources/v2.4/rnasequin_decoychr_2.4.fa"
REF_TRANSCRIPT = "/data/cephfs/punim0614/yupei/deep_sequins/Sequin_resources/v2.4/rnasequin_sequences_2.4.fa"

#python_script_path
PY_SCRIPT_PATH = "/data/cephfs/punim0614/yupei/bulk_nanopore_Ritchie/analysis/script"

#arg
SEARCH_WIN = 20 # window size of candidate splicing site searching
ACCEPT_THRES = 3 # minimum support  # best supported site:
FLANK_SIZE = 20 # junction flanking in each side

PY_validation = PY_SCRIPT_PATH + "validation.py"
#arg
TRIM_SIGNAL = 6
TRIM_MODEL = 2
DTW_ADJ = False
DTW_BAND = 0.4

# minimap2 settings
MAP_CANONICAL = "-ub" 

# command line alias:
PAFTOOLS = "k8 /usr/script/paftools.js "

#temp setting
TRANSCRIPT_ID_LIST = "/home/ubuntu/data/Sequin_resources/transcriptID_list.txt"

# max number of threads
THREAD = 32
MIN_MAPQ = 60

##############################load transcript ID#################################
#trID = []
#with open("Validation/tested_transID.txt", 'r') as f:
#    for line in f:
#        trID.append(line.strip())
##trID = ["R1_41_2"]

