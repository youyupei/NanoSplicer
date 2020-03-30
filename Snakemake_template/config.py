#Suffix
SUFFIX = '_f15_t4_spikeT3_local_guppy'

#Data path
INDEX_FILE = "/home/ubuntu/data/cDNA_Selected_subset/guppy/fast5s.index"
FASTQ_PATH = "/home/ubuntu/data/cDNA_Selected_subset/guppy/pass"
FAST5_DIR = "/home/ubuntu/data/cDNA_Selected_subset/guppy/workspace"
SEQUENCING_SUMMARY = "/home/ubuntu/data/cDNA_Selected_subset/guppy/sequencing_summary.txt"
#Sequins path
TRANSCRIPT_ANNOTATION = "/home/ubuntu/data/Sequin_resources/rnasequin_annotation_2.2.gtf"
REF_TRANSCRIPT = "/home/ubuntu/data/Sequin_resources/rnasequin_sequences_2.2.fa" 
REF_GENOME = "/home/ubuntu/data/Sequin_resources/rnasequin_decoychr_2.2.fa" 

#python_script_path
PY_SCRIPT_PATH = "/home/ubuntu/PhD_proj/guppy_based_pipeline/script/"
PY_Candidates_from_gtf = PY_SCRIPT_PATH + "candidate_from_bed.py"
#arg
SEARCH_WIN = 20 # window size of candidate splicing site searching
ACCEPT_THRES = 3 # minimum support  # best supported site:
FLANK_SIZE = 15 # junction flanking in each side

PY_validation = PY_SCRIPT_PATH + "validation.py"
#arg
TRIM_SIGNAL = 0
TRIM_MODEL = 4
DTW_ADJ = False

# minimap2 settings
MAP_CANONICAL = "-ub" 

# command line alias:
PAFTOOLS = "k8 /usr/script/paftools.js "

#temp setting
TRANSCRIPT_ID_LIST = "/home/ubuntu/PhD_proj/guppy_based_pipeline/script/tested_transID.txt"

# max number of threads
THREAD = 32
MIN_MAPQ = 60

#############################load transcript ID#################################
trID = []
with open(TRANSCRIPT_ID_LIST, 'r') as f:
    for line in f:
        trID.append(line.strip())
