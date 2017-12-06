# Configuration file 

import os
import sys
import glob



# Input data

FASTQ_DIR = "/raid1/lab/ben/ATAC_playground/Fastqs"
METADATA_FILE = "None"
BEDS = {"TSS" : '/shr/Downloaded_data/hg19_data/RefSeq_genes/hg19_snakeATAC_parsed_merged_ANS.tss.bed'}
NARROWPEAKS = {}
BROADPEAKS = {}
PCA_COLOR = "cell_type"
PCA_SHAPE = "None"
CONTRAST_FILE = "None"

# Resources
REFERENCE_FILE = '/shr/genomes/bowtie2/hg19/hg19'
P2_ACTIVATE = '/lab/greenleaf/p2/bin/activate'
P3_ACTIVATE = '/lab/greenleaf/p3/bin/activate'
CHROM_SIZES = '/lab/greenleaf/snakeATAC/resources/hg19/hg19.chrom.sizes'
EFFECTIVE_GENOME_SIZE = 2.7e9 #hg19 = 2.7e9, mm9 = 1.87e9, sacCer3 = 1.2*10**7
#FASTQ_SCREEN_CONF = "/app/fastqscreen/fastq_screen_v0.4.4/fastq_screen.conf"
BLACKLIST = "/raid1/lab/ben/reference_files/blacklists/Hg19_blacklist.bed"

# this directory is added to the system path for the purposes of running snakeATAC
EXE_DIR = '/lab/greenleaf/bin'

PICARD_JAR = '/lab/greenleaf/bin/picard.jar'
SNAKE_DIR = '/mnt/raid1/lab/greenleaf/snakeATAC'
ATAC_TOOLS = SNAKE_DIR + '/atac_tools'

# Parameters
TRIMMING_THREADS = 8
ALIGNING_THREADS = 16
JAVA_GC_CORES = 16
CALLPEAKS_PVAL = 1e-4

# Functions

def make_meta(fastq_dir):
    r1_files = list(map(os.path.abspath,glob.glob(os.path.join(fastq_dir,"*_R1*.f*"))))
    if (len(r1_files) < 1):
        sys.exit("No fastqs with _R1 found.")
    r2_files = [os.path.join(os.path.dirname(r1_file), os.path.basename(r1_file).replace('R1', 'R2')) for r1_file in r1_files]
    if  all([os.path.isfile(r2_file) for r2_file in r2_files]) is False:
        sys.exit("Not all matching _R2 files found.")
    sample_labels = [os.path.basename(r1_file).split("_R1")[0] for r1_file in r1_files]
    meta_prefix = os.path.commonprefix(sample_labels)
    filename = fastq_dir + meta_prefix + "_metadata.txt"
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, 'w') as outfile:
        outfile.write("\t".join(["Name","Read1","Read2"]) + "\n")
        for sample_label, r1_file, r2_file in zip(sample_labels, r1_files, r2_files):
            if len(sample_label) > 30:
                sample_label = sample_label[:20] + "..." + sample_label[-10:]
            outfile.write("\t".join([sample_label, r1_file, r2_file]) + "\n")
    return filename