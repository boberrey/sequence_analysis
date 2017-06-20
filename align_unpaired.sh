#!/usr/bin/env bash
# get_stranded_fastas.sh

### bash script for aligning reads. Output is a coordinate-sorted bam file.
### Use like so:

# ./align_unpaired.sh input.fastq ref_genome output_dir output_prefix

# Ben Ober-Reynolds

# Check that the proper number of parameters were given
if [ "$#" -ne 4 ]; then
	echo "	bash script for getting insert sequence fasta files from paired end reads."
	echo "	Output is a coordinate-sorted bam file."
	echo "	Usage:"
	echo "	align_unpaired.sh input.fastq ref_genome output_dir output_prefix"
	exit
fi

# input parameters
f1=$1
ref_genome=$2
output_dir=$3
output_prefix=$4

# bowtie2 settings:
n_cores="18"

# working filenames:
sam_file="$output_dir$output_prefix.sam"
unsorted_bam_file="$output_dir$output_prefix-unsorted.bam"
coordinate_sorted_bam="$output_dir$output_prefix-c-sorted.bam"
dups_removed_bam="$output_dir$output_prefix-no_dups.bam"
name_sorted_bam="$output_dir$output_prefix.bam"
metrics_file="$output_dir$output_prefix-metrics"


# align
echo ""
echo "Aligning with bowtie2..."
echo ""

# Redirect bowtie2 metrics (normally stderr) to stdout for better wrapper script handling:
bowtie2 -p $n_cores $ref_genome -U $f1 -S $sam_file 2>&1
samtools view -b -S -o $unsorted_bam_file $sam_file 2>&1

# sort and clean up
echo ""
echo "Sorting and removing duplicates..."
echo ""

# sort by coordinate (required for removing duplicates)
# Redirect picard metrics (normally stderr) to stdout for better wrapper script handling:
java -jar -XX:ParallelGCThreads=$n_cores \
-Djava.io.tmpdir=`pwd`/tmp /app/picard/picard-tools-1.130/picard.jar \
SortSam SO=coordinate I=$unsorted_bam_file O=$coordinate_sorted_bam \
VALIDATION_STRINGENCY=SILENT 2>&1

# Remove PCR duplicates:
java -jar -XX:ParallelGCThreads=$n_cores \
-Djava.io.tmpdir=`pwd`/tmp /app/picard/picard-tools-1.130/picard.jar \
MarkDuplicates INPUT=$coordinate_sorted_bam OUTPUT=$dups_removed_bam \
METRICS_FILE=$metrics_file REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT 2>&1

# clean up intermediate files:
rm $sam_file $unsorted_bam_file
