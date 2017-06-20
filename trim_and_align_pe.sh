#!/usr/bin/env bash
# get_stranded_fastas.sh

### bash script for trimming and aligning reads. Output is a coordinate-sorted bam file.
### Use like so:

# ./trim_and_align.sh r1.fastq r2.fastq adapters.fa ref_genome output_dir output_prefix

# Ben Ober-Reynolds

# Check that the proper number of parameters were given
if [ "$#" -ne 6 ]; then
	echo "	bash script for getting insert sequence fasta files from paired end reads."
	echo "	Output is a coordinate-sorted bam file."
	echo "	Usage:"
	echo "	trim_and_align.sh r1.fastq r2.fastq adapters.fa ref_genome output_dir output_prefix"
	exit
fi

# input parameters
r1=$1
r2=$2
adapters=$3
ref_genome=$4
output_dir=$5
output_prefix=$6

# bowtie2 settings:
max_insert_size="5000"
n_cores="18"

# working filenames:
tr1="$output_prefix-trimmed-pair1.fastq"
tr2="$output_prefix-trimmed-pair2.fastq"
trim_logfile="$output_prefix-trimmed.log"
sam_file="$output_dir$output_prefix.sam"
unsorted_bam_file="$output_dir$output_prefix-unsorted.bam"
coordinate_sorted_bam="$output_dir$output_prefix-c-sorted.bam"
dups_removed_bam="$output_dir$output_prefix-no_dups.bam"
name_sorted_bam="$output_dir$output_prefix.bam"
metrics_file="$output_dir$output_prefix-metrics"

# trim adapters
echo ""
echo "Trimming adapters with skewer..."
echo ""

skewer --quiet -y $adapters -o $output_prefix -m pe $1 $2 

# align
echo ""
echo "Aligning with bowtie2..."
echo ""

# Redirect bowtie2 metrics (normally stderr) to stdout for better wrapper script handling:
bowtie2 -X $max_insert_size -p $n_cores $ref_genome -1 $tr1 -2 $tr2 -S $sam_file 2>&1
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

# sort by name (required to make bedpe)
java -jar -XX:ParallelGCThreads=$n_cores \
-Djava.io.tmpdir=`pwd`/tmp /app/picard/picard-tools-1.130/picard.jar \
SortSam SO=queryname I=$dups_removed_bam O=$name_sorted_bam \
VALIDATION_STRINGENCY=SILENT 2>&1

# clean up intermediate files:
rm $tr1 $tr2 $trim_logfile $sam_file $unsorted_bam_file $coordinate_sorted_bam $dups_removed_bam
