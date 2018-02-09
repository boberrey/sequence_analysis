#!/usr/bin/env bash
# align_fasta.sh

### bash script for aligning a fasta file to a reference genome.
### will output a .bam and a .bed file of alignment
### Use like so:

# ./align_fasta.sh seqs.fa ref_genome output_dir output_prefix

# Ben Ober-Reynolds

# Check that the proper number of parameters were given
if [ "$#" -ne 4 ]; then
	echo "	bash script for aligning a fasta file."
	echo "	Output is bam and bed file of original fasta sequences"
	echo "	Usage:"
	echo "	align_fasta.sh seqs.fa ref_genome output_dir output_prefix"
	exit
fi

# input parameters
seqs=$1
ref_genome=$2
output_dir=$3
output_prefix=$4

# bowtie2 settings:
n_cores="8"

# working filenames:
sam_file="$output_dir$output_prefix.sam"
unsorted_bam_file="$output_dir$output_prefix.bam"
unsorted_bed_file="$output_dir$output_prefix.bed"

# align
echo ""
echo "Aligning with bowtie2..."
echo ""

# Redirect bowtie2 metrics (normally stderr) to stdout for better wrapper script handling:
bowtie2 -f -p $n_cores -x $ref_genome -U $seqs -S $sam_file 2>&1

# Make bam file:
samtools view -b -S -o $unsorted_bam_file $sam_file 2>&1

# Make bed file:
bedtools bamtobed -i $unsorted_bam_file > $unsorted_bed_file

# clean up intermediate files:
rm $sam_file
