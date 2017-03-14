#!/usr/bin/env bash
# get_stranded_fastas.sh

### bash script for getting strand-maintained fasta files from paired end reads.
### Use like so:
# ./get_stranded_fastas.sh r1.fastq r2.fastq adapters.fa ref_genome output_prefix

# input files and variables
r1=$1
r2=$2
adapters=$3
ref_genome=$4
output_name=$5
error_log="$output_name-pipe.err"
# trim adapters

echo ""
echo "Trimming adapters with skewer..."
echo ""

skewer -y $adapters -o $output_name -m pe $1 $2

# trimming output
tr1="$output_name-trimmed-pair1.fastq"
tr2="$output_name-trimmed-pair2.fastq"
trim_logfile="$output_name-trimmed.log"

# align
sam_file="$output_name.sam"
unsorted_bam_file="$output_name-unsorted.bam"

# There's got to be a better way to skip lines...
echo ""
echo "Aligning with bowtie2..."
echo ""

bowtie2 -X5000 -p18 $ref_genome -1 $tr1 -2 $tr2 -S $sam_file
samtools view -b -S -o $unsorted_bam_file $sam_file



# sort by chrom position
sorted_bam="$output_name.bam"

echo ""
echo "Sorting by chrom position..."
echo ""

java -jar -Djava.io.tmpdir=`pwd`/tmp /app/picard/picard-tools-1.130/picard.jar \
SortSam SO=coordinate I=$unsorted_bam_file O=$sorted_bam VALIDATION_STRINGENCY=SILENT
samtools index $sorted_bam

# cleanup intermediate files:
rm $tr1 $tr2 $trim_logfile $sam_file $unsorted_bam_file

# convert to mate1-maintained bedpe file:

echo ""
echo "Creating strand-maintained fasta files..."
echo ""
bedpe_file="$output_name.bedpe"
bedtools bamtobed -mate1 -bedpe -i $sorted_bam > $bedpe_file 2> $error_log

rm $sorted_bam "$sorted_bam.bai" $error_log

# 'clean-up' bedpe file by removing any unmapped pairs
cleaned_bedpe="$output_name-cleaned.bedpe"
awk '{if($1 > 0 && $4 > 0)print}' < $bedpe_file > $cleaned_bedpe
rm $bedpe_file

# convert the mate1 bedpe file into an 'insert' bed file
stranded_bed_file="$output_name-insert.bed"
awk 'BEGIN{OFS = "\t"}{if($2 > $5){print $1,$5,$3,$7,$8,$9} else if($2 <= $5){print $1,$2,$6,$7,$8,$9}}' < $cleaned_bedpe > $stranded_bed_file
rm $cleaned_bedpe

# Now get fasta from insert bed file (-s option important for maintaining strandedness!):
fasta_outfile="$output_name.fa"
bedtools getfasta -fi "$ref_genome.fa" -bed $stranded_bed_file -fo $fasta_outfile -s
rm $stranded_bed_file

echo ""
echo "Done."
echo ""