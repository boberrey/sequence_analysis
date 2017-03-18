#!/usr/bin/env bash
# get_stranded_fastas.sh

### bash script for getting strand-maintained fasta files from paired end reads.
### Use like so:

# ./get_stranded_fastas.sh r1.fastq r2.fastq adapters.fa ref_genome output_prefix

# Ben Ober-Reynolds

# input files and variables
r1=$1
r2=$2
adapters=$3
ref_genome=$4
output_name=$5
error_log="$output_name-pipe.err"

# variables and filenames:
tr1="$output_name-trimmed-pair1.fastq"
tr2="$output_name-trimmed-pair2.fastq"
trim_logfile="$output_name-trimmed.log"
sam_file="$output_name.sam"
unsorted_bam_file="$output_name-unsorted.bam"
sorted_bam="$output_name.bam"
bedpe_file="$output_name.bedpe"
bam_bai_file="$sorted_bam.bai"
cleaned_bedpe="$output_name-cleaned.bedpe"
stranded_bed_file="$output_name-insert.bed"
genome_fasta="$ref_genome.fa"
fasta_outfile="$output_name.fa"

# bedpe columns:
r1_chr="\$1"
r1_strt="\$2"
r1_stop="\$3"
r2_chr="\$4"
r2_strt="\$5"
r2_stop="\$6"
name="\$7"
score="\$8"
r1_strand="\$9"

# trim adapters
echo ""
echo "Trimming adapters with skewer..."
echo ""

skewer -y $adapters -o $output_name -m pe $1 $2

# align
echo ""
echo "Aligning with bowtie2..."
echo ""

bowtie2 -X5000 -p18 $ref_genome -1 $tr1 -2 $tr2 -S $sam_file
samtools view -b -S -o $unsorted_bam_file $sam_file

# sort by chrom position
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

bedtools bamtobed -mate1 -bedpe -i $sorted_bam > $bedpe_file 2> $error_log

# Currently, I don't want to keep bam files.
rm $sorted_bam $bam_bai_file

# 'clean-up' bedpe file by removing any unmapped pairs
clean_bedpe_cmd="{if($r1_strt > 0 && $r2_strt > 0) print}"
awk "$clean_bedpe_cmd" < $bedpe_file > $cleaned_bedpe
rm $bedpe_file

# convert the mate1 bedpe file into an 'insert' bed file
make_insert_cmd="BEGIN{OFS = \"\t\"}{
	if($r1_strt > $r2_strt){
		print $r1_chr, $r2_strt, $r1_stop, $name, $score, $r1_strand
	} else if($r1_strt <= $r2_strt){
		print $r1_chr, $r1_strt, $r2_stop, $name, $score, $r1_strand
	}
}"
awk "$make_insert_cmd" < $cleaned_bedpe > $stranded_bed_file

rm $cleaned_bedpe

# Now get fasta from insert bed file (-s option important for maintaining strandedness!):

bedtools getfasta -fi $genome_fasta -bed $stranded_bed_file -fo $fasta_outfile -s -name
rm $stranded_bed_file

echo ""
echo "Done."
echo ""