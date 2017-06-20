#!/usr/bin/env bash
# get_stranded_fastas.sh

### bash script for getting insert sequence fasta files from name-sorted bam file.
### Use like so:

# ./get_insert_bed_and_fasta.sh bam_file.bam ref_genome output_dir output_prefix 

# Ben Ober-Reynolds

# Check that the proper number of parameters were given
if [ "$#" -ne 4 ]; then
	echo "	bash script for getting insert bed and insert fasta files from a name-sorted bam file."
	echo "	Use like so:"
	echo "	get_insert_bed_and_fasta.sh bam_file.bam ref_genome output_dir output_prefix"
	exit
fi

# input parameters:
bam_file=$1
ref_genome=$2
output_dir=$3
output_prefix=$4

# working filenames:
error_log="$output_dir$output_prefix.err"
bedpe_file="$output_dir$output_prefix.bedpe"
cleaned_bedpe="$output_dir$output_prefix-cleaned.bedpe"
insert_bed_file="$output_dir$output_prefix-insert.bed"
genome_fasta="$ref_genome.fa"
fasta_output="$output_dir$output_prefix.fa"

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

# generate mate1-maintained paired-end bedfile:

bedtools bamtobed -mate1 -bedpe -i $bam_file > $bedpe_file 2> $error_log


# "clean-up" bedpe file by removing any unmapped pairs
clean_bedpe_cmd="{if($r1_strt > 0 && $r2_strt > 0) print}"
awk "$clean_bedpe_cmd" < $bedpe_file > $cleaned_bedpe
rm $bedpe_file

# convert the bedpe file into an "insert" bed file
make_insert_cmd="BEGIN{OFS = \"\t\"}{
	if($r1_strt > $r2_strt){
		print $r1_chr, $r2_strt, $r1_stop, $name, $score, $r1_strand
	} else if($r1_strt <= $r2_strt){
		print $r1_chr, $r1_strt, $r2_stop, $name, $score, $r1_strand
	}
}"
awk "$make_insert_cmd" < $cleaned_bedpe > $insert_bed_file
rm $cleaned_bedpe $error_log

# Now get 'strand forced' fasta from insert bed file:

bedtools getfasta -fi $genome_fasta -bed $insert_bed_file -fo $fasta_output -name -s
