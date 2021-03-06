#!/usr/bin/env bash
# get_stranded_fastas.sh

### bash script for getting insert sequence fasta files from name-sorted bam file.
### Use like so:

# ./get_insert_bed_and_fasta.sh bam_file.bam output_dir output_prefix 

# Ben Ober-Reynolds

# Check that the proper number of parameters were given
if [ "$#" -ne 3 ]; then
	echo "	bash script for getting insert bed and insert fasta files from a name-sorted bam file."
	echo "	Use like so:"
	echo "	get_insert_bed.sh bam_file.bam output_dir output_prefix"
	exit
fi

# input parameters:
bam_file=$1
output_dir=$2
output_prefix=$3

# working filenames:
error_log="$output_dir$output_prefix.err"
bedpe_file="$output_dir$output_prefix.bedpe"
cleaned_bedpe="$output_dir$output_prefix-cleaned.bedpe"
insert_bed_file="$output_dir$output_prefix-tmp_ins.bed"
trimmed_bed_file="$output_dir$output_prefix-trimmed.bed"
sorted_bed_file="$output_dir$output_prefix-insert.bed"

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

# insert bed columns:
i_chr="\$1"
i_strt="\$2"
i_stop="\$3"
i_name="\$4"
i_score="\$5"
i_strand="\$6"

# run parameters:
max_len=2000
quality_cutoff=20


### Script ###

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
rm $cleaned_bedpe

# Trim insert bed file to remove overly long or low quality inserts
remove_long_cmd="BEGIN{OFS = \"\t\"}{
	if($i_stop - $i_strt <= $max_len && $i_score >= $quality_cutoff){
		print $i_chr, $i_strt, $i_stop, $i_name, $i_score, $i_strand
	}
}"

awk "$remove_long_cmd" < $insert_bed_file > $trimmed_bed_file
rm $insert_bed_file

# Sorting the bedfile will dramatically speed up later bedtools operations
# This command will sort a bed file first by chromosome, then by start position, then by stop position
sort -k 1,1 -k2,2n -k3,3n $trimmed_bed_file > $sorted_bed_file
rm $trimmed_bed_file