#!/usr/bin/env python
""" 
Count occurances of specific sequences in a set of FASTQ files

 Inputs:
   FASTQ files
   seq table

 Outputs:
 	counts file
   

 Ben Ober-Reynolds, boberrey@stanford.edu
 20180316
 """

import sys
import os
import argparse
import string
import pandas as pd
import numpy as np
import time
import glob
import gzip
import copy
from difflib import SequenceMatcher
from joblib import Parallel, delayed

import time


### Global Vars ###
r1_id = "_R1_"
r2_id = "_R2_"
i1_id = "_I1_"
i2_id = "_I2_"


read1_primers = [
"GGATCCAGGAACGTCTTCCATACAACCTCCTTACTACAT", # Stall sequence
"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  # TruSeqR1 sequence
"GCATGCTAGCGACGTTCCTGGATCC" # pSTAC.3 F_primer sequence
]

# 20171103: We now have multiple R2 primers. I've modified this script to allow for
# multiple R2 primers when running in R1 mode, but maybe the smart thing is just to run
# in R2 mode from now on?
#read2_primer = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG"	# TruSeqR2
#read2_primer = "GATCGTGAGCCGGACCCAGCGTTGAGAAGAGGCAAAG"  # TtEndoR2

read2_primers = [
"AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG",  # TruSeqR2
"GATCGTGAGCCGGACCCAGCGTTGAGAAGAGGCAAAG",  # TtEndoR2
"CGGACGCGGGAAGACAGAATAGAAGTGTACGCCGTCG",  # let7-pShortR2
"CACAGGGAGCGAGGTATACCGAGAGCTGACAGCGTAG"   # let7-pLongR2
"ATAGTCGACACCGCTCTTCCGATCT" # pSTAC.3 R_primer sequence
]

primer_overlap = 15
#max_seq_length = 70
max_seq_length = 50

### MAIN ###

def main():
	start = time.time()
	################ Parse input parameters ################

	#set up command line argument parser
	parser = argparse.ArgumentParser(description='Script for generating a \
		counts file of observed sub-sequences (or regex paterns)')
	group = parser.add_argument_group('required arguments')
	group.add_argument('-fs', '--fastq_substring', required=True,
		help='common substring for fastqs you want to count from')
	group.add_argument('-s', '--seq_file', required=True,
		help='A tab-delimited table where the first column contains the \
		sub-sequences (or regex patterns) to search for.')
	group = parser.add_argument_group('optional arguments for processing data')
	group.add_argument('-od','--output_directory',
		help='output directory for series files with labeled \
		variants (default will use current directory)')
	group.add_argument('-e','--exact', action='store_true',
                        help='Flag to indicate if sub-sequences can be exact matched \
                        i.e. the sequenced insert will be identical to substring \
                        (much faster when there are lots of subsequences to search with).')
	group.add_argument('-p','--paired', action='store_true',
                        help='Flag to indicate if sub-sequences must be read in both reads')
	group.add_argument('-r', '--read_num', type=int, default=1,
		help='which read to use for matching if not using "exact" or "paired" (1 | 2)')
	group.add_argument('-n','--num_cores', type=int, default=1,
                        help='number of cores to use (parallelizes based on fastq group)')

	if not len(sys.argv) > 1:
	    parser.print_help()
	    sys.exit()

	#parse command line arguments
	args = parser.parse_args()
	numCores = args.num_cores

	# If no output directory given, use current directory
	if not args.output_directory:
		args.output_directory = "./"
	output_directory = args.output_directory
	if not os.path.isdir(output_directory):
		print("Error: invalid output directory selection. Exiting...")
		sys.exit()

	# Construct sequence dict:
	print("Reading in sequence file: {}".format(args.seq_file))
	seq_dict = get_seq_dict(args.seq_file, args.read_num)

	# Find fastqs matching prefix:
	print("Finding fastqs with sub-string: {}".format(args.fastq_substring))
	fastq_bundle_list = gather_fastq_bundles(args.fastq_substring)

	# Adjust numCores if necessary
	if numCores > len(fastq_bundle_list):
		numCores = len(fastq_bundle_list)

	# Count subsequences:
	if numCores > 1:
		print("Counting sub-sequences on {} cores...".format(numCores))
		count_dict_list = (Parallel(n_jobs=numCores, verbose=10)\
			(delayed(count_subsequences)(
				seq_dict, fastq_group, args.read_num, exact=args.exact, paired=args.paired) for fastq_group in fastq_bundle_list))
	else:
		print("Counting sub-sequences on a single core")
		count_dict_list = [count_subsequences(
			seq_dict, fastq_group, args.read_num, exact=args.exact, paired=args.paired) for fastq_group in fastq_bundle_list]

	# Format and save results:
	output_df = pd.read_table(args.seq_file, header=None)
	output_df.columns = ['sequence'] + output_df.columns.tolist()[1:]
	for counts_list in count_dict_list:
		colname = counts_list[0]
		new_df = pd.DataFrame.from_dict(counts_list[1], orient='index')
		new_df.columns = [colname]
		output_df = output_df.merge(new_df, how='outer', left_on='sequence', right_index=True)
	output_df.to_csv(args.fastq_substring + "_counts.txt", sep='\t', index=None, header=None)


def get_seq_dict(filename, read_num):
	"""
	Read in a seq table and extract the necessary information for 
	constructing the seq dict:
	Inputs:
		filename (str) - the filename for the variant dict
	Outputs:
		seq_dict (dict) - the seq dict, keyed by sequence, 
		with 0's as values (num seen)
	"""
	seq_dict = {}
	with open(filename, 'r') as f:
		for line in f:
			split_line = line.split('\t')
			seq = split_line[0]
			if read_num == 2:
				seq = rev_comp(seq)

			seq_dict[seq] = 0
	return seq_dict


def gather_fastq_bundles(fastq_substring):
	"""
	Gather all fastqs matching a common basename, then create fastq bundles 
	for each unique group of fastqs.
	Inputs:
		fastq_substring (str) - substring contained by fastqs of interest
	Outputs:
		fastq_bundle_list (list) - list of fastq bundles
	"""
	all_fastq_files = glob.glob("*{}*".format(fastq_substring))
	bundle_headers = set()
	for fastq in all_fastq_files:
		for splitter in [r1_id, r2_id, i1_id, i2_id]:
			if len(fastq.split(splitter)) > 1:
				bundle_headers.add(fastq.split(splitter)[0])
	fastq_bundle_list = []
	for bn in bundle_headers:
		fastq_bundle_list.append(FastqGroup(bn))
	return fastq_bundle_list


def count_subsequences(seq_dict, fastq_group, read_num, exact=True, paired=False):
	"""
	Count occurances of each sub-sequence in the seq dict
	Inputs:
		seq_dict (dict) - the dict keyed by sequence with counts as variables
		fastq_group (FastqGroup) - The FastqGroup object to search through
		read_num (int) - Which read to search through first
		paired (bool) - whether or not subsequences must be identified in both reads
	Outputs:
		counts_list (list) - A list where the first element is the FastqGroup
			basename and the second element is a FastqGroup-specific counts dict
	"""
	new_counts_dict = copy.deepcopy(seq_dict)
	linenum = 0
	# num_reads = fastq_group.get_num_reads() # very slow to calculate
	if exact:
		# Fastest search option when you have many subsequences. 
		# Find insert from paired reads (and a known primer sequence)
		# and perform exact matching in seq dict
		with gzip.open(fastq_group.r1_file, 'r') as r1_f, gzip.open(fastq_group.r2_file, 'r') as r2_f:
			for r1_line, r2_line in zip(r1_f, r2_f):
				linenum += 1
				if linenum % 4 == 2:
					r1_seq = r1_line.decode("utf-8").strip()
					r2_seq = r2_line.decode("utf-8").strip()
					insert = get_insert_seq(r1_seq, r2_seq, read2_primers)
					if insert in new_counts_dict:
						new_counts_dict[insert] += 1
	elif paired:
		with gzip.open(fastq_group.r1_file, 'r') as r1_f, gzip.open(fastq_group.r2_file, 'r') as r2_f:
			for r1_line, r2_line in zip(r1_f, r2_f):
				linenum += 1
				if linenum % 4 == 2:
					r1_seq = r1_line.decode("utf-8").strip()
					r2_seq = r2_line.decode("utf-8").strip()
					check_seq_for_matches(new_counts_dict, r1_seq, r2_seq)
	else:
		if read_num == 1:
			with gzip.open(fastq_group.r1_file, 'r') as r1_f:
				for r1_line in r1_f:
					linenum += 1
					if linenum % 4 == 2:
						r1_seq = r1_line.decode("utf-8").strip()
						check_seq_for_matches(new_counts_dict, r1_seq=r1_seq, r2_seq=None)
		else:
			with gzip.open(fastq_group.r2_file, 'r') as r2_f:
				for r2_line in r2_f:
					linenum += 1
					if linenum % 4 == 2:
						r2_seq = r1_line.decode("utf-8").strip()
						check_seq_for_matches(new_counts_dict, r1_seq=None, r2_seq=r2_seq)
	return [fastq_group.group_name, new_counts_dict]


def check_seq_for_matches(seq_dict, r1_seq=None, r2_seq=None):
	"""
	Check r1_seq and/or r2_seq for any matches in the seq_dict
	Most naive python substring searching option (but makes no
	assumptions about where substring will occur)
	"""
	if r1_seq and r2_seq:
		for seq, count in seq_dict.items():
			if seq in r1_seq:
				if rev_comp(seq) in r2_seq:
					seq_dict[seq] += 1
		return
	if r1_seq:
		for seq, count in seq_dict.items():
			if seq in r1_seq:
				seq_dict[seq] += 1
		return
	if r2_seq:
		for seq, count in seq_dict.items():
			if rev_comp(seq) in r2_seq:
				seq_dict[seq] += 1
		return
	return


def get_merged_insert(r1, r2, min_overlap=15):
	"""
	Get the merged insert between two paired reads.
	Uses difflib SequenceMatcher to find the longest substring

	Warning: fairly slow
	"""
	seq_match = SequenceMatcher(None, r1, rev_comp(r2))
	match = seq_match.find_longest_match(0, len(r1), 0, len(r2))
	if match.size < min_overlap:
		return None
	# If overlap is fully contained within both reads, return insert
	if match.a + match.size < len(r1):
		return r1[match.a: match.a + match.size]
	else:
		# if overlap is not fully contained, return the joined sequence
		return r1 + rev_comp(r2)[match.size:]


def get_insert_seq(readA, readB, primers):
	"""
	Find the insert sequence of two paired reads. If no overlap is found, will
	return false. Requires information about distal primer in readA
	Inputs:
		readA (str) - the read of focus
		readB (str) - the other read
	Outputs:
		insert_seq (str) - the insert sequence or the whole read if no insert
			found
	"""
	primer_seqs = [primer[:primer_overlap] for primer in primers]
	insert_end = -1
	for p in primer_seqs:
		insert_end = readA.find(p)
		if insert_end >= 0:
			break
	# If the primer sequence isn't found, just return the whole read 1
	if insert_end < 0:
		return readA[:max_seq_length]
	insert = readA[:insert_end]
	revB = rev_comp(readB)
	# It seems like there is some difficulty in matching the full length thing...
	if revB.find(insert[1:-1]) >= 0:
		return insert
	else:
		return readA[:max_seq_length]


def rev_comp(seq, complement_dict={'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}, missing='N'):
	# Reverse complement a sequence
	return "".join(complement_dict.get(base, missing) for base in reversed(seq))


def get_first_match(lst, substr):
# Return the first element of the list that contains the substring
	for i in lst:
		if substr in i:
			return i
	return None


class FastqGroup:
	"""
	A class for tracking information about a group of fastqs 
	(i.e. fastqs generated from one index set).
	"""
	
	def __init__(self, fastq_basename):
		self.group_name = fastq_basename
		self.all_matching_files = glob.glob("*{}*".format(fastq_basename))
		self.r1_file = get_first_match(self.all_matching_files, r1_id)
		self.r2_file = get_first_match(self.all_matching_files, r2_id)
		self.i1_file = get_first_match(self.all_matching_files, i1_id)
		self.i2_file = get_first_match(self.all_matching_files, i2_id)


	def __repr__(self):
		# representation of object instance
		return("{}_nr:{}_rl:{}".format(self.group_name, self.reads, self.read_length))

	def get_num_reads(self):
		# return the number of reads in the fastq group
		# (Requires readint through entire file. May be slow)
		reads = 0
		file_to_open = self.all_matching_files[0]
		with gzip.open(file_to_open, 'r') as f:
			for line in f:
				reads += 1
		return int(reads/4)

	def get_read_length(self, indicator=r1_id):
		# return the read length
		file_to_open = get_first_match(self.all_matching_files, indicator)
		with gzip.open(file_to_open, 'r') as f:
			junk_header = f.readline()
			read_length = len(f.readline())
		return read_length





if __name__ == '__main__':
    main()
