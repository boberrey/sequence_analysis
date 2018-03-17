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
from joblib import Parallel, delayed

import time


### Global Vars ###
r1_id = "_R1_"
r2_id = "_R2_"
i1_id = "_I1_"
i2_id = "_I2_"

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
	group.add_argument('-r', '--read_num', type=int, required=True,
		help='which read to use for initial matching (1 | 2)')
	group = parser.add_argument_group('optional arguments for processing data')
	group.add_argument('-od','--output_directory',
		help='output directory for series files with labeled \
		variants (default will use current directory)')
	group.add_argument('-p','--paired', action='store_true',
                        help='Flag to indicate if sub-sequences must be read in both reads')
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

	# Construct variant dict:
	print("Reading in sequence file: {}".format(args.seq_file))
	seq_dict = get_seq_dict(args.seq_file, args.read_num)

	# Find fastqs matching prefix:
	print("Finding fastqs with subtring: {}".format(args.fastq_substring))
	fastq_bundle_list = gather_fastq_bundles(args.fastq_substring)

	print(fastq_bundle_list)

	"""
	if numCores > 1:
		print "Annotating clusters in parallel on {} cores...".format(numCores)
		annotated_cluster_lists = (Parallel(n_jobs=numCores, verbose=10)\
			(delayed(annotate_clusters)(
				args.seq_directory + CPseq, variant_dict, args.read_num) for CPseq in CPseqFiles))
	else:
		print "Annotating clusters on a single core"
		annotated_cluster_lists = [annotate_clusters(
			args.seq_directory + CPseq, variant_dict, args.read_num) for CPseq in CPseqFiles]

	# Combine cluster lists:
	print "Formatting and saving CPannot file..."
	all_annotations = []
	map(all_annotations.extend, annotated_cluster_lists)
	CPannot_df = pd.DataFrame(all_annotations)
	try:
		CPannot_df.columns = ['cluster_ID', 'variant_ID']
	except:
		print "No variants annotated!"
		sys.exit()
	
	# Save the CPannot file as a pickle
	CPannotFilename = "_".join(longestSubstring(CPseqFiles).split("_")[:-1])+".CPannot.pkl"
	print "Creating CPannot.pkl file: {}...".format(CPannotFilename)
	CPannot_df = CPannot_df.set_index("cluster_ID")
	
	CPannot_df.to_pickle(output_directory+CPannotFilename)
	print "Done. {} minutes".format(round((time.time() - start)/60, 2))
	"""





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


def annotate_clusters(CPseq_filename, variant_dict, read_num):
	"""
	Annotate cluster IDs with their appropriate variants
	Inputs:
		CPseq_filename (str) - the CPseq filename
		variant_dict (dict) - the variant dict
	Outputs:
		annotated_clusters (list) - list with annotated clusters
	"""
	annotated_clusters = []
	with open(CPseq_filename, 'r') as f:
		for line in f:
			split_line = line.split('\t')
			clusterID = split_line[clusterID_column]
			read1 = split_line[r1_column][trim_length:]
			read2 = split_line[r2_column][trim_length:]

			# Get the insert sequence from paired reads:
			if read_num == 1:
				insert_seq = get_insert_seq(read1, read2, read2_primer)
			else:
				insert_seq = get_insert_seq(read2, read1, rev_comp(read1_primer))

			# if insert seq not in variant dict, continue to next line
			if not insert_seq in variant_dict:
				continue

			# If still going, it means there is a match, so add that annotation
			annotated_clusters.append([clusterID, variant_dict[insert_seq]])
	return annotated_clusters


def get_insert_seq(readA, readB, primer):
	"""
	Find the insert sequence of two paired reads. If no overlap is found, will
	return false
	Inputs:
		readA (str) - the read of focus
		readB (str) - the other read
	Outputs:
		insert_seq (str) - the insert sequence or the whole read if no insert
			found
	"""
	primer_seq = primer[:primer_overlap]
	insert_end = readA.find(primer_seq)
	# If the primer sequence isn't found, just return the whole read 1
	if insert_end < 0:
		return readA[:max_seq_length]
	insert = readA[:insert_end]
	revB = rev_comp(readB)
	# It seems like there is some difficulty in matching the full length thing...
	if revB.find(insert[1:-1]) > 0:
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
		self.r2_file = get_first_match(self.all_matching_files, i2_id)
		self.reads = 0
		self.read_length = 0
		print(self.group_name)
		with gzip.open(self.r1_file, 'r') as f:
			for line in f:
				self.reads += 1
		self.reads = int(self.reads/4)
		with gzip.open(self.r1_file, 'r') as f:
			junk_header = f.readline()
			self.read_length = len(f.readline())

	def __repr__(self):
		# representation of object instance
		return("{}_nr{}_rl{}".format(self.group_name, self.reads, self.read_length))


if __name__ == '__main__':
    main()
