#!/usr/bin/env python

"""
Calculate enrichment statistics for two sets of fasta files

Inputs:
   two fasta files to compare
   file containing patterns to check

Outputs:
   pickled dictionary of pattern enrichments

Ben Ober-Reynolds
"""


import os
import sys
import re
import time
import argparse
import numpy as np
import pandas as pd
import pickle
from collections import OrderedDict
from Bio import SeqIO
from joblib import Parallel, delayed


def main():

    # set up command line argument parser
    parser = argparse.ArgumentParser(description='Calculate motif densities \
        for a target and a background set of fastas.')
    group = parser.add_argument_group('required arguments:')
    group.add_argument('-fa', '--fasta_files', required=True, nargs='+',
        help='fasta files to be used for analysis')
    group.add_argument('-pf', '--pattern_file', required=True,
        help='file containing patterns to check for. Format: \
        {pattern name}\\t{regex_pattern}')
    group = parser.add_argument_group('optional arguments')
    group.add_argument('-od', '--output_directory', default=".",
        help='output directory for statistics file and figures. \
        Default is current directory')
    group.add_argument('-op', '--output_prefix', default="enrichment",
        help='output prefix for results file and figures')
    group.add_argument('-fd', '--fasta_descriptors', nargs='+', 
        help='The descriptors to use for each fasta file. (If you provide any \
            descriptors, you must provide as many descriptors as files. If you \
            provide none, it will use the original fasta filename.)')
    group.add_argument('-rc', '--reverse_comp', default="y",
        help='also calculate enrichment in reverse complement of each pool \
        [y/n]? Default is y.')
    group.add_argument('-nb', '--num_bootstraps', type=int, default=1000,
        help='number of times to resample pools for enrichment calculation. \
        Default is 1000.')
    group.add_argument('-n', '--num_cores', type=int, default=1,
        help='number of cores to use for bootstrapping.')

    # print help if no arguments provided
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    # parse command line arguments
    args = parser.parse_args()
    numCores = args.num_cores

    # Pre-defined variables, constants, and settings
    input_file_format = 'fasta'
    rev_c_tag = "Rev-Comp"
    output_prefix = time.strftime("%Y%m%d") + "_" + args.output_prefix
    max_resample_size = 100000

    # Do some error checking before running this long script:
    output_dir = args.output_directory
    if not os.path.isdir(output_dir):
        print("Error: invalid output directory. Exiting...")
        sys.exit()
    
    # Read in files:
    # Get the fasta names
    fasta_names = []
    for i in range(len(args.fasta_files)):
        fasta_file = args.fasta_files[i]
        if args.fasta_descriptors:
            fasta_name = args.fasta_descriptors[i]
        else:
            fasta_name = os.path.splitext(os.path.basename(fasta_file))[0]
        fasta_names.append(fasta_name)

    # seq pool dict:
    seq_pool_dict = OrderedDict()
    for i in range(len(fasta_names)):
        fasta_file = args.fasta_files[i]
        fasta_name = fasta_names[i]
        print("Reading in seq pool for {}".format(fasta_name))
        seq_pool_dict[fasta_name] = read_fasta(fasta_file, input_file_format)
    pattern_dict = read_pattern_file(args.pattern_file)

    # Find smallest pool size:
    pool_size = min([len(x) for x in seq_pool_dict.values()])
    pool_size = min([pool_size, max_resample_size])

    # Results dictionary:
    density_result_dict = OrderedDict()
    for pname in pattern_dict.keys():
        density_result_dict[pname] = OrderedDict()

    # compare to reverse complement?
    if args.reverse_comp == 'y':
        for pool_name, seq_pool in seq_pool_dict.items():
            rev_comp_name = pool_name + " " + rev_c_tag
            seq_pool_dict[rev_comp_name] = reverse_comp(seq_pool)

    # calculate motif density for each pattern
    if numCores > 1:
        with Parallel(n_jobs=numCores, verbose=10) as parallel: 
            for pname in pattern_dict.keys():
                for pool_name in seq_pool_dict.keys():
                    densities = []
                    print("Calculating density of pattern '{}' in pool '{}'\
                        ".format(pname, pool_name))
                    densities = parallel(delayed(calc_resampled_motif_density)\
                        (seq_pool_dict[pool_name], pool_size, pattern_dict[pname])
                    for i in range(args.num_bootstraps))
                    density_result_dict[pname][pool_name] = densities
    else:
        for pname in pattern_dict.keys():
            for pool_name in seq_pool_dict.keys():
                densities = []
                print("Calculating density of pattern '{}' in pool '{}'\
                    ".format(pname, pool_name))
                densities = [calc_resampled_motif_density(
                    seq_pool_dict[pool_name], pool_size, pattern_dict[pname])
                    for i in range(args.num_bootstraps)]
                density_result_dict[pname][pool_name] = densities

    # Dump results to pickle for latter replotting
    with open(output_dir + '/' + output_prefix + '.p', 'wb') as f:
        pickle.dump(density_result_dict, f)


def read_fasta(filename, input_file_format):
    """
    Read in a fasta file, and return sequences as a list.
    Input: fasta filename
    Output: sequence array 
    """
    fasta_list = []
    with open(filename, 'r') as f:
        for seq_rec in SeqIO.parse(f, input_file_format):
            seq_rec = seq_rec.upper()
            fasta_list.append(str(seq_rec.seq))
    return np.array(fasta_list)


def read_pattern_file(filename):
    """
    Read in a pattern file. Note that pattern files must be two-column,
    tab-delimited files with the first column being the pattern name, and
    the second column the regular expression defining that pattern.
    """
    pattern_dict = OrderedDict()
    with open(filename, 'r') as f:
        for line in f:
            pname, reg_exp = line.strip().split('\t')
            reg_exp = re.compile(reg_exp)
            pattern_dict[pname] = reg_exp
    return pattern_dict


def reverse_comp(fasta_array):
    """
    Reverse complement a list of sequences
    Input: list of sequences
    Output: reverse complement of same sequence list
    """
    trans_table = str.maketrans('AGCT', 'TCGA')
    rev_list = []
    for seq in fasta_array:
        rev_list.append(seq.translate(trans_table)[::-1])
    return np.array(rev_list)


def calc_resampled_motif_density(seq_array, samp_size, regex):
    """
    Calculate the length-normalized density of a specific regular
    expression pattern in a resampled sequence pool.
    Inputs: list of sequences, number of seqs to draw, regular expression pattern
    Output: length-normalized motif density
    """
    resampled_pool = np.random.choice(seq_array, size=samp_size, replace=True)
    total_seq_space = 0
    patterns_found = 0
    for seq in resampled_pool:
        patterns_found += len(re.findall(regex, seq))
        total_seq_space += len(seq)
    return patterns_found/total_seq_space


if __name__ == '__main__':
    main()
