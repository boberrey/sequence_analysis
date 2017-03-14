#!/usr/bin/env python

"""
Calculate enrichment statistics for two sets of fasta files

Inputs:
   two fasta files to compare
   file containing patterns to check

Outputs:
   pattern enrichment statistics
   enrichment plots

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

# Force matplotlib to use a non-interactive backend
import matplotlib
matplotlib.use('Agg')
# This will prevent it from trying to show plots during code execution

import matplotlib.pyplot as plt
import seaborn as sns
from joblib import Parallel, delayed

def main():
    ################ Parse input parameters ################

    #set up command line argument parser
    parser = argparse.ArgumentParser(description='script for isolating specific clusters from fastq files')
    group = parser.add_argument_group('required arguments:')
    group.add_argument('-fi', '--fasta_of_interest', required=True,
                        help='file containing clusters of interest')
    group.add_argument('-fb', '--background_fasta', required=True,
                        help='file containing background clusters')
    group.add_argument('-pf', '--pattern_file', required=True,
                        help='file containing patterns to check for. Format: {pattern name}\\t{regex_pattern}')
    group = parser.add_argument_group('optional arguments')
    group.add_argument('-od','--output_directory', default=".",
                        help='output directory for statistics file and figures. Default is current directory')
    group.add_argument('-op','--output_prefix', default="enrichment",
                        help='output prefix for results file and figures')
    group.add_argument('-isn','--interesting_seq_name', default="Sequences of Interest",
                        help='The name of the sequence of interest pool. Default is "Sequences of Interest"')
    group.add_argument('-bsn','--background_seq_name', default="Background Sequences",
                        help='The name of the background sequence pool. Default is "Background Sequences"')
    group.add_argument('-rc','--reverse_comp', default="y",
                        help='also calculate enrichment in reverse complement of each pool [y/n]? Default is y.')
    group.add_argument('-nb','--num_bootstraps', type=int, default=1000,
                        help='number of times to resample pools for enrichment calculation. Default is 1000.')
    group.add_argument('-n','--num_cores', type=int, default=1,
                        help='number of cores to use for bootstrapping.')


    # print help if no arguments provided
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    #parse command line arguments
    args = parser.parse_args()
    numCores = args.num_cores

    # Pre-defined variables, constants, and settings
    rev_c_tag = "Rev-Comp"
    color_scheme = ['r', 'b', 'lightcoral', 'lightskyblue']
    plot_col_order = [args.interesting_seq_name, args.background_seq_name]
    plot_ylabels = 'normalized enrichment'
    plot_xlabels = ''
    # orientation is either h or v:
    plot_orientation = 'v'
    output_file_prefix = time.strftime("%Y%m%d") + "_" + args.output_prefix
    resolution = 200
    file_format = 'png'
    sns.set_style('white')

    # Do some error checking before running this long script:
    output_dir = args.output_directory
    if not os.path.isdir(output_dir):
        print("Error: invalid output directory. Exiting...")
        sys.exit()

    # Now run the actual script:

    # Read in files:
    seqs_of_interest = read_fasta(args.fasta_of_interest)
    background_seqs = read_fasta(args.background_fasta)
    pattern_dict = read_pattern_file(args.pattern_file)

    # Find smallest pool size:
    pool_size = min([len(seqs_of_interest), len(background_seqs)])

    # seq pool dict:
    seq_pool_dict = {args.interesting_seq_name: seqs_of_interest, 
    args.background_seq_name: background_seqs}

    # Results dictionary:
    density_result_dict = {}
    for pname in pattern_dict.keys():
        density_result_dict[pname] = {}

    # compare to reverse complement?
    if args.reverse_comp == 'y':
        interesting_seq_rc_name = args.interesting_seq_name + " " + rev_c_tag
        background_seq_rc_name = args.background_seq_name + " " + rev_c_tag
        rc_seqs_of_interest = reverse_comp(seqs_of_interest)
        rc_background_seqs = reverse_comp(background_seqs)
        seq_pool_dict[interesting_seq_rc_name] = rc_seqs_of_interest
        seq_pool_dict[background_seq_rc_name] = rc_background_seqs
        plot_col_order = plot_col_order + [interesting_seq_rc_name, background_seq_rc_name]

    # calculate motif density for each pattern
    if numCores > 1:
        with Parallel(n_jobs = numCores, verbose = 10) as parallel: 
            for pname in pattern_dict.keys():
                for pool_name in seq_pool_dict.keys():
                    densities = []
                    print("Calculating density of pattern '{}' in pool '{}'".format(pname, pool_name))
                    densities = parallel(delayed(calc_resampled_motif_density)\
                        (seq_pool_dict[pool_name], pool_size, pattern_dict[pname])\
                    for i in range(args.num_bootstraps))
                    density_result_dict[pname][pool_name] = densities
    else:
        for pname in pattern_dict.keys():
            for pool_name in seq_pool_dict.keys():
                densities = []
                print("Calculating density of pattern '{}' in pool '{}'".format(pname, pool_name))
                densities = [calc_resampled_motif_density(
                    seq_pool_dict[pool_name], pool_size, pattern_dict[pname])\
                    for i in range(args.num_bootstraps)]
                density_result_dict[pname][pool_name] = densities

    # Dump results to pickle for latter replotting
    with open(output_dir + '/' + output_file_prefix + '.p', 'wb') as f:
        pickle.dump(density_result_dict, f)

    # Generate plots of motif density distributions:
    for pname in density_result_dict.keys():
        plot_distributions(pname, density_result_dict, color_scheme, 
            plot_xlabels, plot_ylabels, plot_orientation, plot_col_order, 
            output_dir, output_file_prefix, resolution, file_format)


def read_fasta(filename):
    """
    Read in a fasta file, and return sequences as a list.
    Input: fasta filename
    Output: sequence array 
    """
    fasta_list = []
    file_len = 0
    with open(filename, 'r') as f:
        file_len = sum(1 for line in f)
    with open(filename, 'r') as f:
        for dublet in range(int(file_len/2)):
            header = f.readline().strip()
            seq = f.readline().strip().upper()
            fasta_list.append(seq)
    return np.array(fasta_list)


def read_pattern_file(filename):
    """
    Read in a pattern file. Note that pattern files must be two-column,
    tab-delimited files with the first column being the pattern name, and
    the second column the regular expression defining that pattern.
    """
    pattern_dict = {}
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


def plot_distributions(pattern_name, results_dict, color_scheme, xlab, ylab, 
    plot_orientation, plot_col_order, output_directory, output_prefix, resolution, 
    file_format):
    """
    Save a plot of the motif density distribution.
    Inputs: motif pattern, bootstrapped results, color_scheme, x label, 
    y label,plot orientation, output directory, output prefix, resolution, file_format
    Output: save a plot of the motif density distributions
    """
    # how many colors are required:
    num_pools = len(results_dict[pattern_name].keys())
    pal = color_scheme[:num_pools]

    # Format data for plotting:
    pattern_df = pd.DataFrame(results_dict[pattern_name])
    melted_frame = pd.melt(pattern_df.T.reset_index(), id_vars='index')

    # Generate plot:
    if plot_orientation == 'h':
        ax = sns.boxplot(x='value', y='index', data=melted_frame, order=plot_col_order,
            orient=plot_orientation, palette=pal)
        ax.set_title(pattern_name)
        ax.set_ylabel(xlab)
        ax.set_xlabel(ylab)
    else:
        ax = sns.boxplot(x='index', y='value', data=melted_frame, order=plot_col_order,
            orient=plot_orientation, palette=pal)
        ax.set_title(pattern_name)
        ax.set_ylabel(ylab)
        ax.set_xlabel(xlab)
    output_file = output_directory + '/' + output_prefix + '_' + \
    '_'.join(pattern_name.split()) + '.' + file_format
    print("Saving {} plot to {}...".format(pattern_name, output_file))
    plt.savefig(output_file, dpi=resolution, format=file_format)

    # close plot
    plt.close()






if __name__ == '__main__':
    main()