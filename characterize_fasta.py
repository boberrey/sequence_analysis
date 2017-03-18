#!/usr/bin/env python

"""
Generate plots of length distribution and base composition for one or more
fasta files.

Inputs:
   fasta file(s)

Outputs:
   plots of fasta characteristics

Ben Ober-Reynolds
"""


import os
import sys
import argparse
import time
from Bio import SeqIO
import pandas as pd
# Force matplotlib to use a non-interactive backend
import matplotlib
matplotlib.use('Agg')
# This will prevent it from trying to show plots during code execution
import matplotlib.pyplot as plt
import seaborn as sns


def main():
    # set up command line argument parser
    parser = argparse.ArgumentParser(description='Characterize fasta files by \
        base composition and sequence length distribution')
    group = parser.add_argument_group('required arguments:')
    group.add_argument('-fa', '--fasta_files', required=True, nargs='+',
        help='file containing list of clusters to select')
    group = parser.add_argument_group('optional arguments')
    group.add_argument('-fd', '--fasta_descriptors', nargs='+', 
        help='The descriptors to use for each fasta file. (If you provide any \
            descriptors, you must provide as many descriptors as files. If you \
            provide none, it will use the original fasta filename.)')
    group.add_argument('-od', '--output_directory', default=".",
        help='output directory for filtered fastq files (default is current \
            directory)')
    group.add_argument('-op', '--output_prefix', default='characterize',
        help='output prefix for fasta plots (default is "characterize")')

    # print help if no arguments provided
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    # parse command line arguments
    args = parser.parse_args()

    # Pre-defined variables, constants, and settings
    input_file_format = 'fasta'
    color_scheme = ['r', 'b', 'lightcoral', 'lightskyblue']

    base_comp_ymax = 0.4
    base_comp_title = 'Base Composition'
    base_comp_ylabel = ''
    base_comp_xlabel = ''

    length_col_header = 'length'
    fasta_col_header = 'fasta'
    seq_len_title = 'Insert Length Distribution'
    seq_len_ylabel = ''
    seq_len_xlabel = 'insert length'
    max_plot_col = 2
    seq_xmin, seq_xmax = [-100, 1200]

    output_prefix_base_comp = time.strftime("%Y%m%d") + "_" + \
        args.output_prefix + "_base_comp"
    output_prefix_seq_len = time.strftime("%Y%m%d") + "_" + \
        args.output_prefix + "_seq_lengths"

    resolution = 200
    file_format = 'png'
    sns.set_style('white')
    
    # set output directory:
    output_dir = args.output_directory.strip('/')
    if not os.path.isdir(output_dir):
        print("Error: invalid output directory selection. Exiting...")
        sys.exit()

    # initialize necessary data structures
    base_comp_dict = {}
    seq_len_dict = {}

    # Get the fasta names
    fasta_names = []
    for i in range(len(args.fasta_files)):
        fasta_file = args.fasta_files[i]
        if args.fasta_descriptors:
            fasta_name = args.fasta_descriptors[i]
        else:
            fasta_name = os.path.splitext(os.path.basename(fasta_file))[0]
        fasta_names.append(fasta_name)
    
    # gather information from fastas
    extract_data_from_fastas(args.fasta_files, fasta_names, input_file_format, 
        base_comp_dict, seq_len_dict)
    
    # Plot results:
    
    plot_base_compositions(base_comp_dict, fasta_names, color_scheme, 
        base_comp_title, base_comp_xlabel, base_comp_ylabel, base_comp_ymax, 
        output_dir, output_prefix_base_comp, resolution, file_format)
    
    if len(fasta_names) < max_plot_col:
        seq_len_col_num = len(fasta_names)
    else:
        seq_len_col_num = max_plot_col

    plot_seq_distributions(seq_len_dict, fasta_names, length_col_header, 
        fasta_col_header, seq_len_col_num, color_scheme, seq_len_title, 
        seq_len_xlabel, seq_len_ylabel, seq_xmin, seq_xmax, output_dir, 
        output_prefix_seq_len, resolution, file_format)
    

def extract_data_from_fastas(fasta_files, fasta_names, input_file_format, 
    base_comp_dict, seq_len_dict):
    """
    Extract the necessary information from a list of fasta files.
    Inputs: list of fasta files, list of fasta names, input file format,
    base composition dict, seq length dict
    Outputs: None (dictionaries modified in place)
    """
    # loop through each fasta file:
    for i in range(len(fasta_files)):
        fasta_file = fasta_files[i]
        fasta_name = fasta_names[i]
        base_comp_dict[fasta_name] = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
        seq_len_dict[fasta_name] = []

        # And then loop through each sequence in each file
        with open(fasta_file, 'r') as f:
            for seq_rec in SeqIO.parse(f, input_file_format):
                seq_rec = seq_rec.upper()
                base_comp_dict[fasta_name]['A'] += seq_rec.seq.count('A')
                base_comp_dict[fasta_name]['G'] += seq_rec.seq.count('G')
                base_comp_dict[fasta_name]['C'] += seq_rec.seq.count('C')
                base_comp_dict[fasta_name]['T'] += seq_rec.seq.count('T')
                seq_len_dict[fasta_name].append(len(seq_rec.seq))
            
        # Now get the proportion of each base:
        total_bases = sum(base_comp_dict[fasta_name].values())
        base_comp_dict[fasta_name] = {key: value / total_bases \
            for key, value in base_comp_dict[fasta_name].items()}


def plot_base_compositions(base_comp_dict, plot_order, color_scheme, title, 
    xlab, ylab, ymax, output_directory, output_prefix, resolution, file_format):
    """
    Plot the base compositions for each input fasta
    Input: base composition dict, plot order, color_scheme, plot title, 
    x label, y label, y max, output filepath, resolution, file format
    Output: none (save a plot)
    """
    # how many colors are required?
    num_fastas = len(base_comp_dict.keys())
    pal = color_scheme[:num_fastas]

    # Format data for plotting:
    base_comp_df = pd.DataFrame(base_comp_dict).reset_index()
    base_comp_melt = pd.melt(base_comp_df, id_vars='index')

    # Set order for plotting:
    base_comp_melt.loc[:,'variable'] = pd.Categorical(
        base_comp_melt.loc[:,'variable'], categories=plot_order)
    
    # Now plot:
    ax = sns.barplot(data=base_comp_melt, x='index', y='value', hue='variable', 
        palette=pal)
    ax.set_ylim(0, ymax)
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.legend(title='')

    # Save figure and close
    output_filepath = output_directory + '/' + output_prefix + '.' + file_format
    print("Saving base composition plot to {}...".format(output_filepath))
    plt.savefig(output_filepath, dpi=resolution, format=file_format)
    plt.close()


def plot_seq_distributions(seq_len_dict, plot_order, length_col_header, 
    fasta_col_header, col_num, color_scheme, title, xlab, ylab, xmin, xmax, 
    output_directory, output_prefix, resolution, file_format):
    """
    Plot the seq length distributions for each input fasta
    Input: sequence length dict, plot order, length column header, 
    fasta column header, column number, color_scheme, plot title, 
    x label, y label, x min, x max, output filepath, resolution, file format
    Output: none (save a plot)
    """
    # how many colors are required?
    num_fastas = len(seq_len_dict.keys())
    pal = color_scheme[:num_fastas]

    # Format data for plotting:
    seq_len_dfs = []
    for key, value, in seq_len_dict.items():
        seq_df = pd.DataFrame({length_col_header: value})
        seq_df[fasta_col_header] = key
        seq_len_dfs.append(seq_df)
    seq_length_df = pd.concat(seq_len_dfs)

    # Set order for plotting:
    seq_length_df.loc[:,fasta_col_header] = pd.Categorical(
        seq_length_df.loc[:, fasta_col_header], categories=plot_order)

    # Now plot:
    g = sns.FacetGrid(seq_length_df, col=fasta_col_header, hue=fasta_col_header, 
        sharex=True, xlim=[xmin, xmax], palette=pal, col_wrap=col_num)
    g.map(sns.kdeplot, length_col_header, shade=True)

    # update axes and labels
    for i in range(num_fastas):
        g.axes[i].set_title(plot_order[i])
    g.set_xlabels(xlab)
    g.set_ylabels(ylab)

    # Save plot:
    output_filepath = output_directory + '/' + output_prefix + '.' + file_format
    print("Saving seq length distribution plot to {}...".format(output_filepath))
    plt.savefig(output_filepath, dpi=resolution, format=file_format)
    plt.close()


if __name__ == '__main__':
    main()
