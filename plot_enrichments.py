#!/usr/bin/env python

"""
Generate plot of motif enrichments from previously generated motif density 
dictionary

Inputs:
   motif density dictionary
   pattern list file (one column file containing patterns to plot)

Outputs:
   enrichment plot

Ben Ober-Reynolds
"""
import os
import sys
import time
import argparse
import pandas as pd
# Force matplotlib to use a non-interactive backend
import matplotlib
matplotlib.use('Agg')
# This will prevent it from trying to show plots during code execution
import matplotlib.pyplot as plt
import seaborn as sns


def main():

    # set up command line argument parser
    parser = argparse.ArgumentParser(description='Plot the previously \
        calculated motif enrichments from a pickled dictionary.')
    group = parser.add_argument_group('required arguments:')
    group.add_argument('-dd', '--densities_dict', required=True,
        help='pickled dictionary containing previously calculated motif \
        densities.')
    group.add_argument('-pl', '--pattern_list', required=True,
        help='file containing patterns from the densities dict that should be \
        plotted. File is a single column of patterns, and order of patterns \
        will determine plot order.')
    group = parser.add_argument_group('optional arguments')
    group.add_argument('-po', '--pool_order',
        help='single column file containing sequence pool order.')
    group.add_argument('-od', '--output_directory', default=".",
        help='output directory for statistics file and figures. Default is \
        current directory')
    group.add_argument('-op', '--output_prefix', default="enrichment",
        help='output prefix for results file and figures')

    # print help if no arguments provided
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()
        
    # parse command line arguments
    args = parser.parse_args()

    # Pre-defined variables, constants, and settings
    color_scheme = ['r', 'b', 'lightcoral', 'lightskyblue']
    col_num_lim = 3
    pattern_col_header = 'pattern'
    plot_ylabel = 'normalized enrichment'
    plot_xlabel = ''
    output_prefix = time.strftime("%Y%m%d") + "_" + args.output_prefix
    resolution = 200
    file_format = 'png'
    sns.set_style('white')

    # Check to make sure output directory is valid:
    output_dir = args.output_directory.strip('/')
    if not os.path.isdir(output_dir):
        print("Error: invalid output directory. Exiting...")
        sys.exit()

    # Read in pattern file:
    pattern_list = read_in_list(args.pattern_list)

    # Set the number of columns for the facet grid
    if len(pattern_list) < col_num_lim:
        col_num = len(pattern_list)
    else:
        col_num = col_num_lim

    # Read in pool file if provided:
    if args.pool_order:
        pool_order = read_in_list(args.pool_order)
    else:
        pool_order = []
    
    # Read in results file:
    density_results_dict = pd.read_pickle(args.densities_dict)

    # generate plot:
    plot_facet_boxplots(pattern_list, density_results_dict, pool_order, 
        pattern_col_header, color_scheme, col_num, plot_xlabel, plot_ylabel, 
        output_dir, output_prefix, resolution, file_format)


def read_in_list(list_file):
    """
    Read in a ordered list from file. List file should be a 
    single column of elements.
    Input: List file
    Output: list
    """
    return_list = []
    with open(list_file, 'r') as f:
        for line in f:
            return_list.append(line.strip())
    return return_list


def plot_facet_boxplots(pattern_list, results_dict, pool_order, pattern_col_header, 
    color_scheme, col_num, xlab, ylab, output_directory, output_prefix, 
    resolution, file_format):
    """
    Generate boxplots in facet grid format from pattern enrichment data.
    Inputs: motif patterns, bootstrapped results, pool order, pattern column header, 
    color_scheme, number of columns in facet plot, x label, y label, 
    output directory, output prefix, resolution, file_format
    Output: saved plot
    """
    # how many colors are required:
    pool_number = max([len(x) for x in [y.keys() 
        for y in results_dict.values()]])
    pal = color_scheme[:pool_number]

    # Format data for plotting
    pattern_dfs = []
    for pattern in pattern_list:
        pat_df = pd.DataFrame(results_dict[pattern])
        pat_df[pattern_col_header] = pattern
        pattern_dfs.append(pat_df)
    full_df = pd.concat(pattern_dfs)
    data_melt = pd.melt(full_df, id_vars=pattern_col_header)

    # Sort by pool order (for some reason passing order to the boxplot mapping
    # doesn't work, so this is a less ideal workaround):
    if pool_order:
        data_melt.loc[:, 'variable'] = pd.Categorical(
            data_melt.loc[:, 'variable'], categories=pool_order)

    # Generate plot:
    g = sns.FacetGrid(data_melt, col=pattern_col_header, sharey=False, 
        sharex=False, col_wrap=col_num, col_order=pattern_list)
    if pool_order:
        g.map(sns.boxplot, pattern_col_header, 'value', 'variable', 
            palette=pal).add_legend(label_order=pool_order)
    else:
        g.map(sns.boxplot, pattern_col_header, 'value', 'variable', 
            palette=pal).add_legend()

    # update axes and labels
    for i in range(len(pattern_list)):
        g.axes[i].set_title(pattern_list[i])
        g.axes[i].set_xticklabels([xlab])
    g.set_xlabels(xlab)
    g.set_ylabels(ylab)

    # save plot
    output_file = output_directory + '/' + output_prefix + '.' + file_format
    print("Saving plot to {}...".format(output_file))
    plt.savefig(output_file, dpi=resolution, format=file_format)
    plt.close()


if __name__ == '__main__':
    main()
