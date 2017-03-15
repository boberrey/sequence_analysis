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
    group.add_argument('-dd', '--densities_dict', required=True,
                        help='pickled dictionary containing previously calculated motif densities.')
    group = parser.add_argument_group('optional arguments')
    group.add_argument('-od','--output_directory', default=".",
                        help='output directory for statistics file and figures. Default is current directory')
    group.add_argument('-op','--output_prefix', default="enrichment",
                        help='output prefix for results file and figures')


    # print help if no arguments provided
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    # Pre-defined variables, constants, and settings
    color_scheme = ['r', 'b', 'lightcoral', 'lightskyblue']
    plot_ylabels = 'normalized enrichment'
    plot_xlabels = ''
    plot_orientation = 'h'
    output_file_prefix = time.strftime("%Y%m%d") + "_" + args.output_prefix
    resolution = 200
    file_format = 'png'
    sns.set_style('white')

    # Do some error checking before running this long script:
    output_dir = args.output_directory
    if not os.path.isdir(output_dir):
        print("Error: invalid output directory. Exiting...")
        sys.exit()

    # Read in results file:
    density_result_dict = pd.read_pickle(args.densities_dict)

    # Generate plots of motif density distributions:
    for pname in density_result_dict.keys():
        plot_distributions(pname, density_result_dict, color_scheme, 
            plot_xlabels, plot_ylabels, plot_orientation, output_dir, 
            output_file_prefix, resolution, file_format)




def plot_distributions(pattern_name, results_dict, color_scheme, xlab, ylab, 
    plot_orientation, output_directory, output_prefix, resolution, file_format):
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
    ax = sns.boxplot(x='index', y='value', data=melted_frame, 
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