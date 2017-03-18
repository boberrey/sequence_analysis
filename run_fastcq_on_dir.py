#!/usr/bin/env python

"""
Run fastqc on a directory of fastqs and pipe 
the summaries to stdout.

Inputs:
   directory containing fastqs to check

Outputs:
   fastqc reports

Ben Ober-Reynolds
"""

import os
import sys
import argparse
import subprocess
from termcolor import colored
from filter_fastqs import find_files_in_directory


def main():
    
    # set up command line argument parser
    parser = argparse.ArgumentParser(description='Script for running fastqc \
        on a directory of fastqs and piping the summaries to stdout.')
    group = parser.add_argument_group('required arguments:')
    group.add_argument('-fd', '--fastq_directory', required=True,
        help='directory containing fastq files')
    group = parser.add_argument_group('optional arguments')
    group.add_argument('-od', '--output_directory',
        help='output directory for fastqc reports (default is original \
            fastq_directory)')
    group.add_argument('-k', '--kmer_size', type=str, default="5",
        help='k-mer length option for fastqc. Default is 5.')

    # print help if no arguments provided
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    # parse command line arguments
    args = parser.parse_args()

    # Pre-defined variables, constants, and settings
    fastq_extension = 'fastq'
    fastqc_folder_ext = '_fastqc/'
    fastqc_command = 'fastqc'
    summary_file_name = 'summary.txt'
    full_summary_filename = 'full_summary.txt'
    pass_flag = 'PASS'
    warn_flag = 'WARN'
    fail_flag = 'FAIL'
    pass_color = 'green'
    warn_color = 'yellow'
    fail_color = 'red'

    # Check input directory
    fastq_dir = args.fastq_directory
    if not os.path.isdir(fastq_dir):
        print("Error: invalid fastq directory selection. Exiting...")
        sys.exit()

    # If no output directory given, use input directory
    output_dir = args.output_directory.strip('/')
    if not output_dir:
        output_dir = fastq_dir
    if not os.path.isdir(output_dir):
        print("Error: invalid output directory selection. Exiting...")
        sys.exit()

    # update full summary filepath
    full_summary_filename = output_dir + '/' + full_summary_filename

    # Gather fastq files:
    print("Finding fastq files in directory {}".format(fastq_dir))
    fastq_list = find_files_in_directory(fastq_dir, 
        extensionList=[fastq_extension])

    # Run fastqc on files:
    print("Running {} on found files...".format(fastqc_command))
    subprocess.call([fastqc_command] + fastq_list + \
        ['-o', output_dir, '-k', args.kmer_size, '-q'])

    # Pipe summaries to stdout:
    for fastq_file in fastq_list:
        # Where the summary file got saved to
        summary_file = output_dir + '/' + \
            os.path.splitext(os.path.basename(fastq_file))[0] + \
            fastqc_folder_ext + summary_file_name
        # Now write to stdout and to file:
        print("Summary of fasqc report for {}".format(fastq_file))
        with open(summary_file, 'r') as infile, \
        open(full_summary_filename, 'a') as outfile:
            outfile.write("summary of {}: \n".format(fastq_file))
            for line in infile:
                if warn_flag in line:
                    color = warn_color
                elif fail_flag in line:
                    color = fail_color
                else:
                    color = pass_color
                print(colored(line.strip(), color))
                outfile.write(line)
            print('')


if __name__ == '__main__':
    main()
