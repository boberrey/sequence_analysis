#!/usr/bin/env python

"""
Convert a directory of paired-end bam files into bed files.

Inputs:
   directory containing bam files
   reference genome

Outputs:
   bed files

Ben Ober-Reynolds
"""

import os
import sys
import argparse
import subprocess
import time
from collections import OrderedDict
from basic_utils import find_files_in_directory

def main():
    # start the timer:
    start_time = time.time()
    # set up command line argument parser
    parser = argparse.ArgumentParser(description='script for generating fasta \
        files from name-sorted paired-end bam files.')
    group = parser.add_argument_group('required arguments:')
    group.add_argument('-bd', '--bam_directory', required=True,
        help='directory containing fastq files')
    group = parser.add_argument_group('optional arguments')
    group.add_argument('-od', '--output_directory',
        help='output directory for filtered fastq files (default is original \
            fastq_directory)')

    # print help if no arguments provided
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    # parse command line arguments
    args = parser.parse_args()

    # Pre-defined variables, constants, and settings
    bam_extension = 'bam'
    get_bed_script = 'get_insert_bed.sh'
    log_file = "alignment.log"
    bytestring_encoding = "UTF-8"

    # Check input directory
    bam_dir = args.bam_directory
    if not os.path.isdir(bam_dir):
        print("Error: invalid bam directory selection. Exiting...")
        sys.exit()

    # If no output directory given, use input directory
    output_dir = args.output_directory
    if not output_dir:
        output_dir = bam_dir
    if not os.path.isdir(output_dir):
        print("Error: invalid output directory selection. Exiting...")
        sys.exit()
    # ensure correct formatting of output_dir
    output_dir = output_dir.strip('/') + '/'

    # update logfile path:
    log_file = output_dir + log_file 
    
    # Gather fastq files:
    print("Finding bam files in directory {}".format(bam_dir))
    bam_list = find_files_in_directory(bam_dir, 
        extensionList=[bam_extension])

    # Report found files:
    print("Found bam files: ")
    for bam_file in bam_list:
        print("\t{}".format(bam_file))
    
    # Run fasta script for each bam file:
    print("Running bed file script: {}".format(get_bed_script))

    for bam_file in bam_list:
        print("Getting insert bed file for {}...".format(bam_file))
        # fasta script:
        # get_insert_fastas.sh bam_file.bam ref_genome output_dir output_prefix
        # subprocess.check_output returns stdout as a string
        output_prefix = os.path.splitext(os.path.basename(bam_file))[0]
        log = subprocess.check_output([get_bed_script, bam_file, 
            output_dir, output_prefix])

        # Convert bytestring to writable string:
        log = log.decode(bytestring_encoding)
        # save stdout information as logfile:
        with open(log_file, 'a') as f:
            f.write("Log for {}:\n".format(output_prefix))
            f.write(log)

    print("Finished generating {} bed files.".format(len(bam_list)))
    print("{} minutes".format(round((time.time() - start_time)/60, 2)))


if __name__ == '__main__':
    main()
