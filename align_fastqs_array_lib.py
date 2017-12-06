#!/usr/bin/env python

"""
Trim and align a directory of fastqs for an array library.
Will output an 'insert bed file', and an 'insert fasta file'

Inputs:
   directory containing fastqs
   adapter list (for adapter trimming)
   reference genome

Outputs:
    'insert bed file'
    'insert fasta file'

Other required scripts:
    trim_and_align_array.sh


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
    parser = argparse.ArgumentParser(description='script for trimming and \
        aligning fastqs for an array experiment, and subsequent generation \
        of useful downstream files.')
    group = parser.add_argument_group('required arguments:')
    group.add_argument('-fd', '--fastq_directory', required=True,
        help='directory containing fastq files')
    group.add_argument('-a1', '--adapter_file1', required=True,
        help='file containing adapters for trimming')
    group.add_argument('-a2', '--adapter_file2', required=True,
        help='file containing adapters for trimming')
    group.add_argument('-g', '--ref_genome', required=True,
        help='reference genome for bowtie2')
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
    fastq_extension = 'fastq'
    bam_extension = 'bam'
    r1_identifier = '_R1_'
    r2_identifier = '_R2_'
    alignment_script = 'trim_and_align_array.sh'
    get_insert_script = 'get_insert_bed_and_fasta.sh'
    alignment_log_file = "_alignment.log"
    insert_log_file = "_inserts.log"
    bytestring_encoding = "UTF-8"

    # Check input directory
    fastq_dir = args.fastq_directory
    if not os.path.isdir(fastq_dir):
        print("Error: invalid fastq directory selection. Exiting...")
        sys.exit()

    # Other required arguments:
    adapter_file1 = args.adapter_file1
    adapter_file2 = args.adapter_file2
    ref_genome = args.ref_genome

    # If no output directory given, use input directory
    output_dir = args.output_directory
    if not output_dir:
        output_dir = fastq_dir
    if not os.path.isdir(output_dir):
        print("Error: invalid output directory selection. Exiting...")
        sys.exit()
    # ensure correct formatting of output_dir
    output_dir = output_dir.strip('/') + '/'
    
    # Gather fastq files:
    print("Finding fastq files in directory {}".format(fastq_dir))
    fastq_list = find_files_in_directory(fastq_dir, 
        extensionList=[fastq_extension])
    
    # Split fastqs based of r1 or r2
    r1_fastqs = []
    r2_fastqs = []
    for fastq_file in fastq_list:
        if r1_identifier in fastq_file:
            r1_fastqs.append(fastq_file)
        if r2_identifier in fastq_file:
            r2_fastqs.append(fastq_file)

    # get paired fastq files:
    paired_fastq_dict = get_paired_fastq_dict(r1_fastqs, r2_fastqs, r1_identifier, r2_identifier)

    # Report found files:
    print("Found fastq file pairs: ")
    for prefix, file_list in paired_fastq_dict.items():
        print("\t{}\t{}".format(file_list[0], file_list[1]))
    
    # Run alignment script for each pair:
    print("Running alignment script: {}".format(alignment_script))
    for prefix, file_list in paired_fastq_dict.items():
        
        r1_fastq, r2_fastq = file_list
        print("Trimming and aligning {} and {}...".format(r1_fastq, r2_fastq))
        # alignment script:
        # trim_and_align.sh r1.fastq r2.fastq adapters.fa ref_genome output_dir output_prefix
        # subprocess.check_output returns stdout as a string
        log = subprocess.check_output([alignment_script, r1_fastq, r2_fastq, 
            adapter_file1, adapter_file2, ref_genome, output_dir, prefix])

        # Convert bytestring to writable string:
        log = log.decode(bytestring_encoding)
        # save stdout information as logfile:
        current_log_file = output_dir + prefix + alignment_log_file

        with open(current_log_file, 'a') as f:
            f.write("Log for {}:\n".format(prefix))
            f.write(log)
    
    print("Finished alignment of {} fastq pairs.".format(len(paired_fastq_dict.keys())))
    print("{} minutes".format(round((time.time() - start_time)/60, 2)))


    # Generate the 'insert' files for later use:
    print("Generating insert files with script: {}".format(get_insert_script))

    # Get all the bam files that were generated in the previous step:
    bam_list = []
    for prefix, file_list in paired_fastq_dict.items():
        bam_list.append(output_dir + prefix + '.' + bam_extension)


    for bam_file in bam_list:
        print("Getting insert bed and fasta for {}...".format(bam_file))
        # insert script:
        # get_insert_bed_and_fasta.sh bam_file.bam ref_genome output_dir output_prefix 
        # subprocess.check_output returns stdout as a string
        output_prefix = os.path.splitext(os.path.basename(bam_file))[0]
        log = subprocess.check_output([get_insert_script, bam_file, ref_genome, 
            output_dir, output_prefix])

        # Convert bytestring to writable string:
        log = log.decode(bytestring_encoding)
        # save stdout information as logfile:
        current_log_file = output_dir + output_prefix + insert_log_file

        with open(current_log_file, 'a') as f:
            f.write("Log for {}:\n".format(output_prefix))
            f.write(log)

    print("Finished generating {} sets of insert files.".format(len(bam_list)))
    print("{} minutes".format(round((time.time() - start_time)/60, 2)))



def get_paired_fastq_dict(r1_fastqs, r2_fastqs, r1_identifier, r2_identifier):
    """
    Generate a dictionary of paired fastq files, keyed by their common prefix.
    Inputs:
        r1_fastqs (list) - list of read 1 fastq files
        r2_fastqs (list) - list of read 2 fastq files
        r1_identifier (str) - identifier that marks a read 1 fastq file
    Outputs:
        paired_fastq_dict (dict) - dictionary of paired fastq files
    """
    paired_fastq_dict = OrderedDict()
    for r1_fastq in r1_fastqs:
        # Using basename ensures that path elements are not part of downstream naming
        prefix = os.path.basename(r1_fastq.split(r1_identifier)[0])
        r2_fastq = ""
        for r2_file in r2_fastqs:
            if prefix == os.path.basename(r2_file.split(r2_identifier)[0]):
                r2_fastq = r2_file
        paired_fastq_dict[prefix] = [r1_fastq, r2_fastq]
    return paired_fastq_dict



if __name__ == '__main__':
    main()
