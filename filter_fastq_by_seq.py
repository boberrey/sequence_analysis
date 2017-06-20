#!/usr/bin/env python

"""
Remove or extract specific clusters from a directory of fastq files based on 
a file of sequences of interest.

This script can be used to demultiplex sequencing runs in a more flexible way
than bcltofastq.

Inputs:
   directory containing fastqs
   file containing sequences of interest

Outputs:
   filtered fastq files (new versions of the original fastqs with either the 
    isolated sequences, or the indicated sequences removed.)

Ben Ober-Reynolds
"""

import os
import sys
import argparse
from joblib import Parallel, delayed


def main():  
    # set up command line argument parser
    parser = argparse.ArgumentParser(description='script for isolating specific \
        clusters from fastq files')
    group = parser.add_argument_group('required arguments:')
    group.add_argument('-fd', '--fastq_directory', required=True,
        help='directory containing fastq files')
    group.add_argument('-sl', '--seq_list', required=True,
        help='file containing sequences to filter. Must be a tab-delimited \
        2-column file, where the first column is the read to inspect \
        (e.g. R1, R2, I1, I2) and the second column is the seqeunce to find.')

    group = parser.add_argument_group('optional arguments')
    group.add_argument('-a', '--action', type=str, default='keep',
        help='Action to take: [keep/remove], default is "keep"')
    group.add_argument('-od', '--output_directory',
        help='output directory for filtered fastq files (default is original \
            fastq_directory)')
    group.add_argument('-op', '--output_prefix',
        help='output prefix for filtered fastq files (default is "filtered_" + \
            filename)')
    group.add_argument('-n', '--num_cores', type=int, default=1,
        help='number of cores to use. Will parallelize based on fastq prefix.')

    # print help if no arguments provided
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    # parse command line arguments
    args = parser.parse_args()
    numCores = args.num_cores

    # Pre-defined variables, constants, and settings
    fastq_extension = 'fastq'

    # Check input directory
    fastq_dir = args.fastq_directory
    if not os.path.isdir(fastq_dir):
        print("Error: invalid fastq directory selection. Exiting...")
        sys.exit()

    # If no output directory given, use input directory
    output_dir = args.output_directory
    if not output_dir:
        output_dir = fastq_dir
    if not os.path.isdir(output_dir):
        print("Error: invalid output directory selection. Exiting...")
        sys.exit()
    output_dir = output_dir.strip('/') + '/'

    # Make sure a valid action option was provided:
    action = args.action
    if action != 'keep' and action != 'remove':
        print("Error: {} is not a valid action; must be either 'keep' or \
            'remove'.".format(action))
        sys.exit()

    # If no prefix is given, use 'filter'
    output_prefix = args.output_prefix
    if not output_prefix:
        output_prefix = 'filter'
    
    # Gather fastq files:
    print("Finding fastq files in directory {}".format(fastq_dir))
    fastq_list = find_files_in_directory(fastq_dir, 
        extensionList=[fastq_extension])

    # Batch fastqs by prefix:
    fastq_batches = batch_fastqs(fastq_list)

    # Read in sequence file:
    seq_dict = read_in_seq_list(args.seq_list)

    
    # loop thorugh fastq batches in parallel or in sequence
    results = []
    if numCores > 1:
        print("Filtering fastq files on {} cores...".format(numCores))
        results = (Parallel(n_jobs=numCores, verbose=10)\
            (delayed(filter_fastq_list)(fastq_list, batch, seq_dict, action, output_prefix, 
            output_dir) for prefix, fastq_list in fastq_batches.items()))
    else:
        results = [filter_fastq_list(fastq_list, batch, seq_dict, action, output_prefix, 
            output_dir) for batch, fastq_list in fastq_batches.items()]
    
    # Report results of filtering:
    # (This is a horrible way to report this info, but oh well.)
    print(results)
    action_line = ''
    if action == 'keep':
        action_line = 'kept'
    else:
        action_line = 'removed'
    for d in results:
        for batch, metrics in d.items():
            print("{} sequences in fastq batch {}".format(metrics['total_seqs'], batch))
            del metrics['total_seqs']
            for key, val in metrics.items():
                print("{}: {} reads {}".format(key, val, action_line))


def find_files_in_directory(dirPath, extensionList=None, 
                            excludedExtensionList=None):
    """
    Locate files in a given directory path. Optionally, desired files are 
    identified as matching one of the extension types provided in 
    'extensionList'
    Input: directory path, list of approved extensions, (list of excluded 
        extensions)
    Output: List of found files 
    """
    def extension_match(filename, extensionList=None):
        # from CPlibs
        if extensionList is not None:
            for currExt in extensionList:
                if filename.lower().endswith(currExt.lower()):
                    return True
        return False

    dirList = os.listdir(dirPath)
    fileList = []
    for currFilename in dirList:
        if (extension_match(currFilename, extensionList) 
        and not extension_match(currFilename, excludedExtensionList)): 
            fileList.append(dirPath+currFilename)
    if len(dirList) == 0:
        print('\tNONE FOUND')
    else:
        for filename in fileList:
            print("found:\t\t{}".format(filename))
        return fileList


def batch_fastqs(fastq_list):
    """
    Batch fastq files by common file prefix.
    Inputs:
        fastq_list (list) - previously identified list of fastq files
    Outputs:
        fastq_batches (dict) - dictionary of fastqs keyed by common prefix
    """
    fastq_batches = {}
    for fastq_file in fastq_list:
        # Defining prefixes this way allows different lanes to be split into
        # different batches for better parallelization
        prefix = "_".join(os.path.basename(fastq_file).split('_')[0:3])
        if prefix in fastq_batches:
            fastq_batches[prefix].append(fastq_file)
        else:
            fastq_batches[prefix] = [fastq_file]
    return fastq_batches


def read_in_seq_list(filename):
    """
    Read in a list of sequences to filter by. Must be a two-column, 
    tab-delimited file, with the first column containing the read to inspect
    (e.g. R1, R2, I1, I2), and the second column containing the sequence of
    interest.
    Inputs:
        filename (str) - filename for sequence file
    Outputs:
        seq_dict (dict) - dictionary of sequences to identify, keyed by read
    """
    seq_dict = {}
    with open(filename, 'r') as f:
        for line in f:
            read_key, seq = line.strip().split()
            if read_key in seq_dict:
                seq_dict[read_key].append(seq)
            else:
                seq_dict[read_key] = [seq]
    return seq_dict


def filter_fastq_list(fastq_list, batch, seq_dict, action, output_prefix, output_dir):
    """
    Filter a batch of fastq files based on sequences in the previously created
    sequence_dict. 
    Inputs:
        fastq_list (list) - a list containing a batch of fastqs 
            (a batch contains each of the reads for a single lane of sequencing)
        seq_dict (dict) - dictionary containing sequences to filter by, keyed 
            by read_id (e.g. R1, R2, I1, I2)
        action (str) - the action to take upon identification of a sequence match 
            (keep or remove)
        output_prefix (str) - the output filename prefix
        output_dir (str) - the output directory path
    Outputs:
        write new files
        results_dict (dict) - some metrics on filtering
    """
    # Example fastq format:
    # @M00653:218:000000000-AYC5G:1:1101:20964:1096 1:N:0:1
    # CNTATAATGATTCTTATTGACCAAAAAGCTGACAATTCACTTATTTTGCTTGACTATTTATTATACTTTCA
    # +
    # C#8BCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGG

    # Construct dict of fastq files, assigning each file to its read ID
    read_keys = ['R1', 'R2', 'I1', 'I2']
    fastq_dict = {}
    for fastq_file in fastq_list:
        file_basename = os.path.basename(fastq_file)
        for key in read_keys:
            if key in file_basename:
                fastq_dict[key] = fastq_file

    # construct output file dict:
    output_file_dict = {}
    for key, filename in fastq_dict.items():
        output_file_dict[key] = output_dir + output_prefix + "_" + os.path.basename(filename)

    # get the total number of lines:
    total_lines = 0
    with open(fastq_dict[read_keys[0]], 'r') as f:
        total_lines = sum(1 for line in f)

    # Every four lines is a new cluster:
    num_clusters = int(total_lines/4)

    # Initialize a 'results_dict' with a single key containin the common prefix
    # for this batch of fastqs
    results_dict = {batch: {'total_seqs': num_clusters}}

    # Loop through all files in chunks of four lines:
    # Note: this currently is a crappy way to do this. Ideally, we would have 
    # a setup where we don't need to provide exactly these sequence files each
    # time. 
    with open(fastq_dict['R1'], 'r') as R1_in, open(output_file_dict['R1'], 'w') as R1_out, \
    open(fastq_dict['R2'], 'r') as R2_in, open(output_file_dict['R2'], 'w') as R2_out, \
    open(fastq_dict['I1'], 'r') as I1_in, open(output_file_dict['I1'], 'w') as I1_out, \
    open(fastq_dict['I2'], 'r') as I2_in, open(output_file_dict['I2'], 'w') as I2_out:
        # Store the file handles in dictionaries according to read key:
        in_buffers = {'R1': R1_in, 'R2': R2_in, 'I1': I1_in, 'I2': I2_in}
        out_buffers = {'R1': R1_out, 'R2': R2_out, 'I1': I1_out, 'I2': I2_out}
        # initialize a 'working cluter' dict that will store the relevant info 
        # for the cluster currently being processed
        working_cluster_dict = {}
        for key in read_keys:
            working_cluster_dict[key] = {'cluster':"", 'seq':"",'spacer':"",'quality':""}

        # now loop through each cluster:
        for chunk in range(num_clusters):
            # has any match been found?
            any_match = False

            # Read in all data:
            for key in read_keys:
                working_cluster_dict[key]['cluster'] = in_buffers[key].readline()
                working_cluster_dict[key]['seq'] = in_buffers[key].readline()
                working_cluster_dict[key]['spacer'] = in_buffers[key].readline()
                working_cluster_dict[key]['quality'] = in_buffers[key].readline()

                # Check for a sequence match:
                if key in seq_dict:
                    for match_seq in seq_dict[key]:
                        if match_seq in working_cluster_dict[key]['seq']:
                            any_match = True
                            if match_seq in results_dict[batch]:
                                results_dict[batch][match_seq] += 1
                            else:
                                results_dict[batch][match_seq] = 1

            # Write data according to action (i.e. if 'keep', then write for matches,
            # if 'remove', write for non-matches)
            if (action == 'keep' and any_match) or (action == 'remove' and not any_match):
                for key in read_keys:
                    out_buffers[key].write(working_cluster_dict[key]['cluster'])
                    out_buffers[key].write(working_cluster_dict[key]['seq'])
                    out_buffers[key].write(working_cluster_dict[key]['spacer'])
                    out_buffers[key].write(working_cluster_dict[key]['quality'])

    return results_dict



if __name__ == '__main__':
    main()
