#!/usr/bin/env python

"""
Pull specific clusters from fastq file

Inputs:
   file containing list of clusters to isolate
   directory of fastq files

Outputs:
   filtered fastq files

Ben Ober-Reynolds
"""

import os
import sys
import argparse
from joblib import Parallel, delayed

def main():
    ################ Parse input parameters ################
    
    #set up command line argument parser
    parser = argparse.ArgumentParser(description='script for isolating specific \
        clusters from fastq files')
    group = parser.add_argument_group('required arguments:')
    group.add_argument('-cl', '--cluster_list', required=True,
                        help='file containing list of clusters to select')
    group.add_argument('-fd', '--fastq_directory', required=True,
                        help='directory containing fastq files')
    group = parser.add_argument_group('optional arguments')
    group.add_argument('-od','--output_directory',
                        help='output directory for filtered fastq files \
                        (default is original fastq_directory)')
    group.add_argument('-op','--output_prefix',
                        help='output prefix for filtered fastq files \
                        (default is cluster list filename)')
    group.add_argument('-n','--num_cores', type=int, default=1,
                        help='number of cores to use (should be same as number \
                            of fastq files)')


    # print help if no arguments provided
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    #parse command line arguments
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

    # If no prefix is given, use cluster list filename
    output_prefix = args.output_prefix
    if not output_prefix:
        output_prefix = os.path.splitext(os.path.basename(args.cluster_list))[0]
    

    # Gather fastq files:
    print("Finding fastq files in directory {}".format(fastq_dir))
    fastq_list = find_files_in_directory(fastq_dir, extensionList=[fastq_extension])

    # Make set of clusters to compare against:
    cluster_set = get_clusters_to_keep(args.cluster_list)
    print("Will keep reads based on cluster IDs in {}".format(args.cluster_list))
    
    # loop thorugh fastq files in parallel or in sequence
    results = []
    if numCores > 1:
        print("Filtering fastq files on {} cores...".format(numCores))
        results = (Parallel(n_jobs=numCores, verbose = 10)\
            (delayed(filter_fastq)(cluster_set, 
            fastq_file, output_prefix, 
            output_dir, fastq_extension)\
             for fastq_file in fastq_list))
    else:
        results = [filter_fastq(cluster_set, 
            fastq_file, output_prefix, 
            output_dir, fastq_extension)\
             for fastq_file in fastq_list]
    
    # Report results of filtering:
    for result in results:
        print("file {} has {} clusters, filtered down from {}".format(
            result[0], result[1], result[2]))



def find_files_in_directory(dirPath, extensionList=None, 
                            excludedExtensionList=None):
    """
    Locate files in a given directory path. Optionally, desired files are 
    identified as matching one of the extension types provided in 
    'extensionList'
    Input: directory path, list of approved extensions, (list of excluded extensions)
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


def get_clusters_to_keep(filename):
    """
    Generate a set of clusters extracted from a provided filename
    Input: filename
    Output: set of clusters
    """
    cluster_set = set()
    with open(filename, 'r') as f:
        for line in f:
            cluster_set.add(line.strip())
    return cluster_set


def filter_fastq(filter_set, fastq_filename, output_prefix, 
    output_dir, fastq_extension):
    """
    filter a fastq file by clusters that exist in the cluster set, then
    save the filtered file as a new file
    (Note: I tried this function using biopython tools first, but it was an 
        order of magnitude slower than writing my own fastq parser. From what
        I could find online, this is likely because the biopython 
        implementations of the parsing and writing functions include 
        significantly more error checking than required for standard
        4-line fastq's)
    Input: filter_set, fastq_filename, fastq_identifier, output_prefix
    Output: saved filtered file
    """
    # Example fastq format:
    # @M00653:218:000000000-AYC5G:1:1101:20964:1096 1:N:0:1
    # CNTATAATGATTCTTATTGACCAAAAAGCTGACAATTCACTTATTTTGCTTGACTATTTATTATACTTTCATCATA
    # +
    # C#8BCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGG
    
    fastq_basename = os.path.splitext(os.path.basename(fastq_filename))[0]
    new_filename = output_dir + output_prefix + '_' + fastq_basename + '.' + \
        fastq_extension

    # get the total number of lines:
    total_lines = 0
    with open(fastq_filename, 'r') as f:
        total_lines = sum(1 for line in f)

    # Every four lines is a new cluster:
    num_clusters = int(total_lines/4)
    cluster_count = 0

    # Loop through file in chunks of four lines:
    with open(fastq_filename, 'r') as infile, open(new_filename, 'w') as outfile:
        for chunk in range(num_clusters):
            cluster = infile.readline()
            seq = infile.readline()
            spacer = infile.readline()
            qual_score = infile.readline()

            # The first character of the cluster_ID in the fastq ('@') is not
            # present in the cluster_ID of the filter set.
            if cluster[1:].split()[0] in filter_set:
                cluster_count += 1
                outfile.write(cluster)
                outfile.write(seq)
                outfile.write(spacer)
                outfile.write(qual_score)
    
    # return the new filename, the starting number of clusters,
    # and the number of clusters kept
    return [new_filename, cluster_count, num_clusters]



if __name__ == '__main__':
    main()