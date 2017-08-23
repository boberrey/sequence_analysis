#!/usr/bin/env python
""" Split a folder of sequenced library fastqs into subdirectories


 Inputs:
   Directory with fastqs 

 Outputs:
   Fastqs split into new directories by library

 Ben Ober-Reynolds, boberrey@stanford.edu 
20170823
 """

import sys
import os
import argparse
from collections import OrderedDict

def main():
    ################ Parse input parameters ################

    #set up command line argument parser
    parser = argparse.ArgumentParser(description='Script for splitting a directory of fastqs into multiple directories')
    group = parser.add_argument_group('required arguments')
    group.add_argument('-fd', '--fastq_dir', required=True,
                        help='directory that holds fastqs to be split (.fastq or .fastq.gz)')

    group = parser.add_argument_group('optional arguments for processing data')
    group.add_argument('-s','--suffix', default="lib",
                        help='suffix for new directories. default = lib')
    group.add_argument('-od','--output_directory', default='./',
                        help='directory in which new directories will be made')
    group.add_argument('-a','--action', default='l',
                        help='what to do with the images (m = move, l = symbolic link). Default is to link.')

    if not len(sys.argv) > 1:
        parser.print_help()
        sys.exit()

    # Pre-defined variables, constants, and settings
    r1_identifier = '_R1_'
    r2_identifier = '_R2_'
    fastq_extensions = ['fastq', 'fastq.gz']

    #parse command line arguments
    args = parser.parse_args()
    if args.action != "m" and args.action != "l":
        print("Error: action must be either 'm' (move) or 'l' (link)!")
        sys.exit()

    fastq_dir = os.path.abspath(args.fastq_dir.strip("/")) + "/"

    # Gather fastq files:
    print("Finding fastq files in directory {}".format(fastq_dir))
    fastq_list = find_files_in_directory(fastq_dir, 
        extensionList=fastq_extensions)
    
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

    # Do the splitting:
    outputPath = args.output_directory
    if not os.path.exists(outputPath):
        print("Error: directory {} does not exist!".format(outputPath))

    for lib_prefix, fastq_files in paired_fastq_dict.items():
        dirname = outputPath + lib_prefix + "_" + args.suffix
        os.mkdir(dirname)
        for fastq in fastq_files:
            if args.action == "m":
                os.rename(fastq, dirname + "/" + os.path.basename(fastq))
            if args.action == "l":
                os.symlink(fastq, dirname + "/" + os.path.basename(fastq))
    print("Files split successfully")


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
            fileList.append(dirPath + currFilename)
    if len(dirList) == 0:
        print('\tNONE FOUND')
    else:
        #for filename in fileList:
        #    print "found:\t\t{}".format(filename)
        return fileList


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

