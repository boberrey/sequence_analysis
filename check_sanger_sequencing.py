#!/usr/bin/env python

"""
Note: Python 3

Check a directory of sanger sequencing outputs against a file of known 
sequences.

Inputs:
   directory containing sanger sequencing outputs in fasta format
   A tab-delimited 2-column text file of known sequences

Outputs:
   prints out results of sanger sequencing alignment

Ben Ober-Reynolds
"""

import os
import sys
import argparse
import pandas as pd
from Bio import pairwise2


# Local alignment parameters
match_score = 1
mismatch_score = 0
gap_open = -0.5
gap_extension = -0.1
match_threshold = 0.9

def main():
    
    # set up command line argument parser
    parser = argparse.ArgumentParser(description='Script for analyzing sanger \
        sanger sequencing results.')
    group = parser.add_argument_group('required arguments:')
    group.add_argument('-sd', '--seq_directory', required=True,
        help='directory containing sanger sequencing output')
    group.add_argument('-mf', '--match_file', required=True,
        help='file containing sequences to match (2 column, tab-delimited format)')
    group = parser.add_argument_group('optional arguments')
    group.add_argument('-ex', '--file_extension', default="seq",
        help='specify the file extension for sanger sequencing results, without the "." \
        (The default is "seq")')
    group.add_argument('-mm', '--max_matches', default=1,
        help='how many potential matches should be reported. If more matches are found, \
        the top scoring matches will be reported.')
    group.add_argument('-v', '--verbose', default=False, action='store_true',
        help='how many potential matches should be reported. If more matches are found, \
        the top scoring matches will be reported.')



    # print help if no arguments provided
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    # parse command line arguments
    args = parser.parse_args()

    # Gather sanger sequencing files
    sanger_files = find_files_in_directory(args.seq_directory, extensionList=[args.file_extension])

    # Read in match sequences
    match_seq_dict = read_in_match_seqs(args.match_file)

    # Loop through fasta files and report on matches
    print("Aligning sequences in {} to sequencing data...".format(args.match_file))

    # Loop through all sanger fasta files
    final_matches = []
    failed_alignment_counter = 0

    for fasta_file in sanger_files:
        sequenced_dict = read_fasta_file(fasta_file)
        # Loop through each sequence found in the fasta file (probably just 1)
        for seq_name, exp_seq in sequenced_dict.items():
            alignments = []
            # Loop through all possible match sequences
            for match_seq, match_name in match_seq_dict.items():
                best_alignment = pairwise2.align.localms(exp_seq, match_seq,
                    match_score, mismatch_score, gap_open, gap_extension)[0]
                # Ignore match if it doesn't meet threshold
                if float(best_alignment[2])/len(match_seq) < match_threshold:
                    continue

                # Match found. Report it (if indicated) and add it to saved alignments
                if args.verbose:
                    match_subseq = exp_seq[best_alignment[3]-5:best_alignment[4]+6]
                    pretty_alignment = pairwise2.align.localms(match_subseq, match_seq, 1, 0, -1, -1)[0]
                    print("Match found for {} to {}".format(match_name, seq_name))
                    print("Scored {} of possible {}".format(pretty_alignment[2], len(match_seq)))
                    print(pairwise2.format_alignment(*pretty_alignment))

                alignments.append([match_name] + [len(match_seq)] + list(best_alignment[2:]))
            # Add all alignments to the final alignment structure.
            # Sort alignments by score and limit to specified number of alignments
            if len(alignments) >= 1:
                alignments.sort(key=lambda x: x[2], reverse=True)
            else:
                alignments = [['NA']*4]
                failed_alignment_counter += 1
            if len(alignments) > args.max_matches:
                alignments = alignments[:args.max_matches]
            for aligned in alignments:
                final_matches.append([seq_name] + aligned)

    # Construct final alignment report
    print("Found potential alignments for {} of {} sequencing files".format(
        len(sanger_files) - failed_alignment_counter, len(sanger_files)))

    alignment_df = pd.DataFrame(final_matches)
    alignment_df.columns = ["seq_name", "match_name", "max_possible_score", "alignment_score", "start", "end"]
    alignment_df.to_csv(args.seq_directory.strip('/') + "/alignment_results.tsv", sep='\t', index=None)




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

def read_in_match_seqs(match_file):
    """
    Read in the sequences to use for matching against found fasta files
    Inputs:
        match_file (str) -- filepath for sequences to match. Match file must
            be tab-delimited with the first column being the sequences name and the
            second column being the actual sequence
    Outputs:
        seq_dict (dict) -- dictionary mapping sequences to seq names
    """
    seq_dict = {}
    with open(match_file, 'r') as f:
        for line in f:
            name, seq = line.strip().split()
            seq = seq.upper()
            seq_dict[seq] = name
            seq_dict[rev_comp(seq)] = "RC_" + name
    return seq_dict


def read_fasta_file(fasta_file, name_indicator=">"):
    """
    Extract a dict mapping names to sequences from a fasta file. 
    (Note: this implementation stores sequences as single strings. Don't use
        this code for very long (> a few kb) sequences!)
    Inputs:
        fasta_file (str) -- fasta filepath
    Outputs:
        fasta_dict (dict) -- fasta dict mapping names to sequences
    """
    fasta_dict = {}
    with open(fasta_file, 'r') as f:
        current_seq = ""
        current_name = ""
        for line in f:
            if line.startswith(name_indicator):
                # New sequence name found. Store current information if necessary
                if current_name != "":
                    fasta_dict[current_name] = current_seq.upper()
                # flush previous name and seq
                current_name = line.strip()[1:]
                current_seq = ""
            else:
                # not a new name -- continue building sequence
                current_seq += line.strip()
        # EOF, so add final sequence if necessary
        if current_name != "":
            fasta_dict[current_name] = current_seq.upper()
    return fasta_dict


def rev_comp(seq, complement_dict={'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}, missing='N'):
    # Reverse complement a sequence
    return "".join(complement_dict.get(base, missing) for base in reversed(seq))


if __name__ == '__main__':
    main()