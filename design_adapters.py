#!/usr/bin/env python

"""
Note: Python 3

Design adapter sequences for use in a sequencing application

Inputs:
   forward primer sequence
   reverse primer sequence
   file containing potential index sequences
   file containing previously used index sequences

Outputs:
   file with newly designed adapter sequences

Ben Ober-Reynolds
"""

import os
import sys
import argparse
import re
import random
import xlrd
import pandas as pd


### GLOBAL VARS ###

# Illumina adapter sequences

C_R2_orientation = "AATGATACGGCGACCACCGAGATCTACAC"
D_R1_orientation = "CAAGCAGAAGACGGCATACGAGAT"

i5_indicator = "5-"
i7_indicator = "7-"

def main():
    
    # set up command line argument parser
    parser = argparse.ArgumentParser(description='Script for designing adapters')
    group = parser.add_argument_group('required arguments:')
    group.add_argument('-fp', '--forward_primer', required=True,
        help='sequence of forward primer')
    group.add_argument('-rp', '--reverse_primer', required=True,
        help='sequence of reverse primer')
    group.add_argument('-bn', '--base_name', required=True,
        help='base name for new adapters')
    group.add_argument('-n', '--adapter_pairs', required=True,
        help='number of adapter pairs to generate')
    group = parser.add_argument_group('optional arguments')
    group.add_argument('-af', '--adapter_file', default="/Users/boberry/Documents/git_clones/sequence_analysis/edittag_barcodes_10nt.xls",
        help='excel file containing adapters to pick from')
    group.add_argument('-bf', '--blacklist_file', default=None,
        help='txt file containing blacklisted adapter sequences')
    group.add_argument('-il', '--index_length', default=8,
        help='desired index length (default is 8)')
    group.add_argument('-ed', '--edit_distance', default=4,
        help='desired edit distance between tags (default is 4)')


    # print help if no arguments provided
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    # parse command line arguments
    args = parser.parse_args()

    # Read in parameters
    fwd_primer = args.forward_primer.upper()
    rev_primer = args.reverse_primer.upper()
    adapter_pairs = int(args.adapter_pairs)
    idx_len = int(args.index_length)
    idx_editdist = int(args.edit_distance)

    if idx_len > 10:
        print("Index lengths cannot exceed 10")
        sys.exit()
    if idx_editdist >= idx_len:
        print("Edit distance must be less than index length")
        sys.exit()

    # Read in all indexes
    index_df_dict = pd.read_excel(args.adapter_file, sheetname=None)

    # Construct tree of indexes for convenient selection of features
    index_search_dict = build_index_tree(index_df_dict)
    
    # Read in blacklist, if provided
    blacklist = set()
    if args.blacklist_file:
        blacklist = set([x.upper() for x in pd.read_table(args.blacklist_file).sequence.tolist()])
        blacklist = blacklist.update(set([rev_comp(x) for x in blacklist]))

    # Generate adapter sequences
    possible_indices = set(index_search_dict[idx_len][idx_editdist]) - blacklist
    try:
        final_indices = set(random.sample(possible_indices, 2*adapter_pairs))
    except ValueError:
        # Not enough indicees
        print("Not enough indicees available with provided settings")
        sys.exit()

    fwd_adapters = format_new_adapters(adapter_pairs, final_indices, 
        C_R2_orientation, fwd_primer, args.base_name, forward=True)

    rev_adapters = format_new_adapters(adapter_pairs, final_indices, 
        D_R1_orientation, rev_primer, args.base_name, forward=False)

    all_adapters = fwd_adapters + rev_adapters

    # Save final results
    adapter_df = pd.DataFrame(all_adapters)
    adapter_df.columns = ["adapter_name", "full_sequence", "NextSeq_index", "MiSeq_index"]

    output_file = "{}_designed_adapters.txt".format(args.base_name)
    print("Adapters successfully generated. Saving file to {}".format(os.getcwd()+"/"+output_file))
    adapter_df.to_csv(output_file, sep='\t', index=None)





def build_index_tree(df_dict):
    """
    Construct search tree of indexes
    Inputs:
        df_dict (dict) -- dictionary of DataFrames from index excel sheet
    Outputs:
        search_dict (dict) -- search tree of indexes
    """
    search_dict = {}
    for key, df in df_dict.items():
        length, editdist = [int(x) for x in re.findall(r'\d+', key)]
        if length in search_dict:
            search_dict[length][editdist] = df.sequence.tolist()
        else:
            search_dict[length] = {editdist: df.sequence.tolist()}
    return search_dict


def format_new_adapters(n, idx_set, adapter, primer, base_name, forward):
    """
    Return components of new adapter
    Inputs:
        n (int) -- number of adapters to generate
        idx_set (set) -- set of possible indices
        adapter (str) -- sequence of base adapter
        primer (str) -- sequence of primer
        base_name (str) -- base name of new adapters
        forward (bool) -- is this the forward adapter (i5)?
    Outputs:
        adapter_list (list) - list of lists of adapter components
    """
    new_adapters = []
    if forward:
        idx_indicator = i5_indicator
    else:
        idx_indicator = i7_indicator
    for i in range(n):
        idx = idx_set.pop()
        full_seq = adapter + idx + primer
        new_name = base_name + idx_indicator + str(i + 1)
        idx_seq = rev_comp(idx)
        alt_idx = idx_seq
        if forward:
            # MiSeq i5 is the same orientation as the C adapter
            alt_idx = idx
        new_adapters.append([new_name, full_seq, idx_seq, alt_idx])
    return new_adapters



def rev_comp(seq, complement_dict={'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}, missing='N'):
    # Reverse complement a sequence
    return "".join(complement_dict.get(base, missing) for base in reversed(seq))


if __name__ == '__main__':
    main()