#!/usr/bin/env python3

"""This script takes the 50k SNP_Utils lookup table and converts
it to fasta format for SNPMeta (https://github.com/MorrellLAB/SNPMeta).
Currently, this script outputs a fasta file using the illumina assay
format (example: [A/G] SNP).

Usage: ./convert_to_fasta.py [snp_utils_lookup_table] > out_filename.fasta

Where:
[snp_utils_lookup_table] is a file output from SNP_Utils and contains two columns:
    1) SNP name, 2) sequence where SNP is formatted as [A/B]

Example line from SNP_Utils lookup table:
BK_01	CCTCGTTTCACCAAATAGCTTCTTTATTGCAGMATYGTCGAGAGCGTCGGAGAGGGCGTGACTGAGCTTGTGCCGGGYGACCATGTMCTCCCGGT[T/G]TTCACCGGCGAGTGCAAGGACTGTGCCCACTGCAAGTCAGAGGAGAGCAACCTTTGTGATCTCCTTAGGATCAATG
"""

import os
import sys

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)


def main(snp_table):
    """Driver function."""
    with open(os.path.expanduser(snp_table), "rt") as file:
        for line in file:
            tmp = line.strip().split('\t')
            # Take the sequence and reformat to fasta format
            # Print sequence descriptor line
            print('>' + tmp[0])
            # Print sequence
            print(tmp[1])


main(sys.argv[1]) # Run the program
