#!/usr/bin/env python3

"""This script was intended to take the barley 50k SNP_Utils lookup table and
convert the snps formatted as [A/B] to IUPAC nucleotide base code and return
a fasta file. This is so we can do an IPK Barley BLAST search to resolve
duplicate/failed SNPs. IPK BLAST does not accept sequences with snps formatted
as [A/B].

Usage: ./snp_to_iupac_fasta.py [snp_utils_lookup_table] > out_filename.fasta

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

# IUPAC nucleotide base code: [Base, Reverse complement]
IUPAC_TABLE = {
    'R': [['A', 'G'], 'Y'],
    'Y': [['C', 'T'], 'R'],
    'S': [['G', 'C'], 'S'],
    'W': [['A', 'T'], 'W'],
    'K': [['T', 'G'], 'M'],
    'M': [['A', 'C'], 'K'],
    'B': [['C', 'G', 'T'], 'V'],
    'D': [['A', 'G', 'T'], 'H'],
    'H': [['A', 'C', 'T'], 'D'],
    'V': [['A', 'C', 'G'], 'B'],
    'N': [['A', 'C', 'G', 'T'], 'N']
}


def main(snp_table):
    """Driver function"""
    with open(os.path.expanduser(snp_table), "rt") as file:
        for line in file:
            tmp = line.strip().split('\t')
            # Take the sequence and split up SNP in format [A/B]
            tmp_seq1 = tmp[1].split(']')
            tmp_seq = tmp_seq1[0].split('[')
            tmp_seq.append(tmp_seq1[1])
            # Convert A/B format to IUPAC code
            snp = tmp_seq[1].split('/')
            for key in IUPAC_TABLE:
                if sorted(snp) == sorted(IUPAC_TABLE[key][0]):
                    iupac_snp = key
            new_seq = ''.join([tmp_seq[0], iupac_snp, tmp_seq[2]])
            print('>' + tmp[0])
            print(new_seq)


main(sys.argv[1]) # Run the program
