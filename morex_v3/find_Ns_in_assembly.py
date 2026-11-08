#!/usr/bin/env python3
"""Find places in the reference genome assembly where there
are stretches of Ns and prints output to stdout in BED format.

Usage: ./find_Ns_in_assembly.py fasta_fp [--with-length] > out_file.bed

Where:
    fasta_fp: the full filepath to the assembly in fasta format.
    --with-length: (optional) adds a fourth column with the length of each N-run.
"""

import sys
import os
import re
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('fasta_fp', help='Full filepath to the assembly in fasta format')
parser.add_argument('--with-length', action='store_true', help='Add a fourth column with the length of each N-run')

args = parser.parse_args()

# Open FASTA, search for masked regions, print in BED3 format (or BED4 if --with-length)
fasta_fp = os.path.expanduser(args.fasta_fp)
with open(fasta_fp, 'r') as handle:
    for record in SeqIO.parse(handle, "fasta"):
        for match in re.finditer('N+', str(record.seq)):
            if args.with_length:
                run_length = match.end() - match.start()
                print(record.id, match.start(), match.end(), run_length, sep='\t')
            else:
                print(record.id, match.start(), match.end(), sep='\t')
