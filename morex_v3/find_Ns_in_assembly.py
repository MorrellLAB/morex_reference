#!/usr/bin/env python3
"""Find places in the reference genome assembly where there
are stretches of Ns and prints output to stdout in BED format.

Usage: ./find_Ns_in_assembly.py [fasta_fp] > out_file.bed

Where:
1) [fasta_fp] is the full filepath to the assembly in fasta format.
"""

import sys
import os
import re
from Bio import SeqIO

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)

# User provided input arguments
fasta_fp = os.path.expanduser(sys.argv[1])

# Open FASTA, search for masked regions, print in BED3 format
with open(fasta_fp, 'r') as handle:
    for record in SeqIO.parse(handle, "fasta"):
        for match in re.finditer('N+', str(record.seq)):
            print(record.id, match.start(), match.end(), sep='\t')
