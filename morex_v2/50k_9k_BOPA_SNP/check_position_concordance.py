#!/usr/bin/env python3

"""This script takes in 2 VCF files and was specifically intended
to check for agreement in physical positions between the BOPA and
9k SNPs that overlap. These SNP positions were identified with
SNP_Utils initially. The duplicate SNPs were manually resolved and
the failed SNPs were programatically resolved.

So, this script serves as a quick check to make sure overlapping SNP
positions agree.

Usage: ./check_position_concordance.py [BOPA_VCF] [9k_VCF]
"""

import os
import sys
import gzip

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)


def read_vcf(filename):
    """Read in a VCF file and store in dictionary."""
    vcf_dict = {}
    if '.gz' in filename:
        with gzip.open(filename, "rt") as file:
            for line in file:
                if line.startswith('#'):
                    continue
                else:
                    tmp = line.strip().split('\t')
                    vcf_dict[tmp[2]] = [tmp]
    else:
        with open(filename, "rt") as file:
            for line in file:
                if line.startswith('#'):
                    continue
                else:
                    tmp = line.strip().split('\t')
                    vcf_dict[tmp[2]] = tmp
    return vcf_dict


def main(VCF_BOPA, VCF_9K):
    """Driver function."""
    # Read in VCF files
    vcf_bopa = read_vcf(os.path.expanduser(VCF_BOPA))
    vcf_9k = read_vcf(os.path.expanduser(VCF_9K))
    # Check that positions match, if not print to stdout
    for key in vcf_9k:
        if key in vcf_bopa:
            if vcf_bopa[key][0] != vcf_9k[key][0] or vcf_bopa[key][1] != vcf_9k[key][1]:
                print('9k\t', '\t'.join(vcf_9k[key]))
                print('bopa\t', '\t'.join(vcf_bopa[key]))


main(sys.argv[1], sys.argv[2]) # Run the program
