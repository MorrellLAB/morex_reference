#!/bin/bash

set -e
set -o pipefail

# This script prints a summary of the total number of SNPs, number of unique SNPs, number of duplicate SNPs, and the number of failed SNPs after running SNP_Utils

# Usage: ./summarize_snp_count.sh [vcf] [failed_log_file]

# User provided input arguments
vcf="$1"
failed_log_file="$2"

#------------------------
# Total number of SNPs
total_snps=$(grep -v "#" ${vcf} | cut -f 3 | wc -l)

# Total number of unique SNPs
total_uniq_snps=$(grep -v "#" ${vcf} | cut -f 3 | sort -Vu | wc -l)

# Total number of duplicates
total_dup=$(grep -v "#" ${vcf} | cut -f 3 | sort -V | uniq -c | sort -n -r | grep -vw "1" | wc -l)

# Total number of failed SNPs
failed_snps=$(wc -l ${failed_log_file} | cut -d' ' -f 1)

# Send summaries to stdout
echo "Total SNPs: ${total_snps}"
echo "Total unique SNPs: ${total_uniq_snps}"
echo "Total duplicate SNPs: ${total_dup}"
echo "Total failed SNPs: ${failed_snps}"
