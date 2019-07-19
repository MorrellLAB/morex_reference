#!/bin/bash

set -e
set -o pipefail

# User provided input arguments
# List of SNPs, one SNP name per line
SNP_LIST=$1
# File extension of SNP_LIST. Ex: .log, .txt
SNP_LIST_EXT=$2
# Full filepath to contextual sequences for SNPs
CONTEXTUAL_FASTA=$3
# Full filepath to our output directory
OUT_DIR=$4

# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR}

# Store SNP list in array
SNP_ARR=($(cat ${SNP_LIST}))

# Extract sequence for SNPs
for i in ${SNP_ARR[@]}
do
    prefix=$(basename ${SNP_LIST} ${SNP_LIST_EXT})
    # Identify line number of SNP in FASTA file
    tmp_line=$(grep -n ${i} ${CONTEXTUAL_FASTA} | cut -d':' -f 1)
    tmp_seq_line=$((${tmp_line}+1))
    # Extract sequence associated with SNP
    sed -n ${tmp_line},${tmp_seq_line}p ${CONTEXTUAL_FASTA} >> ${OUT_DIR}/${prefix}.fasta
done
