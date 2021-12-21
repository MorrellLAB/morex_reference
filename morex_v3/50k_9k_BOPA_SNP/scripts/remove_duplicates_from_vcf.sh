#!/bin/bash

set -e
set -o pipefail

# Dependencies
module load python3/3.8.3_anaconda2020.07_mamba

# User provided input arguments
VCF="$1"
PREFIX="$2"
OUT_DIR="$3"
PSEUDO_TO_PARTS_SCRIPT="$4"
REF_VERSION="$5"

#------------------------
# Check if out dir and out subdir exist, if not make them
mkdir -p ${OUT_DIR} ${OUT_DIR}/duplicate_snps

# Make a list of duplicate SNP names
grep -v "#" ${VCF} | cut -f 3 | sort | uniq -c | sort -n -r | grep -vw "1" | sed -e 's,      ,,' | cut -d' ' -f 2 | sort -V > ${OUT_DIR}/duplicate_snps/${PREFIX}_duplicate_snp_names.txt
# Save duplicates to vcf
grep -f ${OUT_DIR}/duplicate_snps/${PREFIX}_duplicate_snp_names.txt ${VCF} | sort -V -k3,3 > ${OUT_DIR}/duplicate_snps/${PREFIX}_duplicates.vcf
# Make a copy for resolving duplicates
cp ${OUT_DIR}/duplicate_snps/${PREFIX}_duplicates.vcf ${OUT_DIR}/duplicate_snps/${PREFIX}_duplicates_resolved.vcf

# Remove duplicates from main VCF
grep -vf ${OUT_DIR}/duplicate_snps/${PREFIX}_duplicate_snp_names.txt ${VCF} > ${OUT_DIR}/${PREFIX}_noRescuedSNPs.vcf

# Convert pseudomolecular to parts positions
${PSEUDO_TO_PARTS_SCRIPT} --vcf ${OUT_DIR}/${PREFIX}_noRescuedSNPs.vcf ${REF_VERSION} > ${OUT_DIR}/${PREFIX}_noRescuedSNPs_partsRef.vcf
