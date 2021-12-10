#!/bin/bash

set -e
set -o pipefail

# Dependencies
module load python3/3.6.3_anaconda5.0.1
# SNP_Utils directory
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/liux1299/software_development/SNP_Utils
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP/scripts
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts

# User defined input arguments
REF_FASTA="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta"
LOOKUP_TABLE="/panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/lookup_tables/9k_lookup.txt"
GENETIC_MAP="/panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/genetic_maps/BOPA/2013_iSelect_Genetic_Map.txt"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k"
DISTANCE_THRESHOLD="100000"
# The SNP_ARRAY variable will only be used for file naming
SNP_ARRAY="9k"
PERCENT_IDENTITY="95"
OUT_PREFIX="${SNP_ARRAY}_morex_v3_idt${PERCENT_IDENTITY}"

#--------------------
# Go into output directory
cd ${OUT_DIR}

# Configure BLAST options
snp_utils.py CONFIG -d ${REF_FASTA} -k -i ${PERCENT_IDENTITY} -c blast_${OUT_PREFIX}

# Run SNP Utils BLAST
snp_utils.py BLAST -l ${LOOKUP_TABLE} -c blast_${OUT_PREFIX} -b -m ${GENETIC_MAP} -d -t ${DISTANCE_THRESHOLD} -o ${OUT_PREFIX}

# Remove duplicate SNPs from the VCF
# This outputs a VCF with the naming scheme *noRescuedSNPs.vcf
remove_duplicates_from_vcf.sh ${OUT_PREFIX}.vcf ${SNP_ARRAY}_idt${PERCENT_IDENTITY} ${OUT_DIR}
