#!/bin/bash
#PBS -l mem=40gb,nodes=1:ppn=6,walltime=06:00:00
#PBS -m abe
#PBS -M liux1299@umn.edu
#PBS -q mesabi

set -e
set -o pipefail

module load bowtie2/2.3.2

# This script does a quick alignment as a check for the reference genome
# Following what Li did for morex v1: https://github.com/lilei1/9k_BOPA_SNP/blob/master/script/commandlines

# User provided input arguments
FASTA_FILE=/panfs/roc/groups/9/morrellp/shared/Datasets/Genotyping/Contextual_Sequences/BOPA_SNPs/BOPA_SNPs.fasta
DB_NAME=/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules
OUT_PREFIX=BOPA_contextual_morex_v2-bowtie2
OUT_DIR=~/scratch/barley_ref/testing

# Check if our dir exists, if not make it
mkdir -p ${OUT_DIR}

# Go into reference dir
cd ${OUT_DIR}

# Align with bowtie2
bowtie2 -f ${FASTA_FILE} -x ${DB_NAME} -S ${OUT_PREFIX}.sam
