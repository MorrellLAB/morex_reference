#!/bin/bash
#PBS -l mem=22gb,nodes=1:ppn=6,walltime=03:00:00
#PBS -m abe
#PBS -M liux1299@umn.edu
#PBS -q mesabi

set -e
set -o pipefail

module load bwa/0.7.17

# This script does a quick alignment as a check for the reference genome

# User provided input arguments
FASTA_FILE=/panfs/roc/groups/9/morrellp/shared/Datasets/Genotyping/Contextual_Sequences/BOPA_SNPs/BOPA_SNPs.fasta
REF_FILE=/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules.fasta
OUT_PREFIX=BOPA_contextual_morex_v2-bwa-mem
OUT_DIR=~/scratch/barley_ref/testing

# Check if our dir exists, if not make it
mkdir -p ${OUT_DIR}

# Go into reference dir
cd ${OUT_DIR}

# Align with bwa mem
bwa mem ${REF_FILE} ${FASTA_FILE} > ${OUT_DIR}/${OUT_PREFIX}.sam

