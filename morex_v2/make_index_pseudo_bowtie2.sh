#!/bin/bash
#PBS -l mem=48gb,nodes=1:ppn=4,walltime=06:00:00
#PBS -m abe
#PBS -M liux1299@umn.edu
#PBS -q mesabi

set -e
set -o pipefail

module load bowtie2/2.3.2

# User provided input arguments
REF_DIR=/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2
REF_FILENAME=Barley_Morex_V2_pseudomolecules.fasta

# Go into reference dir
cd ${REF_DIR}

# Generate database prefix
PREFIX=$(basename ${REF_FILENAME} .fasta)

# Create bowtie2 index database
bowtie2-build ${REF_DIR}/${REF_FILENAME} ${PREFIX}

# Check the content of the database
bowtie2-inspect -n ${PREFIX}
