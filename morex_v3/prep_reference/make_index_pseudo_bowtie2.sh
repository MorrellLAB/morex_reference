#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=48gb
#SBATCH --tmp=12gb
#SBATCH -t 06:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

module load bowtie2/2.3.4.1

# User provided input arguments
REF="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta"

# Go into reference dir
REF_DIR=$(dirname ${REF})
cd ${REF_DIR}

# Generate database prefix
PREFIX=$(basename ${REF} .fasta)

# Create bowtie2 index database
bowtie2-build ${REF} ${PREFIX}

# Check the content of the database
bowtie2-inspect -n ${PREFIX}
