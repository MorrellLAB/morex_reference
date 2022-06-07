#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=22gb
#SBATCH --tmp=12gb
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

module load bwa/0.7.17

# User provided input arguments
REF="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_plastids_parts.fasta"

#-----------
# Generate directory path only
REF_DIR=$(dirname ${REF})
cd ${REF_DIR}

# Index reference
bwa index ${REF}
