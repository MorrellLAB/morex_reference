#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=22gb
#SBATCH --tmp=12gb
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

module load bwa/0.7.17

# User provided input arguments
REF="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/assembly/HvulgareMorex_702_V3_parts.fasta"
# Hardmasked REF
REF_HM="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/assembly/HvulgareMorex_702_V3.hardmasked_parts.fasta"
# Softmasked REF
REF_SM="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/assembly/HvulgareMorex_702_V3.softmasked_parts.fasta"

#-----------
# Generate directory path only
# In this case, all 3 references above are in the same directory so only
#   need to run this once
REF_DIR=$(dirname ${REF})
cd ${REF_DIR}

# Index reference
bwa index ${REF}
bwa index ${REF_HM}
bwa index ${REF_SM}
