#!/bin/bash
#PBS -l mem=22gb,nodes=1:ppn=4,walltime=03:00:00
#PBS -m abe
#PBS -M liux1299@umn.edu
#PBS -q mesabi

set -e
set -o pipefail

module load bwa/0.7.17

# User provided input arguments
REF_DIR=/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2
REF_FILENAME=Barley_Morex_V2_pseudomolecules.fasta

cd ${REF_DIR}

bwa index ${REF_DIR}/${REF_FILENAME}
