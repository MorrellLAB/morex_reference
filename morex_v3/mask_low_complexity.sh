#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=12gb
#SBATCH --tmp=8gb
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

# Dependencies
module load java/openjdk-8_202
module load python3/3.8.3_anaconda2020.07_mamba

# User provided input arguments
REF_FASTA="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"
# See bbmask documentation for picking entropy value: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmask-guide/
# Experiment with multiple values
ENTROPY1="0.7"
ENTROPY2="0.8"
ENTROPY3="0.9"
# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/entropy_masked"
# Directory that contains suite of BBTools
BBTOOLS_DIR="/panfs/jay/groups/9/morrellp/shared/Software/bbmap"
PYTHON_SCRIPT_DIR="/panfs/jay/groups/9/morrellp/liux1299/GitHub/morex_reference/morex_v3"

#------------------------------
# Define function
function mask_low_complexity() {
    local ref_fasta="$1"
    local entropy_val="$2"
    local out_prefix="$3"
    local out_dir="$4"
    local bbtools_dir="$5"
    local python_script_dir="$6"
    # Mask low complexity regions
    bash ${bbtools_dir}/bbmask.sh in=${ref_fasta} out=${out_dir}/${out_prefix}.entropy_${entropy_val}_masked.fasta entropy=${entropy_val}
    # Masked bases are converted to N
    # Get BED file of regions that are masked
    python ${python_script_dir}/find_Ns_in_assembly.py ${out_dir}/${out_prefix}.entropy_${entropy_val}_masked.fasta > ${out_dir}/${out_prefix}.entropy_${entropy_val}_masked.bed
}

export -f mask_low_complexity

# Make output directory
mkdir -p ${OUT_DIR}

# Get fasta file prefix
if [[ "${REF_FASTA}" == *"fasta" ]]; then
    OUT_PREFIX=$(basename ${REF_FASTA} .fasta)
elif [[ "${REF_FASTA}" == *"fa" ]]; then
    OUT_PREFIX=$(basename ${REF_FASTA} .fa)
fi

# Mask low complexity regions - ENTROPY1
mask_low_complexity ${REF_FASTA} ${ENTROPY1} ${OUT_PREFIX} ${OUT_DIR} ${BBTOOLS_DIR} ${PYTHON_SCRIPT_DIR}

# Mask low complexity regions - ENTROPY2
mask_low_complexity ${REF_FASTA} ${ENTROPY2} ${OUT_PREFIX} ${OUT_DIR} ${BBTOOLS_DIR} ${PYTHON_SCRIPT_DIR}

# Mask low complexity regions - ENTROPY3
mask_low_complexity ${REF_FASTA} ${ENTROPY3} ${OUT_PREFIX} ${OUT_DIR} ${BBTOOLS_DIR} ${PYTHON_SCRIPT_DIR}
