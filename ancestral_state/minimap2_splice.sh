#!/bin/bash -l
#==============================================================================
# SLURM DIRECTIVES
#==============================================================================
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=48g
#SBATCH --tmp=48g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

#==============================================================================
# CONFIGURATION
#==============================================================================
THREADS=10

# Exit immediately if a command exits with a non-zero status
set -e
# Exit if any command in a pipeline fails
set -o pipefail

# Peter L. Morrell - Revised Dec 2025

module purge
module load minimap2/2.30   # Load minimap2 (adjust version as needed)
module load pigz/2.7-gcc-13.1.0-bwi6ag3

# --- 1. Define Genotype Names ---
REF_NAME="Morex_v3"
QUERY_NAME="FB19_011_3_hap1"

# --- 2. Define Source File Paths (Compressed) ---
REF_SOURCE="/projects/standard/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta.gz"
QUERY_SOURCE="/scratch.global/pmorrell/Ancestral/FB19_011_3/FB19_011_3_hap1.fasta"
CDS_FASTA="/scratch.global/pmorrell/Ancestral/Anchorwave/Morex_v3_CDS.fasta"

# --- 3. Define Temporary Output Paths (Fast Local Scratch) ---
REF_FASTA_TMP="${TMPDIR}/${REF_NAME}.fasta"
QUERY_FASTA_TMP="${TMPDIR}/${QUERY_NAME}.fasta"

# --- 4. Define Final Output Directory ---
OUTPUT_DIR="/scratch.global/pmorrell/Ancestral/AnchorWave"
mkdir -p "${OUTPUT_DIR}"

log() {
    local msg="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - ${msg}"
}

# Output directories
REF_DIR="${OUTPUT_DIR}/${REF_NAME}"
QUERY_DIR="${OUTPUT_DIR}/${QUERY_NAME}"

mkdir -p "${REF_DIR}"
mkdir -p "${QUERY_DIR}"

#==============================================================================
# DECOMPRESSION
#==============================================================================

log "Step 0: Preparing genomes in local scratch..."

# Handle QUERY not compressed -
#log "Decompressing ${QUERY_NAME}..."
#pigz -dc -p "${THREADS}" "${QUERY_SOURCE}" > "${QUERY_FASTA_TMP}"
QUERY_FASTA_TMP="${QUERY_SOURCE}"

# Handle REF (Morex_v3)
if [[ "${REF_SOURCE}" == *.gz ]]; then
    log "Decompressing ${REF_NAME} (detected .gz)..."
    pigz -dc -p "${THREADS}" "${REF_SOURCE}" > "${REF_FASTA_TMP}"
else
    log "Copying ${REF_NAME} to local scratch for speed..."
    cp "${REF_SOURCE}" "${REF_FASTA_TMP}"
fi

#==============================================================================
# MINIMAP2 EXECUTION
#==============================================================================

log "Step 1: Mapping CDS to ${REF_NAME} reference using minimap2..."
# -ax splice: cDNA/mRNA alignment to a genome (splice-aware)
# -k14: Use a shorter k-mer for higher sensitivity, recommended for accurate transcripts
minimap2 -x splice -k12 -t "${THREADS}" -a -p 0.4 -N 20 "${REF_FASTA_TMP}" "${CDS_FASTA}" > "${OUTPUT_DIR}/ref.sam"

log "Step 2: Mapping CDS to ${QUERY_NAME}..."
# Use the same parameters for consistency
minimap2 -x splice -k12 -t "${THREADS}" -a -p 0.4 -N 20 "${QUERY_FASTA_TMP}" "${CDS_FASTA}" > "${OUTPUT_DIR}/cds.sam"

# Cleanup Temp Files
rm "${REF_FASTA_TMP}" "${QUERY_FASTA_TMP}"

log "Done! Output files in ${OUTPUT_DIR}"
