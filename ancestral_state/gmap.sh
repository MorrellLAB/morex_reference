#!/bin/bash -l
#==============================================================================
# SLURM DIRECTIVES (Optimized for 10-thread GMAP job)
#==============================================================================
#SBATCH --time=24:00:00
#SBATCH --ntasks=1              # Crucial: Use 1 task for a multithreaded (OpenMP) job
#SBATCH --cpus-per-task=10      # Crucial: Request 10 cores to match THREADS variable
#SBATCH --mem=64g
#SBATCH --tmp=64g
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

# Peter L. Morrell - 28 November 2025 - Falcon Heights, MN

module load gmap/2018-07-04

# File paths
OUTPUT_DIR="/scratch.global/pmorrell/Ancestral/Anchorwave"
REF_FASTA="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta"
QUERY_FASTA="/scratch.global/pmorrell/Ancestral/FB19_011_3/FB19_011_3_hap1.fasta"
CDS_FASTA="/scratch.global/pmorrell/Ancestral/Anchorwave/Morex_v3_CDS.fasta"

mkdir -p "${OUTPUT_DIR}"

log() {
    local msg="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - ${msg}"
}

# Output directories
REF_DIR="${OUTPUT_DIR}/Morex_v3"
QUERY_DIR="${OUTPUT_DIR}/FB19_011_3"

mkdir -p "${REF_DIR}"
mkdir -p "${QUERY_DIR}"

#==============================================================================
# GMAP EXECUTION
#==============================================================================

log "Step 1: Building GMAP database for Morex_v3..."
# gmap_build is kept single-threaded as it does not use a -t flag for this purpose
gmap_build -D "${REF_DIR}" -d Morex_v3 --genomedb=Morex_v3 "${REF_FASTA}"

log "Step 2: Building GMAP database for FB19_011_3..."
# gmap_build is kept single-threaded
gmap_build -D "${QUERY_DIR}" -d FB19_011_3 --genomedb=FB19_011_3 "${QUERY_FASTA}"

log "Step 3: Mapping CDS to Morex_v3 reference using ${THREADS} threads..."
gmap -t "${THREADS}" -A -f samse -d Morex_v3 -D "${REF_DIR}" "${CDS_FASTA}" > "${OUTPUT_DIR}/ref.sam"

log "Step 4: Mapping CDS to FB19_011_3 using ${THREADS} threads..."
gmap -t "${THREADS}" -A -f samse -d FB19_011_3 -D "${QUERY_DIR}" "${CDS_FASTA}" > "${OUTPUT_DIR}/cds.sam"

log "Done! Output files in ${OUTPUT_DIR}"

