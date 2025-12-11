#!/bin/bash -l
#==============================================================================
# SLURM DIRECTIVES
#==============================================================================
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100g
#SBATCH --tmp=100g
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
module load pigz/2.7-gcc-13.1.0-bwi6ag3

ANCHORWAVE="/users/6/pmorrell/Apps/anchorwave/anchorwave"

# File paths - Morex_v2 to Morex_v3 liftover
OUTPUT_DIR="/scratch.global/pmorrell/Ancestral/AnchorWave"
REF_FASTA="/projects/standard/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta"
#QUERY_FASTA_GZ="/scratch.global/pmorrell/Ancestral/FB19_011_3/FB19_011_3_hap1.fasta"
QUERY_FASTA="/scratch.global/pmorrell/Ancestral/FB19_011_3/FB19_011_3_hap1.fasta"
CDS_FASTA="/scratch.global/pmorrell/Ancestral/Anchorwave/Morex_v3_CDS.fasta"
GFF3="/scratch.global/pmorrell/Ancestral/Hv_Morex.pgsb.Jul2020.re_sorted.gff3"
CDS_SAM="${OUTPUT_DIR}/cds.sam"
CDS_REF="${OUTPUT_DIR}/ref.sam"

mkdir -p "${OUTPUT_DIR}"

log() {
    local msg="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - ${msg}"
}

#==============================================================================
# DECOMPRESSION
#==============================================================================

#log "Step 0: Decompressing to local scratch..."
#pigz -dc -p "${THREADS}" "${QUERY_FASTA_GZ}" > "${QUERY_FASTA}"

#==============================================================================
# genoAli EXECUTION
#==============================================================================

log "Step 1: Running AnchorWave genoAli (Morex_v2 -> Morex_v3)"
$ANCHORWAVE genoAli -i $GFF3 -as $CDS_FASTA -r $REF_FASTA -a $CDS_SAM -ar $CDS_REF \
    -s $QUERY_FASTA -n ${OUTPUT_DIR}/anchors -o ${OUTPUT_DIR}/anchorwave.maf \
    -f ${OUTPUT_DIR}/anchorwave.f.maf -IV -m 0 -t $THREADS

# Cleanup
rm "${QUERY_FASTA}"

log "Done! MAF output files in ${OUTPUT_DIR}"
log "Note: Convert anchorwave.maf to chain format for liftOver"
