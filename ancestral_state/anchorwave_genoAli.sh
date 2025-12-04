#!/bin/bash -l
#==============================================================================
# SLURM DIRECTIVES (Optimized for 10-thread GMAP job)
#==============================================================================
#SBATCH --time=24:00:00
#SBATCH --ntasks=1              # Crucial: Use 1 task for a multithreaded (OpenMP) job
#SBATCH --cpus-per-task=10      # Crucial: Request 10 cores to match THREADS variable
#SBATCH --mem=320g
#SBATCH --tmp=320g
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

ANCHORWAVE="/users/6/pmorrell/Apps/anchorwave/anchorwave"

# File paths
OUTPUT_DIR="/scratch.global/pmorrell/Ancestral/Anchorwave"
REF_FASTA="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta"
QUERY_FASTA="/scratch.global/pmorrell/Ancestral/FB19_011_3/FB19_011_3_hap1.fasta"
CDS_FASTA="/scratch.global/pmorrell/Ancestral/Anchorwave/Morex_v3_CDS.fasta"
GFF3="/scratch.global/pmorrell/Ancestral/sorted.gff3"
CDS_SAM="/scratch.global/pmorrell/Ancestral/Anchorwave/cds.sam"
CDS_REF="/scratch.global/pmorrell/Ancestral/Anchorwave/ref.sam"

mkdir -p "${OUTPUT_DIR}"

log() {
    local msg="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - ${msg}"
}


#==============================================================================
# genoAli EXECUTION
#==============================================================================

log "Step 1: Anchorwave align reference to query"
$ANCHORWAVE genoAli -i $GFF3 -as $CDS_FASTA -r $REF_FASTA -a $CDS_SAM -ar $CDS_REF \
    -s $QUERY_FASTA -n anchors -o ${OUTPUT_DIR}/anchorwave.maf \
    -f ${OUTPUT_DIR}/anchorwave.f.maf -IV -m 0 -t $THREADS

log "Done! Output files in ${OUTPUT_DIR}"

