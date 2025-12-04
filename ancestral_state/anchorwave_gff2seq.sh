#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --ntasks=4
#SBATCH --mem=32g
#SBATCH --tmp=32g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

module load pigz

# Peter L. Morrell - 28 November 2025 - Falcon Heights, MN

ANCHORWAVE="/users/6/pmorrell/Apps/anchorwave/anchorwave"
OUTPUT_DIR="/scratch.global/pmorrell/Ancestral/Anchorwave"
REF_FASTA="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta"
GFF3="/scratch.global/pmorrell/Ancestral/sorted.gff3"

mkdir -p "${OUTPUT_DIR}"

log() {
    local msg="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - ${msg}"
}

$ANCHORWAVE gff2seq -i $GFF3 -r $REF_FASTA -o $OUTPUT_DIR/Morex_v3_CDS.fasta
