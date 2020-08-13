#!/bin/bash

set -e
set -o pipefail

# Fix VCF header by adding contig lengths to VCF header
# Many GATK tools require the VCF headers to contain contig lengths

function Usage() {
    echo -e "\
    Usage:  ./evaluate_variant_filtering.sh [VCF] [REF_DICT]
    Where:
        [VCF] is the full filepath to the VCF file we are evaluating
        [REF_DICT] is the full filepath to the reference dictionary
    
    The fixedHeader.vcf file will be output in the same directory as the VCF.
    Note: The Picard Jar file in the script may need to be modified if you are not using MSI.
    " >&2
    exit 1
}

export -f Usage

# If we have less than 2 arguments, display the usage message
if [[ "$#" -lt 2 ]]; then
    Usage
fi

# Dependencies
# Note: Picard Jar filepath will need to be modified if not running on MSI
PICARD_JAR=/home/morrellp/public/Software/picard_ML_2.20.2/picard.jar

# User provided input arguments
VCF="$1"
REF_DICT="$2"

# Generate out dir path
OUT_DIR=$(dirname ${VCF})
# Prepare output file prefix
if [[ ${VCF} == *.vcf.gz ]]; then
    OUT_PREFIX=$(basename ${VCF} .vcf.gz)
else
    # Asssume vcf files ends in .vcf extension
    OUT_PREFIX=$(basename ${VCF} .vcf)
fi

# Fix VCF file header (contigs don’t have lengths in the header)
java -jar ${PICARD_JAR} UpdateVcfSequenceDictionary \
    I=${VCF} \
    O=${OUT_DIR}/${OUT_PREFIX}_fixedHeader.vcf \
    SD=${REF_DICT}
