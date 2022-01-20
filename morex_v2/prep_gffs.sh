#!/bin/bash

set -e
set -o pipefail

# The GFF files need to be sorted otherwise they will cause issues with
#   programs like VeP and ANNOVAR.
# They also need to be converted to parts positions.
# This script serves as a log of commands that were run.

# Dependencies
module load python3/3.8.3_anaconda2020.07_mamba
module load bedops_ML/2.4.38
module load perl/5.28.1
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/liux1299/GitHub/File_Conversions
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/shared/Software/gff3sort

# Directory containing downloaded gene annotations
#   All outputs will also be in this directory
GA_DIR="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/gene_annotation"

# Change into working dir
cd ${GA_DIR}

# Converting pseudomolecules to parts GFF files
# All
Barley_Pseudomolecules_to_Parts.py --gff Barley_Morex_V2_gene_annotation_PGSB.all.gff3 morex_v2 > Barley_Morex_V2_gene_annotation_PGSB.all.parts.gff3
# HC
Barley_Pseudomolecules_to_Parts.py --gff Barley_Morex_V2_gene_annotation_PGSB.HC.gff3 morex_v2 > Barley_Morex_V2_gene_annotation_PGSB.HC.parts.gff3
# LC
Barley_Pseudomolecules_to_Parts.py --gff Barley_Morex_V2_gene_annotation_PGSB.LC.gff3 morex_v2 > Barley_Morex_V2_gene_annotation_PGSB.LC.parts.gff3

# Sorting parts GFF file
# All
gff3sort.pl --precise Barley_Morex_V2_gene_annotation_PGSB.all.parts.gff3 > Barley_Morex_V2_gene_annotation_PGSB.all.parts.sorted.gff3
# HC
gff3sort.pl --precise Barley_Morex_V2_gene_annotation_PGSB.HC.parts.gff3 > Barley_Morex_V2_gene_annotation_PGSB.HC.parts.sorted.gff3
# LC
gff3sort.pl --precise Barley_Morex_V2_gene_annotation_PGSB.LC.parts.gff3 > Barley_Morex_V2_gene_annotation_PGSB.LC.parts.sorted.gff3

# Sorting pseudomolecules GFF file
# All
gff3sort.pl --precise Barley_Morex_V2_gene_annotation_PGSB.all.gff3 > Barley_Morex_V2_gene_annotation_PGSB.all.sorted.gff3
# HC
gff3sort.pl --precise Barley_Morex_V2_gene_annotation_PGSB.HC.gff3 > Barley_Morex_V2_gene_annotation_PGSB.HC.sorted.gff3
# LC
gff3sort.pl --precise Barley_Morex_V2_gene_annotation_PGSB.LC.gff3 > Barley_Morex_V2_gene_annotation_PGSB.LC.sorted.gff3

# GFF3 to BED for parts files
# All
gff2bed < Barley_Morex_V2_gene_annotation_PGSB.all.parts.sorted.gff3 > Barley_Morex_V2_gene_annotation_PGSB.all.parts.sorted.bed
# HC
gff2bed < Barley_Morex_V2_gene_annotation_PGSB.HC.parts.sorted.gff3 > Barley_Morex_V2_gene_annotation_PGSB.HC.parts.sorted.bed
# LC
gff2bed < Barley_Morex_V2_gene_annotation_PGSB.LC.parts.sorted.gff3 > Barley_Morex_V2_gene_annotation_PGSB.LC.parts.sorted.bed

# Cleanup
# Remove unsorted GFF files to avoid future issues
rm Barley_Morex_V2_gene_annotation_PGSB.all.parts.gff3 Barley_Morex_V2_gene_annotation_PGSB.HC.parts.gff3 Barley_Morex_V2_gene_annotation_PGSB.LC.parts.gff3

rm Barley_Morex_V2_gene_annotation_PGSB.all.gff3 Barley_Morex_V2_gene_annotation_PGSB.HC.gff3 Barley_Morex_V2_gene_annotation_PGSB.LC.gff3
