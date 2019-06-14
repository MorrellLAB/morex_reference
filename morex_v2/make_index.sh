#!/bin/bash
#PBS -l mem=50gb,nodes=1:ppn=4,walltime=6:00:00
#PBS -m abe
#PBS -M wyant008@umn.edu
#PBS -q mesabi

cd /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2

module load minimap2_ML/2.17.0

module load bwa/0.7.17

minimap2 -x map-ont -d test.mmi Barley_Morex_V2_pseudomolecules_plastids.fasta

bwa index /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules_parts_plastids.fasta




