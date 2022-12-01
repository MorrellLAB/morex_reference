#!/bin/bash

set -e
set -o pipefail

# User provided input arguments
cent_pos="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/centromere_positions/MorexV3_centromere_positions.tsv"

ref_fai="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta.fai"

out_dir="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3"

out_prefix="parts"

#-----------
# Prepare array
pos_arr=()
for i in $(cat ${cent_pos} | sed -e "s, ,:,")
do
    pos_arr+=(${i})
done

# Generate file header
printf "#CHROM\tSTART\tEND\tNAME\n" > ${out_dir}/${out_prefix}.bed
# Add the split positions
for i in ${pos_arr[@]}
do
    curr_chr=$(echo ${i} | cut -d':' -f 1)
    if [ ${curr_chr} != "chrUn" ]; then
        curr_cent_pos=$(echo ${i} | cut -d':' -f 2)
        curr_chr_end=$(grep ${curr_chr} ${ref_fai} | cut -f 2)
        printf "${curr_chr}\t0\t${curr_cent_pos}\t${curr_chr}_part1\n" >> ${out_dir}/${out_prefix}.bed
        printf "${curr_chr}\t${curr_cent_pos}\t${curr_chr_end}\t${curr_chr}_part2\n" >> ${out_dir}/${out_prefix}.bed
    fi
done

# Add chrUn
curr_chr_end=$(grep chrUn ${ref_fai} | cut -f 2)
printf "chrUn\t0\t${curr_chr_end}\tchrUn\n" >> ${out_dir}/temp_${out_prefix}.bed

# Create parts.bed with new chromosome names in the first column
grep -v "#CHROM" ${out_dir}/temp_${out_prefix}.bed | awk '{ print $4"\t"$2"\t"$3 }' > ${out_dir}/${out_prefix}.bed
# Create version without chrUn
grep -v "chrUn" ${out_dir}/${out_prefix}.bed > ${out_dir}/${out_prefix}.nochrUn.bed
