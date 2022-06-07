#!/bin/bash

set -e
set -o pipefail

# User provided input arguments
cent_pos="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/centromere_positions/MorexV3_centromere_positions.tsv"

ref_fai="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_plastids.fasta.fai"

out_dir="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3"

out_prefix="parts_plastids"

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
    if [ ${curr_chr} != "chrUn" ] || [ ${curr_chr} != "EF115541.1" ] || [ ${curr_chr} != "AP017301.1" ]; then
        curr_cent_pos=$(echo ${i} | cut -d':' -f 2)
        curr_chr_end=$(grep ${curr_chr} ${ref_fai} | cut -f 2)
        printf "${curr_chr}\t0\t${curr_cent_pos}\t${curr_chr}_part1\n" >> ${out_dir}/${out_prefix}.bed
        printf "${curr_chr}\t${curr_cent_pos}\t${curr_chr_end}\t${curr_chr}_part2\n" >> ${out_dir}/${out_prefix}.bed
    fi
done

# Add chrUn
curr_chr_end=$(grep chrUn ${ref_fai} | cut -f 2)
printf "chrUn\t0\t${curr_chr_end}\tchrUn\n" >> ${out_dir}/${out_prefix}.bed

# Add plastids
curr_p1_end=$(grep "EF115541.1" ${ref_fai} | cut -f 2)
printf "EF115541.1\t0\t${curr_p1_end}\tEF115541.1\n" >> ${out_dir}/${out_prefix}.bed
curr_p2_end=$(grep "AP017301.1" ${ref_fai} | cut -f 2)
printf "AP017301.1\t0\t${curr_p2_end}\tAP017301.1\n" >> ${out_dir}/${out_prefix}.bed
