# Prepare reference genome

### Index reference

```bash
sbatch make_index_pseudo_bwa.sh
```

### Break chromosomes into parts by chromosome arm

Create `.fai` index using samtools.

```bash
# In dir: /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3
module load samtools/1.10
# Index fasta file
samtools faidx Barley_MorexV3_pseudomolecules.fasta
```

Make parts bed file "manually" using centromere positions at:

```bash
# In dir: ~/GitHub/morex_reference/morex_v3/prep_reference
./make_parts_bed_file.sh
```

Break up pseudomolecules reference into chromosome arms.

```bash
# In dir: /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3
# Split reference into parts
module load bedtools/2.29.2
# Longer naming scheme. Example: >chr1H_part1::chr1H:0-206486643
# This is due to a change in the -name flag for v2.29.2 compared to v2.28.0
bedtools getfasta -name -fi Barley_MorexV3_pseudomolecules.fasta \
    -fo Barley_MorexV3_pseudomolecules_parts_longSeqID.fasta \
    -bed parts.bed

# Index parts reference
module load samtools/1.10
samtools faidx Barley_MorexV3_pseudomolecules_parts_longSeqID.fasta

# BWA indexing
sbatch make_index_parts_bwa_longSeqID.sh

#--------------------------
# Split reference into parts, using shorter naming scheme.
#   Example: >chr1H_part1
bedtools getfasta -nameOnly -fi Barley_MorexV3_pseudomolecules.fasta \
    -fo Barley_MorexV3_pseudomolecules_parts.fasta \
    -bed parts.bed

# Index parts reference
module load samtools/1.10
samtools faidx Barley_MorexV3_pseudomolecules_parts.fasta

# BWA indexing
sbatch make_index_parts_bwa.sh
```
