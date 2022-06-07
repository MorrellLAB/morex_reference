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
# Load dependencies
module load bedtools/2.29.2
module load samtools/1.10
# Split reference into parts
# Longer naming scheme. Example: >chr1H_part1::chr1H:0-206486643
# This is due to a change in the -name flag for v2.29.2 compared to v2.28.0
bedtools getfasta -name -fi Barley_MorexV3_pseudomolecules.fasta \
    -fo Barley_MorexV3_pseudomolecules_parts_longSeqID.fasta \
    -bed parts.bed
# Index parts reference
samtools faidx Barley_MorexV3_pseudomolecules_parts_longSeqID.fasta
# BWA indexing
cd ~/GitHub/morex_reference/morex_v3/prep_reference
sbatch make_index_parts_bwa_longSeqID.sh

# Split reference into parts, using shorter naming scheme.
#   Example: >chr1H_part1
bedtools getfasta -nameOnly -fi Barley_MorexV3_pseudomolecules.fasta \
    -fo Barley_MorexV3_pseudomolecules_parts.fasta \
    -bed parts.bed
# Index parts reference
samtools faidx Barley_MorexV3_pseudomolecules_parts.fasta
# BWA indexing
cd ~/GitHub/morex_reference/morex_v3/prep_reference
sbatch make_index_parts_bwa.sh
```

#### Split plastids reference into parts

`.fai` index already exists for the `/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_plastids.fasta` file.

Generate `parts.bed` file using centromere positions to split:

```bash
# In dir: ~/GitHub/morex_reference/morex_v3/prep_reference
./make_parts_bed_file-plastids.sh
```

Break up pseudomolecules reference into chromosome arms

```bash
# In dir: /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3
# Load dependencies
module load bedtools/2.29.2
module load samtools/1.10
# Split reference into parts
# Longer naming scheme. Example: >chr1H_part1::chr1H:0-206486643
# This is due to a change in the -name flag for v2.29.2 compared to v2.28.0
bedtools getfasta -name -fi Barley_MorexV3_pseudomolecules_plastids.fasta \
    -fo Barley_MorexV3_pseudomolecules_plastids_parts_longSeqID.fasta \
    -bed parts_plastids.bed
# Index parts reference
samtools faidx Barley_MorexV3_pseudomolecules_plastids_parts_longSeqID.fasta
# BWA indexing
cd ~/GitHub/morex_reference/morex_v3/prep_reference
sbatch make_index_plastids_parts_bwa_longSeqID.sh

# Split reference into parts, using shorter naming scheme
bedtools getfasta -nameOnly -fi Barley_MorexV3_pseudomolecules_plastids.fasta \
    -fo Barley_MorexV3_pseudomolecules_plastids_parts.fasta \
    -bed parts_plastids.bed
# Index parts reference
samtools faidx Barley_MorexV3_pseudomolecules_plastids_parts.fasta
# BWA indexing
cd ~/GitHub/morex_reference/morex_v3/prep_reference
sbatch make_index_plastids_parts_bwa.sh
```

#### Split Phytozome reference assemblies into parts

This includes their softmasked and hardmasked assemblies.

Create `.fai` index using samtools.

```bash
module load samtools/1.10

# In dir: /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/assembly
samtools faidx HvulgareMorex_702_V3.fa
samtools faidx HvulgareMorex_702_V3.hardmasked.fa
samtools faidx HvulgareMorex_702_V3.softmasked.fa
```

We can use the `/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/parts.bed` generated above.

Break up pseudomolecules reference into chromosome parts.

```bash
# Load dependencies
module load bedtools/2.29.2
module load samtools/1.10

# In dir: /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/assembly
# Split reference into parts, using shorter naming scheme
bedtools getfasta -nameOnly -fi HvulgareMorex_702_V3.fa \
    -fo HvulgareMorex_702_V3_parts.fasta \
    -bed /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/parts.bed

# Split hardmasked reference
bedtools getfasta -nameOnly -fi HvulgareMorex_702_V3.hardmasked.fa \
    -fo HvulgareMorex_702_V3.hardmasked_parts.fasta \
    -bed /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/parts.bed

# Split softmasked reference
bedtools getfasta -nameOnly -fi HvulgareMorex_702_V3.softmasked.fa \
    -fo HvulgareMorex_702_V3.softmasked_parts.fasta \
    -bed /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/parts.bed

# Index parts references
samtools faidx HvulgareMorex_702_V3_parts.fasta
samtools faidx HvulgareMorex_702_V3.hardmasked_parts.fasta
samtools faidx HvulgareMorex_702_V3.softmasked_parts.fasta

# BWA indexing
cd ~/GitHub/morex_reference/morex_v3/prep_reference
sbatch make_index_phytozome_parts_bwa.sh
```
