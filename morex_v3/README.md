# Morex v3

This directory contains scripts and submodules for file processing/creation that are necessary for using the reference genome.

Morex v3 reference genome can be accessed through e!DAL - Plant Genomics & Phenomics Research Data Repository: http://dx.doi.org/10.5447/ipk/2021/3

Phytozome V13 has Morex v3 annotations (including repeat annotations) and assembly files available here: https://data.jgi.doe.gov/refine-download/phytozome?organism=HvulgareMorex&expanded=702

---

## Step 1: Prepare Reference

Run scripts in `prep_reference` subdirectory to get reference genome ready for downstream use. See readme in subdirectory for detailed documentation.

## Step 2: Get Morex v3 BOPA, 9K, and 50K positions

Run scripts in `50k_9k_BOPA_SNP` subdirectory. See readme in subdirectory for detailed documentation.

## Step 3: Find stretches of N's in reference and softmasked regions

Get a BED file of places in the reference where there are stretches of N's.

```bash
# In dir: ~/GitHub/morex_reference/morex_v3
module load python3/3.8.3_anaconda2020.07_mamba
# Find stretches of N's in reference genome
# partsRef
./find_Ns_in_assembly.py /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta > /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/stretches_of_Ns/Barley_MorexV3_pseudomolecules_parts_missing.bed

# pseudomolecules
./find_Ns_in_assembly.py /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta > /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/stretches_of_Ns/Barley_MorexV3_pseudomolecules_missing.bed
```

Generate BED file of hardmasked and softmasked regions in the Phytozome 13 Morex v3 reference assembly files.

```bash
# In dir: /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/assembly
module load python3/3.8.3_anaconda2020.07_mamba
# partsRef hardmasked
~/GitHub/morex_reference/morex_v3/find_Ns_in_assembly.py HvulgareMorex_702_V3.hardmasked_parts.fasta > HvulgareMorex_702_V3.hardmasked_parts.bed
# partsRef softmasked
~/GitHub/morex_reference/morex_v3/find_softmasked_in_assembly.py HvulgareMorex_702_V3.softmasked_parts.fasta > HvulgareMorex_702_V3.softmasked_parts.bed

# pseudomolecules hardmasked
~/GitHub/morex_reference/morex_v3/find_Ns_in_assembly.py HvulgareMorex_702_V3.hardmasked.fa > HvulgareMorex_702_V3.hardmasked.bed
# pseudomolecules softmasked
~/GitHub/morex_reference/morex_v3/find_softmasked_in_assembly.py HvulgareMorex_702_V3.softmasked.fa > HvulgareMorex_702_V3.softmasked.bed
```

Mask low complexity regions using [JGI's bbmask](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmask-guide/).

```bash
# In dir: ~/GitHub/morex_reference/morex_v3
sbatch mask_low_complexity.sh
```
