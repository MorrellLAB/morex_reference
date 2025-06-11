# 50k, 9k, and BOPA SNPs Morex v3

This directory contains the scripts used to re-map the BOPA, 9k iSelect, and 50k iSelect SNPs to the Morex v3 reference genome.

### Navigation: Jump to Section

- [Data](#data)
- [Data Preparation](#data-preparation)
- [Methods: BOPA](#methods-bopa)
- [Methods: 9K](#methods-9k)
- [Methods: 50K](#methods-50k)

### Overview

Summary of the number of bopa, 9K, and 50K SNPs that had unique alignments with a percent identify of 95%:

Files below are in the directory:
```bash
/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k
```

| BOPA, 9K, and 50K VCF | Number of SNPs |
| --------------------- | -------------- |
| `50k_idt95_noRescuedSNPs_partsRef.vcf` | 41,813 SNPs |
| `50k_idt90_noRescuedSNPs_partsRef.vcf` | 42,241 SNPs |
| `9k_idt95_noRescuedSNPs_partsRef.vcf` | 7,541 SNPs |
| `9k_idt90_noRescuedSNPs_partsRef.vcf` | 7,622 SNPs |
| `bopa_idt95_noRescuedSNPs_partsRef.vcf` | 2,975 SNPs |
| `bopa_idt90_noRescuedSNPs_partsRef.vcf` | 3,006 SNPs |

---

## Data

The barley 50k iSelect SNP array data is published in [Bayer et al. 2017 Frontiers Plant Science](https://doi.org/10.3389/fpls.2017.01792). The 50k data can be [downloaded here](https://ics.hutton.ac.uk/50k/index.pl).

All contextual fasta sequences for BOPA, 9k, and 50k are located in:

```bash
/panfs/roc/groups/9/morrellp/shared/Datasets/Genotyping/Contextual_Sequences
```

## Data Preparation

To use the BLAST approach, we need to make a BLAST database first. We need to download the taxonomy information otherwise we will get errors when checking our BLAST database.

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k
wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
# Extract files
tar -xzvf taxdb.tar.gz
```

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k
# Using ncbi-blast-2.9.0+
# Pseudomolecules reference
makeblastdb -in ../Barley_MorexV3_pseudomolecules.fasta -dbtype nucl -parse_seqids
# Check integrity of blast database
blastdbcheck -db ../Barley_MorexV3_pseudomolecules.fasta

# Parts reference
makeblastdb -in ../Barley_MorexV3_pseudomolecules_parts.fasta -dbtype nucl -parse_seqids
# Check integrity of blast database
blastdbcheck -db ../Barley_MorexV3_pseudomolecules_parts.fasta
```

---

### Methods: BOPA

Steps to identify the physical positions of the BOPA SNPs relative to the Morex v3 reference genome.

**Step 1:** Run SNP_Utils (https://github.com/mojaveazure/SNP_Utils)

Previously when we ran this process for Morex v2 we found that when using the 95 identity threshold, there were some overlapping SNPs between the BOPA and 9k sets where the BLAST results did not match and resulted in discordant positions. So, here we'll run both 95 and 90 identity thresholds then do our comparisons with the 9k before deciding which final set to use.

We'll also prepare a VCF that removes the duplicate alignments to prepare for downstream steps.

```bash
# In dir: ~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP/scripts
bash run_snp_utils-bopa_idt90.sh
bash run_snp_utils-bopa_idt95.sh
```

Summarizing SNPs:

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k
# BOPA idt95
~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP/scripts/summarize_snp_count.sh bopa_morex_v3_idt95.vcf bopa_morex_v3_idt95_failed.log 
Total SNPs: 2981
Total unique SNPs: 2978
Total duplicate SNPs: 3
Total failed SNPs: 52

# BOPA idt90
~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP/scripts/summarize_snp_count.sh bopa_morex_v3_idt90.vcf bopa_morex_v3_idt90_failed.log
Total SNPs: 3014
Total unique SNPs: 3010
Total duplicate SNPs: 4
Total failed SNPs: 19
```

Add chromosome lengths to header lines so that the VCF works with GATK's Variant Recalibrator.

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k
# BOPA parts ref
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts/add_contig_length_to_header.sh bopa_idt95_noRescuedSNPs_partsRef.vcf ~/Shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.dict
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts/add_contig_length_to_header.sh bopa_idt90_noRescuedSNPs_partsRef.vcf ~/Shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.dict
# BOPA pseudomolecules ref
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts/add_contig_length_to_header.sh bopa_idt95_noRescuedSNPs.vcf ~/Shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.dict
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts/add_contig_length_to_header.sh bopa_idt90_noRescuedSNPs.vcf ~/Shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.dict
```

Rename files so output name is clean:

```bash
mv bopa_idt90_noRescuedSNPs_fixedHeader.vcf bopa_idt90_noRescuedSNPs.vcf
mv bopa_idt90_noRescuedSNPs_partsRef_fixedHeader.vcf bopa_idt90_noRescuedSNPs_partsRef.vcf

mv bopa_idt95_noRescuedSNPs_fixedHeader.vcf bopa_idt95_noRescuedSNPs.vcf
mv bopa_idt95_noRescuedSNPs_partsRef_fixedHeader.vcf bopa_idt95_noRescuedSNPs_partsRef.vcf
```

---

### Methods: 9K

Steps to identify the physical positions of the 9K SNPs relative to the Morex v3 reference genome. We follow the same steps as the BOPA SNPs, so we won't duplicate detailed documentation that is present there.

**Step 1:** Run SNP_Utils (https://github.com/mojaveazure/SNP_Utils)

We'll also prepare a VCF that removes the duplicate alignments to prepare for downstream steps.

```bash
# In dir: ~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP/scripts
bash run_snp_utils-9K_idt95.sh
bash run_snp_utils-9K_idt90.sh
```

Summarizing SNPs:

```bash
# 9k idt90
~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP/scripts/summarize_snp_count.sh 9k_morex_v3_idt90.vcf 9k_morex_v3_idt90_failed.log 
Total SNPs: 7763
Total unique SNPs: 7679
Total duplicate SNPs: 57
Total failed SNPs: 117

# 9k idt95
~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP/scripts/summarize_snp_count.sh 9k_morex_v3_idt95.vcf 9k_morex_v3_idt95_failed.log 
Total SNPs: 7674
Total unique SNPs: 7594
Total duplicate SNPs: 53
Total failed SNPs: 205
```

Add chromosome lengths to header lines so that the VCF works with GATK's Variant Recalibrator.

```bash
# 9K parts ref
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts/add_contig_length_to_header.sh 9k_idt95_noRescuedSNPs_partsRef.vcf ~/Shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.dict
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts/add_contig_length_to_header.sh 9k_idt90_noRescuedSNPs_partsRef.vcf ~/Shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.dict
# 9K pseudomolecules ref
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts/add_contig_length_to_header.sh 9k_idt95_noRescuedSNPs.vcf ~/Shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.dict
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts/add_contig_length_to_header.sh 9k_idt90_noRescuedSNPs.vcf ~/Shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.dict
```

Rename files so output name is clean:

```bash
mv 9k_idt95_noRescuedSNPs_partsRef_fixedHeader.vcf 9k_idt95_noRescuedSNPs_partsRef.vcf
mv 9k_idt95_noRescuedSNPs_fixedHeader.vcf 9k_idt95_noRescuedSNPs.vcf

mv 9k_idt90_noRescuedSNPs_fixedHeader.vcf 9k_idt90_noRescuedSNPs.vcf
mv 9k_idt90_noRescuedSNPs_partsRef_fixedHeader.vcf 9k_idt90_noRescuedSNPs_partsRef.vcf
```

---

### Check overlapping BOPA and 9k iSelect SNPs

Check for BOPA and 9k position concordance. For now, we'll check the `*_noRescuedSNPs_partsRef.vcf` VCF files so we can move forward with some analyses. Later on, one duplicates have been resolved and failed SNPs have been rescued, we'll re-run this process on those SNPs.

**Check 1.** BOPA and 9k used 95idt threshold. Make sure all SNP positions for SNPs that overlap in the two sets are concordant.

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k
module load python3/3.8.3_anaconda2020.07_mamba

~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts/check_position_concordance.py bopa_idt95_noRescuedSNPs_partsRef.vcf 9k_idt95_noRescuedSNPs_partsRef.vcf > temp_discordant_bopa_9k_snps.txt

~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts/check_position_concordance.py bopa_idt95_noRescuedSNPs_partsRef.vcf 9k_idt90_noRescuedSNPs_partsRef.vcf > temp_discordant_bopa_idt95_and_9k_idt90_snps.txt

cat temp_discordant_bopa_idt95_and_9k_idt90_snps.txt
# Output
9k	 chr2H_part2	319699727	11_10092	T	C	.	.	B
bopa	 chr5H_part2	321478839	11_10092	T	C	.	.	B

# BOPA and 50k BOPA SNPs
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts/check_position_concordance.py bopa_idt95_noRescuedSNPs_partsRef.vcf 50k_idt95_noRescuedSNPs_partsRef_cleanNamesBOPA.vcf > temp_discordant_bopa_idt95_and_50k_idt95_snps.txt

~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts/check_position_concordance.py bopa_idt95_noRescuedSNPs_partsRef.vcf 50k_idt90_noRescuedSNPs_partsRef_cleanNamesBOPA.vcf > temp_discordant_bopa_idt95_and_50k_idt90_snps.txt
```

BOPA idt95 vs 9k idt90 have one discordant SNP. Remaining don't have discordant BOPA SNPs. For the discordant SNPs in `temp_discordant_bopa_idt95_and_9k_idt90_snps.txt`, we manually blasted the contextual sequence using IPK's BLAST server to investigate.

---

### Methods: 50K

**Step 1:** Run SNP_Utils (https://github.com/mojaveazure/SNP_Utils)

We'll also prepare a VCF that removes the duplicate alignments to prepare for downstream steps.

```bash
# In dir: ~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP/scripts
bash run_snp_utils-50K_idt95.sh
bash run_snp_utils-50K_idt90.sh
```

Summarizing SNPs:

```bash
# 50k idt90
~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP/scripts/summarize_snp_count.sh 50k_morex_v3_idt90.vcf 50k_morex_v3_idt90_failed.log 
Total SNPs: 44801
Total unique SNPs: 43942
Total duplicate SNPs: 617
Total failed SNPs: 89

# 50k idt95
~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP/scripts/summarize_snp_count.sh 50k_morex_v3_idt95.vcf 50k_morex_v3_idt95_failed.log 
Total SNPs: 44333
Total unique SNPs: 43488
Total duplicate SNPs: 602
Total failed SNPs: 548
```

Add chromosome lengths to header lines so that the VCF works with GATK's Variant Recalibrator.

```bash
# idt90
# 50K parts ref
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts/add_contig_length_to_header.sh 50k_idt90_noRescuedSNPs_partsRef_cleanNamesBOPA.vcf ~/Shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.dict
# 50K pseudomolecules ref
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts/add_contig_length_to_header.sh 50k_idt90_noRescuedSNPs.vcf ~/Shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.dict
```

Rename files so output name is clean:

```bash
mv 50k_idt90_noRescuedSNPs_fixedHeader.vcf 50k_idt90_noRescuedSNPs.vcf
mv 50k_idt90_noRescuedSNPs_partsRef_cleanNamesBOPA_fixedHeader.vcf 50k_idt90_noRescuedSNPs_partsRef.vcf
rm 50k_idt90_noRescuedSNPs_partsRef_cleanNamesBOPA.vcf

mv 50k_idt95_noRescuedSNPs_partsRef_fixedHeader.vcf 50k_idt95_noRescuedSNPs_partsRef.vcf
mv 50k_idt95_noRescuedSNPs_fixedHeader.vcf 50k_idt95_noRescuedSNPs.vcf
```

Add a merged file with all of the SNPs from genotyping platforms

```bash
bcftools merge bopa_idt90_noRescuedSNPs.vcf.gz 9k_idt90_noRescuedSNPs.vcf.gz 50k_idt90_noRescuedSNPs.vcf.gz -Oz -o 50k_9k_bopa_idt90_noRescuedSNPs.vcf.gz
```
Total SNPs: 43833