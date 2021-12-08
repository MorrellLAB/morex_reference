# 50k, 9k, and BOPA SNPs Morex v3

This directory contains the scripts used to re-map the BOPA, 9k iSelect, and 50k iSelect SNPs to the Morex v3 reference genome.

### Navigation: Jump to Section

- [Data](#data)
- [Data Preparation](#data-preparation)
- [Methods: BOPA](#methods-bopa)
- [Methods: 9K](#methods-9k)
- [Methods: 50K](#methods-50k)

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

First, we'll run with **idt95**.

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k
# Prepare dir for intermediate files that can be deleted later
mkdir Intermediates

# Dependencies
module load python3/3.6.3_anaconda5.0.1

# Define shared variables for files
LOOKUP_TABLE="/panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/lookup_tables/BOPA_lookup.txt"
GENETIC_MAP="/panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/genetic_maps/BOPA/2013_iSelect_Genetic_Map.txt"

# Define variable specific to idt95
OUT_PREFIX="BOPA_morex_v3_idt95"

# Run SNP_Utils for idt95
~/software_development/SNP_Utils/snp_utils.py CONFIG -d ../Barley_MorexV3_pseudomolecules.fasta -k -i 95 -c blast_bopa_morex_v3_idt95

# Run SNP Utils BLAST for idt95
~/software_development/SNP_Utils/snp_utils.py BLAST -l ${LOOKUP_TABLE} -c blast_bopa_morex_v3_idt95 -b -m ${GENETIC_MAP} -d -t 100000 -o ${OUT_PREFIX}
# Relevant part of output messages
    Filtering SNPs by hit chromsome/contig
    Using genetic map /panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/genetic_maps/BOPA/2013_iSelect_Genetic_Map.txt
    Filtering SNPs with a minimum distance threshold of 100000
    Filtering SNPs by relative location on the genetic map
    Using genetic map /panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/genetic_maps/BOPA/2013_iSelect_Genetic_Map.txt
    Parsing blast database ../Barley_MorexV3_pseudomolecules.fasta
    Found 8 chromosomes
    Filtering 2978 SNP IDs
    No genetic map for chromosome chr7H
    No genetic map for chromosome chr5H
    No genetic map for chromosome chr2H
    Writing 2981 SNPs to BOPA_morex_v3_idt95.vcf
    Removing masked SNPs that were actually found
    Writing 2 masked SNPs to BOPA_morex_v3_idt95_masked.vcf
    Writing 52 failed SNPs to BOPA_morex_v3_idt95_failed.log
```

Summarizing SNPs for **idt95**:

```bash
# Total SNPs
grep -v "#" BOPA_morex_v3_idt95.vcf | cut -f 3 | wc -l
    2981
# Total unique SNPs
grep -v "#" BOPA_morex_v3_idt95.vcf | cut -f 3 | sort -u | wc -l
    2978
# Total number of duplicates
grep -v "#" BOPA_morex_v3_idt95.vcf | cut -f 3 | sort | uniq -c | sort -n -r | grep -vw "1" | wc -l
    3
```

Next, we'll run with idt90. We'll re-use the dependencies and shared variables defined above.

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k
# Define variable specific to idt90
OUT_PREFIX_idt90="BOPA_morex_v3_idt90"

# Run SNP_Utils for idt90
~/software_development/SNP_Utils/snp_utils.py CONFIG -d ../Barley_MorexV3_pseudomolecules.fasta -k -i 90 -c blast_bopa_morex_v3_idt90

# Run SNP Utils BLAST for idt90
~/software_development/SNP_Utils/snp_utils.py BLAST -l ${LOOKUP_TABLE} -c blast_bopa_morex_v3_idt90 -b -m ${GENETIC_MAP} -d -t 100000 -o ${OUT_PREFIX_idt90}
# Relevant part of output messages
    Filtering SNPs by hit chromsome/contig
    Using genetic map /panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/genetic_maps/BOPA/2013_iSelect_Genetic_Map.txt
    Filtering SNPs with a minimum distance threshold of 100000
    Filtering SNPs by relative location on the genetic map
    Using genetic map /panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/genetic_maps/BOPA/2013_iSelect_Genetic_Map.txt
    Parsing blast database ../Barley_MorexV3_pseudomolecules.fasta
    Found 8 chromosomes
    Filtering 3011 SNP IDs
    No genetic map for chromosome chr7H
    No genetic map for chromosome chr5H
    No genetic map for chromosome chr2H
    Writing 3014 SNPs to BOPA_morex_v3_idt90.vcf
    Removing masked SNPs that were actually found
    Writing 3 masked SNPs to BOPA_morex_v3_idt90_masked.vcf
    Writing 19 failed SNPs to BOPA_morex_v3_idt90_failed.log
```

Summarizing SNPs for **idt90**:

```bash
# Total SNPs
grep -v "#" BOPA_morex_v3_idt90.vcf | cut -f 3 | wc -l
    3014
# Total unique SNPs
grep -v "#" BOPA_morex_v3_idt90.vcf | cut -f 3 | sort -u | wc -l
    3010
# Total number of duplicates
grep -v "#" BOPA_morex_v3_idt90.vcf | cut -f 3 | sort | uniq -c | sort -n -r | grep -vw "1" | wc -l
    4
```

After running `snp_utils.py` for **idt95**, we get:

- 2,978 unique SNPs
- 3 SNPs with duplicates
- 52 failed SNPs
- 2 masked SNPs

For **idt90**, we get:

- 3,010 unique SNPs
- 4 SNPs with duplicates
- 19 failed SNPs

For the ones that failed, we will "manually" BLAST using the IPK server: https://webblast.ipk-gatersleben.de/barley_ibsc/. Previously (for Morex v1), Li had many undergrads help us manually BLAST search the duplicate and failed SNPs. But, this is somewhat more prone to human error and is not possible without a many undergrads. So, we will use the script below to resolve the duplicates and failed SNPs by parsing an HTML file containg the IPK Barley BLAST server results. (*Note:* this is the same script we used when processing positions relative to Morex v2)

**Step 2:** Resolve duplicate SNPs

Let's resolve the duplicate SNPs, since there are only a few, we can do this manually.

List of BOPA SNPs with duplicates:

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k
grep -v "#" BOPA_morex_v3_idt90.vcf | cut -f 3 | sort | uniq -c | sort -n -r | grep -vw "1"
      2 12_31396
      2 12_30861
      2 12_30822
      2 11_20074
```

Make a list of duplicate SNPs and remove duplicates from main VCFs:

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k
# idt90
~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP/remove_duplicates_from_vcf.sh BOPA_morex_v3_idt90.vcf bopa_idt90 ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k

# idt95
~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP/remove_duplicates_from_vcf.sh BOPA_morex_v3_idt95.vcf bopa_idt95 ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k
```

The general steps are as follows: Check if SNP hits a gene by searching in the [SNPMeta output file](https://conservancy.umn.edu/bitstream/handle/11299/181367/Barley_SNP_Annotations.txt?sequence=5&isAllowed=y). This step is intended to reduce noise since we expect all the BOPA and 9k SNPs to only be in genic regions. If there is a gene hit, BLAST search the gene with IPK server, pick the best hit, and write notes in the INFO field of the VCF. If there is no hit, choose the smallest position (smallest chromosome and smallest physical position) and write notes in the INFO field of the VCF file.

We will list ALTCHR first before ALTPOS since that is the format that was already present in the current VCF file we are working with. Example:
> ALTCHR=chr7H,chr7H;ALTPOS=6011918,6045065;B

We will repeat this process for the remaining duplicates. 

**Step 2a:** First search the [SNP meta file for the 9k SNPs](https://conservancy.umn.edu/bitstream/handle/11299/181367/Barley_SNP_Annotations.txt?sequence=5&isAllowed=y) and figure out which gene the SNP is located.

We'll do this programmatically by downloading the SNPMeta output file containing info for SNPs hitting genes and generate a list of protein IDs associated with the SNP.

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k/temp
# Download SNPMeta output file
wget https://conservancy.umn.edu/bitstream/handle/11299/181367/Barley_SNP_Annotations.txt?sequence=5&isAllowed=y
# Cleanup output file name
mv Barley_SNP_Annotations.txt\?sequence\=5 Barley_SNP_Annotations.txt

# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k/duplicate_snps
# Generate list of protein IDs associated with the SNP
grep -f bopa_idt90_duplicate_snp_names.txt ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k/temp/Barley_SNP_Annotations.txt | cut -f 1,4
11_20074	BAJ97892.1
12_30822	BAK06853.1
12_30861	ABV59386.1
12_31396	-

grep -f bopa_idt90_duplicate_snp_names.txt ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k/temp/Barley_SNP_Annotations.txt | cut -f 4 | grep -v "-" > bopa_idt90_duplicate_snps_proteinIDs.txt
```

The purpose for this is to reduce noise since we know that the 9k iSelect genotyping was designed from exome capture samples.

If we find the SNP in the SNP meta file, copy the ProteinID for that SNP and go to NCBI https://www.ncbi.nlm.nih.gov/protein and search that protein. This will take you to a page where you can get the fasta sequence for that ProteinID. Due to changes in how IPK's BLAST server is set up, we'll need to generate a single file containing the FASTA sequences from the NCBI protein database. Copy the FASTA sequence into a the file: `~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k/duplicate_snps/bopa_idt90_duplicate_snp_proteinIDs.fasta`

**Step 2b:** Take the fasta sequence and use the [IPK Barley BLAST server](https://webblast.ipk-gatersleben.de/barley_ibsc/) to do a BLAST search.

Navigate to "Blast barley genomes" > "NBCI BLAST+ blastp". Then select the following:

- Protein query sequence(s): Select `Single dataset` and then upload fasta file `bopa_idt90_duplicate_snp_proteinIDs.fasta`
- Subject database/sequences: Select `Locally installed BLAST database` then under the "Protein BLAST database" section select `Barley all Proteins Morex V3 (Jul 2020)`
- Type of BLAST: `blastp`
- Set expectation value cutoff: `0.001`
- Output format: `Tabular (extended 25 columns)`

Resolved duplicates are in the file `bopa_duplicates_idt90_resolved.vcf` and will be combined with the file `temp_bopa_idt90_noDups.vcf` and failed SNPs (see step 3) later on and sorted.

---

### Methods: 9K

Steps to identify the physical positions of the 9K SNPs relative to the Morex v3 reference genome. We follow the same steps as the BOPA SNPs, so we won't duplicate detailed documentation that is present there.

**Step 1:** Run SNP_Utils (https://github.com/mojaveazure/SNP_Utils)

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k
# Dependencies
module load python3/3.6.3_anaconda5.0.1

# Define shared variables for files
LOOKUP_TABLE="/panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/lookup_tables/9k_lookup.txt"
GENETIC_MAP="/panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/genetic_maps/BOPA/2013_iSelect_Genetic_Map.txt"

# Define variable specific to idt95
OUT_PREFIX="9k_morex_v3_idt95"

# Run SNP_Utils for idt95
~/software_development/SNP_Utils/snp_utils.py CONFIG -d ../Barley_MorexV3_pseudomolecules.fasta -k -i 95 -c blast_${OUT_PREFIX}

# Run SNP Utils BLAST for idt95
~/software_development/SNP_Utils/snp_utils.py BLAST -l ${LOOKUP_TABLE} -c blast_${OUT_PREFIX} -b -m ${GENETIC_MAP} -d -t 100000 -o ${OUT_PREFIX}
# Relevant part of output messages
Using genetic map /panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/genetic_maps/BOPA/2013_iSelect_Genetic_Map.txt
Filtering SNPs with a minimum distance threshold of 100000
Filtering SNPs by relative location on the genetic map
Using genetic map /panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/genetic_maps/BOPA/2013_iSelect_Genetic_Map.txt
Parsing blast database ../Barley_MorexV3_pseudomolecules.fasta
Found 8 chromosomes
Filtering 2978 SNP IDs
No genetic map for chromosome chr7H
No genetic map for chromosome chr5H
No genetic map for chromosome chr2H
Writing 2981 SNPs to 9k_morex_v3_idt95.vcf
Removing masked SNPs that were actually found
Writing 2 masked SNPs to 9k_morex_v3_idt95_masked.vcf
Writing 52 failed SNPs to 9k_morex_v3_idt95_failed.log
```

Summarizing SNPs for **idt95**:

```bash
# Total SNPs
grep -v "#" 9k_morex_v3_idt95.vcf | cut -f 3 | wc -l
    7674
# Total unique SNPs
grep -v "#" 9k_morex_v3_idt95.vcf | cut -f 3 | sort -u | wc -l
    7594
# Total number of duplicates
grep -v "#" 9k_morex_v3_idt95.vcf | cut -f 3 | sort | uniq -c | sort -n -r | grep -vw "1" | wc -l
    53
```

After running `snp_utils.py`, we get:
7,594 SNPs without duplicates
53 SNPs with duplicates
205 failed SNPs

**Step 2:** Resolve duplicate SNPs

Make a list of duplicate SNPs and remove duplicates from main VCFs:

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k
# idt95
~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP/remove_duplicates_from_vcf.sh 9k_morex_v3_idt95.vcf 9k_idt95 ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k

# Check counts
grep -v "#" 9k_idt95_noRescuedSNPs.vcf | wc -l
    7541
```

---

### Check overlapping BOPA and 9k iSelect SNPs

Check for BOPA and 9k position concordance. For now, we'll check the `*_noRescuedSNPs_partsRef.vcf` VCF files so we can move forward with some analyses. Later on, one duplicates have been resolved and failed SNPs have been rescued, we'll re-run this process on those SNPs.

**Check 1.** BOPA and 9k used 95idt threshold. Make sure all SNP positions for SNPs that overlap in the two sets are concordant.

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k
module load python3/3.8.3_anaconda2020.07_mamba(default)

~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts/check_position_concordance.py bopa_idt95_noRescuedSNPs_partsRef.vcf 9k_idt95_noRescuedSNPs_partsRef.vcf > temp_discordant_bopa_9k_snps.txt
```

There are no discordant SNP positions.

Add chromosome lengths to header lines so that the VCF works with GATK's Variant Recalibrator.

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k
# BOPA parts ref
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts/add_contig_length_to_header.sh bopa_idt95_noRescuedSNPs_partsRef.vcf ~/Shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.dict
# BOPA pseudomolecules ref
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts/add_contig_length_to_header.sh bopa_idt95_noRescuedSNPs.vcf ~/Shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.dict

# 9K parts ref
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts/add_contig_length_to_header.sh 9k_idt95_noRescuedSNPs_partsRef.vcf ~/Shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.dict
# 9K pseudomolecules ref
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts/add_contig_length_to_header.sh 9k_idt95_noRescuedSNPs.vcf ~/Shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.dict

# 50K parts ref
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts/add_contig_length_to_header.sh 50k_idt95_noRescuedSNPs_partsRef.vcf ~/Shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.dict
# 50K pseudomolecules ref
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts/add_contig_length_to_header.sh 50k_idt95_noRescuedSNPs.vcf ~/Shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.dict
```

---

### Methods: 50K

**Step 1:** Run SNP_Utils (https://github.com/mojaveazure/SNP_Utils)

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k
# Dependencies
module load python3/3.6.3_anaconda5.0.1

# Define shared variables for files
LOOKUP_TABLE="/panfs/roc/groups/9/morrellp/shared/Datasets/Genotyping/Contextual_Sequences/barley_50k/snp_utils_50k_lookup_table.txt"
GENETIC_MAP="/panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/genetic_maps/50k_iSelect_Genetic_Map.txt"

# Define variable specific to idt95
OUT_PREFIX="50k_morex_v3_idt95"

# Run SNP_Utils for idt95
~/software_development/SNP_Utils/snp_utils.py CONFIG -d ../Barley_MorexV3_pseudomolecules.fasta -k -i 95 -c blast_${OUT_PREFIX}

# Run SNP Utils BLAST for idt95
~/software_development/SNP_Utils/snp_utils.py BLAST -l ${LOOKUP_TABLE} -c blast_${OUT_PREFIX} -b -m ${GENETIC_MAP} -d -t 100000 -o ${OUT_PREFIX}
```

Summarizing SNPs:

```bash
# Total number of SNPs
grep -v "#" 50k_morex_v3_idt95.vcf | cut -f 3 | wc -l
    44333
# Total number of unique SNPs
grep -v "#" 50k_morex_v3_idt95.vcf | cut -f 3 | sort -u | wc -l
    43488
# Total number of duplicates
grep -v "#" 50k_morex_v3_idt95.vcf | cut -f 3 | sort | uniq -c | sort -n -r | grep -vw "1" | wc -l
    602
```

After running `snp_utils.py`, we get:
43,488 Total SNPs (unique)
602 SNPs with duplicates
548 failed SNPs

**Step 2:** Resolve duplicate SNPs

Make a list of duplicate SNPs and remove duplicates from main VCFs:

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k
# idt95
~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP/remove_duplicates_from_vcf.sh 50k_morex_v3_idt95.vcf 50k_idt95 ~/Shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k

# Check counts
grep -v "#" 50k_idt95_noRescuedSNPs.vcf | wc -l
    41813
```
