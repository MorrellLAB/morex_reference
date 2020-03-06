# 50k, 9k, and BOPA SNPs

This directory contains the scripts used to re-map the BOPA, 9k iSelect, and 50k iSelect SNPs to the Morex v2 reference genome.

### Navigation: Jump to Section

- [Data](#data)
- [Data exploration](#data-exploration)
- [Methods: BOPA](#methods-bopa)
- [Methods: 9k iSelect](#methods-9k-iselect)
- [Methods: 50k iSelect](#methods-50k-iselect)

---

## Data

All contextual fasta sequences for BOPA, 9k, and 50k are located in:

```bash
/panfs/roc/groups/9/morrellp/shared/Datasets/Genotyping/Contextual_Sequences
```

## Data exploration

We ran some initial quick alignments to get an idea of the quality of the Morex v2 genome. The scripts used for this are located in the subdirectory `prep_reference`.

```bash
# Index reference
qsub make_index_pseudo_bowtie2.sh
qsub make_index_pseudo_bwa.sh
```

Do a quick alignment with default parameters.

```bash
qsub check_by_aligning_bowtie2_BOPA.sh
qsub check_by_aligning_bwa-mem_BOPA.sh
```

Summary from Bowtie2 alignment.

```bash
# Alignment summary located in the file:
# /panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/morex_v2/check_by_aligning_bowtie2_BOPA.sh.e15261636
3072 reads; of these:
  3072 (100.00%) were unpaired; of these:
    433 (14.10%) aligned 0 times
    2472 (80.47%) aligned exactly 1 time
    167 (5.44%) aligned >1 times
85.90% overall alignment rate
```

Get summary stats from BWA alignment.

```bash
module load samtools/1.9

samtools flagstat BOPA_contextual_morex_v2-bwa-mem.sam
# Output
3675 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
603 + 0 supplementary
0 + 0 duplicates
3666 + 0 mapped (99.76% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

Previously, Li only used Bowtie2 as an initial check but found that a lot of SNPs from their algorithm were not correct when using the physical positions to check. So, she decided to use the BLAST approach.

## Data Preparation

To use the BLAST approach, we need to make a BLAST database first. We need to download the taxonomy information otherwise we will get errors when checking our BLAST database.

```bash
wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
# Extract files
tar -xzvf taxdb.tar.gz
```

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v2
makeblastdb -in Barley_Morex_V2_pseudomolecules.fasta -dbtype nucl -parse_seqids
# Check integrity of blast database
blastdbcheck -db Barley_Morex_V2_pseudomolecules.fasta

# Parts reference
makeblastdb -in Barley_Morex_V2_pseudomolecules_parts_plastids.fasta -dbtype nucl -parse_seqids
# Check integrity of blast database
blastdbcheck -db Barley_Morex_V2_pseudomolecules_parts_plastids.fasta
```

Do a BLAST search for the BOPA SNPs against the reference database.

```bash
blastn -db Barley_Morex_V2_pseudomolecules.fasta -query /panfs/roc/groups/9/morrellp/shared/Datasets/Genotyping/Contextual_Sequences/BOPA_SNPs/BOPA_SNPs.fasta -out /panfs/roc/groups/9/morrellp/shared/Datasets/Blast_Alignments/BOPA_SNPs.out
```

For the genetic map data, do some cleanup.

```bash
# In dir: ~/GitHub/morex_reference/genetic_maps
# Sort by chromosome and SNP name
(head -n 1 GeneticMap_iSelect_9k.txt && tail -n +2 GeneticMap_iSelect_9k.txt | sed -e 's/^X//' | sort -V -k3,3 -k1,1) > 9k_iSelect_Genetic_Map.txt
git rm GeneticMap_iSelect_9k.txt
```

---

### Methods: BOPA

`BOPA_morex_v2_idt95.vcf` stores all of the BOPA SNPs mapped to the Morex v2 reference genome. It includes **xxx** SNPs out of 3,072 contextual sequences because we cannot figure out the position of the remaining **xx** SNPs due to deletions in the reference or the SNPs were not located in the matched region when we BLAST the contextual sequence against the Morex v2 reference genome. Below, I have detailed the steps we ran to produce the VCF files.

**Step 1:** Run SNP_Utils (https://github.com/mojaveazure/SNP_Utils)

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v2
 ~/software_development/SNP_Utils/snp_utils.py CONFIG -d Barley_Morex_V2_pseudomolecules.fasta -k -i 95 -c blast_morex_v2_idt95

# Define variables for files
LOOKUP_TABLE=/panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/lookup_tables/BOPA_lookup.txt
GENETIC_MAP=/panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/genetic_maps/BOPA/2013_iSelect_Genetic_Map.txt
OUT_PREFIX=BOPA_morex_v2_idt95

# Run SNP Utils BLAST
~/software_development/SNP_Utils/snp_utils.py BLAST -l ${LOOKUP_TABLE} -c blast_morex_v2_idt95 -b -m ${GENETIC_MAP} -d -t 100000 -o ${OUT_PREFIX}
```

Summarizing SNPs:

```bash
# Total number of SNPs
grep -v "#" BOPA_morex_v2_idt95.vcf | cut -f 3 | wc -l
    2992

# Total number of unique SNPs
grep -v "#" BOPA_morex_v2_idt95.vcf | cut -f 3 | sort -u | wc -l
    2985

# Total number of duplicates
grep -v "#" BOPA_morex_v2_idt95.vcf | cut -f 3 | sort | uniq -c | sort -n -r | head
      3 12_31124
      3 12_30822
      2 12_30861
      2 12_21135
      2 11_20074
      1 12_31536
      1 12_31535
      1 12_31529
      1 12_31528
      1 12_31527
```

After running `snp_utils.py`, we get:
2,985 SNPs without duplicates
5 SNPs with duplicates
52 failed SNPs

For the ones that failed, we "manually" BLAST using the IPK server: https://webblast.ipk-gatersleben.de/barley_ibsc/. Previously (for Morex v1), Li had many undergrads help us manually BLAST search the duplicate and failed SNPs. But, this is somewhat more prone to human error and is not possible without a many undergrads. So, we will use the script below to resolve the duplicates and failed SNPs by parsing an HTML file containg the IPK Barley BLAST server results.

**Step 2:** Resolve duplicate SNPs

Let's resolve the duplicate SNPs, since there are only 5, we can do this manually.

List of 5 BOPA SNPs with duplicates:
- 12_31124
- 12_30822
- 12_30861
- 12_21135
- 11_20074

SNP 12_31124:

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v2
# Duplicates of this SNP
grep 12_31124 BOPA_morex_v2_idt95.vcf
chr1H    37719516    12_31124    C    T    .    .    B
chrUn    56705296    12_31124    C    T    .    .    B
chrUn    75911901    12_31124    G    A    .    .    B

# In dir: ~/Shared/Datasets/Genotyping/Contextual_Sequences/BOPA_SNPs
# Identify line number of duplicate SNP in FASTA file
grep -n 12_31124 BOPA_SNPs.fasta
5637:>12_31124

# View contextual sequence associated with that SNP
sed -n '5637,5638p' BOPA_SNPs.fasta
>12_31124
CTATTAGAGTCTTGTATATGTATTATATCATAGACAAAGCACACGAAATGATGTCCAGATYATTCTTCTTCTTCATCAGTCCACACGAGAGGTTTAAATTGTATATGTAAATCCAGAATTC
```

Check if SNP hits a gene by searching in the SNPMeta output file `~/GitHub/morex_reference/morex_v2/data/bopa_and_9k_snpmeta_output.txt`. his step is intended to reduce noise since we expect all the BOPA and 9k SNPs to only be in genic regions. If there is a gene hit, BLAST search the gene with IPK server using default settings. If there is a gene hit, BLAST search the gene with IPK server using default settings. (**Note:** See 9k step 2 section for a more automated way of doing this).

For SNP 12_31124, there was no BLAST hit according to SNPMeta. In this case, we pick the smaller physical position on chr1H to use in the VCF and write notes in the INFO field of the VCF file for the other duplicates with the following format. Format from previous BOPA SNPs relative to Morex v1: https://github.com/lilei1/9k_BOPA_SNP/blob/617faed6534ddbf94c287636a068ac4c4f5b25c8/BOPA_9k_vcf_Morex_refv1/sorted_all_BOPA_masked_95idt.vcf

More generally, IPK blast search the SNP fasta sequences and pick the best hit to resolve duplicate. If there is no gene for a SNP (e.g., BK_02) or the blast search has multiple identical/equally good results, we will choose the smaller (left) position and put notes in the Info field.

But we will list ALTCHR first before ALTPOS since that is the format that was already present in the current VCF file we are working with. Example:
> ALTCHR=chr1H;ALTPOS=2525021;B

Example of full line (except in our case we will list ALTCHR first and ALTPOS second:
> chr1H    2514495    11_10460    T    C    .    .    ALTPOS=2525021;ALTCHR=chr1H;B

If there are more than 2 duplicates, we can format the INFO field as follows:
> chr1H    479986322    11_21126    G    C    .    .    ALTPOS=151517613,188397002;ALTCHR=chrUn,chrUn;B

Now, repeat this process to resolve remaining duplicates.

**Step 3:** Rescue failed SNPs

Here are the general rules the script uses for picking SNPs when they fail out of SNP_Utils:

First, we will use the `BOPA_morex_v2_idt95_failed.fasta` file containing failed SNPs only and BLAST using the IPK server with the following parameters:
- Program: blastn
- Database(s): Barley Pseudomolecules Morex v2.0 2019

We will do a basic search. This may take a minute or two. Next, we will save the page as an HTML file (we only want the HTML only and NOT the complete webpage), this file will be given as input to the script `parse_ipk_blast_results.py`.

Let's create a fasta file containing only the failed SNPs.

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v2
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/get_fasta_seq.sh BOPA_morex_v2_idt95_failed.log .log ~/Shared/Datasets/Genotyping/Contextual_Sequences/BOPA_SNPs/BOPA_SNPs.fasta ~/Shared/References/Reference_Sequences/Barley/Morex_v2
```

Use script to resolve failed SNPs.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
./parse_ipk_blast_results.py \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed/BOPA_morex_v2_idt95_failed.fasta \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed/IPK_BLAST_bopa_95idt_failed.html \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed/bopa_95idt_failed_resolved.vcf
```

Manually check resolved SNPs, if reference allele is `-`, this allele needs to be removed since we cannot resolve this.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed
grep -v "-" bopa_95idt_failed_resolved.vcf > bopa_95idt_failed_resolved_noIns.vcf
```

Append rescued failed SNPs to VCF and sort.

```bash
# In dir: /Users/chaochih/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
tail -n +2 duplicates_and_failed/bopa_95idt_failed_resolved_noIns.vcf >> BOPA_morex_v2_idt95.vcf
grep "#" BOPA_morex_v2_idt95.vcf > BOPA_morex_v2_idt95_sorted.vcf
grep -v "#" BOPA_morex_v2_idt95.vcf | sort -V -k1,2 >> BOPA_morex_v2_idt95_sorted.vcf
```

---

### Methods: 9k iSelect

**Step 1:** Run SNP_Utils (https://github.com/mojaveazure/SNP_Utils)

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v2
module load python3/3.6.3_anaconda5.0.1
~/software_development/SNP_Utils/snp_utils.py CONFIG -d Barley_Morex_V2_pseudomolecules.fasta -k -i 90 -e 0.000001 -s 350 -c blast_9k_morex_v2_idt90

# Define variables for files
LOOKUP_TABLE=/panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/lookup_tables/9k_lookup.txt
GENETIC_MAP=/panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/genetic_maps/BOPA/2013_iSelect_Genetic_Map.txt
OUT_PREFIX=9k_morex_v2_idt90

# Run SNP Utils BLAST
~/software_development/SNP_Utils/snp_utils.py BLAST -l ${LOOKUP_TABLE} -c blast_9k_morex_v2_idt90 -b -m ${GENETIC_MAP} -d -t 100000 -o ${OUT_PREFIX}
```

Summarizing SNPs:

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v2/bopa_9k_50k
# Total number of SNPs
grep -v "#" 9k_morex_v2_idt90.vcf | cut -f 3 | wc -l
    7792

# Total number of unique SNPs
grep -v "#" 9k_morex_v2_idt90.vcf | cut -f 3 | sort -u | wc -l
    7687

# Total number of duplicates
# Note: Need to run head -n 65 to see all duplicates, not showing here to save space
grep -v "#" 9k_morex_v2_idt90.vcf | cut -f 3 | sort | uniq -c | sort -n -r | head
     10 SCRI_RS_149262
      9 SCRI_RS_10020
      5 SCRI_RS_69139
      5 SCRI_RS_4728
      5 SCRI_RS_235421
      4 SCRI_RS_167988
      4 SCRI_RS_114487
      3 SCRI_RS_69385
      3 SCRI_RS_69163
      3 SCRI_RS_235806

# Save the duplicate SNP names to a file
grep -v "#" 9k_morex_v2_idt90.vcf | cut -f 3 | sort | uniq -c | sort -n -r | grep -wv 1 | awk '{ print $2 }' | sort -V > 9k_duplicate_snp_names.txt
```

After running `snp_utils.py`, we get:
7,687 SNPs without duplicates
63 SNPs with duplicates
119 failed SNPs

**Step 2:** Resolve duplicate SNPs

Check if SNP hits a gene by searching in the SNPMeta output file `~/GitHub/morex_reference/morex_v2/data/bopa_and_9k_snpmeta_output.txt`. This step is intended to reduce noise since we expect all the BOPA and 9k SNPs to only be in genic regions. If there is a gene hit, BLAST search the gene with IPK server using default settings.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed
# Save duplicate SNPs in snp meta output file
head -n 1 ../../data/bopa_and_9k_snpmeta_output.txt > snpmeta_output_9k_duplicates_only.txt
grep -f 9k_duplicate_snp_names.txt ../../data/bopa_and_9k_snpmeta_output.txt >> snpmeta_output_9k_duplicates_only.txt
# Generate a list of SNP names that hit a gene
grep -f 9k_duplicate_snp_names.txt ../../data/bopa_and_9k_snpmeta_output.txt | cut -f 1 | sort -V > temp_9k_dup_in_snpmeta_output.txt

wc -l temp_9k_dup_in_snpmeta_output.txt 
      62 temp_9k_dup_in_snpmeta_output.txt

wc -l 9k_duplicate_snp_names.txt 
      63 9k_duplicate_snp_names.txt

# Figure out which SNP did not hit a gene
grep -vf temp_9k_dup_in_snpmeta_output.txt 9k_duplicate_snp_names.txt 
    BK_02
```

Generate a list of contextual sequences with the duplicate SNPs for the IPK blast search.

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v2/bopa_9k_50k
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/get_fasta_seq.sh \
    9k_duplicate_snp_names.txt \
    .txt \
    ~/Shared/Datasets/Genotyping/Contextual_Sequences/BOPA_and_9k_SCRI_SNPs.fasta \
    ~/Shared/References/Reference_Sequences/Barley/Morex_v2/bopa_9k_50k

# Output message
BK_02 not in contextual fasta file

# Create a VCF file containing only the duplicate SNPs
grep -f 9k_duplicate_snp_names.txt 9k_morex_v2_idt90.vcf | sort -k3,3 > 9k_duplicates_idt90.vcf
```

IPK blast search the SNP fasta sequences and pick the best hit to resolve duplicate. If there is no gene for a SNP (e.g., BK_02) or the blast search has multiple identical/equally good results, we will choose the smaller (left) position and put notes in the Info field of the file `9k_duplicates_idt90_resolved.vcf`.

**Step 3:** Rescue failed SNPs

Following the same procedure as the BOPA SNPs, we "manually" blast using the IPK server https://webblast.ipk-gatersleben.de/barley_ibsc/.

First, let's create a fasta file containing only the failed SNPs.

```bash
# In dir: ~/Shared/Datasets/Genotyping/Contextual_Sequences
# First, create fasta file containing all BOPA and 9k SCRI SNPs
cp BOPA_SNPs/BOPA_SNPs.fasta BOPA_and_9k_SCRI_SNPs.fasta
cat 9k_snp_SCRI/SCRI_SNPs.fasta >> BOPA_and_9k_SCRI_SNPs.fasta

# Then, extract only the failed SNPs
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/get_fasta_seq.sh \
    9k_morex_v2_idt90_failed.log \
    .log \
    ~/Shared/Datasets/Genotyping/Contextual_Sequences/BOPA_and_9k_SCRI_SNPs.fasta \
    ~/Shared/References/Reference_Sequences/Barley/Morex_v2
```

Then, we will use the `9k_morex_v2_idt90_failed.fasta` file containing failed SNPs only and BLAST using the IPK server with the following parameters:
- Program: blastn
- Database(s): Barley Pseudomolecules Morex v2.0 2019

Use script to resolve failed SNPs.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
./parse_ipk_blast_results.py \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed/9k_morex_v2_idt90_failed.fasta \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed/IPK_BLAST_9k_90idt_failed.html \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed/9k_90idt_failed_resolved.vcf
```

Manually check resolved SNPs, if reference allele is `-`, this allele needs to be removed since we cannot resolve this.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed
grep -v "-" 9k_90idt_failed_resolved.vcf > 9k_90idt_failed_resolved_noIns.vcf

# Save list of snps we could not resolve.
grep "-" 9k_90idt_failed_resolved.vcf > 9k_90idt_ref_deletion.vcf
```

Append duplicate SNPs and rescued failed SNPs to VCF and sort.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
cat duplicates_and_failed/9k_duplicates_idt90_resolved.vcf duplicates_and_failed/9k_90idt_failed_resolved_noIns.vcf | grep -v "#" >> 9k_morex_v2_idt90.vcf

# Save header lines
grep "#" 9k_morex_v2_idt90.vcf > 9k_morex_v2_idt90_sorted.vcf 
# Sort VCF
grep -v "#" 9k_morex_v2_idt90.vcf | sort -k1,2 -V >> 9k_morex_v2_idt90_sorted.vcf
```

Cleanup, move intermediate files in `duplicates_and_failed` directory to subdirectory called `intermediates`.

On MSI, copy final sorted VCF to `/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/bopa_9k_50k`.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
cp *_sorted.vcf /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/bopa_9k_50k
```

---

### Methods: 50k iSelect

**Step 1:** Run SNP_Utils (https://github.com/mojaveazure/SNP_Utils)

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v2
module load python3/3.6.3_anaconda5.0.1
~/software_development/SNP_Utils/snp_utils.py CONFIG -d Barley_Morex_V2_pseudomolecules.fasta -k -i 90 -e 0.000001 -s 350 -c blast_50k_morex_v2_idt90

# In dir: ~/Shared/Datasets/Genotyping/Contextual_Sequences/barley_50k
# Reformat lookup table for SNP Utils (sym linked in Github repo)
cut -f 1,3 barley_50k_snps.txt | tail -n +2 > snp_utils_50k_lookup_table.txt
```

To get the genetic map, I downloaded it from Supplemental Table 3 [Bayer et al. 2017 Frontiers in Plant Science](https://doi.org/10.3389/fpls.2017.01792). Then, converted the `.xlsx` to a `.csv` file and ran the following to format it for SNP utils.

```bash
# In dir: /Users/chaochih/Downloads
awk -F, '{ print $7,$3,$8 }' table3.csv | tail -n +2 > ~/GitHub/morex_reference/genetic_maps/50k_iSelect_Genetic_Map.txt

# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v2
# Define variables for files
LOOKUP_TABLE=/panfs/roc/groups/9/morrellp/shared/Datasets/Genotyping/Contextual_Sequences/barley_50k/snp_utils_50k_lookup_table.txt
GENETIC_MAP=/panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/genetic_maps/50k_iSelect_Genetic_Map.txt
OUT_PREFIX=50k_morex_v2_idt90

# Run SNP Utils BLAST
~/software_development/SNP_Utils/snp_utils.py BLAST -l ${LOOKUP_TABLE} -c blast_50k_morex_v2_idt90 -b -m ${GENETIC_MAP} -d -t 100000 -o ${OUT_PREFIX}
```

Summarizing SNPs:

```bash
# Total number of SNPs
grep -v "#" 50k_morex_v2_idt90.vcf | cut -f 3 | wc -l


# Total number of unique SNPs
grep -v "#" 50k_morex_v2_idt90.vcf | cut -f 3 | sort -u | wc -l


# Total number of duplicates
# Note: Need to run head -n 65 to see all duplicates, not showing here to save space
grep -v "#" 50k_morex_v2_idt90.vcf | cut -f 3 | sort | uniq -c | sort -n -r | head

```

After running `snp_utils.py`, we get:
xxx SNPs without duplicates
xx SNPs with duplicates
xx failed SNPs

**Step 2:** Resolve duplicate SNPs

**Step 3:** Rescue failed SNPs

Following the same procedure as the BOPA SNPs, we "manually" blast using the IPK server https://webblast.ipk-gatersleben.de/barley_ibsc/.

First, let's create a fasta file containing only the failed SNPs.

```bash
# In dir: ~/Shared/Datasets/Genotyping/Contextual_Sequences
# First, create fasta file containing all BOPA and 9k SCRI SNPs
cp BOPA_SNPs/BOPA_SNPs.fasta BOPA_and_9k_SCRI_SNPs.fasta
cat 9k_snp_SCRI/SCRI_SNPs.fasta >> BOPA_and_9k_SCRI_SNPs.fasta

# Then, extract only the failed SNPs
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/get_fasta_seq.sh \
    9k_morex_v2_idt90_failed.log \
    .log \
    ~/Shared/Datasets/Genotyping/Contextual_Sequences/BOPA_and_9k_SCRI_SNPs.fasta \
    ~/Shared/References/Reference_Sequences/Barley/Morex_v2
```

Then, we will use the `9k_morex_v2_idt90_failed.fasta` file containing failed SNPs only and BLAST using the IPK server with the following parameters:
- Program: blastn
- Database(s): Barley Pseudomolecules Morex v2.0 2019

Use script to resolve failed SNPs.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
./parse_ipk_blast_results.py \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed/9k_morex_v2_idt90_failed.fasta \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed/IPK_BLAST_9k_90idt_failed.html \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed/9k_90idt_failed_resolved.vcf
```

