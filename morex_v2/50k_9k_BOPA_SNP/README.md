# 50k, 9k, and BOPA SNPs

This directory contains the scripts used to re-map the BOPA, 9k iSelect, and 50k iSelect SNPs to the Morex v2 reference genome.

### Navigation: Jump to Section

- [Data](#data)
- [Data exploration](#data-exploration)
- [Data Preparation](#data-preparation)
- [Methods: BOPA](#methods-bopa)
- [Methods: 9k iSelect](#methods-9k-iselect)
- [Check overlapping BOPA and 9k iSelect SNPs](#check-overlapping-bopa-and-9k-iselect-snps)
- [Methods: 50k iSelect](#methods-50k-iselect)

---

## Data

The barley 50k iSelect SNP array data is published in [Bayer et al. 2017 Frontiers Plant Science](https://doi.org/10.3389/fpls.2017.01792). The 50k data can be [downloaded here](https://ics.hutton.ac.uk/50k/index.pl).

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
# Using ncbi-blast-2.9.0+
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
module load python3/3.6.3_anaconda5.0.1
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v2/bopa_9k_50k/intermediates
#~/software_development/SNP_Utils/snp_utils.py CONFIG -d ../../Barley_Morex_V2_pseudomolecules.fasta -k -i 95 -c blast_bopa_morex_v2_idt95
~/software_development/SNP_Utils/snp_utils.py CONFIG -d ../../Barley_Morex_V2_pseudomolecules.fasta -k -i 90 -c blast_bopa_morex_v2_idt90

# Define variables for files
LOOKUP_TABLE=/panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/lookup_tables/BOPA_lookup.txt
GENETIC_MAP=/panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/genetic_maps/BOPA/2013_iSelect_Genetic_Map.txt
#OUT_PREFIX=BOPA_morex_v2_idt95
OUT_PREFIX=BOPA_morex_v2_idt90

# Run SNP Utils BLAST
#~/software_development/SNP_Utils/snp_utils.py BLAST -l ${LOOKUP_TABLE} -c blast_bopa_morex_v2_idt95 -b -m ${GENETIC_MAP} -d -t 100000 -o ${OUT_PREFIX}
~/software_development/SNP_Utils/snp_utils.py BLAST -l ${LOOKUP_TABLE} -c blast_bopa_morex_v2_idt90 -b -m ${GENETIC_MAP} -d -t 100000 -o ${OUT_PREFIX}
```

Originally we used a 95 identity threshold but then found that for some overlapping SNPs between the BOPA and 9k sets, the BLAST results did not match and resulted in discordant positions. So, we reduced the identity threshold to 90.

Summarizing SNPs for 90idt:

```bash
# Total SNPs
grep -v "#" BOPA_morex_v2_idt90.vcf | cut -f 3 | wc -l
    3024
# Total unique SNPs
grep -v "#" BOPA_morex_v2_idt90.vcf | cut -f 3 | sort -u | wc -l
    3017
# Total number of duplicates
grep -v "#" BOPA_morex_v2_idt90.vcf | cut -f 3 | sort | uniq -c | sort -n -r | grep -vw "1" | wc -l
    5
```

After running `snp_utils.py` for 95idt, we get:

- 2,985 unique SNPs
- 5 SNPs with duplicates
- 52 failed SNPs

For **90idt**, we get:

- 3,024 unique SNPs
- 5 SNPs with duplicates
- 19 failed SNPs

For the ones that failed, we "manually" BLAST using the IPK server: https://webblast.ipk-gatersleben.de/barley_ibsc/. Previously (for Morex v1), Li had many undergrads help us manually BLAST search the duplicate and failed SNPs. But, this is somewhat more prone to human error and is not possible without a many undergrads. So, we will use the script below to resolve the duplicates and failed SNPs by parsing an HTML file containg the IPK Barley BLAST server results.

**Step 2:** Resolve duplicate SNPs

Let's resolve the duplicate SNPs, since there are only 5, we can do this manually.

List of 5 BOPA SNPs with duplicates:

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
grep -v "#" BOPA_morex_v2_idt90.vcf | cut -f 3 | sort | uniq -c | sort -n -r | grep -vw "1"
   3 12_31124
   3 12_30822
   2 12_30861
   2 12_21135
   2 11_20074
```

Make a list of duplicate SNPs:

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
grep -v "#" BOPA_morex_v2_idt90.vcf | cut -f 3 | sort | uniq -c | sort -n -r | grep -vw "1" | sed -e 's,   ,,' | cut -d' ' -f 2 | sort -V > duplicates_and_failed/bopa_duplicate_snp_names.txt
# Save duplicates to vcf
grep -f duplicates_and_failed/bopa_duplicate_snp_names.txt BOPA_morex_v2_idt90.vcf | sort -V -k3,3 > duplicates_and_failed/bopa_duplicates_idt90.vcf
# Make copy for resolving duplicates
cp bopa_duplicates_idt90.vcf bopa_duplicates_idt90_resolved.vcf
```

Remove duplicates from VCF temporarily:

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
grep -vf duplicates_and_failed/bopa_duplicate_snp_names.txt BOPA_morex_v2_idt90.vcf > temp_bopa_idt90_noDups.vcf
```

Check if SNP hits a gene by searching in the SNPMeta output file `~/GitHub/morex_reference/morex_v2/data/bopa_and_9k_snpmeta_output.txt`. This step is intended to reduce noise since we expect all the BOPA and 9k SNPs to only be in genic regions. If there is a gene hit, BLAST search the gene with IPK server, pick the best hit, and write notes in the INFO field of the VCF. If there is no hit, choose the smallest position (smallest chromosome and smallest physical position) and write notes in the INFO field of the VCF file. For more details, take a look at the file [`Resolve_Duplicate_SNPs_Methods.md`](https://github.com/MorrellLAB/morex_reference/blob/master/morex_v2/50k_9k_BOPA_SNP/Resolve_Duplicate_SNPs_Methods.md).

We will list ALTCHR first before ALTPOS since that is the format that was already present in the current VCF file we are working with. Example:
> ALTCHR=chr1H;ALTPOS=2525021;B

Example of full line (except in our case we will list ALTCHR first and ALTPOS second:
> chr1H    2514495    11_10460    T    C    .    .    ALTPOS=2525021;ALTCHR=chr1H;B

If there are more than 2 duplicates, we can format the INFO field as follows:
> chr1H    479986322    11_21126    G    C    .    .    ALTPOS=151517613,188397002;ALTCHR=chrUn,chrUn;B

Now, repeat this process to resolve remaining duplicates.

Resolved duplicates are in the file `bopa_duplicates_idt90_resolved.vcf` and will be combined with the file `temp_bopa_idt90_noDups.vcf` and failed SNPs (see step 3) at the end and sorted.

**Step 3:** Rescue failed SNPs

*File preparation*. Let's create a fasta file containing only the failed SNPs.

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v2/bopa_9k_50k/intermediates
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/get_fasta_seq.sh BOPA_morex_v2_idt90_failed.log .log ~/Shared/Datasets/Genotyping/Contextual_Sequences/BOPA_SNPs/BOPA_SNPs.fasta ~/Shared/References/Reference_Sequences/Barley/Morex_v2/bopa_9k_50k/intermediates
```

Here are the general rules the script uses for picking SNPs when they fail out of SNP_Utils:

First, we will use the `BOPA_morex_v2_idt90_failed.fasta` file containing failed SNPs only and BLAST using the [IPK server](https://webblast.ipk-gatersleben.de/barley_ibsc/) with the following parameters:

- Program: blastn
- Database(s): Barley Pseudomolecules Morex v2.0 2019

We will do a basic search. This may take a minute or two. Next, we will save the page as an HTML file (we only want the HTML only and NOT the complete webpage), this file will be given as input to the script `parse_ipk_blast_results.py`.

Use script to resolve failed SNPs.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
./parse_ipk_blast_results.py \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed/BOPA_morex_v2_idt90_failed.fasta \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed/BOPA_morex_v2_idt90_failed.html \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed/bopa_90idt_failed_resolved.vcf
```

Manually check resolved SNPs, if reference allele is `-`, this SNP needs to be removed since we cannot resolve this due to a deletion in the reference.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed
# Rename log file
mv unresolved_snps.log bopa_90idt_unresolved_snps.log
# Check resolved snps
# Here are ones we could not resolve
grep -w "-" bopa_90idt_failed_resolved.vcf
chr2H	48775454	11_10084	-	C	.	.	B;Identities=230/246(93%),failed
chr1H	8299812	12_30951	-	T	.	.	B;Identities=119/121(98%),failed

# Remove snps where ref is ins
grep -vw "-" bopa_90idt_failed_resolved.vcf > bopa_90idt_failed_resolved_noRefDel.vcf
```

**Wrapping up and combining files:**

Append resolved duplicate snps and rescued failed SNPs to VCF and sort.

```bash
# In dir: /Users/chaochih/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
tail -n +2 duplicates_and_failed/bopa_90idt_failed_resolved_noRefDel.vcf | cat - duplicates_and_failed/bopa_duplicates_idt90_resolved.vcf >> temp_bopa_idt90_noDups.vcf
grep "#" temp_bopa_idt90_noDups.vcf > BOPA_morex_v2_idt90.vcf
# Sort
grep -v "#" temp_bopa_idt90_noDups.vcf | sort -V -k1,2 >> BOPA_morex_v2_idt90.vcf
```

Convert pseudomolecular positions to parts positions:

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
module load python2/2.7.12_anaconda4.2
~/GitHub/File_Conversions/Pseudomolecules_to_Parts_v2.py --vcf BOPA_morex_v2_idt90.vcf > BOPA_morex_v2_idt90_parts.vcf
```

Make a symbolic link in Shared space on MSI:

```bash
cd ~/Shared/References/Reference_Sequences/Barley/Morex_v2/bopa_9k_50
ln -s /panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/BOPA_morex_v2_idt90.vcf
ln -s /panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/BOPA_morex_v2_idt90_parts.vcf
```

Fix VCF header line by adding contig lengths to vcf header:

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
# Run on MSI, relies on Picard Jar file
# Parts positions
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/add_contig_length_to_header.sh \
  ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/BOPA_morex_v2_idt90_parts.vcf \
  ~/Shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules_parts.dict

# Cleanup the filename
mv BOPA_morex_v2_idt90_parts_fixedHeader.vcf BOPA_morex_v2_idt90_parts.vcf
# Index vcf file
module load gatk/4.1.2
gatk IndexFeatureFile -F BOPA_morex_v2_idt90_parts.vcf

# Pseudomolecular positions
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/add_contig_length_to_header.sh \
  ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/BOPA_morex_v2_idt90.vcf \
  ~/Shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules.dict

# Cleanup the filename
mv BOPA_morex_v2_idt90_fixedHeader.vcf BOPA_morex_v2_idt90.vcf
gatk IndexFeatureFile -F BOPA_morex_v2_idt90.vcf
```

**Troubleshooting SNPs that got left out:**

The following SNPs did not get processed by SNP_Utils for some reason (i.e., these SNPs are not in the VCF file output, duplicates, or failed SNPs log file). So, we will get their positions the same way we resolve the failed SNPs.
- 11_10460
- 12_10571

The FASTA sequences are as follows:

```bash
>11_10460
TGTGCAAGTACAAATACCATTTGTTCATCCATCTATTTTGCAGCAGCTAAACCCATGCAAGGTATTCCTCCAGCARCAGTGCAGCCCTGTGGCAATGTCACAACGTATTGCAAGGTCGCARATGTTGCAACAGAGCAGTTGCCATGTGTTGCAGCAACARTGTTGCCAACAACTGCCGCAAATCCCCGAACAACTCCGCCATGAGGCAGTCCGTGCAATCGTCTACTCTATCGTCCTGCAA
>12_10571
TGGCTGGTGCGGCCAAYGGTGAGAATGGAGAAGGGCCTGGAGCCGCCTGATTTTGACGAGRGCAATTCATATGATTTTCCGTTGTATGCGGTTCACATATCKAATTATCTATGCATTGTTA
```

Now, we will do an IPK Barley BLAST search and use script to resolve these failed SNPs.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts
./parse_ipk_blast_results.py \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/snp_utils_problem_snps/BOPA_snp_utils_problem_snps.fasta \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/snp_utils_problem_snps/IPK_Barley_BLAST_2_BOPA_failed_snpUtils.html \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/snp_utils_problem_snps/BOPA_morex_v2_idt90_snpUtils_problem_snps_resolved.vcf
```

Add the resolved problem SNPs to main vcf file.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/snp_utils_problem_snps
grep -v "#" BOPA_morex_v2_idt90_snpUtils_problem_snps_resolved.vcf >> ../BOPA_morex_v2_idt90.vcf
```

Re-sort VCF file and convert to parts positions too.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
gatk SortVcf -I BOPA_morex_v2_idt90.vcf -O BOPA_morex_v2_idt90_sorted.vcf
# Rename file
mv BOPA_morex_v2_idt90_sorted.vcf BOPA_morex_v2_idt90.vcf
# Index file
gatk IndexFeatureFile --input BOPA_morex_v2_idt90.vcf

# Convert pseudomolecular positions to parts positions
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
module load python2/2.7.12_anaconda4.2
module load gatk_ML/4.1.8
~/GitHub/File_Conversions/Pseudomolecules_to_Parts_v2.py --vcf BOPA_morex_v2_idt90.vcf > BOPA_morex_v2_idt90_parts.vcf
# Index parts vcf
gatk IndexFeatureFile --input BOPA_morex_v2_idt90_parts.vcf
```

---

### Methods: 9k iSelect

**Step 1:** Run SNP_Utils (https://github.com/mojaveazure/SNP_Utils)

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v2/bopa_9k_50k/intermediates
module load python3/3.6.3_anaconda5.0.1
~/software_development/SNP_Utils/snp_utils.py CONFIG -d ../../Barley_Morex_V2_pseudomolecules.fasta -k -i 90 -e 0.000001 -s 350 -c blast_9k_morex_v2_idt90

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
grep -v "#" 9k_morex_v2_idt90.vcf | cut -f 3 | sort | uniq -c | sort -n -r | grep -vw "1" | wc -l
    63
```

After running `snp_utils.py`, we get:
7,687 SNPs without duplicates
63 SNPs with duplicates
119 failed SNPs

**Step 2:** Resolve duplicate SNPs

Make a list of duplicate SNPs:

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
# Save the duplicate SNP names to a file
grep -v "#" 9k_morex_v2_idt90.vcf | cut -f 3 | sort | uniq -c | sort -n -r | grep -wv 1 | awk '{ print $2 }' | sort -V > duplicates_and_failed/9k_duplicate_snp_names.txt
# Save duplicates to vcf
grep -f duplicates_and_failed/9k_duplicate_snp_names.txt 9k_morex_v2_idt90.vcf | sort -k3,3 > duplicates_and_failed/9k_duplicates_idt90.vcf
# Make copy for resolving duplicates
cp duplicates_and_failed/9k_duplicates_idt90.vcf duplicates_and_failed/9k_duplicates_idt90_resolved.vcf
```

Remove duplicates from VCF temporarily:

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
grep -vf duplicates_and_failed/9k_duplicate_snp_names.txt 9k_morex_v2_idt90.vcf > temp_9k_idt90_noDups.vcf
```

We will use the same approach as we did for the BOPA SNPs (above). Check if SNP hits a gene by searching in the SNPMeta output file `~/GitHub/morex_reference/morex_v2/data/bopa_and_9k_snpmeta_output.txt`. This step is intended to reduce noise since we expect all the BOPA and 9k SNPs to only be in genic regions. If there is a gene hit, BLAST search the gene with IPK server, pick the best hit, and write notes in the INFO field of the VCF. If there is no hit, choose the smallest position (smallest chromosome and smallest physical position) and write notes in the INFO field of the VCF file. For more details, take a look at the file [`Resolve_Duplicate_SNPs_Methods.md`](https://github.com/MorrellLAB/morex_reference/blob/master/morex_v2/50k_9k_BOPA_SNP/Resolve_Duplicate_SNPs_Methods.md).

Repeat this process to resolve remaining duplicates.

Resolved duplicates are in the file `9k_duplicates_idt90_resolved.vcf` and will be combined with the file `temp_9k_idt90_noDups.vcf` and failed SNPs (see step 3) at the end and sorted.

**Step 3:** Rescue failed SNPs

Following the same procedure as the BOPA SNPs, we "manually" blast using the IPK server https://webblast.ipk-gatersleben.de/barley_ibsc/.

*File preparation.* Let's create a fasta file containing only the failed SNPs.

```bash
# In dir: ~/Shared/Datasets/Genotyping/Contextual_Sequences
# First, create fasta file containing all BOPA and 9k SCRI SNPs
cp BOPA_SNPs/BOPA_SNPs.fasta BOPA_and_9k_SCRI_SNPs.fasta
cat 9k_snp_SCRI/SCRI_SNPs.fasta >> BOPA_and_9k_SCRI_SNPs.fasta

# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v2/bopa_9k_50k/intermediates
# Then, extract only the failed SNPs
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/get_fasta_seq.sh \
    9k_morex_v2_idt90_failed.log \
    .log \
    ~/Shared/Datasets/Genotyping/Contextual_Sequences/BOPA_and_9k_SCRI_SNPs.fasta \
    ~/Shared/References/Reference_Sequences/Barley/Morex_v2/bopa_9k_50k/intermediates
```

Here are the general rules the script uses for picking SNPs when they fail out of SNP_Utils:

First, we will use the `9k_morex_v2_idt90_failed.fasta` file containing failed SNPs only and BLAST using the [IPK server](https://webblast.ipk-gatersleben.de/barley_ibsc/) with the following parameters:

- Program: blastn
- Database(s): Barley Pseudomolecules Morex v2.0 2019

We will do a basic search. This may take a minute or two. Next, we will save the page as an HTML file (we only want the HTML only and NOT the complete webpage), this file will be given as input to the script `parse_ipk_blast_results.py`.

Use script to resolve failed SNPs.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
./parse_ipk_blast_results.py \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed/9k_morex_v2_idt90_failed.fasta \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed/IPK_BLAST_9k_90idt_failed.html \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed/9k_90idt_failed_resolved.vcf
```

Manually check resolved SNPs, if reference allele is `-`, this SNP needs to be removed since we cannot resolve this due to a deletion in the reference.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed
# Check resolved snps
# Here are ones we could not resolve
grep -w "-" 9k_90idt_failed_resolved.vcf
chr2H	48775454	11_10084	-	C	.	.	B;Identities=230/246(93%),failed
chr1H	8299812	12_30951	-	T	.	.	B;Identities=119/121(98%),failed
chr2H	25095696	SCRI_RS_143250	-	G	.	.	B;Identities=118/121(98%),failed
chr6H	398167515	SCRI_RS_158011	-	T	.	.	B;Identities=119/121(98%),failed
chr2H	18961510	SCRI_RS_185319	-	G	.	.	B;Identities=118/121(98%),failed
chr7H	12800051	SCRI_RS_196319	-	G	.	.	B;Identities=117/121(97%),failed
chr1H	391116331	SCRI_RS_213457	-	A	.	.	B;Identities=61/68(90%),failed
chr4H	1359302	SCRI_RS_235733	-	G	.	.	B;Identities=100/130(77%),failed

# Remove SNPs where ref is ins
grep -v "-" 9k_90idt_failed_resolved.vcf > 9k_90idt_failed_resolved_noRefDel.vcf
```

**Wrapping up and combining files:**

Append duplicate SNPs and rescued failed SNPs to VCF and sort.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
tail -n +2 duplicates_and_failed/9k_90idt_failed_resolved_noRefDel.vcf | cat - duplicates_and_failed/9k_duplicates_idt90_resolved.vcf >> temp_9k_idt90_noDups.vcf
#cat duplicates_and_failed/9k_duplicates_idt90_resolved.vcf duplicates_and_failed/9k_90idt_failed_resolved_noRefDel.vcf | grep -v "#" >> 9k_morex_v2_idt90.vcf
grep "#" temp_9k_idt90_noDups.vcf > 9k_morex_v2_idt90.vcf
# Sort VCF
grep -v "#" temp_9k_idt90_noDups.vcf | sort -V -k1,2 >> 9k_morex_v2_idt90.vcf
# Cleanup
rm temp_9k_idt90_noDups.vcf
```

Cleanup, move intermediate files in `duplicates_and_failed` directory to subdirectory called `intermediates`.

Convert pseudomolecular positions to parts positions:

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
module load python2/2.7.12_anaconda4.2
~/GitHub/File_Conversions/Pseudomolecules_to_Parts_v2.py --vcf 9k_morex_v2_idt90.vcf > 9k_morex_v2_idt90_parts.vcf
```

Make a symbolic link in Shared space on MSI:

```bash
cd ~/Shared/References/Reference_Sequences/Barley/Morex_v2/bopa_9k_50
ln -s /panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/9k_morex_v2_idt90.vcf
ln -s /panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/9k_morex_v2_idt90_parts.vcf
```

Fix VCF header line by adding contig lengths to vcf header:

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
# Run on MSI, relies on Picard Jar file
# Parts positions
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/add_contig_length_to_header.sh \
  ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/9k_morex_v2_idt90_parts.vcf \
  ~/Shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules_parts.dict

# Cleanup the filename
mv 9k_morex_v2_idt90_parts_fixedHeader.vcf 9k_morex_v2_idt90_parts.vcf
# Index vcf file
module load gatk/4.1.2
gatk IndexFeatureFile -F 9k_morex_v2_idt90_parts.vcf

# Pseudomolecular positions
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/add_contig_length_to_header.sh \
  ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/9k_morex_v2_idt90.vcf \
  ~/Shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules.dict

# Cleanup the filename
mv 9k_morex_v2_idt90_fixedHeader.vcf 9k_morex_v2_idt90.vcf
gatk IndexFeatureFile -F 9k_morex_v2_idt90.vcf
```

---

### Check overlapping BOPA and 9k iSelect SNPs

Check for BOPA and 9k position concordance.

**Check 1.** BOPA used 95idt and 9k used 90idt threshold. Make sure all SNP positions for SNPs that overlap in the two sets are concordant.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
./check_position_concordance.py BOPA_morex_v2_idt95_sorted.vcf 9k_morex_v2_idt90_sorted.vcf > temp_discordant_bopa_9k_snps.txt
```

Since there are some discordant SNP positions. To identify what is causing this, I will re-run SNP_Utils on the BOPA set using the same identity threshold as in the 9k set.

```bash
module load python3/3.6.3_anaconda5.0.1
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v2/bopa_9k_50k/intermediates
~/software_development/SNP_Utils/snp_utils.py CONFIG -d ../../Barley_Morex_V2_pseudomolecules.fasta -k -i 90 -c blast_bopa_morex_v2_idt90

# Define variables for files
LOOKUP_TABLE=/panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/lookup_tables/BOPA_lookup.txt
GENETIC_MAP=/panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/genetic_maps/BOPA/2013_iSelect_Genetic_Map.txt
OUT_PREFIX=BOPA_morex_v2_idt90

# Run SNP Utils BLAST
~/software_development/SNP_Utils/snp_utils.py BLAST -l ${LOOKUP_TABLE} -c blast_bopa_morex_v2_idt90 -b -m ${GENETIC_MAP} -d -t 100000 -o ${OUT_PREFIX}
```

**Check2.** BOPA used 90idt and 9k used 90idt threshold.

```bash
./check_position_concordance.py BOPA_morex_v2_idt90.vcf 9k_morex_v2_idt90.vcf
```

Everything looks ok now.

---

### Methods: 50k iSelect

**Step 1:** Run SNP_Utils (https://github.com/mojaveazure/SNP_Utils)

```bash
# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v2/bopa_9k_50k/intermediates
module load python3/3.6.3_anaconda5.0.1
~/software_development/SNP_Utils/snp_utils.py CONFIG -d ../../Barley_Morex_V2_pseudomolecules.fasta -k -i 90 -e 0.000001 -s 350 -c blast_50k_morex_v2_idt90

# In dir: ~/Shared/Datasets/Genotyping/Contextual_Sequences/barley_50k
# Reformat lookup table for SNP Utils (sym linked in Github repo)
cut -f 1,3 barley_50k_snps.txt | tail -n +2 > snp_utils_50k_lookup_table.txt
```

To get the genetic map, I downloaded it from Supplemental Table 3 [Bayer et al. 2017 Frontiers in Plant Science](https://doi.org/10.3389/fpls.2017.01792). Then, converted the `.xlsx` to a `.csv` file and ran the following to format it for SNP utils.

```bash
# In dir: /Users/chaochih/Downloads
awk -F, '{ print $7,$3,$8 }' table3.csv | tail -n +2 > ~/GitHub/morex_reference/genetic_maps/50k_iSelect_Genetic_Map.txt

# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v2/bopa_9k_50k/intermediates
# Define variables for files
LOOKUP_TABLE=/panfs/roc/groups/9/morrellp/shared/Datasets/Genotyping/Contextual_Sequences/barley_50k/snp_utils_50k_lookup_table.txt
LOOKUP_TABLE=
GENETIC_MAP=/panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/genetic_maps/50k_iSelect_Genetic_Map.txt
OUT_PREFIX=50k_morex_v2_idt90

# Run SNP Utils BLAST
~/software_development/SNP_Utils/snp_utils.py BLAST -l ${LOOKUP_TABLE} -c blast_50k_morex_v2_idt90 -b -m ${GENETIC_MAP} -d -t 100000 -o ${OUT_PREFIX}
```

This SNP_Utils BLAST, resulted in the following error.

```bash
# Here are the relevant parts of the error message
No reference allele for JHI-Hv50k-2016-132670:1.0566e-26, have you run Hsp.get_snp_position() yet?
hmm JHI-Hv50k-2016-132670:6.4497e-14
Traceback (most recent call last):
  File "/home/morrellp/liux1299/software_development/SNP_Utils/snp_utils.py", line 285, in <module>
    main()
  File "/home/morrellp/liux1299/software_development/SNP_Utils/snp_utils.py", line 262, in main
    snp_list, no_snps, ref_gen, bconf = blast_based(args, lookup_dict)
  File "/home/morrellp/liux1299/software_development/SNP_Utils/snp_utils.py", line 103, in blast_based
    hsp.add_snp(lookup=lookup)
  File "/home/morrellp/liux1299/.local/lib/python3.6/site-packages/overload.py", line 181, in f
    return callable(*usable_args, **_kw)
  File "/panfs/roc/groups/9/morrellp/liux1299/software_development/SNP_Utils/Objects/blast.py", line 270, in add_snp
    self.add_snp(this_snp=snp.SNP(lookup=lookup, hsp=self))
  File "/home/morrellp/liux1299/.local/lib/python3.6/site-packages/overload.py", line 181, in f
    return callable(*usable_args, **_kw)
  File "/panfs/roc/groups/9/morrellp/liux1299/software_development/SNP_Utils/Objects/snp.py", line 103, in __init__
    self._position = hsp.get_snp_position(lookup=lookup)
  File "/home/morrellp/liux1299/.local/lib/python3.6/site-packages/overload.py", line 181, in f
    return callable(*usable_args, **_kw)
  File "/panfs/roc/groups/9/morrellp/liux1299/software_development/SNP_Utils/Objects/blast.py", line 235, in get_snp_position
    expected=adj
  File "/home/morrellp/liux1299/.local/lib/python3.6/site-packages/overload.py", line 181, in f
    return callable(*usable_args, **_kw)
  File "/panfs/roc/groups/9/morrellp/liux1299/software_development/SNP_Utils/Objects/blast.py", line 198, in get_snp_position
    raise NoSNPError("Cannot find the SNP in " + self) # Error out to avoid infinite loops
TypeError: must be str, not Hsp
```

For now, I will save the problem SNPs and come back to resolve them later on.

```bash
# In dir: ~/Shared/Datasets/Genotyping/Contextual_Sequences/barley_50k
# Save SNP name
grep JHI-Hv50k-2016-132670 snp_utils_50k_lookup_table.txt | cut -f 1 > snp_utils_50k_lookup_table-problem_snp_names.txt
grep JHI-Hv50k-2016-210233 snp_utils_50k_lookup_table.txt | cut -f 1 >> snp_utils_50k_lookup_table-problem_snp_names.txt
grep JHI-Hv50k-2016-451305 snp_utils_50k_lookup_table.txt | cut -f 1 >> snp_utils_50k_lookup_table-problem_snp_names.txt
# Save the SNP Name and sequence into file
grep -f snp_utils_50k_lookup_table-problem_snp_names.txt snp_utils_50k_lookup_table.txt > snp_utils_50k_lookup_table-problem_snps.txt
# Remove problem SNP from lookup table
grep -vf snp_utils_50k_lookup_table-problem_snps.txt snp_utils_50k_lookup_table.txt > snp_utils_50k_lookup_table_filtered.txt

# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v2/bopa_9k_50k/intermediates
# Define variables for files
LOOKUP_TABLE=/panfs/roc/groups/9/morrellp/shared/Datasets/Genotyping/Contextual_Sequences/barley_50k/snp_utils_50k_lookup_table_filtered.txt
GENETIC_MAP=/panfs/roc/groups/9/morrellp/liux1299/GitHub/morex_reference/genetic_maps/50k_iSelect_Genetic_Map.txt
OUT_PREFIX=50k_morex_v2_idt90

# Run SNP Utils BLAST using lookup table without problem snps
~/software_development/SNP_Utils/snp_utils.py BLAST -l ${LOOKUP_TABLE} -c blast_50k_morex_v2_idt90 -b -m ${GENETIC_MAP} -d -t 100000 -o ${OUT_PREFIX}
```

Summarizing SNPs:

```bash
# Total number of SNPs
grep -v "#" 50k_morex_v2_idt90.vcf | cut -f 3 | wc -l
    44964

# Total number of unique SNPs
grep -v "#" 50k_morex_v2_idt90.vcf | cut -f 3 | sort -u | wc -l
    43940

# Total number of duplicates
grep -v "#" 50k_morex_v2_idt90.vcf | cut -f 3 | sort | uniq -c | sort -n -r | grep -vw "1" | wc -l
    582
```

After running `snp_utils.py`, we get:
43,940 Total SNPs (unique)
582 SNPs with duplicates
91 failed SNPs

**Step 2:** Resolve duplicate SNPs

For duplicates, the general procedure we will follow is the same as resolving the 9k duplicate snps (detailed above).

*File preparation:* Pull out only the duplicates to resolve and save to a file.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
# Create list of duplicate SNP names
grep -v "#" 50k_morex_v2_idt90.vcf | cut -f 3 | sort | uniq -c | sort -n -r | grep -vw "1" | sed -e 's,   ,,' -e 's,  ,,' | cut -d' ' -f 2 | sort -V > duplicates_and_failed/50k_duplicate_snp_names.txt
# Create temp vcf file containing only duplicate snp names
grep -wf duplicates_and_failed/50k_duplicate_snp_names.txt 50k_morex_v2_idt90.vcf | sort -V -k1,2 | sort -V -k3,3 > duplicates_and_failed/50k_duplicates_only.vcf
# Make a copy, we will resolve dup snps directly in this copy
cp duplicates_and_failed/50k_duplicates_only.vcf duplicates_and_failed/for_mackenzie/50k_duplicates_resolved.vcf
```

Convert contextual sequence from [A/B] format to IUPAC codes so we can BLAST search the SNP with the [IPK server](https://webblast.ipk-gatersleben.de/barley_ibsc/).

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
./snp_to_iupac_fasta.py ~/GitHub/morex_reference/lookup_tables/snp_utils_50k_lookup_table_filtered.txt > duplicates_and_failed/for_mackenzie/snp_utils_50k_lookup_table_filtered.fasta
```

Use the `morex_reference/morex_v2/data/50k_snpmeta_output.txt` file to reduce noise when resolving duplicate SNPs.

Mackenzie Linane helped resolve duplicate SNPs for the 50k set. The detailed steps and resolved SNPs are in a subdirectory:
- Steps documentation: https://github.com/MorrellLAB/morex_reference/tree/master/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed/for_mackenzie
- Resolved duplicate VCF file: `morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed/for_mackenzie/50k_duplicates_resolved.vcf`

**Step 3:** Rescue failed SNPs

Following the same procedure as the BOPA SNPs, we "manually" blast using the IPK server https://webblast.ipk-gatersleben.de/barley_ibsc/.

*File preparation.* Let's create a fasta file containing only the 50k failed SNPs. But, we need to convert the 50k `[A/B]` formatted SNPs to IUPAC nucleotide base code. The IPK Barley BLAST search does not accept sequences with snps formatted as `[A/B]`.

```bash
# In dir: ~/Shared/Datasets/Genotyping/Contextual_Sequences/barley_50k
module load python3/3.6.3_anaconda5.0.1
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/snp_to_iupac_fasta.py ~/Shared/Datasets/Genotyping/Contextual_Sequences/barley_50k/snp_utils_50k_lookup_table_filtered.txt > 50k_snps_iupac.fasta

# In dir: ~/Shared/References/Reference_Sequences/Barley/Morex_v2/bopa_9k_50k/intermediates
# Then, extract only the failed SNPs
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/get_fasta_seq.sh \
    50k_morex_v2_idt90_failed.log \
    .log \
    ~/Shared/Datasets/Genotyping/Contextual_Sequences/barley_50k/50k_snps_iupac.fasta \
    ~/Shared/References/Reference_Sequences/Barley/Morex_v2/bopa_9k_50k/intermediates
```

Then, we will use the `50k_morex_v2_idt90_failed.fasta` file containing failed SNPs only and BLAST using the IPK server with the following parameters:
- Program: blastn
- Database(s): Barley Pseudomolecules Morex v2.0 2019

We will do a basic search. This may take a minute or two. Next, we will save the page as an HTML file (we only want the HTML only and NOT the complete webpage), this file will be given as input to the script `parse_ipk_blast_results.py`.

Use script to resolve failed SNPs.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
./parse_ipk_blast_results.py \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed/intermediates/50k_morex_v2_idt90_failed.fasta \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed/intermediates/IPK_BLAST_50k_idt90_failed.html \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed/intermediates/50k_idt90_failed_resolved.vcf

# Rename log file
mv unresolved_snps.log 50k_idt90_unresolved_snps.log
```

Manually check resolved SNPs, if reference allele is `-`, this SNP needs to be removed since we cannot resolve this due to a deletion in the reference.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed
# Check resolved snps
# Here are ones we could not resolve
grep -w "-" 50k_idt90_failed_resolved.vcf 
chr2H	25095696	SCRI_RS_143250	-	G	.	.	B;Identities=118/121(98%),failed
chr6H	398167515	SCRI_RS_158011	-	T	.	.	B;Identities=119/121(98%),failed
chr2H	18961510	SCRI_RS_185319	-	G	.	.	B;Identities=118/121(98%),failed
chr7H	12800051	SCRI_RS_196319	-	G	.	.	B;Identities=117/121(97%),failed
chr4H	1359302	SCRI_RS_235733	-	G	.	.	B;Identities=100/130(77%),failed

# Remove SNPs where ref is ins
grep -vw "-" 50k_idt90_failed_resolved.vcf > ../50k_idt90_failed_resolved_noRefDel.vcf
```

**Checks**

Did a few checks and found we somehow missed resolving a few of the duplicate SNPs.

```bash
# Make a list of SNPs we missed in round1 of resolving duplicates
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/duplicates_and_failed $ diff -y 50k_duplicate_snp_names.txt temp_dup_resolved_snp_names.txt | grep "<" | cut -f 1 | sort -V > 50k_duplicate_snp_names_missed_round1.txt
# Pull out duplicates we missed
grep -wf 50k_duplicate_snp_names_missed_round1.txt ../50k_morex_v2_idt90.vcf | sort -k3,3 -k1,1 > for_mackenzie/50k_duplicates_resolved_missed_round1.vcf
# Add these to the duplicates vcf for future reference
grep -wf 50k_duplicate_snp_names_missed_round1.txt ../50k_morex_v2_idt90.vcf | sort -k3,3 -k1,1 >> 50k_duplicates_only.vcf
```

Follow above steps for resolving these duplicates.

**Fixing 3 problem SNPs that returned an error from SNP_Utils.**

The problem SNPs are as follows:

```bash
JHI-Hv50k-2016-132670   AGRAATGTATATAGACATATTTTAGAGTATAGATTCACTGATTTTACTYYGTATGTAGTC[T/C]CTTAGTGAAATCTCTAAAATGACTTATATTTAGGAGGAAGTATTTTGGGATGCATNAAAA
JHI-Hv50k-2016-210233   AATTTTGTACTAGAACTAGTACAAAGTTRAGACAGTTATTTTGGGACAGAGGGRGTATAA[A/G]ATTGGCGAATGAGGACTTGGAATTGGGACTGGAACGGTTGGTTACCTAGAGCAAMAKGAT
JHI-Hv50k-2016-451305   TGTTKGGGAACATAATAGATTATAAAAACTGGATTATATARTCYRGGTGRTTCCAAACAG[A/G]GCCTAGGTAAACATATGTGTGATGCCAATCCATTTGGTTGTTGAGTTGTTTGTCATTTTG
```

We will manually use IPK barley BLAST for these and run the results through the `parse_ipk_blast_results.py` script.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/scripts
./parse_ipk_blast_results.py \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/problem_snps/snp_utils_50k_lookup_table-problem_snps_AB.fasta \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/problem_snps/IPK_Barley_BLAST_3_problem_snps.html \
    ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/problem_snps/50k_problem_snps_resolved.vcf
```

**Wrapping up and combining files:**

Append duplicate SNPs and rescued failed SNPs to VCF and sort.

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
tail -n +2 duplicates_and_failed/50k_idt90_failed_resolved_noRefDel.vcf | cat - duplicates_and_failed/for_mackenzie/50k_duplicates_resolved.vcf duplicates_and_failed/for_mackenzie/50k_duplicates_resolved_missed_round1.vcf > temp_50k_idt90_dups_and_failed_resolved.vcf
# Remove duplicates from original VCF
grep -wvf duplicates_and_failed/50k_duplicate_snp_names.txt 50k_morex_v2_idt90.vcf > temp_50k_idt90_noDups.vcf
grep -wvf duplicates_and_failed/50k_duplicate_snp_names_missed_round1.txt temp_50k_idt90_noDups.vcf > temp_50k_idt90_noDups_clean.vcf
# Combine all files
cat temp_50k_idt90_noDups_clean.vcf temp_50k_idt90_dups_and_failed_resolved.vcf > 50k_morex_v2_idt90_unsorted.vcf

# Fix headers
module load gatk/4.1.2
gatk FixVcfHeader --INPUT 50k_morex_v2_idt90_unsorted.vcf --OUTPUT 50k_morex_v2_idt90_unsorted_fixedHeader.vcf
# Make sure contig lengths are accurate
gatk UpdateVcfSequenceDictionary --INPUT 50k_morex_v2_idt90_unsorted_fixedHeader.vcf --OUTPUT 50k_morex_v2_idt90_unsorted_fixedHeader_lengths.vcf --SEQUENCE_DICTIONARY /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules.dict
# Sort VCF
gatk SortVcf -I 50k_morex_v2_idt90_unsorted_fixedHeader_lengths.vcf -O 50k_morex_v2_idt90_sorted.vcf
# Rename file
git rm 50k_morex_v2_idt90.vcf
mv 50k_morex_v2_idt90_sorted.vcf 50k_morex_v2_idt90.vcf
# Index VCF
gatk IndexFeatureFile --input 50k_morex_v2_idt90.vcf
```

Cleanup intermediate files:

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
rm temp_50k*
rm 50k_morex_v2_idt90_unsorted*
```

Add in resolved snp utils problem snps:

```bash
# In dir: ~/Dropbox/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
cp 50k_morex_v2_idt90.vcf 50k_morex_v2_idt90_unsorted.vcf
tail -n +2 snp_utils_problem_snps/50k_problem_snps_resolved.vcf >> 50k_morex_v2_idt90_unsorted.vcf
# Sort VCF
module load gatk/4.1.2
gatk SortVcf -I 50k_morex_v2_idt90_unsorted.vcf -O 50k_morex_v2_idt90_sorted.vcf
# Rename file
mv 50k_morex_v2_idt90_sorted.vcf 50k_morex_v2_idt90.vcf
# Index VCF
gatk IndexFeatureFile --input 50k_morex_v2_idt90.vcf
```

Convert pseudomolecular positions to parts positions:

```bash
# In dir: ~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP
module load python2/2.7.12_anaconda4.2
module load gatk_ML/4.1.8
~/GitHub/File_Conversions/Pseudomolecules_to_Parts_v2.py --vcf 50k_morex_v2_idt90.vcf > 50k_morex_v2_idt90_parts.vcf
# Index parts vcf
gatk IndexFeatureFile --input 50k_morex_v2_idt90_parts.vcf
```
