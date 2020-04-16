# Identifying the centromeric and pericentromeric regions for Barley Morex v2

The general methods follow the approach previously used for Morex v1: https://github.com/lilei1/9k_BOPA_SNP/tree/master/centromere_barley with a few updated scripts and files.

---

## Files and Information Needed

- Consensus genetic map from **Supplementary Table 6** from [Munoz-Amatriain et al. 2011 The Plant Genome](https://doi.org/10.3835/plantgenome2011.08.0023). Here we store this info in a file called `Supp_Table6-tpg2plantgenome2011080023-sup-0009.csv`.
- Pericentromeric regions (cM) for chr2H to 7H from **Table 6** in [Munoz-Amatriain et al. 2011 The Plant Genome](https://doi.org/10.3835/plantgenome2011.08.0023). Here we store this info in a file called `Table6_pericentromeres_cM.txt`.
- [Morex v2 9k and BOPA SNP positions](https://github.com/MorrellLAB/morex_reference/blob/master/morex_v2/50k_9k_BOPA_SNP/9k_morex_v2_idt90.vcf)
- SNP alternate names in the file `Names_9K_BOPA_SNPs_CrossRef.csv` (which is the first sheet from the [Names_9K_BOPA_SNPs_XREF.xlsx](https://github.com/lilei1/9k_BOPA_SNP/blob/master/centromere_barley/Names_9K_BOPA_SNPs_XREF.xlsx) file).

## Methods

We will run the script `get_centromere_positions.R` that was modified from Ana Poet's 2017 version [here](https://github.com/lilei1/9k_BOPA_SNP/blob/master/script/PositionCentromeres.R). This script will take the intersection between Table 6 and Supp Table 6 [Munoz-Amatriain et al. 2011 The Plant Genome](https://doi.org/10.3835/plantgenome2011.08.0023) and output:

1. Pericentromere physical positions (`pericentromere_physPos.txt`)
2. VCF file with SNPs that are in the pericentromeric region (`9k_snps_pericentromere.vcf`).

Pericentromere positions are also listed in the table below:

| Chr | Start | End |
| --- | ----- | --- |
| chr1H | NA | NA |
| chr2H | 203689667 | 365198431 |
| chr3H | 98960757 | 401158614 |
| chr4H | 54094751 | 507869478 |
| chr5H | 53760064 | 332514770 |
| chr6H | 349655750 | 351670614 |
| chr7H | 199910025 | 436392291 |

We are not able to figure out positions for chr1H. Previously, for Morex v1 Tom found the pericentromere positions (including chr1H) by using Supplemental Table 4.4 from [Mascher et al 2017 Nature](https://www.nature.com/articles/nature22043).

As a check, I compared the inferred positions above to Fig S2 in [Monat et al 2019 Genome Biology](https://doi.org/10.1186/s13059-019-1899-5). All chromosomes have pericentromere positions within the same ballpark as shown in Fig S2, except for chr6H. Based on the figure chr6H should be around 80 to 450 Mb. My inference for chr6H was around 349 to 351 Mb. The positions for chr6H will be taken as tentative positions for now.
