# Identifying the pericentromeric regions for Barley Morex v3

The general methods follow the approaches previously used for Morex v1 (https://github.com/lilei1/9k_BOPA_SNP/tree/master/centromere_barley) and Morex v2 (https://github.com/MorrellLAB/morex_reference/tree/master/morex_v2/centromere).


---

## Files and Information Needed

- Consensus genetic map from **Supplementary Table 6** from [Munoz-Amatriain et al. 2011 The Plant Genome](https://doi.org/10.3835/plantgenome2011.08.0023). Here we store this info in a file called `Supp_Table6-tpg2plantgenome2011080023-sup-0009.csv`.
- Pericentromeric regions (cM) for chr2H to 7H from **Table 6** in [Munoz-Amatriain et al. 2011 The Plant Genome](https://doi.org/10.3835/plantgenome2011.08.0023). Here we store this info in a file called `Table6_pericentromeres_cM.txt`.
- [Morex v2 9k and BOPA SNP positions](https://github.com/MorrellLAB/morex_reference/blob/master/morex_v2/50k_9k_BOPA_SNP/9k_morex_v2_idt90.vcf)
- SNP alternate names in the file `Names_9K_BOPA_SNPs_CrossRef.csv` (which is the first sheet from the [Names_9K_BOPA_SNPs_XREF.xlsx](https://github.com/lilei1/9k_BOPA_SNP/blob/master/centromere_barley/Names_9K_BOPA_SNPs_XREF.xlsx) file).
- Morex v3 centromere positions were downloaded from: [M. Mascher (2020-10-13): Pseudomolecules and annotation of the third version of the reference genome sequence assembly of barley cv. Morex [Morex V3]. DOI:10.5447/ipk/2021/3](https://doi.ipk-gatersleben.de/DOI/b2f47dfb-47ff-4114-89ae-bad8dcc515a1/7eb2707b-d447-425c-be7a-fe3f1fae67cb/2)

## Methods

We will run the script `get_pericentromere_positions.R` that was modified from Ana Poet's 2017 version [here](https://github.com/lilei1/9k_BOPA_SNP/blob/master/script/PositionCentromeres.R). This script will take the intersection between Table 6 and Supp Table 6 [Munoz-Amatriain et al. 2011 The Plant Genome](https://doi.org/10.3835/plantgenome2011.08.0023) and output:

1. Pericentromere physical positions (`pericentromere_physPos.txt`)
2. VCF file with SNPs that are in the pericentromeric region (`9k_snps_pericentromere.vcf`).

Files that contain the genetic map positions haven't changed and are located in the morex_v2 subdirectory.

Run the script to get Morex v3 pericentromere positions:

```bash
# In dir: ~/GitHub/morex_reference/morex_v3/pericentromere
module load R/4.0.4
./get_pericentromere_positions.R
```

Pericentromere positions are also listed in the table below:

**Note:** chr6H pericentromere positions doesn't make sense, it doesn't overlap the centromere. The approach used here doesn't seem to work well for chr6H. chr6H will need to be determined another way.

| Chr | Start | End |
| --- | ----- | --- |
| chr1H | NA | NA |
| chr2H | 201406208 | 359516223 |
| chr3H | 97202374 | 397494503 |
| chr4H | 52116492 | 496024940 |
| chr5H | 54152588 | 328984081 |
| chr6H | ~~343875085~~ | ~~345879621~~ |
| chr7H | 192266334 | 434396272 |

Centromere positions:

| Chr | Centromere |
| chr1H | 206486643 |
| chr2H | 301293086 |
| chr3H | 267852507 |
| chr4H | 276149121 |
| chr5H | 204878572 |
| chr6H | 256319444 |
| chr7H | 328847192 |

Currently, we are not able to figure out positions for chr1H. Previously, for Morex v1 Tom found the pericentromere positions (including chr1H) by using Supplemental Table 4.4 from [Mascher et al 2017 Nature](https://www.nature.com/articles/nature22043). We may consider a similar approach to Casale et al. 2021 in defining chr1H pericentromere positions but we'll need to explore this further first.
