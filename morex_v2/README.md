# Morex v2

This directory contains scripts and submodules for file processing/creation that are necessary for using the reference genome.

Preprint: [Monat C, Padmarasu S, Lux T, Wicker T, Gundlach H, Himmelbach A, Ens J, Li C, Muehlbauer GJ et al. 2019. TRITEX: chromosome-scale sequence assembly of Triticeae genomes with open-source tools. BioRxiv 631648.](https://www.biorxiv.org/content/10.1101/631648v1)

Morex v2 reference genome can be accessed through e!DAL - Plant Genomics & Phenomics Research Data Repository: https://doi.org/10.5447/IPK/2019/8

---

Peter L. Morrell
6 June 2019
Falcon Heights, MN

## Source for reference
- The reference genome and gene model GFF files are at the doi below:
https://doi.org/10.5447/IPK/2019/8

- The Morex v2 reference is named:
Barley_Morex_V2_pseudomolecules.fasta.gz

## Also adding plastid genomes to Morex v2
- The chloroplast genome from Morex is at the link below:
https://www.ncbi.nlm.nih.gov/nuccore/EF115541.1
- The mitochondrial genome from Haruna Nijo (there is no published version from Morex) is at the link below:
https://www.ncbi.nlm.nih.gov/nuccore/AP017301.1/
- Appended these files to the reference in a new file using `gzip -c plastid.fasta >> Barley_Morex_V2_pseudomolecules_plastids.fasta.gz`
- minimap2 -d was used to create an index; required a qsubbed job: https://github.com/pmorrell/Utilities/blob/master/genome_indexing.sh
- Morex v2 comes with a GFF (below), will need to convert plastid genome annotations to GFF

Skylar R. Wyant
7 June 2019
Falcon Heights, MN

## Breaking chromosomes into parts by chromosome arm
- Created .fai index using `samtools faidx` (samtools 1.7)
- Downloaded the centromere positions from the link below:
https://doi.ipk-gatersleben.de/DOI/83e8e186-dc4b-47f7-a820-28ad37cb176b/7ffd0036-13e6-49be-89fd-27c0e77b0ae0/0/
- Made bed file by hand using centromere positions from above and chromosome lengths from .fai index (parts_with_plastids.bed) 
- Broke up reference with plastids into chromosome arms using the following command:
`bedtools getfasta -name -fi Barley_Morex_V2_pseudomolecules_plastids.fasta -fo Barley_Morex_V2_pseudomolecules_parts_plastids.fasta -bed parts_with_plastids.bed`

