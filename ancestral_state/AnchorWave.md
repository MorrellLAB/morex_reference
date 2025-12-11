# Anchorwave for barley genome alignments

## The first step is to create a CDS fasta file that Anchorwave likes!

* Use a custom sort order: gene first, then mRNA, then others

```bash
cd /scratch.global/pmorrell/Ancestral
GFF3="/projects/standard/morrellp/pmorrell/Shared/References/Reference_Sequences/Barley/Morex_v3/gene_annotation/Hv_Morex.pgsb.Jul2020.sorted.gff3.gz"
pigz -dc $GFF3 | awk -F'\t' 'BEGIN{OFS="\t"}
  /^#/ {print; next}
  {
    # Add sort key: 1=gene, 2=mRNA, 3=everything else
    if ($3 == "gene") key = 1
    else if ($3 == "mRNA") key = 2
    else key = 3
    print key "\t" $0
  }' | sort -t$'\t' -k2,2 -k5,5n -k1,1n | cut -f2- > Hv_Morex.pgsb.Jul2020.re_sorted.gff3
```

### Running the CDS extraction with bash script
* This will create a CDS fasta file, which only needs to be run once for Morex_v3
Use [anchorwave_gff2seq.sh](https://github.com/MorrellLAB/morex_reference/blob/master/ancestral_state/anchorwave_gff2seq.sh)


## Setup the query genome we want to map with AnchorWave
* [FB19_011_3](https://galaxy-web.ipk-gatersleben.de/libraries/folders/F6505e875ddb66fd2/dataset/1b02120528b47029) is a high-quality assembly of _Hordeum bulbosum_, the genome is diploid 
* Chromosomes are named, chr1H_1, chr1H_2 for the two haplotypes
* Testing with only the first haplotype, so cutting that out with awk
* For AnchorWave, we'll need chromosomes names to match in the two assemblies, so chr1H, chr2H in both, for example

```bash
module load pigz/2.7-gcc-13.1.0-bwi6ag3 

cd /scratch.global/pmorrell/Ancestral/FB19_011_3

awk '
/^>chr[0-9]+H_1$/ {
    printing=1
}
/^>/ && !/^>chr[0-9]+H_1$/ {
    printing=0
}
printing {
    print
}
' <(pigz -cd FB19_011_3.final.fasta.gz.fasta.gz) > FB19_011_3_hap1.fasta

[//]: Remove the suffix so the chromosome names will match 
sed -e 's/_1//g' FB19_011_3_hap1.fasta >temp.fasta
mv temp.fasta FB19_011_3_hap1.fasta
```

* A similar setup can be used for haplotype 2 of this assembly or any other barley assembly

## minimap2 to place CDS fasta positions on the reference and query genome
* minimap2 can apparently struggle with mapping of short exons, but the results from [minimap2_splice.sh](https://github.com/MorrellLAB/morex_reference/blob/master/ancestral_state/minimap2_splice.sh)
* also tried gmap, see [gmaph.sh](https://github.com/MorrellLAB/morex_reference/blob/master/ancestral_state/gmap.sh)

## `AnchorWave` alignment, using our CDS anchors
* The current settings test for inversions within chromosomes, but don't compare chromosomes, so we won't find translocations
Use [anchorwave_genoAli.sh](https://github.com/MorrellLAB/morex_reference/blob/master/ancestral_state/anchorwave_genoAli.sh)

## Converting the MAF directly to BED input files for ESTSFS_to_ancestral.py
* We can use the alignment AnchorWave produces to query the nucleotide state at positions in Morex_v3
* This can feed [ESTSFS_to_ancestral.py](https://github.com/MorrellLAB/File_Conversions/blob/master/ESTSFS_to_ancestral.py)
* Currently [MAF_to_nucleotide.py](https://github.com/MorrellLAB/morex_reference/blob/master/ancestral_state/MAF_to_nucleotide.py) will output "N" for missing and "-" for deletion
* That behavior may need to change to match required ESTSFS_to_ancestral.py input
