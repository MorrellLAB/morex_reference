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

### Running the CDS extraction with my script
* This will create a CDS fasta file, which only needs to be run once!
Use `anchorwave_gff2seq.sh`


## Setup the query genome we want to map with AnchorWave
* [FB19_011_3](https://galaxy-web.ipk-gatersleben.de/libraries/folders/F6505e875ddb66fd2/dataset/1b02120528b47029) is a high-quality assembly of _Hordeum bulbosum_, the genome is diploid 
* Chromosomes are named, chr1H_1, chr1H_2 for the two haplotypes
* Testing with only the first haplotype, so cutting that out with awk
* Also need to remove suffixes from chromosome names so they 

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

## `gmap` to place CDS fasta positions on the reference and query genome
* minimap2 can apparently struggle with mapping of short exons
Use `gmap.sh`


## `AnchorWave` alignment, using our CDS anchors
Use `anchorwave_gff2seq.sh`

