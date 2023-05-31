# Creating a data set for setup of ancestral state information

## Plan is to create est-sfs input
* Chaochih created a VCF of genotypes of the 284 samples from Fang et al. 2014 with Morex_v3 positions
* This is a great test data set!
* The VCF does not have any of the tags normally needed to parse of VCF, so I'm going to add them!


```bash
module load bcftools/1.10.2
module load htslib/1.9
module load bedtools/2.29.2
cd /panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/vcf_morex_v3
[//]: Adding tags, particularly AC and AN
bcftools +fill-tags WBDC_284_Bopa1and2_morex_v3.vcf.gz -Oz -o WBDC_284_Bopa1and2_morex_v3_tags.vcf.gz -- -t all


[//]: Extract positions from VCF
bcftools view -r chr5H:358793050-362657993 WBDC_284_Bopa1and2_morex_v3_tags.vcf.gz
[//]: Only 12_11151 & 12_30747 show up and 12_30747 is not polymorphic

[//]: Extract all positions to get a BED file to use with ancestral state files
bcftools query -f '%CHROM %POS %ID %REF %ALT\n' WBDC_284_Bopa1and2_morex_v3_tags.vcf.gz | awk -v OFS='\t' '{print $1, $2-1, $2, $3, $4, $5}' >BOPA_Morex_v3_positions.bed 
mv BOPA_Morex_v3_positions.bed ~/Workshop/Ancestral_State/ 

[//]: Covert positions to parts of chromosomes to use with ancestral state files
~/Apps/TKono/File_Conversions/Barley_Pseudomolecules_to_Parts.py --bed ~/Workshop/Ancestral_State/BOPA_Morex_v3_positions.bed morex_v3 >~/Workshop/Ancestral_State/BOPA_Morex_v3_positions_parts.bed            
[//]: The SNP names are lost in the conversion. Can add them back by backward conversion!


bedtools getfasta -fi /panfs/jay/groups/9/morrellp/shared/Datasets/Outgroups/morex_v3_outgroups_partsRef/bulbosum_A12_0.03.fa -bed BOPA_Morex_v3_positions_parts.bed -bedOut >BOPA_Morex_v3_positions_parts_bulbosum.bed

bedtools getfasta -fi /panfs/jay/groups/9/morrellp/shared/Datasets/Outgroups/morex_v3_outgroups_partsRef/murinum_BCC2017_0.03.fa -bed BOPA_Morex_v3_positions_parts.bed -bedOut >BOPA_Morex_v3_positions_murinum.bed

 bedtools getfasta -fi /panfs/jay/groups/9/morrellp/shared/Datasets/Outgroups/morex_v3_outgroups_partsRef/pubiflorum_BCC2028_0.05.fa -bed BOPA_Morex_v3_positions_parts.bed -bedOut >BOPA_Morex_v3_positions_pubiflorum.bed

cd /panfs/jay/groups/9/morrellp/pmorrell/Workshop/Ancestral_State 

cp WBDC_284_Bopa1and2_morex_v3_tags.vcf.gz ~/Workshop/Ancestral_State/ 
cd ~/Workshop/Ancestral_State/
gunzip WBDC_284_Bopa1and2_morex_v3_tags.vcf.gz

~/Apps/TKono/File_Conversions/Barley_Pseudomolecules_to_Parts.py --vcf WBDC_284_Bopa1and2_morex_v3_tags.vcf morex_v3 >WBDC_284_Bopa1and2_morex_v3_tags_parts.vcf

```

* The files needed for ancestral state inference are here:
`/panfs/jay/groups/9/morrellp/pmorrell/Workshop/Ancestral_State`
WBDC_284_Bopa1and2_morex_v3_tags_parts.vcf 
BOPA_Morex_v3_positions_bulbosum.bed
BOPA_Morex_v3_positions_murinum.bed
BOPA_Morex_v3_positions_pubiflorum.bed
* For the VCF, we have a few ways to get the non-reference allele count, simplest is perhaps: AC (Allele counts in genotypes)
* The reference allele count is then | AN - AC | (where AN = Total number of alleles in called genotypes)
* For our purposes, we might want a flag for "haploid" that divides those values by 2 because our sample is highly homozygous
* So for the 1st SNP in this dataset, the reference is G, alternate is A, AC = 90, AN = 568
* So | (568/2) - (90/2) | = 239
* The genotypes would be "45 0 239 0" represented in the order A C G T for counts of alleles 

* BOPA SNP 12_11151: 
`chr5H_part2	153932289	12_11151	G	C`
"C" is the minor allele and "G" is present in the reference and all three outgroups, so "G" has to ancestral