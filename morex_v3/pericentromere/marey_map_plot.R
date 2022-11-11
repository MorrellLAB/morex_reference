#!/usr/bin/env Rscript

# Dependencies
library(gtools)
library(data.table)
library(tidyverse)

### User provided input arguments
# File containing the physical position of the pericentromere
pcent_phys_fp <- "~/GitHub/morex_reference/morex_v3/pericentromere/pericentromere_physPos.txt"
centromere_fp <- "~/GitHub/morex_reference/morex_v3/pericentromere/MorexV3_centromere_positions.tsv"
supp_table6_fp <- "~/GitHub/morex_reference/morex_v2/centromere/Supp_Table6-tpg2plantgenome2011080023-sup-0009.csv"
# Cross reference table between POPA (Munoz et al) and BOPA SNPs
crossref_fp <- "~/GitHub/morex_reference/morex_v2/centromere/Names_9K_BOPA_SNPs_CrossRef.csv"
# Vcf file with BOPA SNP position in new reference
vcf_fp <- "~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP/9k_idt95_noRescuedSNPs.vcf"
# Where do we want to output our files?
out_dir <- "~/GitHub/morex_reference/morex_v3/pericentromere"

#---------------------
# Read in data
centromere <- read.table(centromere_fp, header = F)
pericentromere <- read.table(pcent_phys_fp, header = F)
consensus_map <- read.csv(supp_table6_fp, skip = 2)
crossref <- read.csv(crossref_fp)
vcf <- read.table(vcf_fp)
# Make column names more meaningful
colnames(vcf) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
colnames(centromere) <- c("chr", "pos")
colnames(pericentromere) <- c("chr", "startpos", "endpos")

## STEP 1: Change SNP names in consensus map to BOPA names
# Find which SNPs have BOPA names
# 2994 markers in consensus_map df
# 7,864 markers in crossref df
Shared <- intersect(consensus_map$Marker, crossref$POPA)
# Subset data frames
crossref_sh <- crossref[crossref$POPA %in% Shared, ]
# Merge consensus map with cross reference SNP names
tmp_cmap_bopac_names <- merge(x = consensus_map, y = crossref_sh,
                              by.x = "Marker",
                              by.y = "POPA")
# Add physical positions from vcf
tmp_cmap_vcf <- merge(x = tmp_cmap_bopac_names, y = vcf,
                      by.x = "BOPA_C",
                      by.y = "ID")
# Select relevant columns for cleaner df
cmap_vcf_unsorted <- data.frame(chr = tmp_cmap_vcf$CHROM, pos = tmp_cmap_vcf$POS, id = tmp_cmap_vcf$BOPA_C, marker = tmp_cmap_vcf$Marker, cM_pos = tmp_cmap_vcf$cM)
cmap_vcf <- cmap_vcf_unsorted[order(cmap_vcf_unsorted$chr, cmap_vcf_unsorted$pos), ]

## STEP 2: Generate marey map plots
ggplot() +
  geom_rect(data=pericentromere, 
            aes(xmin=startpos, xmax=endpos, ymin=-Inf, ymax=Inf),
            fill = "grey", alpha=0.25) +
  geom_vline(data=centromere, aes(xintercept=pos)) +
  geom_point(data=cmap_vcf, aes(pos, cM_pos), alpha = 0.25) +
  facet_wrap(~chr) +
  theme_bw() +
  labs(x = "Physical Position", y = "Genetic distances (cM)")

# Additional exploration: Output a df with the necessary columns for MareyMapOnline:
# https://lbbe-shiny.univ-lyon1.fr/MareyMapOnline/
out_df <- data.frame(set = "Hordeum vulgare", map = cmap_vcf$chr, mkr = cmap_vcf$id,
                     phys = as.numeric(cmap_vcf$pos), gen = cmap_vcf$cM_pos)
out_fp <- paste0(out_dir, "/", "bopa_9k_snps_physPos_cM.txt")
write.table(x = out_df, file = out_fp, sep = " ", row.names = FALSE)
