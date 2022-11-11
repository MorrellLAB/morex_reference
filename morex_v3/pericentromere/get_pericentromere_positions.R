#!/usr/bin/env Rscript

# This script is based on Ana Poet's approach to finding the pericentromeric positions.
# See: https://github.com/lilei1/9k_BOPA_SNP/blob/master/script/PositionCentromeres.R
# Morex v3

# Dependencies
library(gtools)
library(data.table)
library(tidyverse)

### User provided arguments
# List of SNPs in centromeric region. Taken from supplemental table in Munoz et al 2011.
pericentromere_fp <- "~/GitHub/morex_reference/morex_v2/centromere/Table6_pericentromeres_cM.txt"
supp_table6_fp <- "~/GitHub/morex_reference/morex_v2/centromere/Supp_Table6-tpg2plantgenome2011080023-sup-0009.csv"
# Cross reference table between POPA (Munoz et al) and BOPA SNPs
crossref_fp <- "~/GitHub/morex_reference/morex_v2/centromere/Names_9K_BOPA_SNPs_CrossRef.csv"
# Vcf file with BOPA SNP position in new reference
vcf_fp <- "~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP/9k_idt95_noRescuedSNPs.vcf"
# Where do we want to output our files?
out_dir <- "~/GitHub/morex_reference/morex_v3/pericentromere"

#-----------
# Read in data
pericentromere <- read.table(pericentromere_fp, header = F)
consensus_map <- read.csv(supp_table6_fp, skip = 2)
crossref <- read.csv(crossref_fp)
vcf <- read.table(vcf_fp)
# Make column names more meaningful
colnames(vcf) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

# STEP 0: Prepare data frames

# Subset VCF so that we only use SNP positions that had unique BLAST hits.
# Some SNPs hit in multiple places non-uniquely and can't be used, we'll remove these.
# Duplicate SNP hits all contain the pattern "ALTCHR" in the info field
vcf.sub <- vcf[!vcf$INFO %like% "ALTCHR", ]
# Add "chr" to beginning of string to match vcf
consensus_map$Chromosome.1 <- gsub(pattern="^", replacement="chr", x=consensus_map$Chromosome.1)
# Pull out markers in the pericentromere region
pcentromere <- data.frame(Marker = character(),
                         Chr = character(),
                         cM = numeric())
for (i in pericentromere$V1) {
    if (i == "chr1H") {
        next # "NA" for chr1H, so we can't figure out pericentromere here
    } else {
        start_pos <- pericentromere[pericentromere$V1 == i, 2]
        end_pos <- pericentromere[pericentromere$V1 == i, 3]
        if (start_pos == end_pos) {
            temp_chr <- consensus_map[consensus_map$Chromosome.1 == i, ]
            temp_pcent <- temp_chr[grep(start_pos, temp_chr$cM), ]
        } else {
            temp_pcent <- consensus_map[consensus_map$Chromosome.1 == i & consensus_map$cM >= start_pos & consensus_map$cM <= end_pos, ]
        }
        # Pull out only relevant columns
        current_chr <- data.frame(Marker = temp_pcent$Marker,
                                  Chr = temp_pcent$Chromosome.1,
                                  cM = temp_pcent$cM)
        # Combine with main dataframe
        pcentromere <- rbind(pcentromere, current_chr)
    }
}

# STEP 1: change centromere SNP names to BOPA names

# Find which SNPs have BOPA names
# 129 markers in pcentromere df
# 7,864 markers in crossref df
Shared <- intersect(pcentromere$Marker, crossref$POPA) # 182 remaining in Shared
# Subset data frames
cen_sh <- pcentromere[pcentromere$Marker %in% Shared, ]
crossref_sh <- crossref[crossref$POPA %in% Shared, ]
# Order cross ref SNPs as in CENTROMERE
crossref_sh_or <- crossref_sh[match(cen_sh$Marker, crossref_sh$POPA), ]
# Merge to get all info together
pcent_bopa <- merge(x = cen_sh, y = crossref_sh_or,
                    by.x = "Marker",
                    by.y = "POPA")

# STEP 2: Find position of BOPA SNPs in the new reference genome

vcf_pcent <- vcf.sub[vcf.sub$ID %in% pcent_bopa$BOPA_C | vcf.sub$ID %in% pcent_bopa$BOPA, ] # 118 are present
# Sort VCF by chromosome then by position
vcf_pcent_ordered <- vcf_pcent[order(vcf_pcent$CHROM, vcf_pcent$POS), ]

# Find max and min of pericentromere positions
pcent_table <- vcf_pcent_ordered %>%
    group_by(CHROM) %>%
    summarise(min_pos = min(POS),
              max_pos = max(POS))

# "NA" for chr1H positions
final_table <- rbind(tibble(CHROM = "chr1H", min_pos = "NA", max_pos = "NA"), pcent_table)

# STEP 3: Save to output files

# Save table to file
write.table(final_table,
            paste(out_dir, "pericentromere_physPos.txt", sep = "/"),
            quote=F,
            row.names=F,
            col.names=F,
            sep="\t")

# Also save vcf containing SNPs in the pericentromeric region identified here
write.table(vcf_pcent_ordered,
            paste(out_dir, "9k_snps_pericentromere.vcf", sep = "/"),
            quote=F,
            row.names=F,
            col.names=F,
            sep="\t")
