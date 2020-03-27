# Temporary subdirectory for resolving duplicate SNPs

We will reorganize files here later.

### Some useful files

- SNPMeta output file for 9k set: https://conservancy.umn.edu/bitstream/handle/11299/181367/Barley_SNP_Annotations.txt?sequence=5&isAllowed=y

### Methods

To resolve duplicates:

**Step 1:** First search the [SNP meta file for the 9k SNPs](https://conservancy.umn.edu/bitstream/handle/11299/181367/Barley_SNP_Annotations.txt?sequence=5&isAllowed=y) and figure out which gene the SNP is located.

The purpose for this is to reduce noise since we know that the 9k iSelect genotyping was designed from exome capture samples.

If you find the SNP in the SNP meta file, copy the ProteinID for that SNP and go to NCBI https://www.ncbi.nlm.nih.gov/protein and search that protein. This will take you to a page where you can get the fasta sequence for that ProteinID.

**Step 2:** Take the fasta sequence and use the [IPK Barley BLAST server](https://webblast.ipk-gatersleben.de/barley_ibsc/) to do a BLAST search.

If you found the fasta sequence from the NCBI protein database, use the following options:

- Program: `blastp`
- Database(s): `Barley AA (HC and LC) Morex v2.0 2019`

**Note:** Directly pasting the fasta sequence from NCBI may introduce weird spaces, make sure you remove these before searching. I also noticed through experimenting between blastn and blastp that blastp seems to be a little more accurate than blastn.

Now, click `Basic search`.

**Step 3:** Pick the best BLAST hit.

If we can't find the gene that the SNP hit or the gene BLAST search results have multiple identical results, choose the left most position (smallest chromosome and smallest physical position) and put notes for the other identical hits in the last field of the VCF file.
