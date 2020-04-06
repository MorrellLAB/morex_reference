# SNPMeta output files

[`bopa_and_9k_snpmeta_output.txt`](http://conservancy.umn.edu/bitstream/handle/11299/181367/Barley_SNP_Annotations.txt?sequence=5&isAllowed=y) file was from https://conservancy.umn.edu/handle/11299/181367.

`50k_snpmeta_output.txt` file was generating using the following commands:

Run [SNPMeta](https://github.com/MorrellLAB/SNPMeta) to generate metadata for each SNP so we can reduce noise when resolving duplicate SNPs.

```bash
# In dir: ~/Shared/Datasets/Genotyping/Contextual_Sequences/barley_50k
# Prepare file for SNPMeta
module load python3/3.7.1_anaconda
~/GitHub/morex_reference/morex_v2/50k_9k_BOPA_SNP/convert_to_fasta.py snp_utils_50k_lookup_table_filtered.txt > snpmeta_50k_snps.fasta

# Run SNPMeta on duplicates only first
qsub run_dups_snpmeta_50k_snps.job

# Then, run SNPMeta on remaining SNPs
qsub run_snpmeta_50k_snps.job
```

Re-submit `run_snpmeta_50k_snps.job` script until SNPMeta has pulled all annotations. The script is written to check which snps are already completed and will only run those that have not been run yet. So, we don't need to modify the script.
