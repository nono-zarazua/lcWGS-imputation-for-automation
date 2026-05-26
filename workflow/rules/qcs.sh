bcftools stats -s - quilt.down0x.mspbwt.genome.vcf.gz > sample_stats.txt

# 1. Create the file and write the correct headers
echo -e "Sample_ID\tnRefHom\tnHomAlt(1/1)\tHeterozygous(0/1)\tRatio(Het/Hom-Alt)" > clean_ratios.tsv

# # 2. Extract the data, calculate the ratio, and append it to the file
grep "^PSC" sample_stats.txt | awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%.2f\n", $3, $4, $5, $6, $6/$5}' >> clean_ratios.tsv

# Create a temporary dataset with fixed IDs
# CHR:POS:REF:ALT for VCF ID col
plink2 --vcf quilt.down0x.mspbwt.genome.vcf.gz \
  --set-all-var-ids '@:#:$r:$a' \
  --make-pgen \
  --out temporary_named_batch \
  --autosome

# Extract variants 
plink2 --pfile temporary_named_batch \
  --extract /home/ec2-user/workdir/project-quilt-workdir/data/pruned/qc_pruning.prune.in \
  --make-pgen \
  --out clean_batch_for_pca

# Clean data set
plink2 --pfile clean_batch_for_pca --rm-dup exclude-all --make-pgen --out clean_batch_dedup

# PCA
plink2 --pfile clean_batch_dedup --read-freq /home/ec2-user/workdir/project-quilt-workdir/data/pruned/ref_panel_freqs.afreq --pca --out clean_pca

# Kinship
plink2 --pfile clean_batch_dedup --make-king-table --out clean_kinship
