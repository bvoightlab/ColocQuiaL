# make_significant_pairs_files.R
# By Brian Chen
# Takes eQTL Catalogue datasets, generates significant pairs files, then creates tabix index files
library(vroom)
library(qvalue)

# Command line argument: all.tsv.gz file
args = commandArgs(trailingOnly = TRUE)
SNP_file = args[1]

print(paste("Reading file:", SNP_file))
all_SNPs <- vroom(SNP_file)

# Calculate q-values
print("Calculating q-values")
qvalue <- qvalue(all_SNPs$pvalue, fdr.level = 0.05, pfdr = TRUE)

all_SNPs$qvalue = qvalue$qvalues
all_SNPs$significant = qvalue$significant

# Filter for signficant SNPs
print("Filtering for significant SNPs")
all_sig_SNPs = all_SNPs[all_SNPs$significant,]

# Write to file
out_file = gsub("all.tsv.gz", "signif_variant_gene_pairs.tsv.gz", SNP_file)
vroom_write(all_sig_SNPs, out_file)

# Create tabix file using create_tabix_eqtl_catalogue.sh script
system(paste0("bsub -o create_tabix_eqtl_catalogue.%J-%I.out  -e create_tabix_eqtl_catalogue.%J-%I.err -N 'bash create_tabix_eqtl_catalogue.sh ", out_file,"'"))

