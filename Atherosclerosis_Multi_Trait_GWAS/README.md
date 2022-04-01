# eQTL_colocalizer

This repository contains the code needed to generate dependency files and run the eQTL colocalizer pipeline. Given a GWAS signal of your choosing, this pipeline will run COLOC on your signal and all the GTEx v8 single-tissue eQTLs.

This directory contains the version of this code used in Bellomo et al. (https://www.medrxiv.org/content/10.1101/2021.05.21.21257493v1).

**Preparing the pipeline:**
- Download files:
  - Download GRCh38 dbSNP BED files from here: https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/BED/
  - Download all SNPs tissue specific GTEx v8 files from the GTEx website: https://www.gtexportal.org/home/datasets
- Generate dependency files:
  - Generate tabix index files for the GTEx v8 significant pair files and all pairs files
  - Run dbsnp_hash_table_maker.py to create hash tables from the GRCh38 dbSNP BED files.
  - Make sure you also have GTEx_Tissue_Summary_with_filenames.csv and update the paths to these files in eqtl_colocalizer.R

**Running the pipeline:**
- Create an analysis directory, and add eqtl_colocalizer.R and the eQTL_config.R file modified from eQTL_config_template.R to correspond to the GWAS signal on which you would like to perform eQTL colocalization analysis. 
- Then, execute the eqtl_colocalizer.R in that directory (eg Rscript ./eqtl_colocalizer.R)
