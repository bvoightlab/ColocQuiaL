# ColocQuiaL

This repository contains the code needed to generate dependency files and run ColocQuiaL. Given GWAS signals of your choosing, this pipeline will run COLOC on your signals and all eQTL or sQTL summary statistics that you have available to you. For example, this pipeline will work with the single-tissue eQTL or sQTL datasets available from GTEx or eQTL Catalogue.

## Preparing the pipeline
###  Download files
- Download GRCh38 dbSNP BED files from here: https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/BED/

- Download tissue specific eQTL or sQTL dataset. For example, all SNPs tissue specific GTEx v8 files can be downloaded from the GTEx website at https://www.gtexportal.org/home/datasets and additional QTL datasets from eQTL Catalogue can be downloaded at http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/.
 
### Generate dependency files
- Run ``dbsnp_hash_table_maker.py`` to create hash tables from the GRCh38 dbSNP BED files in the directory containing the files.
  ```
  python dbsnp_hash_table_maker.py
  ```

- LD and recombination rate reference files: need to have merge plink files (.bed, .bim, .fam) files, a filie of which samples to use, and recombination rates to generate the regional association.
  - Recombination rate files can be downloaded from here: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates

#### For GTEx datasets
- Generate tabix index files for the significant pair files and all pairs files. Examples of how to generate the files for GTEx datasets are provided in ``create_tabix_gtex_eqtl_allpairs.sh``, ``create_tabix_gtex_eqtl_sigpairs.sh``, ``create_tabix_gtex_sqtl_allpairs.sh``, and ``create_tabix_gtex_sqtl_sigpairs.sh``.
  ```
  bash create_tabix_gtex_eqtl_sigpairs.sh Adipose_Subcutaneous.v8.signif_variant_gene_pairs.txt
  ```

- Create a tissue summary CSV file (e.g., ``GTEx_v8_Tissue_Summary_with_filenames.csv`` and ``GTEx_v8_sQTL_Tissue_Summary_with_filenames.csv``) containing the tissue names, sample sizes, and file names of the tabix files corresponding to your QTL dataset.

#### For eQTL Catalogue datasets
- Create significant pair files from the downloaded all pairs files. An example of how to create significant pair files is provided in ``make_significant_pairs_files.R``
  ```
  Rscript make_significant_pairs_files.R GENCORD_ge_fibroblast.all.tsv.gz
  ```

- Generate tabix index files for the significant pair files. The tabix index files for the all pair files are provided by eQTL Catalogue. An example of how to generate the files for eQTL Catalogue datasets is provided in ``create_tabix_eqtl_catalogue.sh``. Note that the example ``make_significant_pairs_files.R`` from the previous step does this step automatically.

- For each QTL dataset, create a tissue summary CSV file containing the tissue names, sample sizes, and file names of the tabix files corresponding to the dataset (Similar to ``GTEx_v8_Tissue_Summary_with_filenames.csv`` and ``GTEx_v8_sQTL_Tissue_Summary_with_filenames.csv``).

- Create a configuration file for each downloaded QTL dataset based on ``setup_config.R``.

#### For other datasets
Assuming the file types are similar to the GTEx or eQTL Catalogue datasets, the overall setup structure should be similar: 
1. If necessary, create significant pair files from all pairs files. Note that this step was not required for GTEx as significant pair files are provided. 
2. Generate tabix index files for the significant pair files and all pairs files. 
3. For each QTL dataset, create a tissue summary CSV file (e.g., ``GTEx_v8_Tissue_Summary_with_filenames.csv`` and ``GTEx_v8_sQTL_Tissue_Summary_with_filenames.csv``) containing the tissue names, sample sizes, and file names of the tabix files corresponding to the dataset.
4. Create a configuration file for each downloaded QTL dataset based on ``setup_config.R``.

### Modify path to files
- In ``setup_config.sh``, update ``colocquial_dir`` to the path to the directory containing ``colocquial.R``, ``qtl_coloc_template.bsub``, and ``QTL_config_template.R``. 

- Make sure all paths in ``setup_config.sh`` and ``setup_config.R`` files are correct. 

## Running the pipeline
Create an analysis directory, and add ``colocquial_wrapper.sh`` and a ``qtl_config.sh`` file modified to correspond to the GWAS signals on which you would like to perform eQTL or sQTL colocalization analysis.

- The provided example ``qtl_config.sh`` file provides further detail on what is required in this file and what the expected entries are for user defined parameters.

- In ``qtl_config.sh``, make sure to set ``setup_config_sh`` and ``setup_config_R`` to the configuration files of the dataset you want to use. 

Then, execute ``colocquial_wrapper.sh`` in that directory.
```
bash colocquial_wrapper.sh
```

### Running the pipeline on a single locus
Create an analysis directory, and add ``colocquial.R`` and a ``QTL_config.R`` file modified from ``QTL_config_template.R`` to correspond to the GWAS signal on which you would like to perform eQTL or sQTL colocalization analysis.

Import the necessary tools.
```
module load R/3.6.3
module load plink/1.90Beta4.5
module load tabix
module load liftOver
```

Then, execute ``colocquial.R`` in that directory.
```
Rscript colocquial.R
```

## Additional Notes
- For entries in configuration files that are not pertinent to your available data, set to empty string.

- The pipeline is currently set up to work with Pulit et al. BMI data. Directions on how to download the Pulit et al. data and example ``qtl_config.sh`` and lead SNP files are provided in the ``test_data/`` directory.

- The version of this code used in Bellomo et al. (https://doi.org/10.3389/fgene.2021.787545) is in the ```Atherosclerosis_Multi_Trait_GWAS/``` directory.

- Previous versions of ColocQuiaL are located here: https://github.com/bychen9/eQTL_colocalizer 
