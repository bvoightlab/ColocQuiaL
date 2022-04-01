#this R config file allows you to set the paths to all the dependency file you need to run the qtl colocalizer pipeline
# this config file is required for running the colocalizer pipeline on one locus

#provide the path to plink reference files to be used for plink commands
plink_bfile = "/project/voight_selscan/ksiewert/CardioMetaAnalysis/LDL_CHD_Bivar/LDClump/PlinkFilesOnlyRs/mergedBed"
#provide the path to the plink ped files for the list of individuals you wish to use in your LD reference panel
plink_keep = "/project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/AD_bivariate_scan_code/EUR.final.plink"

#provide path to GRCh37 to GRCh38 variant hash table directory
hash_table_dir = "/project/voight_GWAS/bychen9/human_9606_b151_GRCh38p7/BED/hash_tables_2/"

#provide path to eqtl tissue table
eQTL_tissue_table = "/project/voight_GWAS/bychen9/eQTL_colocalizer/GTEx_v8_Tissue_Summary_with_filenames.csv"
#provide path to sqtl tissue table
sQTL_tissue_table = "/project/voight_GWAS/bychen9/eQTL_colocalizer/GTEx_v8_sQTL_Tissue_Summary_with_filenames.csv"

#provide path to significant eqtl data tabix directory NOTE: don't include last slash
eQTL_sig_qtl_tabix_dir = "/project/voight_datasets/GTEx_v8/eQTL/GTEx_Analysis_v8_eQTL_tabix"
#column number of column in significant eqtl files that contains geneID
eQTL_sig_geneID_col = 7 

#provide path to the all eqtl data tabix directory
eQTL_all_qtl_tabix_dir = "/project/voight_datasets_01/GTEx_v8/TissueSpecific_tabix/"
#provide the header to the eqtl data
eQTL_all_header = c("chrom_b38", "chromStart_b38", "chromEnd_b38", "eGeneID", "A1_eqtl", "A2_eqtl", "build", "tss_distance", "ma_samples", "ma_count", "maf", "pvalue_eQTL", "slope", "slope_se")
#name of column in header representing geneID
eQTL_all_geneID = "eGeneID" 
#name of column in header representing chromosome
eQTL_all_chrom = "chrom_b38"
#name of column in header representing end coordinate 
eQTL_all_chromEnd = "chromEnd_b38"
#name of column in header representing p-value
eQTL_all_pvalue = "pvalue_eQTL"

#provide path to significant sqtl data tabix directory NOTE: don't include last slash
sQTL_sig_qtl_tabix_dir = "/project/voight_datasets/GTEx_v8/sQTL/GTEx_Analysis_v8_sQTL_tabix"
#column number of column in significant sqtl files that contains geneID
sQTL_sig_geneID_col = 12

#provide path to the all sqtl data tabix directory
sQTL_all_qtl_tabix_dir = "/project/voight_viz/bychen9/sQTL_tabix/"
#provide the header to the sqtl data
sQTL_all_header = c("chrom_b38", "chromStart_b38", "chromEnd_b38", "eGeneID", "intron_chr", "intron_bp_first", "intron_bp_end", "intron_clu", "A1_sqtl", "A2_sqtl", "build", "tss_distance", "ma_samples", "ma_count", "maf", "pvalue_sQTL", "slope", "slope_se")
#name of column in header representing geneID
sQTL_all_geneID = "eGeneID" 
#name of column in header representing chromosome
sQTL_all_chrom = "chrom_b38"
#name of column in header representing end coordinate 
sQTL_all_chromEnd = "chromEnd_b38"
#name of column in header representing p-value
sQTL_all_pvalue = "pvalue_sQTL"
#intron information
sQTL_all_intron_chr = "intron_chr" 
sQTL_all_intron_bp_first = "intron_bp_first"
sQTL_all_intron_bp_end = "intron_bp_end"
sQTL_all_intron_clu = "intron_clu"


#path to GRCh37 to GRCh38 liftOver chain file
liftOver_chain = "/appl/liftOver-20180423/chains/hg19ToHg38.over.chain"
#path to recombination rate directory (expected to be 1K genome chromosome recombination rate files or in the same format as these files: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/
recomb_rate_data="/project/voight_GWAS/wbone/recomb_data_dir/CEU/CEU"
