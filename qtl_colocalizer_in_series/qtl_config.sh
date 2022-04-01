#ID for user's trait of interest. (Can be any string)
trait="BRCAN"
#path to the input files
traitFilePath="/project/voight_datasets/GWAS/18_brcan/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics_jun20_withrsIDs.txt"
#column IDs from trait file
trait_A1col="Effect.Gwas"
trait_A2col="Baseline.Gwas"
trait_SNPcol="rsid"
trait_CHRcol="chr.iCOGs"
trait_BPcol="Position.iCOGs"
trait_Pcol="P.value.Gwas" 
trait_Ncol="NumCalled.iCOGs" 
trait_MAFcol="Freq.Gwas"
#trait info not in the input file
#traitType is set either to "cc" or "quant"
traitType="cc"
#This is the proportion of samples that are cases in a case control GWAS, if you are using a quantitative trait this should be set to "". traitProp = cases / case + controls
traitProp=0.536671 #look this up
#path to lead SNPs file, set as empty string to use plink to identify lead SNPs from the GWAS summary stat file
leadSNPsFilePath=""
#qtlType is set to "eqtl" or "sqtl"
qtlType="eqtl"
#colc window size in bps, this is only required when no leadSNPsFilePath is set. If you are providing a leadSNPsFile leave the empty string field for window
window=""
#plink parameters: Only need to edit these if you wish to change from the default values of 0.0000001, 1000, and 0.2
clumpP1=""
clumpKB=""
clumpR2=""

#provide path to the setup_config.sh (this is the config file with all of the file paths to dependency files such as plink reference files, GTEx eQTL data, etc.
setup_config_sh="/project/voight_GWAS/wbone/qtl_colocalizer/eQTL_colocalizer/setup_config.sh"
#provide the path to plink refernce files to be used for plink commands
#plink_bfile="/project/voight_selscan/ksiewert/CardioMetaAnalysis/LDL_CHD_Bivar/LDClump/PlinkFilesOnlyRs/mergedBed"
#provide the path to the plink ped files for the list of individuals you wish to use in your LD reference panel
#plink_keep="/project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/AD_bivariate_scan_code/EUR.final.plink"
#provide path to GRCh37 to GRCh38 variant hash table directory
#hash_table="/project/voight_GWAS/bychen9/human_9606_b151_GRCh38p7/BED/hash_tables_2"
#provide path to qtl tissue table
#tissue_table="/project/voight_GWAS/bychen9/eqtl_coloc/GTEx_v8_Tissue_Summary_with_filenames.csv"
#provide path to significant qtl data tabix directory
#sig_qtl_tabix_dir="/project/voight_datasets/GTEx_v8/sQTL/GTEx_Analysis_v8_sQTL_tabix"
#provide path to the all qtl data tabix directory
#all_qtl_tabix_dir="/project/voight_datasets_01/GTEx_v8/TissueSpecific_tabix/"
#path to GRCh37 to GRCh38 liftOver chain file
#lifOver_chain="/appl/liftOver-20180423/chains/hg19ToHg38.over.chain"

