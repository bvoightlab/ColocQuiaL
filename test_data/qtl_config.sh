#ID for user's trait of interest. (Can be any string)
trait="BMI"
#path to the input GWAS summary stats files
traitFilePath="/project/ritchie21/personal/wbone/ColocQuiaL_test/test_runs_dir/BMI_ColocQuiaL_test/BMI_Pulit_data_rsID_ready.txt"
#column IDs from trait file
trait_A1col="Tested_Allele"
trait_A2col="Other_Allele"
trait_SNPcol="SNP"
trait_CHRcol="CHR"
trait_BPcol="POS"
trait_Pcol="P" 
trait_Ncol="N" 
trait_MAFcol="Freq_Tested_Allele"

#trait info not in the input file
#traitType is set either to "cc" (case-control) or "quant" (quantitative) depending on the kind of trait the GWAS is.
traitType="quant"

#This is the proportion of samples that are cases in a case control GWAS, if you are using a quantitative trait this should be set to "". traitProp = cases / case + controls
traitProp=""

#path to lead SNPs file. This is a file with information on the lead SNP and the start and stop position for the colocalization analysis at each locus. Set as empty string to use plink to identify lead SNPs from the GWAS summary stat file
#Format of this file is: rsID\tCHR\tPOS\tSTART\tSTOP, with each lead SNP on it's own line)
leadSNPsFilePath="/project/ritchie21/personal/wbone/ColocQuiaL_test/test_runs_dir/BMI_ColocQuiaL_test/GIPR_rsID_20210929_first_10.txt"

#reference genome build of SNPs in lead SNPs file, either "hg19" or "hg38"
build="hg19"

#qtlType is set to "eqtl" or "sqtl"
qtlType="eqtl"

#colc window size in bps, this is only required when no leadSNPsFilePath is set. If you are providing a leadSNPsFile leave the empty string field for window
window=""

#plink parameters for LD clumping Only need to edit these if you wish to change from the default values of 0.0000001, 1000, and 0.2
clumpP1=""
clumpKB=""
clumpR2=""

#provide path to the setup_config.sh (this is the config file with all of the file paths to dependency files such as plink reference files, GTEx eQTL data, etc.
setup_config_sh="/project/ritchie21/personal/wbone/ColocQuiaL_test/ColocQuiaL_test_config_files/setup_config.sh"

#provide path to the desired setup_config.R file
setup_config_R="/project/ritchie21/personal/wbone/ColocQuiaL_test/ColocQuiaL_test_config_files/setup_config.R"

