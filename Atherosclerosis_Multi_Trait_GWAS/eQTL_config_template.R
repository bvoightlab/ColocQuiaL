#ID for user's trait of interest. (Can be any string)
trait = "T2D"
#path to the input files
traitFilePath = "/project/voight_GWAS/wbone/CHD_T2D_2020_bivariate_scan/T2D_Vujkovic_et_al_dbSNP150_onlySNPs_data.txt"
#column IDs from trait file
trait_A1col = "EA"
trait_A2col = "NEA"
trait_SNPcol = "name"
trait_CHRcol = "CHR"
trait_BPcol = "POS"
trait_Pcol = "P"
trait_Ncol = "N"
trait_MAFcol = "EAF"
#trait info not in the input file
#traitType is set either to "cc" or "quant"
traitType = "cc"
#This is the proportion of samples that are cases in a case control GWAS, if you are using a quantitative trait this should be set to "". traitProp = cases / case + controls
traitProp = 0.162369 #look this up
#locus information for running coloc. Currently this assumes these genomic positions to be from hg19
chrom = CHROMOSOME
colocStart = STARTBP
colocStop = STOPBP
lead_SNP = "SNPNUMBER" 

