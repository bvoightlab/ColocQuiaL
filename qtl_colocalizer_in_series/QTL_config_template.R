#ID for user's trait of interest. (Can be any string)
trait = "TRAITNAME" 
#path to the input files
traitFilePath = "TRAITPATH" 
#column IDs from trait file
trait_A1col = "A1COL" 
trait_A2col = "A2COL" 
trait_SNPcol = "SNPCOL" 
trait_CHRcol = "CHRCOL" 
trait_BPcol = "BPCOL" 
trait_Pcol = "PCOL" 
trait_Ncol = "NCOL" 
trait_MAFcol = "MAFCOL" 
#trait info not in the input file
#traitType is set either to "cc" or "quant"
traitType = "TRAITTYPE" 
#This is the proportion of samples that are cases in a case control GWAS, if you are using a quantitative trait this should be set to "". traitProp = cases / case + controls
traitProp = TRAITPROP #look this up
#locus information for running coloc. Currently this assumes these genomic positions to be from hg19
chrom = CHROMOSOME
colocStart = STARTBP
colocStop = STOPBP
lead_SNP = "SNPNUMBER" 
#"eqtl" or "sqtl"
qtlType = "QTL" 
