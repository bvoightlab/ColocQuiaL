#eqtl_colocalizer.R
#This script performs colocalization between a trait of the user's interest and the GTEx single trait eQLTs from
# the GTEx online GUI. It will perform colocalization for each gene-tissue pair in the GTEx csv file downloaded
# from the online GUI. NOTE it currently assumes you wish to use The GTEx_v8 data available on the voight lab LPC
# single trait eQLTs from at /project/voight_datasets_01/GTEx_v8/TissueSpecific_2.

#This script outputs to the directory it is run from and expects a "eQTL_config.R" file to be in the directory 

library(coloc)
library(data.table)
library(rjson)
library(tidyverse)
library(ieugwasr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(glue)
library(vroom)

ld_extract <- function(variants, bfile, plink_bin) {
  # Make textfile
  shell <- ifelse(Sys.info()["sysname"] == "Mac", "cmd", "sh")
  fn <- tempfile()
  
  data.frame(variants) %>% 
    vroom::vroom_write(fn, col_names = FALSE)
  
  fun2 <- paste0(
    shQuote(plink_bin, type = shell),
    " --bfile ", shQuote(bfile, type = shell),
    " --extract ", shQuote(fn, type = shell),
    " --r inter-chr ",
    " --out ", shQuote(fn, type = shell)
  )
  system(fun2)
  
  print("hello")
  res <- data.table::fread(paste0(fn, ".ld")) %>% 
    dplyr::select(SNP_A, SNP_B, R)
  
  system(paste0("rm ", fn, ".ld"))

  print(head(res))
 
  return(bind_rows(res, res %>% rename(SNP_A = "SNP_B", SNP_B = "SNP_A")))
}


gg_regional_association_plink <- function(df, lead_snps = NULL, rsid = rsid, chromosome = chromosome, position = position, p_value = p_value, p_value_threshold = 1, clump_kb = 1000, clump_r2 = 1, plot_distance = 500000, bfile, plink_bin, plot_title = NULL, plot_subtitle = NULL, n_row = 2) {
  df <- df %>%
    dplyr::select(rsid = {{ rsid }}, chromosome = {{ chromosome }}, position = {{ position }}, p_value = {{ p_value }}) %>%
    mutate_if(is.factor, as.character)

  if (!is.null(lead_snps)) {
    indep_snps <- df %>%
      dplyr::select(rsid = {{ rsid }}, pval = {{ p_value }}) %>%
      filter(rsid %in% lead_snps)
  } else {
    indep_snps <- df %>%
      dplyr::select(rsid = {{ rsid }}, pval = {{ p_value }}) %>%
      filter(pval < p_value_threshold) %>%
      ld_clump(bfile = bfile, plink_bin = plink_bin, clump_kb = 500, clump_r2 = 0.2)
  }

  locus_snps <- df %>%
    filter(rsid %in% indep_snps$rsid) %>%
    dplyr::select(chromosome, position, lead_rsid = rsid) %>%
    pmap_dfr(function(chromosome_filter = first, position_filter = second, lead_rsid = third) {
      df %>%
        filter(chromosome == chromosome_filter & between(position, position_filter - plot_distance / 2, position_filter + plot_distance / 2)) %>%
        mutate(lead_position = position_filter) %>%
        mutate(lead_rsid = lead_rsid)
    }) %>%
    mutate(lead_marker = glue::glue("{chromosome}:{lead_position}")) %>%
    group_by(lead_marker) %>%
    mutate(lead_p_value = min(p_value)) %>%
    ungroup() %>%
    mutate(label = paste0(lead_marker, "\n", lead_rsid)) %>%
    mutate(label = fct_reorder(label, lead_p_value))


  locus_snps_ld <- locus_snps %>%
    group_nest(lead_rsid, keep = TRUE) %>%
    ungroup() %>%
    mutate(r2_matrix = map(data, function(x) {
      ld_extract(x$rsid, bfile = bfile, plink_bin = plink_bin) %>%
        filter(SNP_A %in% indep_snps$rsid) %>%
        mutate(r2 = abs(R)^2) %>%
        dplyr::select(-R)
    })) %>%
    dplyr::select(r2_matrix) %>%
    unnest(r2_matrix) %>%
    mutate(color_code = as.character(cut(r2, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("blue4", "blue", "darkgreen", "orange", "red"), include.lowest = TRUE))) %>%
    bind_rows(tibble(
      lead_rsid = unique(locus_snps$lead_rsid),
      rsid = unique(locus_snps$lead_rsid),
      r2 = rep(1, length(unique(locus_snps$lead_rsid)))
    )) %>%
    mutate(color_code = case_when(
      rsid == lead_rsid ~ "purple",
      TRUE ~ color_code
    )) %>%
    mutate(color_code = fct_relevel(color_code, "purple", "red", "orange", "darkgreen", "blue", "blue4")) %>%
    mutate(lead = lead_rsid == rsid)
 
  plot <- locus_snps %>%
    left_join(locus_snps_ld) %>%
    mutate(r2 = ifelse(is.na(r2) == TRUE, 0.1, r2)) %>%
    mutate(lead = ifelse(is.na(lead) == TRUE, FALSE, lead)) %>%
    mutate(color_code = as.character(cut(r2, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("blue4", "blue", "darkgreen", "orange", "red"), include.lowest = TRUE))) %>%
    mutate(color_code = ifelse(lead == TRUE, "purple", color_code)) %>%
    ggplot(aes(position / 1000000, -log10(p_value), fill = color_code, size = lead, alpha = lead, shape = lead)) +
    geom_point() +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
    scale_fill_identity(parse(text = "r^2"), guide = "legend", labels = c("Lead SNP", "0.8 - 1", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0 - 0.2"), breaks = c("purple", "red", "orange", "darkgreen", "blue", "blue4")) +
    scale_size_manual(values = c(3, 8), guide = FALSE) +
    scale_shape_manual(values = c(21, 23), guide = FALSE) +
    scale_alpha_manual(values = c(0.8, 1), guide = FALSE) +
    scale_x_continuous(n.breaks = 3) +
    guides(fill = guide_legend(override.aes = list(shape = 22, size = 8))) +
    facet_wrap(~label, scales = "free", nrow = n_row) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Position (Mb)",
      y = expression(-log[10]("p-value"))
    ) +
    theme_light(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.title.align = 0.5,
      legend.key = element_rect(size = 6, fill = "white", colour = NA)
    )

  return(plot)
}


#Read in the arguments from the config file
source("eQTL_config.R")

#Print all of the config file settings to screen or the stnd out file
print("trait")
print(trait)
print("trait_file_header_info")
print(traitFilePath)
print(trait_A1col)
print(trait_A2col)
print(trait_SNPcol)
print(trait_CHRcol)
print(trait_BPcol)
print(trait_Pcol)
print(trait_Ncol)
print(trait_MAFcol)
print("traitType:")
print(traitType)
print(traitProp)
print("Locus Info:")
print(chrom)
print(colocStart)
print(colocStop)
print(lead_SNP)

#check if the lead_SNP is in the plink BIM file

SNP_str <- paste('"\t',lead_SNP,'\t"',sep="")
SNP_grep_str <- paste('grep','-P', SNP_str, "/project/voight_GWAS/wbone/NGWAMA/data_maf0.01_rs_ref/data_maf0.01_rs_ref.bim", ">", "leadSNP_test_file.txt" , sep=" " )

system(SNP_grep_str)

fileInfo <- file.info("leadSNP_test_file.txt")

if(fileInfo$size == 0){

    SNPin1KG = FALSE

} else{

    SNPin1KG = TRUE

}

#read in the tissue specific eQLT summary file with the file names added
trait_region = fread(file=traitFilePath, sep="\t", header=TRUE)
print("trait input file successfully loaded")

#add "trait" to the Allele fields, the MAF field, and the SNP field to avoid confusion when we compare to the eQTL alleles later 
trait_A1col_str = paste(trait_A1col,"_trait",sep="")
trait_A2col_str = paste(trait_A2col,"_trait",sep="")
trait_SNPcol_str = paste(trait_SNPcol,"_trait",sep="")
trait_MAFcol_str = paste(trait_MAFcol,"_trait",sep="")

colnames(trait_region)[colnames(trait_region)== trait_A1col] <- trait_A1col_str
colnames(trait_region)[colnames(trait_region)== trait_A2col] <- trait_A2col_str
colnames(trait_region)[colnames(trait_region)== trait_SNPcol] <- trait_SNPcol_str
colnames(trait_region)[colnames(trait_region)== trait_MAFcol] <- trait_MAFcol_str

#remove any alphabetical characters from the chromosome column
trait_region[[trait_CHRcol]] <- as.integer(gsub('[a-zA-Z]', '', trait_region[[trait_CHRcol]]))

str(trait_region)

#grab the snps that are within the start stop and on the correct chromosome from the trait file
trait_region = trait_region[trait_region[[trait_CHRcol]] == chrom & trait_region[[trait_BPcol]] >= colocStart & trait_region[[trait_BPcol]] <= colocStop,]

head(trait_region,3)

#grab rs numbers if present
trait_region_rs = trait_region[[trait_SNPcol_str]]
head(trait_region_rs, 3)

#use hash tables to find chromosome positions
hash_table_file = paste0("/project/voight_GWAS/bychen9/human_9606_b151_GRCh38p7/BED/hash_tables_2/chr_", chrom, "_snp151_hash_table.json")
print(hash_table_file)

hash_table <- fromJSON(file = hash_table_file) 
head(hash_table, 3)
print("hash table loaded")

trait_chrom_pos = hash_table[trait_region_rs]

#remove empty elements 
trait_chrom_pos = trait_chrom_pos[!sapply(trait_chrom_pos,is.null)]
head(trait_chrom_pos, 3)

tissueTable = read.table(file="/project/voight_GWAS/bychen9/eqtl_coloc/GTEx_Tissue_Summary_with_filenames.csv", sep=",", header=TRUE)

#Convert to strings from factors
tissueTable$Tissue = as.character(tissueTable$Tissue)
tissueTable$Filename = as.character(tissueTable$Filename)

#Format Tissue names
tissueTable_tissue_noSpace = gsub("\\(","",tissueTable$Tissue)
tissueTable_tissue_noSpace = gsub("\\)","",tissueTable_tissue_noSpace)
tissueTable_tissue_noSpace = gsub("[[:space:]]","_",tissueTable_tissue_noSpace)
tissueTable_tissue_noSpace = gsub("-_", "", tissueTable_tissue_noSpace)
tissueTable$Tissue = tissueTable_tissue_noSpace

#create csv file from significant pair files
signif_pair_files <- list.files(path="/project/voight_datasets/GTEx_v8/eQTL/GTEx_Analysis_v8_eQTL_tabix", pattern="*tab.gz$", full.names=TRUE, recursive=FALSE)
lead_SNP_pos = hash_table[[lead_SNP]]
#convert format for tabix
lead_SNP_pos_tabix = paste0(gsub("_",":",lead_SNP_pos), "-", gsub("^.*?_","",lead_SNP_pos))

lapply(signif_pair_files, function(file) {
	system(paste("tabix", file, lead_SNP_pos_tabix, ">", paste0(lead_SNP, "_temp.csv")))    

	#add tissue name
	tissue_name = gsub(".v8.signif_variant_gene_pairs.tab.gz", "", gsub("/project/voight_datasets/GTEx_v8/eQTL/GTEx_Analysis_v8_eQTL_tabix/", "", file))
	system(paste0("sed -i \"s/$/\t", tissue_name, "/\" ", lead_SNP, "_temp.csv"))
	system(paste0("cat ", lead_SNP, "_temp.csv >> ", lead_SNP, ".csv"))
})

#read in the csv file of eGene-Tissue pairs 
eGenes = tryCatch({
		read.table(file=paste0(lead_SNP, ".csv"), sep="\t", stringsAsFactors=FALSE)
	}, error = function(err){
		print(paste(lead_SNP, "is not a significant eQTL in any tissue of GTEx version 8."))
		quit(status=0)
	})

#loop through the eGene-Tissue pairs in eGenes and prep running COLOC
for(i in 1:nrow(eGenes)){
    geneID <- eGenes[i,7]
    geneID_noDOT <-  gsub("\\..*","", geneID)

    geneSymbol <- tryCatch({

        bitr(geneID_noDOT, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")$SYMBOL

    }, error = function(err){

        print("ENSEMBL ID could not be converted to HGNC Symbol")
        print(paste("geneSymbol will be set to the ENSEMBL ID",geneID_noDOT))

        return(geneID_noDOT)

    })

    tissue <- eGenes[i,18]
    
    print(geneID)
    print(geneSymbol)

	if (length(geneSymbol) > 1) {
		geneSymbol <- geneSymbol[1]
		print(paste("Using first gene symbol: ", geneSymbol))
	}

    print(tissue)

    # find the all pair file that contains the tissue of interest
	
    tissueLine <- tissueTable[tissueTable$Tissue == tissue,] 
    allpair_filename <- tissueLine$Filename

    #if the Filename field is NA then skip this eGene-Tissue pair
    if (is.na(allpair_filename)){
        print(tissue)
        print("This tissue is not available in the all pairs files currently")
        next
    }
    
    tabix_allpair_path = paste0("/project/voight_datasets_01/GTEx_v8/TissueSpecific_tabix/", gsub("txt", "tab", allpair_filename))
    eqtl_N <- tissueLine$NumberRNASeqandGTSamples

    #parentheses are causing issues too
    tissue_noSpace = gsub("\\(","",tissue)
    tissue_noSpace = gsub("\\)","",tissue_noSpace)
    # The tissue names have any whitespace in them and we want to use these in the output file names so replace " " with "_"
    tissue_noSpace = gsub("[[:space:]]","_",tissue_noSpace)
    
    # make a eGene-Tissue and trait prefix for file names
	out_prefix = paste(geneSymbol,geneID,tissue_noSpace,trait,sep="_")

	#run liftover on colocStart and colocStop
	print("running liftOver")
	bed_liftover = data.frame("chr" = c(paste0("chr", chrom), paste0("chr", chrom)), "bp1" = c(colocStart - 1, colocStop - 1), "bp2" = c(colocStart, colocStop)) 
	write.table(bed_liftover,file=paste0("temp_hg19.bed"),sep="\t",quote = FALSE,row.names=FALSE,col.names=FALSE)

	system("liftOver temp_hg19.bed /appl/liftOver-20180423/chains/hg19ToHg38.over.chain temp_hg38.bed temp_hg19.unmapped -bedPlus=3 -tab")

	hg38_positions = as.data.frame(read.table("temp_hg38.bed", header=FALSE, sep="\t")) 

    print("Grabing the all pairs data")    
    #Use tabix to grab data
	eGeneTissueInputFile = paste(geneSymbol,tissue_noSpace,chrom,colocStart,colocStop,".txt", sep="_" )
    system(paste0("tabix ", tabix_allpair_path, " chr", chrom, ":", hg38_positions[1,3], "-", hg38_positions[2,3], " > ", eGeneTissueInputFile))    

    print("reading the all pairs data into R")
    #read the file we just generated from the grep command into R
    eGeneTissueInput = fread(file = eGeneTissueInputFile, header=FALSE)

	#check if tabix result is empty
	print(dim(eGeneTissueInput))
	if (dim(eGeneTissueInput)[1] == 0) {
		print("Did not find gene in all pairs data")
		quit(status=0)
	}	
	
    #remove the eGeneTissueInputFile after has been read into R to save disk space
    system(paste0("rm ", eGeneTissueInputFile, sep=""))

    header = c("chrom_b38", "chromStart_b38", "chromEnd_b38", "eGeneID", "A1_eqtl", "A2_eqtl", "build", "tss_distance", "ma_samples", "ma_count", "maf", "pvalue_eQTL", "slope", "slope_se")

    colnames(eGeneTissueInput) = header

    print("Filtering on the geneID")
    #filter for the geneID of interest
    eGeneTissue_region = eGeneTissueInput[eGeneTissueInput$eGeneID == geneID,]

    if (nrow(eGeneTissue_region) == 0) {
      print("Warning: There was not an exact match on Ensembl ID. Likely this is due to a GTEX version  update.")
      #make a string that removes the everything after in the geneID
      noDecimalGeneID = gsub("\\..*","",geneID)
  
      #recreated eGeneTissue_region
      #eGeneTissue_region = eGeneTissueInput[eGeneTissueInput$chromEnd >= colocStart & eGeneTissueInput$chromEnd <= colocStop,]
	  eGeneTissue_region = eGeneTissueInput
  
      #grep the simplified Ensembl ID
      possible_Ensembl_gene_lines <- eGeneTissue_region[grepl(noDecimalGeneID, eGeneTissue_region$eGeneID),]
  
      #check to make sure there is just one other Ensembl ID
      possible_Ensembl_genes <- unique(possible_Ensembl_gene_lines$eGeneID)
  
      if(length(possible_Ensembl_genes) == 1){
        print("Found a unique Ensembl ID so this analysis will continue using the Ensembl ID:")
        print(possible_Ensembl_genes)
    
        #this will be a single Ensembl ID string
        geneID = possible_Ensembl_genes
        eGeneTissue_region = eGeneTissue_region[eGeneTissue_region$eGeneID == geneID,]

      } else {
        print("The Ensembl ID from your GTEx csv was not able to be reliably mapped to an Enseml ID in the GTEx database, so this gene will be skipped:")
        print(geneID)
        next
      }    
    }


    print("adding rs numbers to the eQTL data")
    #add rs genegene,,numbers to the eGeneTissue_region DF

	#make chromosome_position column for merging
	eGeneTissue_region$chromosome_position <- paste(eGeneTissue_region$chrom_b38,eGeneTissue_region$chromEnd_b38,sep="_")

    #create data frame with rs numbers associated with chromosome_position and add to eGeneTissue region
    uniqID_DF = as.data.frame(t(as.data.frame(trait_chrom_pos)))
    uniqID_DF$SNP <- rownames(uniqID_DF)
    colnames(uniqID_DF) <- c("chromosome_position", "SNP")

    eGeneTissue_region = merge(eGeneTissue_region, uniqID_DF, by = "chromosome_position")  

    print("merging the trait and eqtl data on unique ID")
    #merge the trait and eGeneTissue region DFs on rs numbers
    colocInputFile = merge(eGeneTissue_region, trait_region, by.x="SNP", by.y=trait_SNPcol_str)

    #remove any NAs
    colocInputFile = colocInputFile[complete.cases(colocInputFile), ]

    #write colocInputFile to file for making locus zoom plots
    colocInputFile_outputStr = paste(out_prefix,"coloc_input_data.txt",sep="_")
    write.table(colocInputFile, file= colocInputFile_outputStr, sep="\t", row.names=FALSE, quote=FALSE)

    print("Running coloc")
    #run coloc
    if (traitType == "cc"){

        coloc_results <- coloc.abf(dataset1=list(pvalues=colocInputFile[[trait_Pcol]], N=colocInputFile[[trait_Ncol]], type=traitType, s=traitProp), dataset2=list(pvalues=colocInputFile$pvalue_eQTL, N=eqtl_N, type="quant"),MAF=colocInputFile[[trait_MAFcol_str]])

    } else {

        coloc_results <- coloc.abf(dataset1=list(pvalues=colocInputFile[[trait_Pcol]], N=colocInputFile[[trait_Ncol]], type=traitType), dataset2=list(pvalues=colocInputFile$pvalue_eQTL, N=eqtl_N, type="quant"),MAF=colocInputFile[[trait_MAFcol_str]])

    }

    #prepare useful outputs
    coloc_results_summary = coloc_results$summary
    coloc_results_full = coloc_results$results

    #calculate pp4 / pp3 + pp4
    PP3andPP4 = coloc_results_summary[5] + coloc_results_summary[6]

    pp4_conditional = coloc_results_summary[6] / PP3andPP4

    pp4_conditional = coloc_results_summary[6] / PP3andPP4

    #prep coloc output strings
    coloc_results_summary_outputStr = paste(out_prefix,"coloc_results_summary.txt",sep="_")
    coloc_results_full_outputStr = paste(out_prefix,"coloc_results_full.txt",sep="_")
    coloc_results_pp4_cond_outputStr = paste(out_prefix,"coloc_results_pp4_cond.txt",sep="_")

    #write to file
    write.table(coloc_results_summary, file=coloc_results_summary_outputStr, sep="\t", row.names=TRUE, quote=FALSE)
    write.table(coloc_results_full, file=coloc_results_full_outputStr, sep="\t", row.names=FALSE, quote=FALSE)
    write.table(pp4_conditional, file=coloc_results_pp4_cond_outputStr, sep="\t", row.names=FALSE, quote=FALSE)

	#generate regional association plot

    #find lead SNP in 1KG if needed
    if (SNPin1KG == TRUE){

        print(lead_SNP)

    } else{

        print("lead SNP is not in 1KG so we need to find a different SNP for making the RA plots")

        find_new_lead_SNP_df <- colocInputFile %>% dplyr::select(SNP, chrom_b38, chromEnd_b38, pvalue_eQTL, P) %>% dplyr::select(rsid = SNP, chromosome = chrom_b38, position = chromEnd_b38, pval = P)
        find_new_lead_SNP_df$chromosome = as.integer(gsub('[a-zA-Z]', '', find_new_lead_SNP_df$chromosome)) 

        ld_clump_df <- find_new_lead_SNP_df %>%
        ld_clump(bfile = "/project/voight_GWAS/wbone/NGWAMA/data_maf0.01_rs_ref/data_maf0.01_rs_ref", plink_bin = "plink", clump_kb = 500, clump_r2 = 0.2)

        lead_SNP <- as.character(ld_clump_df[which.min(ld_clump_df$pval),]$rsid)

}


	leadSNP_DF = colocInputFile#[colocInputFile$SNP == lead_SNP,]
	leadSNP_DF$chrom_b38 = as.integer(gsub('[a-zA-Z]', '', leadSNP_DF$chrom_b38)) 
	leadSNP_DF = leadSNP_DF %>% dplyr::select(SNP, chrom_b38, chromEnd_b38, pvalue_eQTL, P)

	if (i == 1) {
		#plot for trait	
		trait_leadSNP_DF = leadSNP_DF %>% dplyr::select(rsid = SNP, chromosome = chrom_b38, position = chromEnd_b38, p_value = P)
		trait_plot_title = paste(lead_SNP, trait)
		RA_plot <- gg_regional_association_plink(trait_leadSNP_DF, p_value_threshold = 5e-8, lead_snps = lead_SNP, bfile = "/project/voight_GWAS/wbone/NGWAMA/data_maf0.01_rs_ref/data_maf0.01_rs_ref", plink_bin = "plink", plot_title = paste(trait_plot_title, "Regional Association Plot"), plot_subtitle = expression(""))
		pdf(file = paste0(lead_SNP, "_", trait,".pdf"), paper = 'USr', width = 15, height = 20)  
		print(RA_plot)
		dev.off()
	}

	eQTL_leadSNP_DF = leadSNP_DF %>% dplyr::select(rsid = SNP, chromosome = chrom_b38, position = chromEnd_b38, p_value = pvalue_eQTL)
	eQTL_plot_title = paste(lead_SNP, geneSymbol, tissue)
	RA_plot <- gg_regional_association_plink(eQTL_leadSNP_DF, p_value_threshold = 5e-8, lead_snps = lead_SNP, bfile = "/project/voight_GWAS/wbone/NGWAMA/data_maf0.01_rs_ref/data_maf0.01_rs_ref", plink_bin = "plink", plot_title = paste(eQTL_plot_title, "Regional Association Plot"), plot_subtitle = expression("GTEx v8"))

	pdf(file = paste0(lead_SNP, "_", geneSymbol, "_", tissue,".pdf"), paper = 'USr', width = 15, height = 20)  
	print(RA_plot)
	dev.off()
}
