#colocquial.R
#This script performs colocalization between a trait of the user's interest and the GTEx single trait eQLTs or sQTLs from
# the GTEx online GUI. It will perform colocalization for each gene-tissue pair in the GTEx csv file downloaded
# from the online GUI. NOTE it currently assumes you wish to use The GTEx_v8 data available on the voight lab LPC
# single trait eQLTs from at /project/voight_datasets_01/GTEx_v8/TissueSpecific_2.

#This script outputs to the directory it is run from and expects a "QTL_config.R" file to be in the directory 

library(coloc)
library(data.table)
library(rjson)
library(tidyverse)
library(ieugwasr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(glue)
library(vroom)
#library(BSgenome.Hsapiens.UCSC.hg19)
library(Homo.sapiens)
library(ggbio)
library(egg)

qtl_vars <- c(
  "eqtl" = c(
    "sig_qtl_tabix_dir" = eQTL_sig_qtl_tabix_dir,
    "sig_geneID_col" = eQTL_sig_geneID_col,
    "all_qtl_tabix_dir" = eQTL_all_qtl_tabix_dir,
    "all_header" = eQTL_all_header,
    "all_geneID" = eQTL_all_geneID,
    "all_chrom" = eQTL_all_chrom,
    "all_chromEnd" = eQTL_all_chromEnd,
    "all_pvalue" = eQTL_all_pvalue,
    "QTL_tissue_table" = eQTL_tissue_table),
  "sqtl" = c(
    "sig_qtl_tabix_dir" = sQTL_sig_qtl_tabix_dir,
    "sig_geneID_col" = sQTL_sig_geneID_col,
    "all_qtl_tabix_dir" = sQTL_all_qtl_tabix_dir,
    "all_header" = sQTL_all_header,
    "all_geneID" = sQTL_all_geneID,
    "all_chrom" = sQTL_all_chrom,
    "all_chromEnd" = sQTL_all_chromEnd,
    "all_pvalue" = sQTL_all_pvalue,
    "QTL_tissue_table" = sQTL_tissue_table,
  ),
)


ld_extract <- function(variants, bfile, plink_bin) {
  # Make textfile
  shell <- ifelse(Sys.info()["sysname"] == "Mac", "cmd", "sh")
  fn <- tempfile()
  
  data.frame(variants) %>% 
    vroom::vroom_write(fn, col_names = FALSE)

  #generate a plink extract command  
  fun2 <- paste0(
    shQuote(plink_bin, type = shell),
    " --bfile ", shQuote(bfile, type = shell),
    " --extract ", shQuote(fn, type = shell),
    " --r inter-chr ",
    " --out ", shQuote(fn, type = shell)
  )
  system(fun2)
  
  res <- data.table::fread(paste0(fn, ".ld")) %>% 
    dplyr::select(SNP_A, SNP_B, R)
  
  system(paste0("rm ", fn, ".ld"))

  print(head(res))
 
  return(bind_rows(res, res %>% rename(SNP_A = "SNP_B", SNP_B = "SNP_A")))
}


gg_regional_association_plink <- function(df, lead_snps = NULL, rsid = rsid, chromosome = chromosome, position = position, p_value = p_value, p_value_threshold = 0.0000001, clump_kb = 1000, clump_r2 = 0.2, plot_distance = 500000, bfile, plink_bin, plot_title = NULL, plot_subtitle = NULL, n_row = 2, region_recomb = region_recomb) {
  df <- df %>%
    dplyr::select(rsid = {{ rsid }}, chromosome = {{ chromosome }}, position = {{ position }}, p_value = {{ p_value }}) %>%
    mutate_if(is.factor, as.character)

  #identify the lead SNP
  if (!is.null(lead_snps)) {
    indep_snps <- df %>%
      dplyr::select(rsid = {{ rsid }}, pval = {{ p_value }}) %>%
      filter(rsid %in% lead_snps)
  } else {
    indep_snps <- df %>%
      dplyr::select(rsid = {{ rsid }}, pval = {{ p_value }}) %>%
      filter(pval < p_value_threshold) %>%
      ld_clump(bfile = bfile, plink_bin = plink_bin, clump_kb = 1000, clump_r2 = 0.2)
  }

  head(df)
  #based on the lead SNP get the other SNPs at the locus
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

  head(locus_snps)
  #add LD data to the SNPs df
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
      SNP_A = unique(locus_snps$lead_rsid),
      SNP_B = unique(locus_snps$lead_rsid),
      r2 = 1
    )) %>%
    mutate(color_code = case_when(
      SNP_A == SNP_B ~ "purple",
      TRUE ~ color_code
    )) %>%
    mutate(color_code = fct_relevel(color_code, "purple", "red", "orange", "darkgreen", "blue", "blue4")) %>%
    mutate(lead = SNP_A == SNP_B)
 
    #rename the df to be more user friendly
    names(locus_snps_ld) <- c("lead_rsid","rsid", "r2", "color_code", "lead")

    RA_plot_data <- locus_snps %>%
    left_join(locus_snps_ld) %>%
    mutate(r2 = ifelse(is.na(r2) == TRUE, 0.1, r2)) %>%
    mutate(lead = ifelse(is.na(lead) == TRUE, FALSE, lead)) %>%
    mutate(color_code = as.character(cut(r2, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("blue4", "blue", "darkgreen", "orange", "red"), include.lowest = TRUE))) %>%
    mutate(color_code = ifelse(lead == TRUE, "purple", color_code))

    #find the max of the recomb rate for scaling
    maxlogP <- max(-log10(RA_plot_data$p_value))

    #scale recomb for plotting with p-values
    region_recomb$`Rate(cM/Mb)` <- region_recomb$`Rate(cM/Mb)` * (maxlogP/200)

    rescale <- 1/(maxlogP/200)

    plot <- ggplot() +
  geom_point(data = RA_plot_data, mapping = aes(x = position / 1000000, y = -log10(p_value), fill = color_code, size = lead, alpha = lead, shape = lead)) +
  geom_line(data = region_recomb, mapping = aes(x= `Position(bp)`/ 1000000 , y = `Rate(cM/Mb)`)) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
  scale_fill_identity(parse(text = "r^2"), guide = "legend", labels = c("Lead SNP", "0.8 - 1", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0 - 0.2"), breaks = c("purple", "red", "orange", "darkgreen", "blue", "blue4")) +    scale_size_manual(values = c(3, 8), guide = FALSE) +
    scale_shape_manual(values = c(21, 23), guide = FALSE) +
    scale_alpha_manual(values = c(0.8, 1), guide = FALSE) +
    scale_x_continuous(n.breaks = 3) +
    scale_y_continuous(name= expression(-log[10]("p-value")), sec.axis = sec_axis(~. * rescale, name = "Recomb Rate cM/Mb"))+
    guides(fill = guide_legend(override.aes = list(shape = 22, size = 5))) +
    facet_wrap(~label, scales = "free", nrow = n_row) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Position (Mb)"
    ) +   
    theme_light(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.title.align = 0.5,
      legend.key = element_rect(size = 6, fill = "white", colour = NA)
    )

  return(plot)
}

#generate the gene tracks for the RA plots
ggbio_genetrack <- function(chrom_str, BPStart, BPStop){

    gene_region <- GRanges(
seqnames = Rle(c(chrom_str), c(1)),
ranges = IRanges(BPStart:BPStop))

    plot <- autoplot(Homo.sapiens, which = gene_region) + xlim(BPStart,BPStop) + scale_x_continuous(expand=c(0,0))

    #convert from ggbio to ggplot object
    plot <- plot@ggplot

    return(plot)

}

#Print all of the config file settings to screen or the stnd out file
print_config_settings <-function(){
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

  print("QTL type:")
  print(qtlType)
}

#NK
prep_coloc_input_file <- function(qtl_type) {
    sprintf("merging the trait and %s on unique ID",qtl_type) 
    colocInputFile = merge(eGeneTissue_region, trait_region, by.x="SNP", by.y=trait_SNPcol_str) # nolint
    colocInputFile = colocInputFile[complete.cases(colocInputFile), ]

    if (0 %in% colocInputFile[[trait_Pcol]]){
    print("WARNING: THERE ARE SNPS WITH P-VALUES OF 0 AT THIS LOCUS. These SNPs have been removed for the Colocalization anlysis and may lead to unusual regional association plots") # nolint
    
    #remove SNPs who's trait P-value is 0 # nolint
    colocInputFile = colocInputFile[colocInputFile[[trait_Pcol]] != 0,]
  }

  return(colocInputFile)
}

#NK
get_coloc_results <- function(colocInputFile, qtl_all_pvalue) {
if (traitType == "cc") {
    return(coloc.abf(dataset1=list(pvalues=colocInputFile[[trait_Pcol]], N=colocInputFile[[trait_Ncol]], type=traitType, s=traitProp), dataset2=list(pvalues=colocInputFile[[qtl_all_pvalue]], N=qtl_N, type="quant"),MAF=colocInputFile[[trait_MAFcol_str]]))
  } else {
    return(oloc.abf(dataset1=list(pvalues=colocInputFile[[trait_Pcol]], N=colocInputFile[[trait_Ncol]], type=traitType), dataset2=list(pvalues=colocInputFile[[qtl_all_pvalue]], N=qtl_N, type="quant"),MAF=colocInputFile[[trait_MAFcol_str]]))
  }
}

#NK
write_coloc_input_file <- function(colocInputFile) {
#write colocInputFile to file for making locus zoom plots
  colocInputFile_outputStr = paste(out_prefix,"coloc_input_data.txt",sep="_")
  write.table(colocInputFile, file= colocInputFile_outputStr, sep="\t", row.names=FALSE, quote=FALSE)
}

write_table_to_file <- function(coloc_results) {
  #prepare useful outputs
  coloc_results_summary = coloc_results$summary
  coloc_results_full = coloc_results$results
  
  #calculate pp4 / pp3 + pp4
  PP3andPP4 = coloc_results_summary[5] + coloc_results_summary[6]
  pp4_conditional = coloc_results_summary[6] / PP3andPP4
  
  #prep coloc output strings
  coloc_results_summary_outputStr = paste(out_prefix,"coloc_results_summary.txt",sep="_")
  coloc_results_full_outputStr = paste(out_prefix,"coloc_results_full.txt",sep="_")
  coloc_results_pp4_cond_outputStr = paste(out_prefix,"coloc_results_pp4_cond.txt",sep="_")
  
  #write to file
  write.table(coloc_results_summary, file=coloc_results_summary_outputStr, sep="\t", row.names=TRUE, quote=FALSE)
  write.table(coloc_results_full, file=coloc_results_full_outputStr, sep="\t", row.names=FALSE, quote=FALSE)
  write.table(pp4_conditional, file=coloc_results_pp4_cond_outputStr, sep="\t", row.names=FALSE, quote=FALSE)
}

#NK
find_lead_snp_in_ld <- function(colocInputFile, qtl_all_chrom, qtl_all_pvalue) {
    
    print("lead SNP is not in the provided LD reference, so we need to find a different SNP for making the RA plots")
    
    find_new_lead_SNP_df <- colocInputFile %>% dplyr::select(SNP, all_of(qtl_all_chrom), all_of(trait_BPcol), all_of(qtl_all_pvalue), all_of(trait_Pcol)) %>% dplyr::select(rsid = SNP, chromosome = all_of(qtl_all_chrom), position = all_of(trait_BPcol), pval = all_of(trait_Pcol))
    find_new_lead_SNP_df$chromosome = as.integer(gsub('[a-zA-Z]', '', find_new_lead_SNP_df$chromosome)) 
    
    ld_clump_df <- find_new_lead_SNP_df %>%
    ld_clump(bfile = plink_bfile, plink_bin = "plink", clump_kb = clump_KB , clump_r2 = clump_R2)
    
    #This was left in commented out in the sQTL code
    #ld_clump(bfile = "/project/voight_GWAS/wbone/NGWAMA/data_maf0.01_rs_ref/data_maf0.01_rs_ref", plink_bin = "plink", clump_kb = 500, clump_r2 = 0.2)

    lead_SNP <- as.character(ld_clump_df[which.min(ld_clump_df$pval),]$rsid)
    return (lead_SNP)
}

#NK
create_lead_snp_df <- function(colocInputFile, qtl_all_chrom, qtl_all_pvalue) {
  leadSNP_DF = colocInputFile#[colocInputFile$SNP == lead_SNP,]
  leadSNP_DF[[eQTL_all_chrom]] = as.integer(gsub('[a-zA-Z]', '', leadSNP_DF[[eQTL_all_chrom]])) 
  leadSNP_DF = leadSNP_DF %>% dplyr::select(SNP, all_of(eQTL_all_chrom), all_of(trait_BPcol), all_of(eQTL_all_pvalue), all_of(trait_Pcol))
  
  return (leadSNP_DF)
}

#NK
plot_trait_data_first_iteration <- function(colocInputFile, qtl_all_chrom) {

  #plot trait data if it is the first time going through the loop	
    trait_leadSNP_DF = leadSNP_DF %>% dplyr::select(rsid = SNP, chromosome = all_of(qtl_all_chrom), position = all_of(trait_BPcol), p_value = all_of(trait_Pcol))
    trait_plot_title = paste(lead_SNP, trait)
    
    RA_plot <- gg_regional_association_plink(trait_leadSNP_DF, p_value_threshold = clump_P1, lead_snps = lead_SNP, bfile = plink_bfile, plink_bin = "plink", plot_distance = bps_in_region, plot_title = paste(trait_plot_title, "Regional Association Plot"), plot_subtitle = expression(""), region_recomb = region_recomb)
    
    #make a gene track plot
    gene_track_plot <- ggbio_genetrack(chrom_str, colocStart, colocStop)
    
    #combine the RA plot and the gene track plot
    RA_plot <- ggarrange(RA_plot, gene_track_plot, widths=c(1,1),heights=c(5,3))
    
    return (RA_plot)
}

#NK
gene_track_plot <- function(qtl_all_chrom, qtl_all_pvalue) {

  eQTL_leadSNP_DF = leadSNP_DF %>% dplyr::select(rsid = SNP, chromosome = all_of(qtl_all_chrom), position = all_of(trait_BPcol), p_value = all_of(qtl_all_pvalue))
  eQTL_plot_title = paste(lead_SNP, geneSymbol, tissue)
  RA_plot <- gg_regional_association_plink(eQTL_leadSNP_DF, p_value_threshold = clump_P1, lead_snps = lead_SNP, bfile = plink_bfile, plink_bin = "plink", plot_distance = bps_in_region, plot_title = paste(eQTL_plot_title, "Regional Association Plot"), plot_subtitle = expression("GTEx v8"), region_recomb = region_recomb)
  
  #make a gene track plot
  gene_track_plot <- ggbio_genetrack(chrom_str, colocStart, colocStop)
  
  #combine the RA plot and the gene track plot
  RA_plot <- ggarrange(RA_plot, gene_track_plot, widths=c(1,1),heights=c(5,3))
  return (RA_plot)
}

#NK
eqtl_colocalization <- function() {
  #NK
  colocInputFile = prep_coloc_input_file()
  write_coloc_input_file(colocInputFile)
  
  #NK - common to both
  print("Running coloc")
  #run coloc
  coloc_results = get_coloc_results(colocInputFile, eQTL_all_pvalue)
  
  write_table_to_file(coloc_results)
  
  #generate regional association plot
  
  #find lead SNP in LD reference if needed
  if (SNPinLDref == TRUE) {

    print(lead_SNP)

  } else {
    
    print("lead SNP is not in the provided LD reference, so we need to find a different SNP for making the RA plots")
    
    lead_SNP <- find_lead_snp_in_ld(colocInputFile, eQTL_all_chrom, eQTL_all_pvalue)
  }  
  
  leadSNP_DF = create_lead_snp_df(colocInputFile, eQTL_all_chrom, eQTL_all_pvalue)
  
  if (i == 1) {
    RA_plot <- plot_trait_data_first_iteration(colocInputFile, eQTL_all_chrom)
    pdf(file = paste0(lead_SNP, "_", trait,".pdf"), paper = 'USr', width = 15, height = 20)  
    print(RA_plot)
    dev.off()
  }
  
  RA_plot <- gene_track_plot(colocInputFile, eQTL_all_chrom, eQTL_all_pvalue)
  
  pdf(file = paste0(lead_SNP, "_", geneSymbol, "_", tissue,".pdf"), paper = 'USr', width = 15, height = 20)  
  print(RA_plot)
  dev.off()
}

#NK
sqtl_colocalization <- function(){

  #NK
  colocInputMasterFile = prep_coloc_input_file(qtl_type)       
  
  #UNIQUE TO SQTL
  {
  #combine intron columns
  colocInputMasterFile$intronID = paste(colocInputMasterFile[[sQTL_all_intron_chr]], colocInputMasterFile[[sQTL_all_intron_bp_first]], colocInputMasterFile[[sQTL_all_intron_bp_end]], colocInputMasterFile[[sQTL_all_intron_clu]], colocInputMasterFile[[sQTL_all_geneID]], sep=":")
  
  print("finding unique introns")
  #find all unique introns
  uniqueIntrons = unique(colocInputMasterFile$intronID)
  }
  
  #loop through unique introns
  for(j in 1:length(uniqueIntrons)) {
    print(uniqueIntrons[j])
    intronID = uniqueIntrons[j]
    
    #grep intron lines from colocInputMasterFile
    intron_lines = grepl(intronID, colocInputMasterFile$intronID)
    colocInputFile = unique(colocInputMasterFile[intron_lines,])
    
    # make a eGene-Tissue and trait prefix for file names
    out_prefix = paste(geneSymbol,intronID,tissue_noSpace,trait,sep="_")
    
    write_coloc_input_file(colocInputFile)

    
    #NK - common to both
    print("Running coloc")
    #run coloc
    coloc_results = get_coloc_results(colocInputFile, eQTL_all_pvalue)

    write_table_to_file(coloc_results)
    
    #generate regional association plot
    
    #find lead SNP in LD reference if needed
    if (SNPinLDref == TRUE) {

      print(lead_SNP)

    } else {
    
      print("lead SNP is not in the provided LD reference, so we need to find a different SNP for making the RA plots")
    
      lead_SNP <- find_lead_snp_in_ld(colocInputFile, sQTL_all_chrom, sQTL_all_pvalue)
    }
    
    leadSNP_DF = create_lead_snp_df(colocInputFile = , sQTL_all_chrom, sQTL_all_pvalue)

    if (i == 1 && j == 1) {
      RA_plot <- plot_trait_data_first_iteration(colocInputFile, sQTL_all_chrom)
      pdf(file = paste0(lead_SNP, "_", trait,".pdf"), paper = 'USr', width = 15, height = 20)  
      print(RA_plot)      
      dev.off()
    }
    
    sQTL_leadSNP_DF = leadSNP_DF %>% dplyr::select(rsid = SNP, chromosome = all_of(sQTL_all_chrom), position = all_of(trait_BPcol), p_value = all_of(sQTL_all_pvalue))
    sQTL_plot_title = paste(lead_SNP, geneSymbol, tissue)
    RA_plot <- gg_regional_association_plink(sQTL_leadSNP_DF, p_value_threshold = clump_P1, lead_snps = lead_SNP, bfile = plink_bfile, plink_bin = "plink", plot_distance = bps_in_region, plot_title = paste(sQTL_plot_title, "Regional Association Plot"), plot_subtitle = expression("GTEx v8"), region_recomb = region_recomb)
    
    #make a gene track plot
    gene_track_plot <- ggbio_genetrack(chrom_str, colocStart, colocStop)
    
    #combine the RA plot and the gene track plot
    RA_plot <- ggarrange(RA_plot, gene_track_plot, widths=c(1,1),heights=c(5,3))
    
    pdf(file = paste0(lead_SNP, "_", geneSymbol, "_", intronID, "_", tissue,".pdf"), paper = 'USr', width = 15, height = 20)
    print(RA_plot)
    dev.off()
    #NK - end common to both
  }
}

#tissue table function - NK
handle_tissue_table <- function(tissueTable){
  for (i in 1:nrow(tissueTable)) {
    sigpair_filename = tissueTable$sigPairsTabixFilename[i]
    if (is.na(sigpair_filename)) {
      print(paste(tissueTable$Tissue[i], "is not available in the significant pairs files"))
      next
    }
    
    file = paste0(sig_qtl_tabix_dir, "/", sigpair_filename) 
    
    system(paste("tabix", file, lead_SNP_pos_tabix_with_chr, ">", paste0(lead_SNP, "_temp.csv"))) 
    system(paste("tabix", file, lead_SNP_pos_tabix_without_chr, ">>", paste0(lead_SNP, "_temp.csv"))) 
    
    system(paste0("sed -i \"s/$/\t", tissueTable$Tissue[i], "/\" ", lead_SNP, "_temp.csv"))
    system(paste0("cat ", lead_SNP, "_temp.csv >> ", lead_SNP, ".csv"))
  }
}

#MAIN WORKFLOW STARTS HERE
#Read in the arguments from the config file
source("QTL_config.R")
source(setup_config_R)

print_config_settings()

# Set up QTL variables for either eQTL or sQTL
if (qtlType == "eqtl" | qtlType == "sqtl") {
    sig_qtl_tabix_dir = qtl_vars[qtlType]["sig_qtl_tabix_dir"]
    sig_geneID_col = qtl_vars[qtlType]["sig_geneID_col"]
    all_qtl_tabix_dir = qtl_vars[qtlType]["all_qtl_tabix_dir"]
    all_header = qtl_vars[qtlType]["all_header"]
    all_geneID = qtl_vars[qtlType]["all_geneID"]
    all_chrom = qtl_vars[qtlType]["all_chrom"]
    all_chromEnd = qtl_vars[qtlType]["all_chromEnd"]
    all_pvalue = qtl_vars[qtlType]["all_pvalue"]
    QTL_tissue_table = qtl_vars[qtlType]["QTL_tissue_table"]
} else {
    print("ERROR: Please specify qtlType: \"eqtl\" or \"sqtl\"")
    quit()
}


#some libraries require the "chr" at the beginning of the chromosome
chrom_str <- paste0("chr",chrom, sep="")

bps_in_region = colocStop - colocStart

#grab the recombination rate data for this region from 1KG recomb file
recomb_file_path <- paste(recomb_rate_data,"-",chrom,"-final.txt.gz", sep="")
chr_recomb <- vroom(recomb_file_path)
region_recomb <- chr_recomb[chr_recomb$`Position(bp)` >= colocStart & chr_recomb$`Position(bp)` <= colocStop, ]

#check if the lead_SNP is in the plink BIM file
SNP_str <- paste('"\t',lead_SNP,'\t"',sep="")
bim_file = paste(plink_bfile,".bim", sep ="")
SNP_grep_str <- paste('grep','-P', SNP_str, bim_file, ">", "leadSNP_test_file.txt" , sep=" " )

system(SNP_grep_str)

fileInfo <- file.info("leadSNP_test_file.txt")
SNPinLDref = (fileInfo$size != 0)


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
hash_table_file = paste0(hash_table_dir, "chr_", chrom, "_snp151_hash_table.json")
print(hash_table_file)

hash_table <- fromJSON(file = hash_table_file) 
head(hash_table, 3)
print("hash table loaded")

trait_chrom_pos = hash_table[trait_region_rs]

#remove empty elements 
trait_chrom_pos = trait_chrom_pos[!sapply(trait_chrom_pos,is.null)]
head(trait_chrom_pos, 3)

#read in tissue table from setup_config.R
tissueTable = read.table(file=QTL_tissue_table, sep=",", header=TRUE)

#Convert to strings from factors
tissueTable$Tissue = as.character(tissueTable$Tissue)
tissueTable$allPairsTabixFilename = as.character(tissueTable$allPairsTabixFilename)
tissueTable$sigPairsTabixFilename = as.character(tissueTable$sigPairsTabixFilename)

#Format Tissue names
tissueTable_tissue_noSpace = gsub("\\(","",tissueTable$Tissue)
tissueTable_tissue_noSpace = gsub("\\)","",tissueTable_tissue_noSpace)
tissueTable_tissue_noSpace = gsub("[[:space:]]","_",tissueTable_tissue_noSpace)
tissueTable_tissue_noSpace = gsub("-_", "", tissueTable_tissue_noSpace)
tissueTable$Tissue = tissueTable_tissue_noSpace

#create csv file from significant pair files
lead_SNP_pos = hash_table[[lead_SNP]]
#convert format for tabix
lead_SNP_pos_tabix_with_chr = paste0(gsub("_",":",lead_SNP_pos), "-", gsub("^.*?_","",lead_SNP_pos))
lead_SNP_pos_tabix_without_chr = paste0(gsub("chr", "", gsub("_",":",lead_SNP_pos)), "-", gsub("^.*?_","",lead_SNP_pos))

handle_tissue_table(tissueTable)

#read in the csv file of eGene-Tissue pairs 
eGenes = tryCatch({
    read.table(file=paste0(lead_SNP, ".csv"), sep="\t", stringsAsFactors=FALSE)
}, error = function(err) {
    print(paste(lead_SNP, "is not a significant", qtlType, "in any tissue in the", qtlType, "dataset."))
    quit(status=0)
})

#NK - call from before "adding rs ..."
get_gene_tissue_region <- function(){
  if(qtl_type == "eqtl"){
    return (eGeneTissue_region %>% dplyr::select(all_of(eQTL_all_chrom), all_of(eQTL_all_chromEnd), all_of(eQTL_all_geneID), all_of(eQTL_all_pvalue)))
  } else if(qtl_type == "sqtl"){
    return (eGeneTissue_region %>% dplyr::select(all_of(sQTL_all_chrom), all_of(sQTL_all_chromEnd), all_of(sQTL_all_geneID), all_of(sQTL_all_pvalue), all_of(sQTL_all_intron_chr), all_of(sQTL_all_intron_bp_first), all_of(sQTL_all_intron_bp_end), all_of(sQTL_all_intron_clu)))
  }
}

#NK
run_liftover <- function(){
  repeat {
    print("running liftOver")
    bed_liftover = data.frame("chr" = c(paste0("chr", chrom), paste0("chr", chrom)), "bp1" = c(colocStart - 1, colocStop - 1), "bp2" = c(colocStart, colocStop)) 
    write.table(bed_liftover,file=paste0("temp_hg19.bed"),sep="\t",quote = FALSE,row.names=FALSE,col.names=FALSE)
    
    #generate liftOver command
    liftOver_command = paste( "liftOver temp_hg19.bed", liftOver_chain, "temp_hg38.bed temp_hg19.unmapped -bedPlus=3 -tab", sep=" ")
    system(liftOver_command)
    
    hg38_positions = as.data.frame(read.table("temp_hg38.bed", header=FALSE, sep="\t")) 
    
    if (!is.na(hg38_positions[2,3])) {
      break
    } else {
      unmapped = fread(file="temp_hg19.unmapped", sep='\t', header=FALSE)
      if (unmapped[1,3] == colocStart) {
        print(paste("Could not map to hg38:", colocStart))
        colocStart = colocStart + 5000
      } else if (unmapped[1,3] == colocStop) {
        print(paste("Could not map to hg38:", colocStop))
        colocStop = colocStop - 5000
      }
      print(paste("Trying new region:", colocStart, "-", colocStop)) 
      
      if (colocStart >= colocStop) {
        print("Could not perform liftover")
        quit(status=0)
      }
    }
  }
}

#NK
format_tissue <- function(tissue){
  #parentheses are causing issues too
  tissue_noSpace = gsub("\\(","",tissue)
  tissue_noSpace = gsub("\\)","",tissue_noSpace)
  # The tissue names have any whitespace in them and we want to use these in the output file names so replace " " with "_"
  return(gsub("[[:space:]]","_",tissue_noSpace))
}

#NK
validate_build <-function(build){
  if (build == "hg38") {
    system(paste0("tabix ", tabix_allpair_path, " ", chrom, ":", colocStart, "-", colocStop, " >> ", eGeneTissueInputFile))
    system(paste0("tabix ", tabix_allpair_path, " chr", chrom, ":", colocStart, "-", colocStop, " >> ", eGeneTissueInputFile))
  } else if (build == "hg19") {
    system(paste0("tabix ", tabix_allpair_path, " ", chrom, ":", hg38_positions[1,3], "-", hg38_positions[2,3], " >> ", eGeneTissueInputFile))    
    system(paste0("tabix ", tabix_allpair_path, " chr", chrom, ":", hg38_positions[1,3], "-", hg38_positions[2,3], " >> ", eGeneTissueInputFile))    
  } else {
    print("ERROR: Please specify build: \"hg19\" or \"hg38\"")
    quit()
  }
}

#loop through the eGene-Tissue pairs in eGenes and prep running COLOC
for(i in 1:nrow(eGenes)){
    geneID <- eGenes[i, sig_geneID_col]
    
    geneID_noDOT <-  gsub("\\..*","", geneID)

    geneSymbol <- tryCatch({

        bitr(geneID_noDOT, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")$SYMBOL

    }, error = function(err){

        print("ENSEMBL ID could not be converted to HGNC Symbol")
        print(paste("geneSymbol will be set to the ENSEMBL ID",geneID_noDOT))

        return(geneID_noDOT)

    })

    
    print(geneID)
    print(geneSymbol)

    if (length(geneSymbol) > 1) {
        geneSymbol <- geneSymbol[1]
        print(paste("Using first gene symbol: ", geneSymbol))
    }

    tissue <- eGenes[i, length(eGenes)]
    print(tissue)

    # find the all pair file that contains the tissue of interest
    
    tissueLine <- tissueTable[tissueTable$Tissue == tissue,] 
    allpair_filename <- tissueLine$allPairsTabixFilename

    #if the Filename field is NA then skip this eGene-Tissue pair
    if (is.na(allpair_filename)){
        print(tissue)
        print("This tissue is not available in the all pairs files currently")
        next
    }
    
    tabix_allpair_path = paste0(all_qtl_tabix_dir, allpair_filename)
    qtl_N <- tissueLine$NumberRNASeqandGTSamples

    tissue_noSpace = format_tissue(tissue)
    
    # make a eGene-Tissue and trait prefix for file names
    if (qtlType == "eqtl") {
        out_prefix = paste(geneSymbol,geneID,tissue_noSpace,trait,sep="_")
    }	

    #run liftover on colocStart and colocStop if in HG19
    if (build == "hg19") {	
        run_liftover()
    }

    print("Grabbing the all pairs data")    
    #Use tabix to grab data, try both with and without "chr"
    eGeneTissueInputFile = paste(geneSymbol,tissue_noSpace,chrom,colocStart,colocStop,".txt", sep="_")
    validate_build(build)
    
    
    print("reading the all pairs data into R")
    #read the file we just generated from the grep command into R
    eGeneTissueInput = fread(file = eGeneTissueInputFile, header=FALSE)

    #check if tabix result is empty
    print(dim(eGeneTissueInput))
    if (dim(eGeneTissueInput)[1] == 0) {
        print("Did not find gene in all pairs data")
        if (qtlType == "eqtl") {
            quit(status=0)
        } else if (qtlType == "sqtl") {
            next
        }
    }	
    
    #remove the eGeneTissueInputFile after has been read into R to save disk space
    system(paste0("rm ", eGeneTissueInputFile, sep=""))

    colnames(eGeneTissueInput) = all_header

    print("Filtering on the geneID")
    #filter for the geneID of interest
    eGeneTissue_region = eGeneTissueInput[eGeneTissueInput[[all_geneID]] == geneID,]

    if (nrow(eGeneTissue_region) == 0) {
        print("Warning: There was not an exact match on Ensembl ID. Likely this is due to a GTEX version  update.")
          #make a string that removes the everything after in the geneID
          noDecimalGeneID = gsub("\\..*","",geneID)
  
          #recreated eGeneTissue_region
          eGeneTissue_region = eGeneTissueInput

          #grep the simplified Ensembl ID
        possible_Ensembl_gene_lines <- eGeneTissue_region[grepl(noDecimalGeneID, eGeneTissue_region[[all_geneID]]),]

        #check to make sure there is just one other Ensembl ID
        possible_Ensembl_genes <- unique(possible_Ensembl_gene_lines[[all_geneID]])
  
        if(length(possible_Ensembl_genes) == 1){
            print("Found a unique Ensembl ID so this analysis will continue using the Ensembl ID:")
            print(possible_Ensembl_genes)
    
            #this will be a single Ensembl ID string
            geneID = possible_Ensembl_genes
            eGeneTissue_region = eGeneTissue_region[eGeneTissue_region[[all_geneID]] == geneID,]
          } else {
            print("The Ensembl ID from your GTEx csv was not able to be reliably mapped to an Enseml ID in the GTEx database, so this gene will be skipped:")
            print(geneID)
            next
          }    
    }

    eGeneTissue_region <- get_gene_tissue_region()
    
    print("adding rs numbers to the QTL data")
    #add rs genegene,,numbers to the eGeneTissue_region DF

    #Add "chr" prefix to chromosome number if not present
    eGeneTissue_region[[all_chrom]] = paste0("chr", sub("chr", "", eGeneTissue_region[[all_chrom]]))

    #make chromosome_position column for merging
    eGeneTissue_region$chromosome_position <- paste(eGeneTissue_region[[all_chrom]],eGeneTissue_region[[all_chromEnd]],sep="_")

    #create data frame with rs numbers associated with chromosome_position and add to eGeneTissue region
    uniqID_DF = as.data.frame(t(as.data.frame(trait_chrom_pos)))
    uniqID_DF$SNP <- rownames(uniqID_DF)
    colnames(uniqID_DF) <- c("chromosome_position", "SNP")

    eGeneTissue_region = merge(eGeneTissue_region, uniqID_DF, by = "chromosome_position")  

    ################################ eQTL colocalization and RA Plots ################################ 
    if (qtlType == "eqtl") {
      eqtl_colocalization()
    ################################ sQTL colocalization and RA Plots ################################  
    } else if (qtlType == "sqtl") {
        sqtl_colocalization()
    }
}