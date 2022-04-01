#!/bin/bash
# Runs colocquial.R 
#This version is designed to use an LSF job submission to parallelize the coloc jobs
# Input files needed in directory: qtl_config.sh (modified) 

module load plink/1.90Beta4.5
module load R/3.6.3
module load bedtools2

echo "modules loaded"

#source the config file for this run (parameters such as GWAS file location and column names) and the setup config (paths to dependency files for the pipeline)
source qtl_config.sh
source $setup_config_sh

#check for changes to the default plink parameters
if [ "$clumpP1" == "" ]
then
    clumpP1=0.0000001
fi

if [ "$clumpKB" == "" ]
then
    clumpKB=1000
fi

if [ "$clumpR2" == "" ]
then
    clumpR2=0.2
fi

if [ "$leadSNPsFilePath" == "" ]
then
	echo "No lead SNPs file specified. Running PLINK to find significant SNPs"

    #calculate half of the window, so we can generate a window centered around the lead SNP
    radius=`echo $((window/2))`

	#Perform LD-clumping from GWAS summary statstics to  get clumps/loci of significant SNPs
    plink --noweb --bfile $plink_bfile  --keep $plink_keep --clump-p1 $clumpP1 --clump-r2 $clumpR2  --clump-kb $clumpKB --clump $traitFilePath --clump-snp-field $trait_SNPcol --clump-field $trait_Pcol --out allChrMergedClumped

	awk '{OFS="\t"} (NR>1 && NF>0) {print "chr"$1,$4,$4+1,$3}' allChrMergedClumped.clumped | sort -k1,1 -k2,2n > allChrMergedClumped.bed

	#create file with rsid, CHR, BP, Start, Stop as columns
	awk -v var="$radius" '{print $4"\t"$1"\t"$2"\t"$2-var"\t"$2+var}' allChrMergedClumped.bed | sed 's/chr//' > $trait"_lead_SNPs.txt"

	leadSNPsFilePath=$trait"_lead_SNPs.txt"
fi

#create bsub template for trait
sed "s/TRAITNAME/$trait/" $colocquial_dir/qtl_coloc_template.bsub > $trait"_template.bsub" 

#Add fields from qtl_config.sh to QTL_config_template.R for trait
sed "s/TRAITNAME/$trait/" $colocquial_dir/QTL_config_template.R | sed "s|TRAITPATH|$traitFilePath|" | sed "s/A1COL/$trait_A1col/" | sed "s/A2COL/$trait_A2col/" | sed "s/SNPCOL/$trait_SNPcol/" | sed "s/CHRCOL/$trait_CHRcol/" | sed "s/BPCOL/$trait_BPcol/" | sed "s/PCOL/$trait_Pcol/" | sed "s/NCOL/$trait_Ncol/" | sed "s/MAFCOL/$trait_MAFcol/" | sed "s/TRAITTYPE/$traitType/" | sed "s/TRAITPROP/$traitProp/" | sed "s/BUILD/$build/" | sed "s/QTLTYPE/$qtlType/" | sed "s/CLUMPP1/$clumpP1/" | sed "s/CLUMPKB/$clumpKB/" | sed "s/CLUMPR2/$clumpR2/" | sed "s:SETUPCONFIGR:$setup_config_R:" > $trait"_QTL_config_template.R"

#for each lead SNP
cat $leadSNPsFilePath | while read line
do 
	echo $line

	#parse the SNP
	SNP=`echo $line | cut -f 1 -d " "`
	echo $SNP

	#parse chr start stop for coloc
	CHR=`echo $line | cut -f 2 -d " "`
	echo $CHR

	Start=`echo $line | cut -f 4 -d " "`
	echo $Start

	Stop=`echo $line | cut -f 5 -d " "`
	echo $Stop

	#make the analysis directory
	mkdir $SNP

	#cd into the dir
	cd $SNP

	#copy the colocquial.R file into the dir
	cp $colocquial_dir/colocquial.R ./

	#add an edited version of the bsub file to the dir
	sed "s/SNPNUMBER/$SNP/" ../$trait"_template.bsub" > ./$SNP".bsub"

	bsubfile=$SNP".bsub"
	echo $bsubfile

	#add an edited version of the QTL_config.R file to the dir (edit gtex file path, chr, start, and stop with sed)
	sed "s/SNPNUMBER/$SNP/" ../$trait"_QTL_config_template.R" | sed "s/CHROMOSOME/$CHR/" | sed "s/STARTBP/$Start/" | sed "s/STOPBP/$Stop/" > ./QTL_config.R

	#fire off the bsub job for this lead SNP
	bsub < $bsubfile -q $bsub_queue -R "rusage[mem=16GB]" -M 16G

	#cd back into the main directory to go to the next SNP 
	cd ..

	echo
	echo
    
done

echo "all lead SNP jobs have been submitted"

#replace "TRAITNAME" with the $trait in out and err file names for summary file bsub
sed "s/TRAITNAME/$trait/" $colocquial_dir/summarize_results.bsub | sed "s|COLOCQUIAL_DIR|$colocquial_dir|" > ./summarize_results.bsub

#run bsub to collect all of the COLOC results into 1 file
bsub < summarize_results.bsub -q $bsub_queue -R "rusage[mem=16GB]" -M 16G

