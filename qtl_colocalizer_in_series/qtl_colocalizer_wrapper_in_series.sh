#!/bin/bash
# Runs qtl_colocalizer.R on a set of lead SNPs in series
# This version is intended for user who do not have cluster computing resourcesa
# Input files needed in directory: qtl_config.sh (modified), qtl_coloc_template.bsub, QTL_config_template.R, qtl_colocalizer.R

module load plink/1.90Beta4.5
module load R/3.6.3
module load bedtools2
module load tabix
module load liftOver

echo "modules loaded"

source qtl_config.sh
source $setup_config_sh

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
#sed "s/TRAITNAME/$trait/" qtl_coloc_template.bsub > $trait"_template.bsub" 

#Add fields from qtl_config.sh to QTL_config_template.R for trait
sed "s/TRAITNAME/$trait/" QTL_config_template.R | sed "s|TRAITPATH|$traitFilePath|" | sed "s/A1COL/$trait_A1col/" | sed "s/A2COL/$trait_A2col/" | sed "s/SNPCOL/$trait_SNPcol/" | sed "s/CHRCOL/$trait_CHRcol/" | sed "s/BPCOL/$trait_BPcol/" | sed "s/PCOL/$trait_Pcol/" | sed "s/NCOL/$trait_Ncol/" | sed "s/MAFCOL/$trait_MAFcol/" | sed "s/TRAITTYPE/$traitType/" | sed "s/TRAITPROP/$traitProp/" | sed "s/QTL/$qtlType/" > $trait"_QTL_config_template.R" 

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

	#copy the qtl_colocalizer.R file into the dir
	cp ../qtl_colocalizer.R ./

	#add an edited version of the bsub file to the dir
	#sed "s/SNPNUMBER/$SNP/" ../$trait"_template.bsub" > ./$SNP".bsub"

	#bsubfile=$SNP".bsub"
	#echo $bsubfile

	#add an edited version of the QTL_config.R file to the dir (edit gtex file path, chr, start, and stop with sed)
	sed "s/SNPNUMBER/$SNP/" ../$trait"_QTL_config_template.R" | sed "s/CHROMOSOME/$CHR/" | sed "s/STARTBP/$Start/" | sed "s/STOPBP/$Stop/" > ./QTL_config.R

	#Run qtl_colocalizer.R for this lead SNP
    Rscript ./qtl_colocalizer.R

    
	#cd back into the main directory to go to the next SNP 
	cd ..

	echo
	echo
    
done

echo "all lead SNP COLOC Analyses are complete!"

echo "Collecting the results into a single file"

#run the summary results code
./summarize_qtl_results.sh

