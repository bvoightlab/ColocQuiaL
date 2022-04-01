#!/bin/bash
#Code to collect the QTL colocalization results after they have all been generated.

#set variables for date and time
today=`date +%F`
now=`date +%T`

pwd

#load config file
source qtl_config.sh
source $setup_config_sh

#write header to the summary results file
printf "SNP\tGene\tGeneID-Tissue\tTrait\tPPID\tPP\n" > $trait"_"$qtlType"_coloc_results_all_summary_"$today"_"$now".txt" 

#go through each lead SNP directory
for dir in ./rs*/
do
    #echo $dir
    dirname=`echo $dir | tr "/" "\t" | cut -f 2`
    echo $dirname

    cd  $dir

    #for each GTEx coloc output in the directory
    for colocOut in ./*coloc_results_summary.txt
    do
        #echo $colocOut
        colocfilename=`basename $colocOut`
        #echo $colocfilename

        #save a "trait" string for formatting
        trait_str="_"$trait"_"

        #grab the information we need from the locus coloc results file
        sed "s/^/$colocfilename /" $colocOut | sed "s/^/$dirname /" | sed "s,$trait_str, $trait," | sed "s/coloc_results_summary.txt//" | sed "s/_ENSG/ ENSG/"| tr " " "\t" | tail -n+3 >> ../$trait"_"$qtlType"_coloc_results_all_summary_"$today"_"$now".txt"

    done

    cd ..

done

#if using sQTL data you need to parse out the gene name
if [ "$qtlType" == "sqtl" ]
then

     sed -i "s/_/\\t/" $trait"_"$qtlType"_coloc_results_all_summary_"$today"_"$now".txt"

fi

#run Rscript to generate a file with PP3, PP4, and PP4/(PP3 + PP4) for each lead SNP-Gene-Tissue result.
Rscript $colocquial_dir"/summarize_qtl_coloc_PP3_PP4_results.R" $trait"_"$qtlType"_coloc_results_all_summary_"$today"_"$now".txt"

echo "Your summary QTL file is ready!" 




