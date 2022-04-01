# create_tabix_gtex_sqtl_allpairs.sh
# by Brian Chen
# This script takes GTEx v8 sQTL all pair files and outputs index file for tabix
# $1 is path to tissue file

module load tabix
module load vcftools

tissueName=$(echo "${1##*/}" | cut -f 1,2 -d '.')
echo $tissueName

zcat $1 | sed '1D' | tr '_' '\t' | tr ':' '\t' | awk 'BEGIN {OFS = "\t"} {t1 = $1; t2 = $2; t3 = $3; t4 = $4; t5 = $5; $1 = $7; $2 = $8 - 1; $3 = $8; $4 = $6; $5 = t1; $6 = t2; $7 = t3; $8 = t4"_"t5; print;}' | vcf-sort -t /scratch -p 8 | bgzip -c > $tissueName".tab.gz"

echo $tissueName".tab.gz created"

#create index file
tabix -p bed $tissueName".tab.gz" 

