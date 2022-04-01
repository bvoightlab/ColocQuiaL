# create_tabix_gtex_sqtl_sigpairs.sh
# by Brian Chen
# This script takes GTEx v8 sQTL significant pair files and outputs index file for tabix
# $1 is path to tissue file

module load tabix
module load vcftools

tissueName=$(echo "${1##*/}" | cut -f 1,2 -d '.')
echo $tissueName

cat $1 | sed '1D' | tr '_' '\t' | tr ':' '\t' | awk 'BEGIN {OFS = "\t"} {$2 = ($2 - 1)"\t"$2; print;}' | vcf-sort -p 8 | bgzip -c > $tissueName".signif_variant_gene_pairs.tab.gz"

echo $tissueName".signif_variant_gene_pairs.tab.gz created"

#create index file
tabix -p bed $tissueName".signif_variant_gene_pairs.tab.gz" 
