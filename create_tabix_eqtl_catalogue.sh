# create_tabix_eqtl_catalogue.sh
# by Brian Chen
# This script creates tabix index files for significant pairs files generated from eQTL Catalogue datasets
# used in make_significant_pairs_files.R
# $1 is path to tissue file

module load tabix
module load vcftools

# Remove .tsv.gz
tissueName=$(echo $1 | cut -f 1,2 -d '.')
echo $tissueName

# Open file | Remove first row (header) | Create BED format | Sort by position
zcat $1 | sed '1D' | awk 'BEGIN {OFS = "\t"} {$1 = "chr"$2"\t"($3 - 1)"\t"$3"\t"$1; print;}' | bgzip -c > $tissueName".tab.gz"

echo $tissueName".tab.gz created"

# create index file
tabix -p bed $tissueName".tab.gz"

