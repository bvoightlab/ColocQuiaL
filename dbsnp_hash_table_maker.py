#!/usr/bin/python2.7
#dbsnp_hash_table_maker.py

#uses pandas commands to take in downloaded BED files and make hash tables with rs numbers as keys and GRCh38 chromosome positions as values
#designed to take bed_chr_#.bed.gz files and output json hash table files in same directory

import pandas as pd
import json

def main():
    chr_list = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "AltOnly", "MT", "Multi", "NotOn", "PAR", "Un", "X", "Y"] 

    for chrom in chr_list:

        print chrom
        
        inputfileStr = "bed_chr_" + chrom + ".bed.gz"
        
        #load the zipped BED file into python as a pandas DF unless empty
        try:
            dbsnp_data = pd.read_table(inputfileStr, low_memory = False, sep = "\t", header = None, skiprows = 1, compression = 'gzip')
        except pd.io.common.EmptyDataError:
            continue

        #add a header
        dbsnp_data.columns = ["chrom", "chromStart", "chromEnd", "SNP", "score", "strand"]
        
        #conver chromEnd to a string
        dbsnp_data['chromEnd'] = dbsnp_data['chromEnd'].apply(str)

        #add column with chromosome and position
        dbsnp_data['chrom_chromEnd'] = dbsnp_data[['chrom', 'chromEnd']].apply(lambda x: '_'.join(x), axis = 1)

        #create hash table
        rsToChrPosition = pd.Series(dbsnp_data.chrom_chromEnd.values, index = dbsnp_data.SNP).to_dict()

        #save to json file
        jsonHashTable = json.dumps(rsToChrPosition)        
        jsonfile = open("chr_" + chrom + "_snp151_hash_table.json", "w")
        jsonfile.write(jsonHashTable)
        jsonfile.close()

if __name__ == "__main__":
    main()
