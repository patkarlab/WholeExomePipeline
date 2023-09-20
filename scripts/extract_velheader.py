#! /usr/bin/env python3
import pandas as pd
import sys

inputFile = sys.argv[1]
outputFile = sys.argv[2]

data = pd.read_csv(inputFile,sep = '\t')



data[['CHROM','POS']] = data['Location'].str.split(':', expand=True)
data[['Start','End']] = data['POS'].str.split('-', expand=True)

extractedData = data[[ 'CHROM','Start', 'End',
                       'SYMBOL',
                       'AF',
                       'MAX_AF_POPS',
                       'VARIANT_CLASS',
                       'CLIN_SIG',
                       'Feature',
                       'Feature_type',
					   'HGVSc',
					   'HGVSp',
                       'Consequence',
                       'cDNA_position',
                       'CDS_position',
                       'Protein_position',
                       'Amino_acids',
                       'Codons',
                       'ENSP',
                       'ALLELE_NUM'
                      ]]

extractedData.to_csv(outputFile, sep = ',', header=True, index=False)
