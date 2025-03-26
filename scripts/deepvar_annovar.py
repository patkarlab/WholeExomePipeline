#! /usr/bin/env python3
import pandas as pd
import sys

input_1 =  sys.argv[1]		#23RSEQ420-Exome.hg19_multianno.csv
input_2 =  sys.argv[2]		#23RSEQ420-Exome_deepvariant.txt
outputFile = sys.argv[3]	#_deepvar_vep_annonvar.txt

data1 = pd.read_csv(input_1,sep = ',')
data2 = pd.read_csv(input_2,sep = '\t')

extractedData = data1[['Chr', 'Start', 'Ref','Alt', 'Func.refGene', 'Gene.refGene', 'ExonicFunc.refGene', 'PopFreqMax', 'InterVar_automated']]
extractedData = extractedData.rename(columns = {'Chr':'CHROM', 'Start':'POS', 'Ref':'REF', 'Alt':'ALT'})

merge = data2.merge(extractedData, how = 'inner', left_on = ['CHROM','POS','REF','ALT'], right_on = ['CHROM','POS','REF','ALT'])
merge.to_csv(outputFile, sep = '\t', header=True, index=False)