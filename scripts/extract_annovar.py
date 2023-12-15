#! /usr/bin/env python3
import pandas as pd
import sys


input_1 =  sys.argv[1]  #23RSEQ420-Exome.somaticseq.hg19_multianno.csv
input_2 =  sys.argv[2]  #23RSEQ420-Exome_somaticseq.vep.txt
outputFile = sys.argv[3] #somatics_vep_annonvar.txt

data1 = pd.read_csv(input_1,sep = ',')
data2 = pd.read_csv(input_2,sep = '\t')
extractedData = data1[['Chr', 'Start', 'Ref','Alt', 'Func.refGene', 'Gene.refGene', 'ExonicFunc.refGene', 'PopFreqMax', 'InterVar_automated']]

extractedData.rename(columns = {'Chr':'CHROM', 'Ref':'REF', 'Alt':'ALT'}, inplace = True)
#print(extractedData)

merge = data2.merge(extractedData, how = 'inner', left_on = ['CHROM','Start','REF','ALT'], right_on = ['CHROM','Start','REF','ALT'])
#print(merge)

merge.to_csv(outputFile, sep = '\t', header=True, index=False)

