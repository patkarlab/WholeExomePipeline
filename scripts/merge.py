#! /usr/bin/env python3

import pandas as pd
import os
import sys

sample = sys.argv[1]
outfile = sys.argv[2]
cava_path = sys.argv[3]
mosdepth_path = sys.argv[4]
mosdepthCov_path = sys.argv[5]
mosCOV50_path = sys.argv[6]
pindel_path = sys.argv[7]
somaticseq_path = sys.argv[8]
pharma = sys.argv[9]
deepvariant = sys.argv[10]

csvfilenames = [ ]

cava_df = pd.read_csv(cava_path,sep =',')
mosdepth_df = pd.read_csv(mosdepth_path, sep='\t')
mosdepthCov_df = pd.read_csv(mosdepthCov_path, sep='\t', names=['Chrom', 'Start', 'End', 'Median'])
mosCOV50_df = pd.read_csv(mosCOV50_path, sep='\t')
pindel_df = pd.read_csv(pindel_path, sep='\t')
somaticseq_df = pd.read_csv(somaticseq_path, sep='\t')
pharma_df = pd.read_csv(pharma, sep='\t')
deepvar_df = pd.read_csv(deepvariant, sep = '\t')

# Create a new Excel writer
writer = pd.ExcelWriter(outfile)

# Write data frames to separate sheets
cava_df.to_excel(writer, sheet_name='cava', index=False)
mosdepth_df.to_excel(writer, sheet_name='mosdepth', index=False)
mosdepthCov_df.to_excel(writer, sheet_name='mosdepthCov', index=False)
mosCOV50_df.to_excel(writer, sheet_name='mosCOV50', index=False)
pindel_df.to_excel(writer, sheet_name='pindel', index=False)
somaticseq_df.to_excel(writer, sheet_name='somaticseq', index=False)
pharma_df.to_excel(writer, sheet_name='pharma_marker', index=False)
deepvar_df.to_excel(writer, sheet_name='deepvariant', index=False) 

# Save the Excel file
writer.save()
