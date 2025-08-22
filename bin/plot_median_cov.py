#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

fastq_file = sys.argv[1]      # File for FASTQ_TO_BAM
dragen_file = sys.argv[2]     # File for Dragen_BAM

df_fastq = pd.read_csv(fastq_file, sep=r"\s+", header=None, names=["Sample", "BEDTOOLS"])
df_dragen = pd.read_csv(dragen_file, sep=r"\s+", header=None, names=["Sample", "COVERVIEW"])

df = pd.merge(df_fastq, df_dragen, on="Sample")

plt.figure(figsize=(8, 6))
plt.scatter(df['BEDTOOLS'], df['COVERVIEW'], color='blue', edgecolor='k', s=80, label='Samples')

max_val = max(df['BEDTOOLS'].max(), df['COVERVIEW'].max())
plt.plot([0, max_val], [0, max_val], 'r--', label='y = x')

slope, intercept = np.polyfit(df['BEDTOOLS'], df['COVERVIEW'], 1)
fit_line = slope * df['BEDTOOLS'] + intercept
plt.plot(df['BEDTOOLS'], fit_line, color='green', label=f'Fit: y = {slope:.2f}x + {intercept:.2f}')

# add_line_x = [ x for x in range(-40, int(df['BEDTOOLS'].min())) ]
# add_line_y =list()
# for values in add_line_x:
# 	add_line_y.append(slope * values + intercept)
# # print (add_line_x, add_line_y)
# plt.plot(add_line_x, add_line_y, color='green')

corr = df['COVERVIEW'].corr(df['BEDTOOLS'])

plt.axis('square')
plt.xlabel('BEDTOOLS Median Coverage')
plt.ylabel('COVERVIEW Median Coverage')
plt.title(f'Scatter Plot - Coverview\nPearson r = {corr:.3f}')
plt.legend()

plt.tight_layout()
plt.savefig("BedtoolsVsCoverview", dpi=300)
# print("Scatter plot saved to scatter_plot_coverview.png")