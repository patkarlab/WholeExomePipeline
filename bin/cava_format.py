#!/usr/bin/env python3

import pandas as pd
import sys
import tempfile

args = sys.argv

file1 = args[1]
outfile = args[2]

# Determine expected number of columns from the header
with open(file1) as f:
    header = f.readline().rstrip("\n").split("\t")
expected_cols = len(header)

# Create a temporary corrected file
with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp:
    with open(file1) as fin:
        for i, line in enumerate(fin):
            fields = line.rstrip("\n").split("\t")

            # Keep header unchanged
            if i == 0:
                tmp.write(line)
                continue

            # Merge any extra fields into the last column
            if len(fields) > expected_cols:
                fields = fields[:expected_cols-1] + [",".join(fields[expected_cols-1:])]

            tmp.write("\t".join(fields) + "\n")

    temp_file = tmp.name

# Read corrected file
df1 = pd.read_csv(temp_file, sep="\t")

# Replace '.' with '-1'
df1.replace(to_replace='.', value='-1', inplace=True)

# Write output
df1.to_csv(outfile, index=False)
