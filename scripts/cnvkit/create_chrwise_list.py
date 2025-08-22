#! /usr/bin/env python3

import sys
import re

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} <input_bed> <output_file>")
    sys.exit(1)

bed_file = sys.argv[1]
output_file = sys.argv[2]

# Function for natural sorting of chromosome names
def natural_sort_key(chrom):
    # Extract numeric part if present
    match = re.match(r"(chr)?(\d+|X|Y|M|MT)$", chrom)
    if match:
        val = match.group(2)
        if val.isdigit():
            return (0, int(val))  # numeric chromosomes first
        elif val in ["X", "Y"]:
            return (1, val)
        else:
            return (2, val)
    return (3, chrom)  # non-standard chroms last

# Read BED file
chrom_dict = {}
label_dict = {}

with open(bed_file) as f:
    for line in f:
        if not line.strip():
            continue
        chrom, start, end, label = line.strip().split("\t")
        
        # Keep only the last part after semicolon in label
        label = label.split(";")[-1]

        key = f"{chrom}:{start}-{end}"
        if chrom not in chrom_dict:
            chrom_dict[chrom] = []
            label_dict[chrom] = []
        chrom_dict[chrom].append(key)
        label_dict[chrom].append(label)

# Sort chromosomes naturally
sorted_chroms = sorted(chrom_dict.keys(), key=natural_sort_key)

# Write output with max 70 entries per line
with open(output_file, "w") as out:
    for chrom in sorted_chroms:
        coords = chrom_dict[chrom]
        labels = label_dict[chrom]
        for i in range(0, len(coords), 70):
            chunk_coords = coords[i:i+70]
            chunk_labels = labels[i:i+70]
            out.write(",".join(chunk_coords) + " " + ",".join(chunk_labels) + "\n")

print(f"Output saved to {output_file}")
