#!/usr/bin/env python3
import sys
import pandas as pd

def vcf_to_tsv(vcf_file, tsv_file):
    records = []
    info_fields = set()
    format_fields = set()

    with open(vcf_file) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                sample_name = header[-1] if len(header) > 9 else None
                continue

            fields = line.strip().split("\t")
            chrom, pos, vid, ref, alt, qual, filt, info = fields[:8]
            fmt = fields[8] if len(fields) > 8 else None
            sample = fields[9] if len(fields) > 9 else None

            entry = {
                "CHROM": chrom,
                "POS": pos,
                "ID": vid,
                "REF": ref,
                "ALT": alt,
                "QUAL": qual,
                "FILTER": filt,
            }

            # Parse INFO fields
            for kv in info.split(";"):
                if "=" in kv:
                    key, value = kv.split("=", 1)
                else:
                    key, value = kv, True
                entry[key] = value
                info_fields.add(key)

            # Parse FORMAT/sample fields (skip GT)
            if fmt and sample:
                fmt_tags = fmt.split(":")
                sample_values = sample.split(":")
                for tag, val in zip(fmt_tags, sample_values):
                    if tag == "GT":  # skip genotype field
                        continue
                    entry[tag] = val
                    format_fields.add(tag)

            records.append(entry)

    # Build DataFrame
    df = pd.DataFrame(records)

    # Ensure all INFO and FORMAT fields are present
    for key in sorted(info_fields | format_fields):
        if key not in df.columns:
            df[key] = "."

    # Write TSV
    df.to_csv(tsv_file, sep="\t", index=False)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.stderr.write(f"Usage: {sys.argv[0]} input.vcf output.tsv\n")
        sys.exit(1)
    vcf_to_tsv(sys.argv[1], sys.argv[2])
