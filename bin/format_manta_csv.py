#!/usr/bin/env python3

import sys
import gzip

def open_vcf(vcf_file):
    if vcf_file.endswith(".gz"):
        return gzip.open(vcf_file, "rt")
    else:
        return open(vcf_file, "r")

def parse_info(info_str):
    info_dict = {}
    for item in info_str.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            info_dict[key] = val
        else:
            info_dict[item] = True
    return info_dict

def extract_alt_coordinates(alt):
    # Example: GGGGG[chr14:106641561[
    import re
    match = re.search(r'(chr[\w]+):(\d+)', alt)
    if match:
        return match.group(1), match.group(2)
    return "NA", "NA"

def parse_format(format_col, sample_col):
    keys = format_col.split(":")
    values = sample_col.split(":")

    fmt = dict(zip(keys, values))

    def get_alt_count(field):
        if field not in fmt:
            return 0
        val = fmt[field]
        if "," in val:
            return int(val.split(",")[1])  # ALT count
        return 0

    pr_alt = get_alt_count("PR")
    sr_alt = get_alt_count("SR")

    return pr_alt, sr_alt

def main(vcf_file, out_file):

    with open_vcf(vcf_file) as f, open(out_file, "w") as out:

        header = [
            "CHROM", "POS", "ID", "SVTYPE",
            "CHR2", "POS2",
            "PR_normal", "SR_normal",
            "PR_tumor", "SR_tumor",
            "SOMATICSCORE", "FILTER"
        ]
        out.write("\t".join(header) + "\n")

        for line in f:
            if line.startswith("#"):
                continue

            cols = line.strip().split("\t")

            chrom = cols[0]
            pos = cols[1]
            vid = cols[2]
            ref = cols[3]
            alt = cols[4]
            filt = cols[6]
            info = cols[7]
            fmt = cols[8]

            normal = cols[9]
            tumor = cols[10]

            info_dict = parse_info(info)

            svtype = info_dict.get("SVTYPE", "NA")
            som_score = info_dict.get("SOMATICSCORE", "NA")

            chr2, pos2 = extract_alt_coordinates(alt)

            pr_n, sr_n = parse_format(fmt, normal)
            pr_t, sr_t = parse_format(fmt, tumor)

            out.write("\t".join(map(str, [
                chrom, pos, vid, svtype,
                chr2, pos2,
                pr_n, sr_n,
                pr_t, sr_t,
                som_score, filt
            ])) + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python newcsv.py input.vcf[.gz] output.tsv")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
