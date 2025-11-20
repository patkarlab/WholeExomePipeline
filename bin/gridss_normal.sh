#!/usr/bin/bash

# Usage:
# ./run_gridss_somatic.sh <ref_fasta> <tumor_bam> <normal_bam> <output_prefix>

ref_fasta=$1
tumor_bam=$2
normal_bam=$3
output_prefix=$4

export R_LIBS="/home/tuhina/R/x86_64-pc-linux-gnu-library/3.6:$R_LIBS"
export R_LIBS="/home/vishram/R/x86_64-pc-linux-gnu-library/3.6:$R_LIBS"
export R_LIBS="/home/arpit/lib/R/library:$R_LIBS"


export PATH="/home/programs/samtools-1.15.1:$PATH"
export PATH="/home/programs/gridss:$PATH"
gridss_path="/home/programs/gridss"

# Output filenames
output_vcf="${output_prefix}.gridss.vcf"
somatic_vcf="${output_prefix}.somatic.gridss.vcf"

# (Optional) Exclude problematic regions such as centromeres or chrM
# sed 's/^/chr/g' ENCFF001TDO.bed | grep -v 'chrM' > exclude_list.bed

# Step 1: Call structural variants (jointly on tumor and normal)
echo "Running GRIDSS SV calling on tumor + normal..."
gridss -r ${ref_fasta} \
    -j ${gridss_path}/gridss-2.13.2-gridss-jar-with-dependencies.jar \
    -t 24 \
    -o ${output_vcf} \
    -b ${gridss_path}/exclude_list.bed \
    ${tumor_bam} ${normal_bam}

# Step 2: Filter for somatic events (tumor-specific)
# GRIDSS provides the SomaticFilter command for this
echo "Filtering somatic SVs..."
/home/arpit/bin/Rscript /home/programs/gridss/gridss_somatic_filter \
    -r BSgenome.Hsapiens.UCSC.hg19 \
    -i ${output_vcf} \
    -o ${somatic_vcf} \
    -n 2 \
    -t 1 \
    --scriptdir ${gridss_path}

bgzip -d ${somatic_vcf}.bgz

# Step 3: Annotate structural variants (optional)
echo "Annotating variants..."
gridss_annotate_vcf_repeatmasker \
    -j ${gridss_path}/gridss-2.13.2-gridss-jar-with-dependencies.jar \
    -t 24 \
    -o RM_${somatic_vcf} \
    ${somatic_vcf}

echo "Done! Somatic SV VCF: RM_${somatic_vcf}"
