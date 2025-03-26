#!/usr/bin/bash

BIN_VERSION="1.6.1"
YOUR_INPUT_DIR=$1
YOUR_OUTPUT_DIR=$2
YOUR_OUTPUT_VCF=$3
YOUR_BAM=$4
genome=$5
bedfile=$6

genome_name=$(basename ${genome})
genome_location=$(echo ${genome} | awk 'BEGIN{OFS=FS="/"} {$NF=""; print $0}')
bedfile_location=$(dirname ${bedfile})
bedfile_name=$(basename ${bedfile})
#echo "${genome_name} ${genome_location}"

#docker run \
#  -v ${YOUR_INPUT_DIR}:"/input" \
#  -v ${YOUR_OUTPUT_DIR}:"/output" \
#  -v ${genome_location}:"/genome_dir" \
#  google/deepvariant:"${BIN_VERSION}" \
#  /opt/deepvariant/bin/run_deepvariant \
#  --model_type=WES #**Replace this string with exactly one of the following [WGS,WES,PACBIO,ONT_R104,HYBRID_PACBIO_ILLUMINA]**
#  --ref=/genome_dir/${genome_name} \
#  --reads=/input/${YOUR_BAM} \
#  --output_vcf=/output/${YOUR_OUTPUT_VCF}
#  --output_gvcf=/output/YOUR_OUTPUT_GVCF \
#  --num_shards=16 #\ **This will use all your cores to run make_examples. Feel free to change.**
#  --logging_dir=/output/logs  # **Optional. This saves the log output for each stage separately.
#  --haploid_contigs="chrX,chrY" \ **Optional. Heterozygous variants in these contigs will be re-genotyped as the most likely of reference or homozygous alternates. For a sample with karyotype XY, it should be set to "chrX,chrY" for GRCh38 and "X,Y" for GRCh37. For a sample with karyotype XX, this should not be used.
#  --par_regions_bed="/input/GRCh3X_par.bed" \ **Optional. If --haploid_contigs is set, then this can be used to provide PAR regions to be excluded from genotype adjustment. Download links to this files are available in this page.
 # --dry_run=false **Default is false. If set to true, commands will be printed out but not executed.

docker run -v ${YOUR_INPUT_DIR}:"/input" -v ${YOUR_OUTPUT_DIR}:"/output" -v ${genome_location}:"/genome_dir" -v ${bedfile_location}:"/bedfile"  google/deepvariant:"${BIN_VERSION}" /opt/deepvariant/bin/run_deepvariant --model_type=WES --ref=/genome_dir/${genome_name} --regions=/bedfile/${bedfile_name} --reads=/input/${YOUR_BAM} --output_vcf=/output/${YOUR_OUTPUT_VCF} --num_shards=16
