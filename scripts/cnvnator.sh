#! /usr/bin/bash

source activate ROOT

final_bam=$1
outdir=$2
sample=$3
program_dir="/home/programs/CNVnator/CNVnator"

${program_dir}/cnvnator -root out.root -tree ${final_bam}
${program_dir}/cnvnator -root out.root -his 1000 -d ${program_dir}/chr_fasta/
${program_dir}/cnvnator -root out.root -stat 1000
${program_dir}/cnvnator -root out.root -partition 1000
${program_dir}/cnvnator -root out.root -call 1000 > ${outdir}/${sample}_Cnvnator_calls.tsv

conda deactivate 
