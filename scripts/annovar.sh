#!/usr/bin/bash
# This script will take a vcf as input and run annovar on it

annovar_path=$1
vcf_infile=$2
Sample=$3

perl ${annovar_path}/convert2annovar.pl -format vcf4 ${vcf_infile} --outfile ${Sample}.avinput -withfreq --includeinfo -allsample
perl ${annovar_path}/table_annovar.pl ${Sample}.avinput --out ${Sample}.annovar --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${annovar_path}/humandb/ --xreffile ${annovar_path}/example/gene_fullxref.txt