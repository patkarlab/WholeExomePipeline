#! /usr/bin/bash

samplesheet=$1

for i in `cat ${samplesheet}` 
do 
	#java -jar /home/programs/picard/build/libs/picard.jar CollectHsMetrics I=/home/diagnostics/pipelines/WholeExomePipeline/Final_Output/$i/$i".final.bam" O=/home/diagnostics/pipelines/WholeExomePipeline/Final_Output/$i/$i".hsmetrics.txt" R=/home/reference_genomes/hg19_broad/hg19_all.fasta BAIT_INTERVALS=/home/diagnostics/pipelines/WholeExomePipeline/bedfiles/Exome_hg19_sortd.interval_list TARGET_INTERVALS=/home/diagnostics/pipelines/WholeExomePipeline/bedfiles/Exome_hg19_sortd.interval_list VALIDATION_STRINGENCY=LENIENT

	echo -ne $i'\t'; grep -v '#' /home/diagnostics/pipelines/WholeExomePipeline/Final_Output/$i/$i.hsmetrics.txt | awk 'BEGIN{FS="\t"; OFS="\t"}NR==3{ print $7,$8}'

done
