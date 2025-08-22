#!/usr/bin/env nextflow

process HSMETRICS {
    tag "${Sample}"
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*_hsmetrics.txt'
	input:
		tuple val(Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*_hsmetrics.txt")
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} CollectHsMetrics I= ${finalBam} O= ${Sample}_hsmetrics.txt BAIT_INTERVALS= ${params.bedfile}.interval_list TARGET_INTERVALS= ${params.bedfile}.interval_list R= ${params.genome} VALIDATION_STRINGENCY=LENIENT
	${params.hsmetrics_all} ${params.output}/hsmetrics.tsv ${Sample} ${Sample}_hsmetrics.txt
	"""
}