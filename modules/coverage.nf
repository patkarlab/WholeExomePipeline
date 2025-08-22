#!/usr/bin/env nextflow

process COVERAGE_MOSDEPTH {
   tag "${Sample}"
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*'
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*")
		tuple val (Sample), file ("${Sample}_cov.mosdepth.summary.txt"), file("${Sample}_cov.regions.bed"), file("${Sample}_median50"), emit: mosdepth_out
	script:
	"""
	${params.mosdepth_script} ${finalBam} ${Sample}_cov ${params.bedfile}.bed
	${params.extract_COV50_script_path} ${Sample}_cov.regions.bed ${Sample}_median50
	sleep 2s
	"""
}

process COVERAGE_BEDTOOLS {
	tag "${Sample}"
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*'
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("${Sample}_cov.regions_bedtools.bed")
	script:
	"""
	${params.bedtools} bamtobed -i ${finalBam} > ${Sample}_temp.bed
	${params.bedtools} coverage -counts -a ${params.bedfile}.bed -b ${Sample}_temp.bed > "${Sample}_cov.regions_bedtools.bed"
	"""
}

process COVERAGE_COVERVIEW {
	tag "${Sample}"
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: "*.coverview_regions.csv"
	input:
		tuple val (Sample), file(finalBam), file(finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.coverview_regions.csv")
	script:
	"""
	${params.coverview_path}/coverview -i ${finalBam} -b ${params.bedfile}.bed -c ${params.coverview_path}/config/config.txt -o ${Sample}.coverview
	python3 ${params.coverview_script_path} ${Sample}.coverview_regions.txt ${Sample}.coverview_regions.csv
	"""
}