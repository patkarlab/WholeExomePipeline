#!/usr/bin/env nextflow

process DEEPVARIANT_GCNV {
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*.vcf'
	input:
		tuple val (Sample), file (DeepVarOut), file (DeepVarVcf), file(genotypedIntervals), file (genotypedSegments), file(copyRatios)

	output:
		tuple val (Sample), file("${Sample}_deepvar_gCNV_merged.vcf")

	script:
	"""
	bgzip ${DeepVarVcf} 
	bcftools index ${DeepVarVcf}.gz
	bcftools index ${genotypedSegments}
	bcftools concat -a ${genotypedSegments} ${DeepVarVcf}.gz -o ${Sample}_deepvar_gCNV_merged.vcf	
	"""
}

process CAVA {
	input:
		tuple val(Sample), file(somaticseqVepAnnovar), file (somaticVcf)
	output:
		tuple val(Sample), file ("*.cava.csv")
	script:
	"""
	${params.cava_path}/cava -c ${params.cava_path}/config_v2.txt -t 10 -i ${somaticVcf} -o ${Sample}.somaticseq
	python3 ${params.cava_script_path} ${Sample}.somaticseq.txt ${Sample}.cava.csv
	"""
}

process MERGE_CSV {
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*.xlsx'
	input:
		tuple val (Sample), file (pindelVep), file (somaticCsv), file (somaticVcf), file (cavaCsv), file (cov_regions)
	output:
		tuple val (Sample), file("${Sample}.xlsx")
	script:
	"""
	merge_csv.py ${Sample} ${Sample}.xlsx  ${cavaCsv} ${cov_regions} ${pindelVep} ${somaticCsv}
	sleep 1s
	"""
}

process MERGE_CSV_WES {
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*.xlsx'
	input:
		tuple val (Sample), file (pindelVep), file (somaticseqVep), file (somaticVcf), file (DeepVarOut), file (DeepVarVcf), file (cavaCsv), file (mosdepth_summary), file (cov_regions), file(median50)
	output:
		val Sample
	script:
	"""
	python3 ${params.pharma_marker_script} $PWD/Final_Output/${Sample}/${Sample}_somaticseq.vep_annonvar.txt ${params.pharma_input_xlxs} Pharma.tsv
	merge_csv_wes.py ${Sample} ${params.output}/${Sample}/${Sample}.xlsx ${cavaCsv} ${mosdepth_summary} ${cov_regions} ${median50} ${pindelVep} ${somaticseqVep} Pharma.tsv ${DeepVarOut}
	sleep 1s
	"""
}
