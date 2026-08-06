
process EXTRACT_COV50 {
	label 'process_medium'
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*'
	input:
		tuple val (Sample), file (mosdepth_summary), file (cov_regions)
	output:
		tuple val (Sample), file ("${Sample}_median50.bed")
	script:
	"""
    extract_COV50.py ${cov_regions} ${Sample}_median50.bed
	"""
}