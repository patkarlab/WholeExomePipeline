process FORMAT_COVERVIEW {
	tag "${Sample}"
	label 'process_low'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: "*.coverview_regions.csv"
	input:
		tuple val (Sample), file(coverview_regions_txt)
	output:
		tuple val (Sample), file ("${Sample}.coverview_regions.csv")
	script:
	"""
	format_coverview.py ${coverview_regions_txt} ${Sample}.coverview_regions.csv
	"""
	stub:
	"""
	touch ${Sample}.coverview_regions.csv
	"""
}
