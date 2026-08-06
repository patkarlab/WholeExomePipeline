
process COVERAGE_MOSDEPTH {
	label 'process_medium'
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*'
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai)
		path (bedfile)
	output:
		tuple val (Sample), file ("${Sample}_cov.mosdepth.summary.txt"), file("${Sample}_cov.regions.bed")
	script:
	"""
	mosdepth -n -t ${task.cpus} --fast-mode -m --by ${bedfile} ${Sample}_cov ${finalBam}
	gunzip ${Sample}_cov.regions.bed.gz
	"""
}