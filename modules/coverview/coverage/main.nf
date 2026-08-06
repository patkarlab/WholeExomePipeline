process COVERAGE_COVERVIEW {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val (Sample), file(bam), file(bamBai)
		path (bedfile)
		path (coverview_config)
	output:
		tuple val (Sample), file ("${Sample}.coverview_regions.txt")
	script:
	"""
	coverview -i ${bam} -b ${bedfile} -c ${coverview_config} -o ${Sample}.coverview
	"""
	stub:
	"""
	touch ${Sample}.coverview_regions.txt
	"""
}
