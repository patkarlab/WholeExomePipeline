
process COVERAGE_BEDTOOLS {
	tag "${Sample}"
	label 'process_high'
	input:
		tuple val (Sample), file(bam), file(bamBai)
		path(bedfile)
	output:
		tuple val (Sample), file ("${Sample}_cov.regions_bedtools.bed")
	script:
	"""
	bedtools coverage -counts -a ${bedfile} -b ${bam} > ${Sample}_cov.regions_bedtools.bed
	"""
	stub:
	"""
	touch ${Sample}_cov.regions_bedtools.bed
	"""
}
