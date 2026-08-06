process FREEBAYES {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val (Sample), file(bam), file (bamBai)
		path (GenFile)
		path (GenDir)
		path (bedfile)
	output:
		tuple val (Sample), file ("${Sample}.freebayes.vcf")
	script:
	"""
	freebayes -f ${GenFile} -b ${bam} -t ${bedfile} > ${Sample}.freebayes.vcf
	"""
	stub:
	"""
	touch ${Sample}.freebayes.vcf
	"""
}
