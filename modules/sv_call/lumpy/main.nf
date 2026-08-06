process LUMPY {
	tag "${Sample}"
	label "process_low"
	input:
		tuple val (Sample), file(bqsrBam), file (bqsrBamBai), file (splitters_bam), file(discordants_bam)
	output:
		tuple val (Sample), file("${Sample}.vcf")
	script:
	"""
	lumpyexpress -B ${bqsrBam} -S ${splitters_bam} -D ${discordants_bam} -o ${Sample}.vcf
	"""
	stub:
	"""
	touch ${Sample}.vcf
	"""
}