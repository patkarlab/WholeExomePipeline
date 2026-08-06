process LUMPY_PREPROCESS {
	tag "${Sample}"
	label "process_medium"
	input:
		tuple val (Sample), file(bqsrBam), file (bqsrBamBai)
	output:
		tuple val (Sample), file("${Sample}.splitters.bam"), file("${Sample}.discordants.bam")
	script:
	"""
	samtools view -@ ${task.cpus} -b -F 1294 ${bqsrBam} >  ${Sample}.discordants.unsorted.bam
	samtools view -@ ${task.cpus} -h ${bqsrBam} | extractSplitReads_BwaMem.py -i stdin | samtools view -@ ${task.cpus} -Sb - > ${Sample}.splitters.unsorted.bam
	samtools sort -@ ${task.cpus} ${Sample}.discordants.unsorted.bam > ${Sample}.discordants.bam
	samtools sort -@ ${task.cpus} ${Sample}.splitters.unsorted.bam > ${Sample}.splitters.bam
	"""
	stub:
	"""
	touch ${Sample}.splitters.bam ${Sample}.discordants.bam
	"""
}