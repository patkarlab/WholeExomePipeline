process LOFREQ {
	tag "${Sample}"
    label 'process_medium'
	input:
		tuple val (Sample), file (bam), file (bamBai)
		path (GenFile)
		path (GenDir)
		file (bedfile)
	output:
		tuple val (Sample), file ("${Sample}.lofreq.filtered.vcf")
	script:
	"""
	lofreq viterbi -f ${GenFile} -o ${Sample}.lofreq.pre.bam ${bam}
	samtools sort -@ ${task.cpus} ${Sample}.lofreq.pre.bam > ${Sample}.lofreq.bam
	lofreq call -b dynamic -C 50 -a 0.00005 -q 30 -Q 30 -m 50 -f ${GenFile} -l ${bedfile} -o ${Sample}.lofreq.vcf ${Sample}.lofreq.bam
	lofreq filter -a 0.005 -i ${Sample}.lofreq.vcf -o ${Sample}.lofreq.filtered.vcf
	"""
	stub:
	"""
	touch ${Sample}.lofreq.filtered.vcf
	"""
}