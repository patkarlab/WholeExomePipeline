process ABRA_BAM {
	label 'process_medium'
	tag "${Sample}"
	input:
		tuple val(Sample), file (bam), file(bamBai)
		path (bedfile)
		path (GenFile)
		path (GenDir)
	output:
		tuple val(Sample), file ("${Sample}.abra.bam")
	script:
	"""
	java -Xmx${task.memory.toGiga()}g -jar /opt/biotools/abra.jar --in ${bam} --out ${Sample}.abra.bam --ref ${GenFile} --threads ${task.cpus} --targets ${bedfile} --working temp
	"""
	stub:
	"""
	touch ${Sample}.abra.bam
	"""
}
