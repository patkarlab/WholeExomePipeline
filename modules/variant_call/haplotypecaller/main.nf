process HAPLOTYPECALLER {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val (Sample), file (bam), file (bamBai)
		path (GenFile)
		path (GenDir)
		path (bedfile)
		path (known_SNPs)
		path (known_SNPs_index)
	output:
		tuple val (Sample), file ("${Sample}.haplotypecaller.vcf")
	script:
	"""
	mv ${GenFile}.dict ${GenFile.simpleName}.dict
	java -Xmx${task.memory.toGiga()}g -jar /usr/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${GenFile} -I ${bam} -L ${bedfile} -o ${Sample}.haplotypecaller.vcf --dbsnp ${known_SNPs} 
	"""
	stub:
	"""
	touch ${Sample}.haplotypecaller.vcf
	"""
}