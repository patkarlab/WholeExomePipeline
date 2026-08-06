process SVABA {
	tag "${Sample}"
	label "process_medium"
	publishDir "${params.outdir}/${Sample}/svaba/", mode: 'copy', pattern: '*svaba*'
	input:
		tuple val (Sample), file(bqsrBam), file (bqsrBamBai)
		path (GenFile)
		path (GenDir)
		path (known_SNPs)
		path (known_SNPs_index)
		path (control_bqsr_bam)
		path (control_bqsr_bamBai)
	output:
		tuple val (Sample), file ("${Sample}.svaba.somatic.sv.vcf"), file ("${Sample}.svaba.germline.sv.vcf"), file ("${Sample}.svaba.somatic.indel.vcf"), file ("${Sample}.svaba.germline.indel.vcf")
	script:
	"""
	svaba run -t ${bqsrBam} -n ${control_bqsr_bam} -G ${GenFile} -p ${task.cpus} -D ${known_SNPs} -a ${Sample}
	"""
	stub:
	"""
	touch ${Sample}.svaba.somatic.sv.vcf ${Sample}.svaba.germline.sv.vcf ${Sample}.svaba.somatic.indel.vcf ${Sample}.svaba.germline.indel.vcf
	"""
}