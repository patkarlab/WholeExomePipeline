process GRIDSS {
	tag "${Sample}"
	label "process_high"
	publishDir "${params.outdir}/${Sample}/translocatn/gridss", mode: 'copy', pattern: '*vcf'
	input:
		tuple val (Sample), file(bqsrBam), file (bqsrBamBai)
		path (GenFile)
		path (GenDir)
		path (control_bqsr_bam)
		path (control_bqsr_bamBai)
		path (gridss_exclude_list)
	output:
		tuple val (Sample), file ("${Sample}.gridss.vcf"), file ("${Sample}.gridss.somatic.vcf.bgz"), file ("RM_${Sample}.gridss.somatic.vcf")
	script:
	"""
	gridss -r ${GenFile} \
	-j /opt/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar \
	-t ${task.cpus} \
	-o ${Sample}.gridss.vcf \
	-b ${gridss_exclude_list} \
	${bqsrBam} ${control_bqsr_bam}

	cp /opt/gridss/gridss_somatic_filter .
	sed -i 's/\r\$//' gridss_somatic_filter
	chmod +x gridss_somatic_filter

	./gridss_somatic_filter \
	--scriptdir /opt/gridss \
	-i ${Sample}.gridss.vcf \
	-o ${Sample}.gridss.somatic.vcf \
	-n 2 \
	-t 1 \

	gridss_annotate_vcf_repeatmasker \
	-j /opt/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar \
	-t ${task.cpus} \
	-o RM_${Sample}.gridss.somatic.vcf \
	${Sample}.gridss.somatic.vcf.bgz
	"""
	stub:
	"""
	touch ${Sample}.gridss.vcf ${Sample}.gridss.somatic.vcf.bgz RM_${Sample}.gridss.somatic.vcf
	"""
}