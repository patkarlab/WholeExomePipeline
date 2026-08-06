process PINDEL_FLT3 {
	tag "${Sample}"
	label 'process_inter'
	input:
		tuple val (Sample), file(bam), file (bamBai)
		path (GenFile)
		path (GenDir)
	output:
		tuple val (Sample), file ("${Sample}_pindel_flt3_SI.vcf")
	script:
	"""
	export BAM_2_PINDEL_ADAPT=./opt/pindel-0.2.5b8/Adaptor.pm
	printf '%s\t%s\t%s' ${bam} "300" ${Sample} > ./config.txt
	pindel -f ${GenFile} -i ./config.txt -T ${task.cpus} -c chr13 -o ${Sample}_pindel
	pindel2vcf -r ${GenFile} -P ${Sample}_pindel -R hg19 -d 07102019 -v ${Sample}_pindel_flt3_SI.vcf
 	"""
	stub:
	"""
	touch ${Sample}_pindel_flt3_SI.vcf
	"""
}