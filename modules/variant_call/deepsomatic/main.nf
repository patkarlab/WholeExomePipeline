process DEEPSOMATIC {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val (Sample), file(bam), file (bamBai)
		path (control_final_bam)
		path (control_final_bamBai)
		path (GenFile)
		path (GenDir)
		path (bedfile)
	output:
		tuple val(Sample), file("${Sample}_DS.vcf")
	script:
	"""
	run_deepsomatic --model_type=WGS --ref=${GenFile} --reads_tumor=${bam} --reads_normal=${control_final_bam} --output_vcf=${Sample}_DS.vcf --num_shards=${task.cpus} --regions=${bedfile}
	"""
	stub:
	"""
	touch ${Sample}_DS.vcf
	"""
}
