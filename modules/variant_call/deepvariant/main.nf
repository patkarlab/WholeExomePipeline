process DEEPVARIANT {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val (Sample), file (bam), file (bamBai)
		path (GenFile)
		path (GenDir)
        file (bedfile)
	output:
		tuple val(Sample), file ("${Sample}_deepvar_filt.vcf")
	script:
	"""
	/opt/deepvariant/bin/run_deepvariant --model_type=WES --ref=${GenFile} --regions=${bedfile} --reads=${bam} --output_vcf=${Sample}.deepvar.vcf --num_shards=${task.cpus}
	bcftools view -f PASS ${Sample}.deepvar.vcf > ${Sample}_deepvar_filt.vcf
	"""
}
