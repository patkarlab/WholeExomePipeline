process DELLY {
	tag "${Sample}"
	label "process_low"
	publishDir "${params.outdir}/${Sample}/translocatn/delly", mode: 'copy', pattern: '*vcf'
	input:
		tuple val (Sample), file(bqsrBam), file (bqsrBamBai)
		path (GenFile)
		path (GenDir)
		path (control_bqsr_bam)
		path (control_bqsr_bamBai)
	output:
		tuple val (Sample), file ("${Sample}_delly.vcf"), file ("${Sample}_delly_somatic.vcf")
	script:
	"""
	tumor_name=`basename ${bqsrBam} .old_final.bam`
	control_name=`basename ${control_bqsr_bam} .old_final.bam`

	echo -e "\${tumor_name}\ttumor" > samples.tsv
	echo -e "\${control_name}\tcontrol" >> samples.tsv

	delly call \
		-g ${GenFile} \
		${control_bqsr_bam} ${bqsrBam} \
		-o ${Sample}_delly.bcf
	
	bcftools view ${Sample}_delly.bcf -Ov -o ${Sample}_delly.vcf

	
	delly filter \
		-f somatic \
		-s samples.tsv \
		${Sample}_delly.bcf \
	| bcftools view -Ov -o ${Sample}_delly_somatic.vcf
	"""
	stub:
	"""
	touch ${Sample}_delly.vcf ${Sample}_delly_somatic.vcf
	"""
}
