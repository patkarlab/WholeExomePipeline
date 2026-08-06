process SOMATICSEQ_CONCAT { 
	tag "${Sample}"
	label 'process_low'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.vcf'
	input:
		tuple val(Sample), file(somaticseq_snv_vcf), file(somaticseq_indel_vcf)
	output:
		tuple val (Sample), path("${Sample}.somaticseq.vcf")
	script:
	"""
	bgzip -@ ${task.cpus} -c ${somaticseq_snv_vcf} > ${Sample}_somaticseq_snv.vcf.gz
	bcftools index --threads ${task.cpus} -t ${Sample}_somaticseq_snv.vcf.gz

	bgzip -@ ${task.cpus} -c ${somaticseq_indel_vcf} > ${Sample}_somaticseq_indel.vcf.gz
	bcftools index --threads ${task.cpus} -t ${Sample}_somaticseq_indel.vcf.gz
	
	bcftools concat --threads ${task.cpus} -a ${Sample}_somaticseq_snv.vcf.gz ${Sample}_somaticseq_indel.vcf.gz -o ${Sample}.somaticseq.vcf

	sed -i 's/##INFO=<ID=VLK012,Number=6,Type=Integer,Description="Calling decision of the 6 algorithms: VarScan2, LoFreq, Strelka, SnvCaller_0, SnvCaller_1, SnvCaller_2">/##INFO=<ID=VLSFPH,Number=6,Type=String,Description="Calling decision of the 6 algorithms:  VarScan2, LoFreq, Strelka, Freebayes, Platypus, Haplotypecaller">/g' ${Sample}.somaticseq.vcf
	sed -i 's/VLK012/VLSFPH/g' ${Sample}.somaticseq.vcf
	"""
	stub:
	"""
	touch ${Sample}.somaticseq.vcf
	"""
}

process SOMATICSEQ_LYMPHOMA_CONCAT { 
	tag "${Sample}"
	label 'process_low'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.vcf'
	input:
		tuple val(Sample), file(somaticseq_snv_vcf), file(somaticseq_indel_vcf)
	output:
		tuple val (Sample), path("${Sample}.somaticseq.vcf")
	script:
	"""
	bgzip -@ ${task.cpus} -c ${somaticseq_snv_vcf} > ${Sample}_somaticseq_snv.vcf.gz
	bcftools index --threads ${task.cpus} -t ${Sample}_somaticseq_snv.vcf.gz

	bgzip -@ ${task.cpus} -c ${somaticseq_indel_vcf} > ${Sample}_somaticseq_indel.vcf.gz
	bcftools index --threads ${task.cpus} -t ${Sample}_somaticseq_indel.vcf.gz
	
	bcftools concat --threads ${task.cpus} -a ${Sample}_somaticseq_snv.vcf.gz ${Sample}_somaticseq_indel.vcf.gz -o ${Sample}.somaticseq.vcf

	sed -i 's/##INFO=<ID=VLK0123,Number=7,Type=Integer,Description="Calling decision of the 7 algorithms: VarScan2, LoFreq, Strelka, SnvCaller_0, SnvCaller_1, SnvCaller_2, SnvCaller_2">/##INFO=<ID=VLSFPHD,Number=7,Type=String,Description="Calling decision of the 7 algorithms:  VarScan2, LoFreq, Strelka, Freebayes, Platypus, Haplotypecaller, Deepsomatic">/g' ${Sample}.somaticseq.vcf
	sed -i 's/VLK0123/VLSFPHD/g' ${Sample}.somaticseq.vcf
	"""
	stub:
	"""
	touch ${Sample}.somaticseq.vcf
	"""
}