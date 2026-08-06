process MANTA {
	tag "${Sample}"
	label "process_low"
	errorStrategy 'ignore'
	publishDir "${params.outdir}/${Sample}/translocatn/manta", mode: 'copy', pattern: '*vcf'
	publishDir "${params.outdir}/${Sample}/translocatn/manta", mode: 'copy', pattern: '*tsv'
	input:
		tuple val (Sample), file(bqsrBam), file (bqsrBamBai)
		path (GenFile)
		path (GenDir)
		path (control_bqsr_bam)
		path (control_bqsr_bamBai)
	output:
		tuple val (Sample), file ("${Sample}_manta_candidateSV.vcf"), file ("${Sample}_manta_diploidSV.vcf"), file ("${Sample}_manta_somaticSV.vcf")
	script:
	"""
	/opt/manta/bin/configManta.py \
		--normalBam ${control_bqsr_bam} \
		--tumorBam ${bqsrBam} \
		--referenceFasta ${GenFile} --runDir ./ \
		--exome

	./runWorkflow.py -j ${task.cpus}

	gunzip -c ./results/variants/candidateSV.vcf.gz > ${Sample}_manta_candidateSV.vcf
	gunzip -c ./results/variants/diploidSV.vcf.gz > ${Sample}_manta_diploidSV.vcf
	gunzip -c ./results/variants/somaticSV.vcf.gz > ${Sample}_manta_somaticSV.vcf
	"""
	stub:
	"""
	touch ${Sample}_manta_candidateSV.vcf ${Sample}_manta_diploidSV.vcf ${Sample}_manta_somaticSV.vcf
	"""
}