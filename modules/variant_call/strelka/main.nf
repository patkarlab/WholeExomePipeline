process STRELKA {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val (Sample), file(bam), file (bamBai)
		path (GenFile)
		path (GenDir)
		path (bedfile)
		path (bedfile_zipped)
	output:
		tuple val (Sample), file ("${Sample}.strelka.vcf")
	script:
	"""
	/opt/strelka/bin/configureStrelkaGermlineWorkflow.py --bam ${bam} --referenceFasta ${GenFile} --callRegions  ${bedfile}.gz --targeted --exome --runDir ./
	./runWorkflow.py -m local -j ${task.cpus}
	gunzip -f ./results/variants/variants.vcf.gz
	mv ./results/variants/variants.vcf ./${Sample}.strelka.vcf
	"""
	stub:
	"""
	touch ${Sample}.strelka.vcf
	"""
}
