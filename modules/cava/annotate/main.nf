process CAVA {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val(Sample), file(somaticseqVcf)
		path(cava_config)
		path(GenFile)
		path(GenDir)
		path(bedfile)
		path(bedfile_zipped)
		path(SNPs)
		path(SNPs_index)

	output:
		tuple val(Sample), file("${Sample}.somaticseq.txt")

	script:
	"""
	cp ${cava_config} cava.config
	sed -i "s|^@target *=.*|@target = ${bedfile}.gz|g" cava.config

	cava.py -c cava.config -t ${task.cpus} -i ${somaticseqVcf} -o ${Sample}.somaticseq
	"""
	stub:
	"""
	touch ${Sample}.somaticseq.txt
	"""
}