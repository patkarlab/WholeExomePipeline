process CAVA {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val(Sample), file (somaticseqVcf)
		path (cava_config)
		path (GenFile)
		path (GenDir)
		path (bedfile)
		path (SNPs)
		path (SNPs_index)
		path (ensembl_db)
		path (ensembl_db_index)        
	output:
		tuple val(Sample), file ("${Sample}.somaticseq.txt")
	script:
	"""
	cava.py -c ${cava_config} -t ${task.cpus} -i ${somaticseqVcf} -o ${Sample}.somaticseq
	"""
	stub:
	"""
	touch ${Sample}.somaticseq.txt
	"""
}
