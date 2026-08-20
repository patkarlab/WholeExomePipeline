process FORMAT_ANNOVAR {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val (Sample), file (somaticseq_multianno)
	output:
		tuple val (Sample), file ("${Sample}.somaticseq.csv")
	script:
	"""
	format_somaticseq_annovar.py ${somaticseq_multianno} ${Sample}.somaticseq.csv
	"""
	stub:
	"""
	touch ${Sample}.somaticseq.csv
	"""
}
