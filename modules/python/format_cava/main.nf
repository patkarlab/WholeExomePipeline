process FORMAT_CAVA {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val (Sample), path(cava_out)
	output:
		tuple val (Sample), path("${Sample}.cava.csv")
	script:
	"""
	cava_format.py ${cava_out} ${Sample}.cava.csv
	"""
	stub:
	"""
	touch ${Sample}.cava.csv
	"""
}
