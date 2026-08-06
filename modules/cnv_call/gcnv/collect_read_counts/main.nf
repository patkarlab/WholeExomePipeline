process COLLECT_READ_COUNTS {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val(Sample), file(finalBam), file(finalBamBai)
		path (GenFile)
		path (GenDir)
		path (preprocessed_intervals)
	output:
		tuple val(Sample), path("${Sample}.tsv"), path("${GenFile.simpleName}.dict")
	script:
	"""
	export HOME="\$PWD"
	mkdir -p "\$PWD/.theano"
	export THEANO_FLAGS="base_compiledir=\$PWD/.theano"

	mv ${GenFile}.dict ${GenFile.simpleName}.dict

	gatk --java-options "-Xmx${task.memory.toGiga()}g" \
		CollectReadCounts \
		-L ${preprocessed_intervals} \
		-R ${GenFile} \
		-imr OVERLAPPING_ONLY \
		-I ${finalBam} \
		--format TSV \
		-O ${Sample}.tsv
	"""
}