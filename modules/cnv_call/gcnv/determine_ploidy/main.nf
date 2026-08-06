process DETERMINE_PLOIDY {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val(Sample), path(counts), path(dict)
		path (ploidy_model)
	output:
		tuple val(Sample), path("${Sample}-ploidy-calls"), path(counts), path(dict)

	script:
	"""
	export HOME="\$PWD"
	mkdir -p "\$PWD/.theano"
	export THEANO_FLAGS="base_compiledir=\$PWD/.theano"

	gatk --java-options "-Xmx${task.memory.toGiga()}g" \
		DetermineGermlineContigPloidy \
		--model ${ploidy_model} \
		-I ${counts} \
		-O . \
		--output-prefix ${Sample}-ploidy \
		--verbosity DEBUG
	"""
}