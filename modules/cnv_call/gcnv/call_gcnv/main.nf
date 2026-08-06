process CALL_GCNV {
	tag "${Sample}"
	label 'process_high'
	input:
		tuple val(Sample), path(ploidy_calls), path(counts), path(dict)
		path (cohort_model)
	output:
		tuple val(Sample), path("${Sample}_gCNV"), path(ploidy_calls), path(dict)
	script:
	"""
	export HOME="\$PWD"
	mkdir -p "\$PWD/.theano"
	export THEANO_FLAGS="base_compiledir=\$PWD/.theano"

	export OMP_NUM_THREADS=${task.cpus}
	export MKL_NUM_THREADS=${task.cpus}

	gatk --java-options "-Xmx${task.memory.toGiga()}g" \
		GermlineCNVCaller \
		--run-mode CASE \
		-I ${counts} \
		--contig-ploidy-calls ${ploidy_calls} \
		--model ${cohort_model} \
		--output ${Sample}_gCNV \
		--output-prefix ${Sample} \
		--verbosity DEBUG
	"""
}