process POSTPROCESS_GCNV {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val(Sample), path(gcnv_dir), path(ploidy_calls), path(dict)
		path (cohort_model)
	output:
		tuple val(Sample), path("${Sample}-genotyped-intervals.vcf.gz"), path("${Sample}-genotyped-segments.vcf.gz"), path("${Sample}-denoised-copy-ratios.tsv")
	script:
	"""
	export HOME="\$PWD"
	mkdir -p "\$PWD/.theano"
	export THEANO_FLAGS="base_compiledir=\$PWD/.theano"

	gatk --java-options "-Xmx${task.memory.toGiga()}g" \
		PostprocessGermlineCNVCalls \
		--model-shard-path ${cohort_model} \
		--calls-shard-path ${Sample}_gCNV/${Sample}-calls \
		--allosomal-contig chrX \
		--allosomal-contig chrY \
		--contig-ploidy-calls ${ploidy_calls} \
		--sample-index 0 \
		--output-genotyped-intervals ${Sample}-genotyped-intervals.vcf.gz \
		--output-genotyped-segments ${Sample}-genotyped-segments.vcf.gz \
		--sequence-dictionary ${dict} \
		--output-denoised-copy-ratios ${Sample}-denoised-copy-ratios.tsv
	"""
}