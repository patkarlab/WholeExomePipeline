process GCNV {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai)
		path (GenFile)
		path (GenDir)
		path (preprocessed_intervals)
		path (ploidy_model)
		path (cohort_model)
	output:
		tuple val (Sample), file("*-genotyped-intervals.vcf.gz"), file("*-genotyped-segments.vcf.gz"), file ("*-denoised-copy-ratios.tsv")
	
	script:
	"""
	export HOME="$PWD"
	mkdir -p "$PWD/.theano"
	export THEANO_FLAGS="base_compiledir=$PWD/.theano"
	
	mv ${GenFile}.dict ${GenFile.simpleName}.dict
	# CollectReadCounts
	gatk --java-options "-Xmx${task.memory.toGiga()}g" CollectReadCounts -L ${preprocessed_intervals} -R ${GenFile} -imr OVERLAPPING_ONLY -I ${finalBam} --format TSV -O ${Sample}.tsv
	
	# DetermineGermlineContigPloidy_case
	gatk --java-options "-Xmx${task.memory.toGiga()}g" DetermineGermlineContigPloidy --model ${ploidy_model} -I ${Sample}.tsv -O . --output-prefix ${Sample}-ploidy --verbosity DEBUG
	
	# GermlineCNVCaller_case
	gatk --java-options "-Xmx${task.memory.toGiga()}g" GermlineCNVCaller --run-mode CASE -I ${Sample}.tsv --contig-ploidy-calls ${Sample}-ploidy-calls --model ${cohort_model} --output ${Sample}_gCNV --output-prefix ${Sample} --verbosity DEBUG
	
	# PostprocessGermlineCNVCalls
	gatk --java-options "-Xmx${task.memory.toGiga()}g" PostprocessGermlineCNVCalls --model-shard-path ${cohort_model} --calls-shard-path ${Sample}_gCNV/${Sample}-calls --allosomal-contig chrX --allosomal-contig chrY --contig-ploidy-calls ${Sample}-ploidy-calls --sample-index 0 --output-genotyped-intervals ${Sample}-genotyped-intervals.vcf.gz --output-genotyped-segments ${Sample}-genotyped-segments.vcf.gz --sequence-dictionary ${GenFile.simpleName}.dict --output-denoised-copy-ratios ${Sample}-denoised-copy-ratios.tsv
	"""
}