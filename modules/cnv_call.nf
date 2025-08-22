#!/usr/bin/env nextflow

process IFCNV {
	conda '/home/miniconda3/envs/new_base/'
	input:
		val Sample
	output:
		//tuple val (Sample), file ("*.html"), file ("*.tsv")
		val (Sample)
	script:
	"""
	${params.links} $PWD/Final_Output/ ${params.input}
	mkdir ifCNV
	${params.ifcnv} ./ ${params.bedfile}.bed ifCNV

	# Making ifCNV's output directory for each sample
	for i in `cat ${params.input}`
	do
		if [ ! -d $PWD/Final_Output/\${i}/ifCNV ]; then
			mkdir $PWD/Final_Output/\${i}/ifCNV
		fi		
	done

	# Copying output of ifCNV to respective samples
	if [ -f ifCNV/ifCNV.tsv ]; then
		for i in `awk 'NR>1{print \$3}' ifCNV/ifCNV.tsv | awk 'BEGIN{FS="."}{print \$1}' | sort | uniq`	
		do
			cp ifCNV/\${i}.*.html $PWD/Final_Output/\${i}/ifCNV/	
		done
	fi
	"""
}

process GCNV {
	conda '/home/diagnostics/.conda/envs/gatk'
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	
	output:
		tuple val (Sample), file("*-genotyped-intervals.vcf.gz"), file("*-genotyped-segments.vcf.gz"), file ("*-denoised-copy-ratios.tsv")
	
	script:
	"""
	# CollectReadCounts
	${params.gatk_path} CollectReadCounts -L ${params.preprocessed_intervals} -R ${params.genome} -imr OVERLAPPING_ONLY -I ${finalBam} --format TSV -O ${Sample}.tsv
	
	# DetermineGermlineContigPloidy_case
	${params.gatk_path} DetermineGermlineContigPloidy --model ${params.ploidy_model} -I ${Sample}.tsv -O . --output-prefix ${Sample}-ploidy --verbosity DEBUG
	
	# GermlineCNVCaller_case
	${params.gatk_path} GermlineCNVCaller --run-mode CASE -I ${Sample}.tsv --contig-ploidy-calls ${Sample}-ploidy-calls --model ${params.cohort_model} --output ${Sample}_gCNV --output-prefix ${Sample} --verbosity DEBUG
	
	# PostprocessGermlineCNVCalls
	${params.gatk_path} PostprocessGermlineCNVCalls --model-shard-path ${params.cohort_model} --calls-shard-path ${Sample}_gCNV/${Sample}-calls --allosomal-contig chrX --allosomal-contig chrY --contig-ploidy-calls ${Sample}-ploidy-calls --sample-index 0 --output-genotyped-intervals ${Sample}-genotyped-intervals.vcf.gz --output-genotyped-segments ${Sample}-genotyped-segments.vcf.gz --sequence-dictionary /home/reference_genomes/hg19_broad/hg19_all.dict --output-denoised-copy-ratios ${Sample}-denoised-copy-ratios.tsv
	"""
}