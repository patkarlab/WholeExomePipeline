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

process CNVKIT {
	tag "${Sample}"
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*gene_scatter.pdf'
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*final-scatter.png'
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*final-diagram.pdf'
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*final-scatter_delIG.png'
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*final-diagram_delIG.pdf'
	input:
		tuple val (Sample), file(finalBam), file(finalBamBai)
		file (gene_scatter)
		file (cnvkit_reference)
		file (cnvkit_reference_delIG)
	output:
		tuple val (Sample), file ("*gene_scatter.pdf"), file ("*final-scatter.png"), file ("*final-diagram.pdf"), file ("*final-scatter_delIG.png"), file ("*final-diagram_delIG.pdf")
	script:
	"""
	[ -d results ] || mkdir results
	[ -d results_delIG ] || mkdir results_delIG

	cnvkit.sh ${finalBam} ${cnvkit_reference} results/
	cp results/${Sample}_final-diagram.pdf ${Sample}_final-diagram.pdf
	cp results/${Sample}_final-scatter.png ${Sample}_final-scatter.png

	cnvkit.sh ${finalBam} ${cnvkit_reference_delIG} results_delIG/
	cp results_delIG/${Sample}_final-diagram.pdf ${Sample}_final-diagram_delIG.pdf
	cp results_delIG/${Sample}_final-scatter.png ${Sample}_final-scatter_delIG.png

	custom_scatter_chrwise.py ${gene_scatter} results/${Sample}_final.cnr results/${Sample}_final.cns ${Sample}_chr_
	"""
}
