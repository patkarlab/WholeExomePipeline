process CNVKIT {
	tag "${Sample}"
	label 'process_medium'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy'
	input:
		tuple val (Sample), file (bam), file (bamBai)
		file (cnvkit_ref)
	output:
		tuple val (Sample), file ("results/${Sample}.final.cnr"), file ("results/${Sample}.final.cns"), emit: cnvkit_files
		tuple val (Sample), file ("results/${Sample}.final-scatter.png"), file ("results/${Sample}.final-diagram.pdf"), emit: cnvkit_plots
	script:
	"""
	[ -d results ] || mkdir results

	cnvkit.py batch ${bam} -r ${cnvkit_ref} -m hybrid -p ${task.cpus} --drop-low-coverage --output-dir ./results --diagram --scatter

	"""
	stub:
	"""
	mkdir -p results

	touch \
		results/${Sample}.final.cnr \
		results/${Sample}.final.cns \
		results/${Sample}.final-scatter.png \
		results/${Sample}.final-diagram.pdf \
	"""
}
