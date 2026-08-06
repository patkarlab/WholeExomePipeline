process CNVKIT {
	tag "${Sample}"
	label 'process_medium'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy'
	input:
		tuple val (Sample), file (bam), file (bamBai)
		file (cnvkit_ref)
		file (cnvkit_ref_delIG)
	output:
		tuple val (Sample), file ("results/${Sample}.final.cnr"), file ("results/${Sample}.final.cns"), emit: cnvkit_files
		tuple val (Sample), file ("results/${Sample}.final-scatter.png"), file ("results/${Sample}.final-diagram.pdf"), file("results_delIG/${Sample}.final-scatter_delIG.png"), file("results_delIG/${Sample}.final-diagram_delIG.pdf"), emit: cnvkit_plots
	script:
	"""
	[ -d results ] || mkdir results
	[ -d results_delIG ] || mkdir results_delIG

	cnvkit.py batch ${bam} -r ${cnvkit_ref} -m hybrid -p ${task.cpus} --drop-low-coverage --output-dir ./results --diagram --scatter
	cnvkit.py batch ${bam} -r ${cnvkit_ref_delIG} -m hybrid -p ${task.cpus} --drop-low-coverage --output-dir ./results_delIG --diagram --scatter

	mv results_delIG/${Sample}.final-scatter.png results_delIG/${Sample}.final-scatter_delIG.png
	mv results_delIG/${Sample}.final-diagram.pdf results_delIG/${Sample}.final-diagram_delIG.pdf
	"""
	stub:
	"""
	mkdir -p results results_delIG

	touch \
		results/${Sample}.final.cnr \
		results/${Sample}.final.cns \
		results/${Sample}.final-scatter.png \
		results/${Sample}.final-diagram.pdf \
		results_delIG/${Sample}.final-scatter_delIG.png \
		results_delIG/${Sample}.final-diagram_delIG.pdf
	"""
}
