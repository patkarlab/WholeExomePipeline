process PLOT_CNVKIT {
	tag "${Sample}"
	label 'process_low'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*gene_scatter.pdf'
	input:
		tuple val (Sample), file (final_cnr), file (final_cns)
		file (gene_scatter_list)
	output:
		tuple val (Sample), file ("${Sample}_chr_gene_scatter.pdf")
	script:
	"""
	custom_scatter_chrwise.py ${gene_scatter_list} ${final_cnr} ${final_cns} ${Sample}_chr_
	"""
	stub:
	"""
	touch ${Sample}_chr_gene_scatter.pdf
	"""
}
