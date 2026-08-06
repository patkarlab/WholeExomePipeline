process MERGE_CSV_WES {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.xlsx'
	input:
		tuple val (Sample), file (pindelVep), file (somaticseqVep), file (deepVarOut), file (cavaCsv), file (mosdepth_summary), file (cov_regions), file(median50)
		path (pharma_input)
	output:
		val Sample
	script:
	"""
	pharamwes.py ${somaticseqVep} ${pharma_input} Pharma.tsv
	merge_csv_wes.py ${Sample} ${Sample}.xlsx ${cavaCsv} ${mosdepth_summary} ${cov_regions} ${median50} ${pindelVep} ${somaticseqVep} Pharma.tsv ${deepVarOut}
	"""
	stub:
	"""
	touch ${Sample}.xlsx
	"""
}

process MERGE_CSV_LYMPHOMA {
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.xlsx'
	input:
		tuple val (Sample), file (pindelVep), file (somaticCsv), file (cavaCsv), file (cov_regions)
	output:
		tuple val (Sample), file("${Sample}.xlsx")
	script:
	"""
	merge_csv_lymphoma.py ${Sample} ${Sample}.xlsx  ${cavaCsv} ${cov_regions} ${pindelVep} ${somaticCsv}
	"""
	stub:
	"""
	touch ${Sample}.xlsx
	"""
}
