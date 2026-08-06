process HSMETRICS {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val(Sample), path(bam), path(bai)
		path (bedfile)
		path (GenFile)
		path (GenDir)
	output:
		path("${Sample}_hsmetrics.txt")
	script:
	"""
	gatk BedToIntervalList I=${bedfile} O=${bedfile}_sortd.interval_list SD=${GenFile}.dict
	gatk CollectHsMetrics I=${bam} O=${Sample}_hsmetrics.txt BAIT_INTERVALS=${bedfile}_sortd.interval_list TARGET_INTERVALS=${bedfile}_sortd.interval_list R= ${GenFile} VALIDATION_STRINGENCY=LENIENT
	"""
	stub:
	"""
	touch ${Sample}_hsmetrics.txt
	"""
}

process HSMETRICS_COLLECT {
	label 'process_low'
	publishDir "${params.outdir}/", mode: 'copy'
	input:
		file (hsmetrics)
	output:
		file("hsmetrics.txt")
	script:
	"""
	echo -e "Sample name\tOn target\tOff target" > hsmetrics.txt
	for i in ${hsmetrics}
	do
		samp_name=\$(basename -s _hsmetrics.txt \${i})
		grep -v '#' \${i} | awk -v name=\${samp_name} 'BEGIN{FS="\t"; OFS="\t"}NR==3{ print name,\$7,\$8}' >> hsmetrics.txt
	done
	"""
	stub:
	"""
	touch hsmetrics.txt
	"""
}
