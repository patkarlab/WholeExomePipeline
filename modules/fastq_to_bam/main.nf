#!/usr/bin/nextflow
process TRIM {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val(Sample), file(read1), file(read2)
		path (adapters)
	output:
		tuple val(Sample), file("${Sample}_trim_R1.fastq"), file("${Sample}_trim_R2.fastq")
	script:
	"""
	fastp -i ${read1} -I ${read2} -o ${Sample}_trim_R1.fastq -O ${Sample}_trim_R2.fastq --adapter_fasta ${adapters} -w ${task.cpus}
	"""
	stub:
	"""
	touch ${Sample}_trim_R1.fastq ${Sample}_trim_R2.fastq
	"""
}

process MAPBAM {
	tag "${Sample}"
	label 'process_high'
	input:
		tuple val(Sample), file(trim1), file(trim2)
		path (GenFile)
		path (GenDir)
	output:
		tuple val(Sample), file ("${Sample}.bam")
	script:
	"""
	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:200" \
	-M -t ${task.cpus} ${GenFile} ${trim1} ${trim2} | samtools sort -@ ${task.cpus} -o ${Sample}.bam -
	"""
	stub:
	"""
	touch ${Sample}.bam
	"""
}

process MARK_DUPS {
	tag "${Sample}"
	label 'process_high'
	input:
		tuple val(Sample), file(sortd_bam)
	output:
		tuple val(Sample), file("${Sample}_markdups.bam"), file("${Sample}_marked_dup_metrics.txt")
	script:
	"""
	picard -Xmx${task.memory.toGiga()}g \
		MarkDuplicates \
		I=${sortd_bam} \
		O=${Sample}_markdups.bam \
		M=${Sample}_marked_dup_metrics.txt \
		TMP_DIR=${Sample}_tmp
	"""
	stub:
	"""
	touch ${Sample}_markdups.bam ${Sample}_marked_dup_metrics.txt
	"""
}

process BQSR {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val(Sample), file(markdups_bam), file(markdups_metrics)
		path (GenFile)
		path (GenDir)
		path (SNPS)
		path (SNPS_index)
		path (INDELS)
		path (INDELS_index)
	output:
		tuple val(Sample), file("${Sample}_recal.table")
	script:
	"""
	mv ${GenFile}.dict ${GenFile.simpleName}.dict
	gatk --java-options "-Xmx${task.memory.toGiga()}g" \
		BaseRecalibrator \
		-I ${markdups_bam} \
		-R ${GenFile} \
		--known-sites ${SNPS} \
		--known-sites ${INDELS} \
		--bqsr-baq-gap-open-penalty 30.0 \
		-O ${Sample}_recal.table
	"""
	stub:
	"""
	touch ${Sample}_recal.table
	"""
}


process APPLY_BQSR {
	label 'process_medium'
	tag "${Sample}"
	input:
		tuple val(Sample), file(markdups_bam), file(markdups_metrics), file(recal_table)
		path (GenFile)
		path (GenDir)
	output:
		tuple val(Sample), file("${Sample}_bqsr.bam"), file("${Sample}_bqsr.bam.bai")
	script:
	"""
	mv ${GenFile}.dict ${GenFile.simpleName}.dict
	gatk --java-options "-Xmx${task.memory.toGiga()}g" \
		ApplyBQSR \
		-R ${GenFile} \
		-I ${markdups_bam} \
		--bqsr-recal-file ${recal_table} \
		-O ${Sample}_bqsr.bam

	mv ${Sample}_bqsr.bai ${Sample}_bqsr.bam.bai
	"""
	stub:
	"""
	touch ${Sample}_bqsr.bam ${Sample}_bqsr.bam.bai
	"""
}


process ALIGNMENT_METRICS {
	label 'process_low'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_alignment_summary_metrics.txt'
	tag "${Sample}"
	input:
		tuple val(Sample), file(bqsr_bam), file(bqsr_bam_bai)
		path (GenFile)
		path (GenDir)
	output:
		tuple val(Sample), file("${Sample}_alignment_summary_metrics.txt")
	script:
	"""
	picard -Xmx${task.memory.toGiga()}g \
		CollectAlignmentSummaryMetrics \
		R=${GenFile} \
		I=${bqsr_bam} \
		O=${Sample}_alignment_summary_metrics.txt
	"""
	stub:
	"""
	touch ${Sample}_alignment_summary_metrics.txt
	"""
}

process INSERT_SIZE_METRICS {
	label 'process_low'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_insert_size_metrics.txt'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_insert_size_metrics.pdf'
	tag "${Sample}"
	input:
		tuple val(Sample), file(bqsr_bam), file(bqsr_bam_bai)
	output:
		tuple val(Sample), file("${Sample}_insert_size_metrics.txt"), file("${Sample}_insert_size_metrics.pdf")
	script:
	"""
	picard -Xmx${task.memory.toGiga()}g \
		CollectInsertSizeMetrics \
		I=${bqsr_bam} \
		O=${Sample}_insert_size_metrics.txt \
		H=${Sample}_insert_size_metrics.pdf \
		HISTOGRAM_WIDTH=500 \
		TMP_DIR=${Sample}_tmp
	"""
	stub:
	"""
	touch ${Sample}_insert_size_metrics.txt ${Sample}_insert_size_metrics.pdf
	"""
}
