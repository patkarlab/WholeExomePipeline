process ABRA_SORT {
	label 'process_medium'
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.final.bam'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.final.bam.bai'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.old_final.bam'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.old_final.bam.bai'
	input:
		tuple val(Sample), file (bam), file(bamBai), file(abra_bam)
	output:
		tuple val(Sample), file ("*.final.bam"), file ("*.final.bam.bai"), emit: final_bam
		tuple val(Sample), file ("*.old_final.bam"), file ("*.old_final.bam.bai"), emit: old_bam
	script:
	"""
	samtools sort -@ ${task.cpus} ${bam} > ${Sample}.old_final.bam
	samtools index ${Sample}.old_final.bam > ${Sample}.old_final.bam.bai
	samtools sort -@ ${task.cpus} ${abra_bam} > ${Sample}.final.bam
	samtools index ${Sample}.final.bam > ${Sample}.final.bam.bai
	"""
	stub:
	"""
	touch ${Sample}.old_final.bam ${Sample}.old_final.bam.bai ${Sample}.final.bam ${Sample}.final.bam.bai
	"""
}
