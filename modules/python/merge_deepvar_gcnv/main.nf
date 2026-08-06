
process MERGE_DEEPVAR_GCNV {
    tag "${Sample}"
    label 'process_low'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.vcf'
	input:
		tuple val (Sample), file (DeepVarVcf), file(genotypedIntervals), file (genotypedSegments), file(copyRatios)

	output:
		tuple val (Sample), file("${Sample}_deepvar_gCNV_merged.vcf")

	script:
	"""
	bgzip ${DeepVarVcf} 
	bcftools index --threads ${task.cpus} ${DeepVarVcf}.gz
	bcftools index --threads ${task.cpus} ${genotypedSegments}
	bcftools concat --threads ${task.cpus} -a ${genotypedSegments} ${DeepVarVcf}.gz -o ${Sample}_deepvar_gCNV_merged.vcf
	"""
}