process FORMAT_PINDEL {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val (Sample), path(Vcf), path(Vep_out)
	output:
		tuple val (Sample), path("${Sample}_pindel.vep.txt")
	script:
	"""
	extract_pindelSI.py ${Vcf} ${Sample}_extractedPindelSI.txt
	extract_pindel.py  ${Vep_out} ${Sample}_extractedPindelVep.txt
	mergepindel.py ${Sample}_extractedPindelSI.txt ${Sample}_extractedPindelVep.txt ${Sample}_pindel.vep.txt
	"""
	stub:
	"""
	touch ${Sample}_pindel.vep.txt
	"""
}
