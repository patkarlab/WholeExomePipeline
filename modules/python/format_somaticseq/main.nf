process FORMAT_SOMATICSEQ {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val (Sample), path(Vcf), path(Vep_out), path(Annovar_out)
	output:
		tuple val (Sample), path("${Sample}_somaticseq.vep_annonvar.txt")
	script:
	"""
	extract_somatic.py ${Vcf} ${Sample}.extractedSomaticseq.txt
	extract_velheader.py ${Vep_out} ${Sample}.extractedvepdelheaders.txt
	mergeSomaticVep.py ${Sample}.extractedSomaticseq.txt ${Sample}.extractedvepdelheaders.txt ${Sample}_somaticseq.vep.txt
	sed -i 's/SYMBOL/Gene/g' ${Sample}_somaticseq.vep.txt
	extract_annovar.py ${Annovar_out} ${Sample}_somaticseq.vep.txt ${Sample}_somaticseq.vep_annonvar.txt
	"""
	stub:
	"""
	touch ${Sample}_somaticseq.vep_annonvar.txt
	"""
}
