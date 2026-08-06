process FORMAT_DEEPVARIANT {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val (Sample), path(Vcf), path(Vep_out), path(Annovar_out)
	output:
		tuple val (Sample), path("${Sample}_deepvar_vep_annonvar.txt")
	script:
	"""
	extract_vaf.py ${Vcf} ${Sample}.extracted.csv
	extract_vepdata.py ${Vep_out} ${Sample}.extractedvepdelheaders.csv
	mergeDeepVariantVep.py ${Sample}.extracted.csv ${Sample}.extractedvepdelheaders.csv ${Sample}_deepvariant.txt
	deepvar_annovar.py ${Annovar_out} ${Sample}_deepvariant.txt ${Sample}_deepvar_vep_annonvar.txt
	"""
}
