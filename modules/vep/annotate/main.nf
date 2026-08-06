process VEP {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val (Sample), path(Vcf)
	output:
		tuple val (Sample), file("${Sample}_vep_delheaders.txt")
	script:
	"""
	vep -i ${Vcf} --cache -o ${Sample}_vep.txt --offline --tab --force_overwrite --symbol --protein --af --max_af  --no_check_alleles --sift b --variant_class --canonical --allele_number --hgvs --shift_hgvs 1 --af_1kg --af_gnomadg --pubmed
	
	grep -v "##" ${Sample}_vep.txt > ${Sample}_vep_delheaders.txt
	"""
	stub:
	"""
	touch ${Sample}_vep_delheaders.txt
	"""
}

process VEP_ANALYSIS {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val (Sample), path(Vcf)
	output:
		tuple val (Sample), file("${Sample}_vep_delheaders.txt")
	script:
	"""
	awk 'length(\$0) <= 100000 || /^#/' ${Vcf} > pindel_filtered.vcf
	vep -i pindel_filtered.vcf --cache -o ${Sample}_vep.txt --offline --tab --force_overwrite --symbol --protein --af --max_af  --no_check_alleles --sift b --variant_class --canonical --allele_number --hgvs --shift_hgvs 1 --af_1kg --af_gnomadg --pubmed
	filter_vep -i ${Sample}_vep.txt -o ${Sample}_filtered.txt --filter "(CANONICAL is YES) and (AF < 0.01 or not AF)" --force_overwrite
	grep -v "##" ${Sample}_filtered.txt > ${Sample}_vep_delheaders.txt
	"""
	stub:
	"""
	touch ${Sample}_vep_delheaders.txt
	"""
}