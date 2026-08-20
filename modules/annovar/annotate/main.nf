process ANNOVAR {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val (Sample), path(Vcf)
		val(variant_caller)
	output:
		tuple val (Sample), path("${Sample}_${variant_caller}.annovar.hg19_multianno.txt")
	script:
	"""
	convert2annovar.pl -format vcf4 ${Vcf} --outfile ${Sample}_${variant_caller}.avinput -withfreq --includeinfo -allsample
	table_annovar.pl ${Sample}_${variant_caller}.avinput --out ${Sample}_${variant_caller}.annovar --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' \
	--thread ${task.cpus} /databases/humandb -xreffile /databases/gene_fullxref.txt

	touch ${Sample}_${variant_caller}.annovar.hg19_multianno.txt
	"""
	stub:
	"""
	touch ${Sample}_${variant_caller}.annovar.hg19_multianno.txt
	"""
}

process ANNOVAR_LYMPHOMA {
	tag "${Sample}"
	label 'process_inter'
	input:
		tuple val (Sample), path(Vcf)
		val(variant_caller)
	output:
		tuple val (Sample), path("${Sample}_${variant_caller}.out.hg38_multianno.csv")
	script:
	"""
	convert2annovar.pl -format vcf4 ${Vcf} --outfile ${Sample}_${variant_caller}.avinput --withzyg --includeinfo
	table_annovar.pl ${Sample}_${variant_caller}.avinput --out ${Sample}_${variant_caller}.out --remove --protocol refGene,cytoBand,avsnp151,intervar_20180118,1000g2015aug_all,cosmic70,clinvar_20250721,gnomad211_exome --operation g,r,f,f,f,f,f,f --buildver hg38 --nastring '-1' --otherinfo --csvout \
	--thread ${task.cpus} /databases/humandb -xreffile /databases/gene_fullxref.txt

	touch ${Sample}_${variant_caller}.out.hg38_multianno.csv
	"""
	stub:
	"""
	touch ${Sample}_${variant_caller}.out.hg38_multianno.csv
	"""
}


