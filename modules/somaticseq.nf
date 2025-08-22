#!/usr/bin/env nextflow

process SOMATICSEQ {
	label 'process_medium'
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*_somaticseq.vep_annonvar.txt'
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '${Sample}.somaticseq.vcf'
	input:
		tuple val (Sample), file (lofreqVcf), file (varscanVcf), file (platypusVcf), file (strelkaVcf), file (haplotypecallerVcf), file (freebayesVcf), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*_somaticseq.vep_annonvar.txt"), file("${Sample}.somaticseq.vcf")
	script:
	"""
	${params.vcf_sorter_path} ${freebayesVcf} ${Sample}.freebayes.sorted.vcf
	${params.vcf_sorter_path} ${platypusVcf} ${Sample}.platypus.sorted.vcf
	${params.vcf_sorter_path} ${haplotypecallerVcf} ${Sample}.haplotypecaller.sorted.vcf

	python3 ${params.splitvcf_path} -infile ${Sample}.platypus.sorted.vcf -snv ${Sample}_platypus_cnvs.vcf -indel ${Sample}_platypus_indels.vcf
	python3 ${params.splitvcf_path} -infile ${Sample}.freebayes.sorted.vcf -snv ${Sample}_freebayes_cnvs.vcf -indel ${Sample}_freebayes_indels.vcf
	python3 ${params.splitvcf_path} -infile ${Sample}.haplotypecaller.sorted.vcf -snv ${Sample}_haplotypecaller_cnvs.vcf -indel ${Sample}_haplotypecaller_indels.vcf

	${params.vcf_sorter_path} ${Sample}_platypus_cnvs.vcf ${Sample}_platypus_cnvs_sort.vcf
	${params.vcf_sorter_path} ${Sample}_platypus_indels.vcf ${Sample}_platypus_indels_sort.vcf
	${params.vcf_sorter_path} ${Sample}_freebayes_cnvs.vcf ${Sample}_freebayes_cnvs_sort.vcf
	${params.vcf_sorter_path} ${Sample}_freebayes_indels.vcf ${Sample}_freebayes_indels_sort.vcf
	${params.vcf_sorter_path} ${Sample}_haplotypecaller_cnvs.vcf ${Sample}_haplotypecaller_cnvs_sort.vcf
	${params.vcf_sorter_path} ${Sample}_haplotypecaller_indels.vcf ${Sample}_haplotypecaller_indels_sort.vcf

	somaticseq_parallel.py --output-directory ${Sample}.somaticseq --genome-reference ${params.genome} --inclusion-region ${params.bedfile}.bed --threads 25 --pass-threshold 0 --lowqual-threshold 0 --algorithm xgboost -minMQ 0 -minBQ 0 -mincaller 0 --dbsnp-vcf /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file ${finalBam} --varscan-vcf ${varscanVcf} --lofreq-vcf ${lofreqVcf} --strelka-vcf ${strelkaVcf} --sample-name ${Sample} --arbitrary-snvs ${Sample}_freebayes_cnvs_sort.vcf ${Sample}_platypus_cnvs_sort.vcf ${Sample}_haplotypecaller_cnvs_sort.vcf --arbitrary-indels ${Sample}_freebayes_indels_sort.vcf ${Sample}_platypus_indels_sort.vcf ${Sample}_haplotypecaller_indels_sort.vcf

	${params.vcf_sorter_path} ${Sample}.somaticseq/Consensus.sSNV.vcf ${Sample}.somaticseq/somaticseq_snv.vcf
	bgzip -c ${Sample}.somaticseq/somaticseq_snv.vcf > ${Sample}.somaticseq/somaticseq_snv.vcf.gz
	${params.bcftools_path} index -t ${Sample}.somaticseq/somaticseq_snv.vcf.gz

	${params.vcf_sorter_path} ${Sample}.somaticseq/Consensus.sINDEL.vcf ${Sample}.somaticseq/somaticseq_indel.vcf
	bgzip -c ${Sample}.somaticseq/somaticseq_indel.vcf > ${Sample}.somaticseq/somaticseq_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.somaticseq/somaticseq_indel.vcf.gz

	${params.bcftools_path} concat -a ${Sample}.somaticseq/somaticseq_snv.vcf.gz ${Sample}.somaticseq/somaticseq_indel.vcf.gz -o ${Sample}.somaticseq.vcf

	sed -i 's/##INFO=<ID=VLK012,Number=6,Type=Integer,Description="Calling decision of the 6 algorithms: VarScan2, LoFreq, Strelka, SnvCaller_0, SnvCaller_1, SnvCaller_2">/##INFO=<ID=VLSFPH,Number=6,Type=String,Description="Calling decision of the 6 algorithms:  VarScan2, LoFreq, Strelka, Freebayes, Platypus, Haplotypecaller">/g' ${Sample}.somaticseq.vcf
	sed -i 's/VLK012/VLSFPH/g' ${Sample}.somaticseq.vcf
	
	# to extract vaf,af,alt and ref count
	${params.extract_somatic_script_path} ${Sample}.somaticseq.vcf ${Sample}.extractedSomaticseq.txt
	#adding vep
	${params.vep_script_path} ${Sample}.somaticseq.vcf ${Sample}
	${params.extract_velheader_script_path} ${Sample}_vep_delheaders.txt ${Sample}.extractedvepdelheaders.txt
	
	# for merging extracted somaticsseq and vepheaders
	${params.mergeSomaticvep_script_path} ${Sample}.extractedSomaticseq.txt ${Sample}.extractedvepdelheaders.txt ${Sample}_somaticseq.vep.txt
	sed -i 's/SYMBOL/Gene/g' ${Sample}_somaticseq.vep.txt
	
	#Annotating ${Sample}.somaticseq.vcf using annovar
	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}.somaticseq.vcf --outfile ${Sample}.somaticseq.avinput -withfreq --includeinfo -allsample
	
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.somaticseq.avinput --out ${Sample}.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt
	
	#extracting columns  Func.refGene,Gene.refGene,ExonicFunc.refGene,PopFreqMax,InterVar_automated from somaticseq.hg19_multianno.csv and adding them to somaticseq.vep.txt
	python3 ${params.extract_annovar} ${Sample}.somaticseq.hg19_multianno.csv ${Sample}_somaticseq.vep.txt ${Sample}_somaticseq.vep_annonvar.txt
	"""
} 