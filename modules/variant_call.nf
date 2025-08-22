#!/usr/bin/env nextflow

process FREEBAYES {
	label 'process_low'
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.freebayes.vcf")
	script:
	"""
	${params.freebayes_path} -f ${params.genome} -b ${finalBam} > ${Sample}.freebayes.vcf
	"""
}

process HAPLOTYPECALLER {
    label 'process_medium'
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.haplotypecaller.vcf")
	script:
	"""
	${params.java_path}/java -Xmx10G -jar ${params.GATK38_path} -T HaplotypeCaller -R ${params.genome} -I ${finalBam} -o ${Sample}.haplotypecaller.vcf --dbsnp ${params.site2} 
	"""
}

process STRELKA {
    label 'process_low'
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.strelka.vcf")
	script:
	"""
	${params.strelka_path}/configureStrelkaGermlineWorkflow.py --bam ${finalBam} --referenceFasta ${params.genome} --targeted --exome --runDir ./
	./runWorkflow.py -m local -j $task.cpus
	gunzip -f ./results/variants/variants.vcf.gz
	mv ./results/variants/variants.vcf ${Sample}.strelka.vcf
	"""
}

process PLATYPUS {
    label 'process_low'
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file("*.platypus.vcf")
	script:
	"""
	 python2.7 ${params.platypus_path} callVariants --bamFiles=${finalBam[0]} --refFile=${params.genome} --output=${Sample}.platypus.vcf --nCPU=$task.cpus --minFlank=0 --filterDuplicates=0 --maxVariants=6 --minReads=6
	"""
}

process VARSCAN {
	label 'process_medium'
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.varscan.vcf")
	script:
	"""
	${params.samtools} mpileup -f ${params.genome} ${finalBam} > ${Sample}.mpileup
	${params.java_path}/java -jar ${params.varscan_path} mpileup2snp ${Sample}.mpileup --min-coverage 10 --min-reads2 5 --min-avg-qual 15 --min-var-freq 0.003 --p-value 1e-4 --output-vcf 1 > ${Sample}.varscan_snp.vcf
	${params.java_path}/java -jar ${params.varscan_path} mpileup2indel ${Sample}.mpileup --min-coverage 10 --min-reads2 5 --min-avg-qual 15 --min-var-freq 0.003 --p-value 1e-4 --output-vcf 1 > ${Sample}.varscan_indel.vcf
	bgzip -c ${Sample}.varscan_snp.vcf > ${Sample}.varscan_snp.vcf.gz
	bgzip -c ${Sample}.varscan_indel.vcf > ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_snp.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} concat -a ${Sample}.varscan_snp.vcf.gz ${Sample}.varscan_indel.vcf.gz -o ${Sample}.varscan.vcf
	"""
}

process DEEPVARIANT {
	label 'process_low'
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*.vcf'
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("${Sample}_deepvar_vep_annonvar.txt"), file("${Sample}_deepvar_filt.vcf")
	script:
	"""
	genome_path=`realpath ${params.genome}`
	bam_path=`realpath ${finalBam} | awk 'BEGIN{OFS=FS="/"} {\$NF=""; print \$0}'`
	pwd=`realpath ./`
	${params.deepvariant} \${bam_path} \${pwd} ${Sample}_deepvar.vcf ${finalBam} \${genome_path} ${params.bedfile}.bed
	${params.bcftools_path} view -f PASS ${Sample}_deepvar.vcf > ${Sample}_deepvar_filt.vcf
	# Combining data from vcf and vep annotation data
	${params.vep_wrapper} ${Sample} ${Sample}_deepvar_filt.vcf ${Sample}_deepvariant.txt

	#Annotating ${Sample}.somaticseq.vcf using annovar
	${params.annovar_wrapper} ${params.annovarLatest_path} ${Sample}_deepvar_filt.vcf ${Sample}

	#extracting columns  Func.refGene,Gene.refGene,ExonicFunc.refGene,PopFreqMax,InterVar_automated from somaticseq.hg19_multianno.csv and adding them to somaticseq.vep.txt
	python3 ${params.deepvar_annovar} ${Sample}.annovar.hg19_multianno.csv ${Sample}_deepvariant.txt ${Sample}_deepvar_vep_annonvar.txt
	"""
}

process LOFREQ {
	label 'process_low'
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file("*.lofreq.filtered.vcf")
	script:
	"""
	${params.lofreq_path} viterbi -f ${params.genome} -o ${Sample}.lofreq.pre.bam ${oldfinalBam}
	${params.samtools} sort ${Sample}.lofreq.pre.bam > ${Sample}.lofreq.bam
	${params.lofreq_path} call -b dynamic -C 50 -a 0.00005 -q 30 -Q 30 -m 50 -f ${params.genome} -l ${params.bedfile}.bed -o ${Sample}.lofreq.vcf ${Sample}.lofreq.bam
	${params.lofreq_path} filter -a 0.005 -i ${Sample}.lofreq.vcf -o ${Sample}.lofreq.filtered.vcf
	"""
}

process PINDEL {
	label 'process_low'
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file("*_pindel.vep.txt")
	script:
	"""
	export BAM_2_PINDEL_ADAPT=${params.pindel}/Adaptor.pm
	sh ${params.pindel_config_script} -s ${Sample} -b ${finalBam} -c config.txt
	${params.pindel}pindel -f ${params.genome} -i config.txt -c chr13 -o ${Sample}_pindel
	${params.pindel}pindel2vcf -r ${params.genome} -P ${Sample}_pindel -R hg19 -d 07102019 -v ${Sample}_pindel_SI.vcf
	
	#extracting required columns from ${Sample}_pindel_SI.vcf
	${params.extract_pindelSI._script_path} ${Sample}_pindel_SI.vcf ${Sample}extractedPindelSI.txt
	
	#using vep
	${params.vep_script_path} ${Sample}_pindel_SI.vcf ${Sample}

	# extracting required columns from ${Sample}_vep_delheaders.txt
	${params.extract_pindel_script_path}  ${Sample}_vep_delheaders.txt ${Sample}extractedPindelVep.txt
	#merge extracted data
	${params.mergepindel_script_path} ${Sample}extractedPindelSI.txt ${Sample}extractedPindelVep.txt ${Sample}_pindel.vep.txt
	sed -i 's/SYMBOL/Gene/g' ${Sample}_pindel.vep.txt
	sed -i 's/Existing_variation/ID/g' ${Sample}_pindel.vep.txt
	"""
}


process SOMATICSEQ {
	label 'process_low'
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