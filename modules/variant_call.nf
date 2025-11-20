#!/usr/bin/env nextflow

process FREEBAYES {
	label 'process_low'
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
		path (genome_dir)
		path (genome_fasta)
		file (bedfile)
	output:
		tuple val (Sample), file ("*.freebayes.vcf")
	script:
	"""
	${params.freebayes_path} -f ${genome_dir}/${genome_fasta} -b ${finalBam} -t ${bedfile} > ${Sample}.freebayes.vcf
	"""
}

process HAPLOTYPECALLER {
    label 'process_low'
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
		path (genome_dir)
		path (genome_fasta)
		file (bedfile)
		file (dbsnp)
	output:
		tuple val (Sample), file ("*.haplotypecaller.vcf")
	script:
	"""
	${params.java_path}/java -Xmx10G -jar ${params.GATK38_path} -T HaplotypeCaller -R ${genome_dir}/${genome_fasta} -I ${finalBam} -L ${bedfile} -o ${Sample}.haplotypecaller.vcf --dbsnp ${dbsnp} 
	"""
}

process STRELKA {
    label 'process_low'
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
		path (genome_dir)
		path (genome_fasta)
	output:
		tuple val (Sample), file ("*.strelka.vcf")
	script:
	"""
	${params.strelka_path}/configureStrelkaGermlineWorkflow.py --bam ${finalBam} --referenceFasta ${genome_dir}/${genome_fasta} --callRegions ${params.bedfile}.bed.gz --targeted --exome --runDir ./
	./runWorkflow.py -m local -j ${task.cpus}
	gunzip -f ./results/variants/variants.vcf.gz
	mv ./results/variants/variants.vcf ${Sample}.strelka.vcf
	"""
}

process PLATYPUS {
    label 'process_low'
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
		path (genome_dir)
		path (genome_fasta)
		file (bed_regions)
	output:
		tuple val (Sample), file("*.platypus.vcf")
	script:
	"""
	 python2.7 ${params.platypus_path} callVariants --bamFiles=${finalBam[0]} --refFile=${genome_dir}/${genome_fasta} --output=${Sample}.platypus.vcf --nCPU=${task.cpus} --minFlank=0 --filterDuplicates=0 --maxVariants=6 --minReads=6 --regions=${bed_regions}
	"""
}

process VARSCAN {
	label 'process_medium'
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
		path (genome_dir)
		path (genome_fasta)
		file (bedfile)
	output:
		tuple val (Sample), file ("*.varscan.vcf")
	script:
	"""
	${params.samtools} mpileup -l ${bedfile} -f ${genome_dir}/${genome_fasta} ${finalBam} > ${Sample}.mpileup
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
		path (genome_dir)
		path (genome_fasta)
		file (bedfile)
	output:
		tuple val (Sample), file("*.lofreq.filtered.vcf")
	script:
	"""
	${params.lofreq_path} viterbi -f ${params.genome} -o ${Sample}.lofreq.pre.bam ${oldfinalBam}
	${params.samtools} sort ${Sample}.lofreq.pre.bam > ${Sample}.lofreq.bam
	${params.lofreq_path} call -b dynamic -C 50 -a 0.00005 -q 30 -Q 30 -m 50 -f ${params.genome} -l ${bedfile} -o ${Sample}.lofreq.vcf ${Sample}.lofreq.bam
	${params.lofreq_path} filter -a 0.005 -i ${Sample}.lofreq.vcf -o ${Sample}.lofreq.filtered.vcf
	"""
}

process PINDEL {
	label 'process_low'
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
		path (genome_dir)
		path (genome_fasta)
	output:
		tuple val (Sample), file("*_pindel.vep.txt")
	script:
	"""
	export BAM_2_PINDEL_ADAPT=${params.pindel}/Adaptor.pm
	make_config.sh -s ${Sample} -b ${finalBam} -c config.txt
	${params.pindel}pindel -f ${genome_dir}/${genome_fasta} -i config.txt -c chr13 -o ${Sample}_pindel
	${params.pindel}pindel2vcf -r ${genome_dir}/${genome_fasta} -P ${Sample}_pindel -R hg19 -d 07102019 -v ${Sample}_pindel_SI.vcf
	
	#extracting required columns from ${Sample}_pindel_SI.vcf
	extract_pindelSI.py ${Sample}_pindel_SI.vcf ${Sample}extractedPindelSI.txt
	
	#using vep
	vep_analysis.sh ${Sample}_pindel_SI.vcf ${Sample}

	# extracting required columns from ${Sample}_vep_delheaders.txt
	extract_pindel.py  ${Sample}_vep_delheaders.txt ${Sample}extractedPindelVep.txt

	#merge extracted data
	mergepindel.py ${Sample}extractedPindelSI.txt ${Sample}extractedPindelVep.txt ${Sample}_pindel.vep.txt
	sed -i 's/SYMBOL/Gene/g' ${Sample}_pindel.vep.txt
	sed -i 's/Existing_variation/ID/g' ${Sample}_pindel.vep.txt
	"""
}

process DEEPSOMATIC {
	tag "${Sample}"
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val(Sample), file("*_DS.vcf")
	script:
	""" 
	mkdir output
	outpath=`realpath output`
	bam_path=`realpath ${finalBam} | awk 'BEGIN{OFS=FS="/"} {\$NF=""; print \$0}'`
	vcf_output=${Sample}_DS.vcf
	control_bam_path=`realpath /home/diagnostics/pipelines/WholeExomePipeline/scripts/cnvkit_lymphoma_pon_exonwise/ARPIT-LYMPHOMA.final.bam`
	echo \$bam_path \$outpath \$vcf_output ${finalBam} \${control_bam_path} ${params.genome} ${params.bedfile}.bed
	deepsomatic.sh \$bam_path \$outpath \$vcf_output ${finalBam} \${control_bam_path} ${params.genome} ${params.bedfile}.bed
	"""
}
