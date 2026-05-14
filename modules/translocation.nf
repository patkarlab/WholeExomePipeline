#!/usr/bin/env nextflow

process SVABA {
	tag "${Sample}"
	label "process_low"
	publishDir "${params.output}/${Sample}/svaba/", mode: 'copy', pattern: '*svaba*'
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("${Sample}_svaba.sortd.bam"), file ("${Sample}_svaba.svaba.sv.vcf")
	script:
	"""
	${params.svaba_path} run -t ${finalBam} -G ${params.genome} -p $task.cpus -k ${params.trans_bedfile}.bed -D ${params.site2} -a ${Sample}_svaba
	${params.samtools} sort -@ $task.cpus ${Sample}_svaba.contigs.bam -o ${Sample}_svaba.sortd.bam
	${params.samtools} index -@ $task.cpus ${Sample}_svaba.sortd.bam
	"""
}

process SVABA_LYMPHOMA {
	tag "${Sample}"
	label "process_low"
	publishDir "${params.output}/${Sample}/translocatn/svaba", mode: 'copy', pattern: '*svaba*'
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
		file (normal_bamfile)
		file (normal_bamBaifile)
	output:
		tuple val (Sample), file ("*svaba.somatic.sv.vcf"), file ("*svaba.somatic.indel.vcf"), file("${Sample}.svaba.somatic.sv.tsv"), file("${Sample}.svaba.somatic.indel.tsv"), file("${Sample}.svaba.germline.sv.tsv"), file("${Sample}.svaba.germline.indel.tsv")
	script:
	"""
	${params.svaba_path} run -t ${oldfinalBam} -n ${normal_bamfile} -G ${params.genome} -p $task.cpus -D ${params.site2} -a ${Sample}_svaba
	format_csv.py ${Sample}_svaba.svaba.somatic.sv.vcf ${Sample}.svaba.somatic.sv.tsv
	format_csv.py ${Sample}_svaba.svaba.germline.sv.vcf ${Sample}.svaba.germline.sv.tsv
	format_csv.py ${Sample}_svaba.svaba.somatic.indel.vcf ${Sample}.svaba.somatic.indel.tsv
	format_csv.py ${Sample}_svaba.svaba.germline.indel.vcf ${Sample}.svaba.germline.indel.tsv	
	"""
}

// process LUMPY {
// 	tag "${Sample}"
// 	label "process_low"
// 	input:
// 		tuple val(Sample), file(read1), file(read2)
// 	output:
// 		tuple val (Sample), file("${Sample}.vcf")
// 	script:
// 	"""
// 	speedseq align -R "@RG\tID:id\tSM:${Sample}\tLB:lib" ${params.genome} ${read1} ${read2} -t $task.cpus -o ${Sample} 
// 	lumpyexpress -B ${read1}.bam -S ${read1}.splitters.bam -D ${read1}.discordants.bam -o ${Sample}.vcf
// 	"""
// }

process LUMPY {
	tag "${Sample}"
	label "process_low"
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file("${Sample}.vcf")
	script:
	"""
	samtools view -@ ${task.cpus} -b -F 1294 ${oldfinalBam} >  ${Sample}.discordants.unsorted.bam
	samtools view -@ ${task.cpus} -h ${oldfinalBam} | extractSplitReads_BwaMem.py -i stdin | samtools view -@ $task.cpus -Sb - > ${Sample}.splitters.unsorted.bam
	samtools sort -@ ${task.cpus} ${Sample}.discordants.unsorted.bam > ${Sample}.discordants.bam
	samtools sort -@ ${task.cpus} ${Sample}.splitters.unsorted.bam > ${Sample}.splitters.bam
	lumpyexpress -B ${oldfinalBam} -S ${Sample}.splitters.bam -D ${Sample}.discordants.bam -o ${Sample}.vcf
	"""
}

process GRIDSS {
	tag "${Sample}"
	label "process_low"
	publishDir "${params.output}/${Sample}/translocatn/gridss", mode: 'copy', pattern: '*vcf'
	publishDir "${params.output}/${Sample}/translocatn/gridss", mode: 'copy', pattern: '*tsv'
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
		file (normal_bamfile)
		file (normal_bamBaifile)
	output:
		tuple val (Sample), file ("*")
	script:
	"""
	gridss_normal.sh ${params.genome} ${oldfinalBam} ${normal_bamfile} ${Sample}
	format_csv.py ${Sample}.gridss.vcf ${Sample}.gridss.tsv
	format_csv.py ${Sample}.somatic.gridss.vcf ${Sample}.somatic.gridss.tsv
	format_csv.py RM_${Sample}.somatic.gridss.vcf RM_${Sample}.somatic.gridss.tsv                                                     
	"""
}

process DELLY {
	tag "${Sample}"
	label "process_low"
	publishDir "${params.output}/${Sample}/translocatn/delly", mode: 'copy', pattern: '*vcf'
	publishDir "${params.output}/${Sample}/translocatn/delly", mode: 'copy', pattern: '*tsv'
	input:
		tuple val (Sample), file(finalBam), file(finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
		file (normal_bamfile)
		file (normal_bamBaifile)
	output:
		tuple val (Sample), file ("${Sample}_delly_somatic.vcf"), file ("${Sample}_delly_somatic.tsv"), file ("${Sample}_delly.vcf"), file("${Sample}_delly.tsv")
	script:
	"""
	tumor_name=`basename ${oldfinalBam} .old_final.bam`
	control_name=`basename ${normal_bamfile} _final.bam`

	echo -e "\${tumor_name}\ttumor" > samples.tsv
	echo -e "\${control_name}\tcontrol" >> samples.tsv

	${params.delly} call \
		-g ${params.genome} \
		${normal_bamfile} ${oldfinalBam} \
		-o ${Sample}_delly.bcf
	
	bcftools view ${Sample}_delly.bcf -Ov -o ${Sample}_delly.vcf

	
	${params.delly} filter \
		-f somatic \
		-s samples.tsv \
		${Sample}_delly.bcf \
	| bcftools view -Ov -o ${Sample}_delly_somatic.vcf

	format_csv.py ${Sample}_delly.vcf ${Sample}_delly.tsv
	format_csv.py ${Sample}_delly_somatic.vcf ${Sample}_delly_somatic.tsv
	"""
}

process MANTA {
	tag "${Sample}"
	label "process_low"
	errorStrategy 'ignore'
	publishDir "${params.output}/${Sample}/translocatn/manta", mode: 'copy', pattern: '*vcf'
	publishDir "${params.output}/${Sample}/translocatn/manta", mode: 'copy', pattern: '*tsv'
	input:
		tuple val (Sample), file(finalBam), file(finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
		file (normal_bamfile)
		file (normal_bamBaifile)
	output:
		tuple val (Sample), file ("${Sample}_manta_candidateSV.vcf"), file ("${Sample}_manta_diploidSV.vcf"), file ("${Sample}_manta_somaticSV.vcf"), file ("${Sample}_manta_somaticSV.tsv")
	script:
	"""
	${params.manta_path}/configManta.py \
		--normalBam ${normal_bamfile} \
		--tumorBam ${oldfinalBam} \
		--referenceFasta ${params.genome} --runDir ./ \
		--exome

	./runWorkflow.py -j ${task.cpus}

	gunzip -c ./results/variants/candidateSV.vcf.gz > ${Sample}_manta_candidateSV.vcf

	gunzip -c ./results/variants/diploidSV.vcf.gz > ${Sample}_manta_diploidSV.vcf

	gunzip -c ./results/variants/somaticSV.vcf.gz > ${Sample}_manta_somaticSV.vcf
	format_manta_csv.py ${Sample}_manta_somaticSV.vcf ${Sample}_manta_somaticSV.tsv

	"""	
}

process TRANSLOCATION {
	tag "${Sample}"
	publishDir "${params.output}/${Sample}/translocatn", mode: 'copy', pattern: '*.tsv'	
	input:
		tuple val (Sample), file (svaba_bam), file (svaba_vcf), file (lumpy_vcf)
	output:
		tuple val (Sample), file ("*.tsv")
	script:
	"""
	svaba_lumpy_common.py ${svaba_vcf} ${lumpy_vcf} ${Sample}_common
	svaba_file=\$(basename ${svaba_vcf} .vcf)
	lumpy_file=\$(basename ${lumpy_vcf} .vcf)
	cat \${svaba_file}.tsv \${lumpy_file}.tsv ${Sample}_common > ${Sample}_translocatns.tsv 
	"""
}

process TRANSLOCATION_LYMPHOMA {
	tag "${Sample}"
	publishDir "${params.output}/${Sample}/translocatn", mode: 'copy', pattern: '*_svaba_lumpy.tsv'	
	input:
		tuple val (Sample), file (svaba_somatic_sv_vcf), file (svaba_somatic_indel_vcf), file (svaba_somatic_sv_tsv), file (svaba_somatic_indel_tsv), file (svaba_germline_sv_tsv), file (svaba_germline_indel_tsv), file (lumpy_vcf)
	output:
		tuple val (Sample), file ("${Sample}_svaba_lumpy.tsv")
	script:
	"""
	svaba_lumpy_common.py ${svaba_somatic_sv_vcf} ${lumpy_vcf} ${Sample}_common
	svaba_file=\$(basename ${svaba_somatic_sv_vcf} .vcf)
	lumpy_file=\$(basename ${lumpy_vcf} .vcf)
	cat \${svaba_file}.tsv \${lumpy_file}.tsv ${Sample}_common > ${Sample}_svaba_lumpy.tsv 
	"""
}

