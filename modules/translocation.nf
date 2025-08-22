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
	samtools view -@ $task.cpus -b -F 1294 ${finalBam} >  ${Sample}.discordants.unsorted.bam
	samtools view -@ $task.cpus -h ${finalBam} | extractSplitReads_BwaMem.py -i stdin | samtools view -@ $task.cpus -Sb - > ${Sample}.splitters.unsorted.bam
	samtools sort -@ $task.cpus ${Sample}.discordants.unsorted.bam > ${Sample}.discordants.bam
	samtools sort -@ $task.cpus ${Sample}.splitters.unsorted.bam > ${Sample}.splitters.bam
	lumpyexpress -B ${finalBam} -S ${Sample}.splitters.bam -D ${Sample}.discordants.bam -o ${Sample}.vcf
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
