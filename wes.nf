#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=
Sample list: ${params.input}
BED file: ${params.bedfile}.bed
Sequences in:${params.sequences}
"""
process fastqc{
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*'
	input:
		val Sample
	output:
		path "*"
	"""
	${params.fastqc} -o ./ -f fastq ${params.sequences}/${Sample}_*R1_*.fastq.gz ${params.sequences}/${Sample}_*R2_*.fastq.gz
	"""
}

process trimming_trimmomatic { 
	input:
		val Sample
	output:
		tuple val (Sample), file("*1P.fq.gz"), file("*2P.fq.gz")
	script:
	"""
	${params.trimmomatic_path}trimmomatic PE \
	${params.sequences}/${Sample}_*R1_*.fastq.gz ${params.sequences}/${Sample}_*R2_*.fastq.gz -baseout ${Sample}.fq.gz \
	ILLUMINACLIP:${params.adaptors}:2:30:10:2:keepBothReads \
	LEADING:3 SLIDINGWINDOW:4:15 MINLEN:40
	sleep 5s
	"""
}

process pair_assembly_pear {	
	input:
		tuple val (Sample), file(paired_forward), file(paired_reverse)
	output:
		tuple val (Sample), file("*.assembled.fastq")
	script:
	"""
	${params.pear_path} -f ${paired_forward} -r ${paired_reverse} -o ${Sample} -n 53 -j 25	
	"""
}

process mapping_reads {
	input:
		tuple val (Sample), file(pairAssembled)
	output: 
		tuple val (Sample), file ("*.sam")	
	script:
	"""
	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:200" -M -t 20 ${params.genome} ${pairAssembled} > ${Sample}.sam 
	"""
}

process mapping_both_reads {
	input:
		tuple val (Sample), file(paired_forward), file(paired_reverse)
	output:
		tuple val (Sample), file ("*.sam")
	script:
	"""
	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:200" -M -t 20 ${params.genome} ${paired_forward} ${paired_reverse} > ${Sample}.sam
	"""
}

process sam_conversion {
	input:
		tuple val (Sample), file (samfile)
	output:
		tuple val (Sample), file ("*.sorted.bam"), file("*.sorted.bam.bai")
	script:
	"""
	${params.samtools} view -bT ${params.genome} ${samfile} > ${Sample}.bam 
	${params.samtools} sort ${Sample}.bam > ${Sample}.sorted.bam
	${params.samtools} index ${Sample}.sorted.bam > ${Sample}.sorted.bam.bai
	"""
}

process mark_duplicates {
	input:
		tuple val (Sample), file (sorted_bam), file (sorted_bam_index)
	output:
		tuple val (Sample), file ("*.bam"), file ("*.bam.bai"),  file ("*.txt")
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} MarkDuplicates \
	I=${sorted_bam} \
	O=${Sample}_sorted_marked.bam \
	M=${Sample}_picard.info.txt
	REMOVE_DUPLICATES=true 
	${params.samtools} index ${Sample}_sorted_marked.bam > ${Sample}_sorted_marked.bam.bai
	"""
}

process RealignerTargetCreator {
	input:
		tuple val (Sample), file (bamFile), file (bamBai), file (pi_card_info)
	output:
		tuple val (Sample), file ("*.intervals")
	script:
	"""
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T RealignerTargetCreator -R ${params.genome} -nt 10 -I ${bamFile} --known ${params.site1} -o ${Sample}_target.intervals
	"""
}

process IndelRealigner {
	input:
		tuple val (Sample), file (target_interval), file (bamFile), file (bamBai), file (pi_card_info)
	output:
		tuple val (Sample), file ("*.realigned.bam")
	script:
	"""
	echo ${Sample} ${target_interval} ${bamFile}
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T IndelRealigner -R ${params.genome} -I ${bamFile} -known ${params.site1} --targetIntervals ${target_interval} -o ${Sample}.realigned.bam
	"""
}

process BaseRecalibrator {
	input:
		tuple val (Sample), file (realignedBam)
	output:
		tuple val(Sample), file ("*.recal_data.table")
	script:
	"""
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T BaseRecalibrator -R ${params.genome} -I ${realignedBam} -knownSites ${params.site2} -knownSites ${params.site3} -maxCycle 600 -o ${Sample}.recal_data.table
	"""
}

process PrintReads {
	input:
		tuple val (Sample), file (realigned_Bam), file (recal_data_table)
	output:
		tuple val (Sample), file ("*.aligned.recalibrated.bam")
	script:
	"""
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T PrintReads -R ${params.genome} -I ${realigned_Bam} --BQSR ${recal_data_table} -o ${Sample}.aligned.recalibrated.bam
	"""
}

process generatefinalbam {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.final.bam*'
	input:
		tuple val (Sample), file (aligned_recalibrated_bam)
	output:
		tuple val(Sample), file ("*.final.bam"), file ("*.final.bam.bai"), file ("*.old_final.bam"), file ("*.old_final.bam.bai")
	script:
	"""
	${params.java_path}/java -Xmx16G -jar ${params.abra2_path}/abra2-2.23.jar --in ${aligned_recalibrated_bam} --out ${Sample}.abra.bam --ref ${params.genome} --threads 8 --tmpdir ./ > abra.log
	${params.samtools} sort ${aligned_recalibrated_bam} > ${Sample}.old_final.bam
	${params.samtools} index ${Sample}.old_final.bam > ${Sample}.old_final.bam.bai
	${params.samtools} sort ${Sample}.abra.bam > ${Sample}.final.bam
	${params.samtools} index ${Sample}.final.bam > ${Sample}.final.bam.bai
	"""
}

process generatefinalbamin {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.final.bam*'
	input:
		val (Sample)
	output:
		tuple val(Sample), file ("*.final.bam"), file ("*.final.bam.bai"), file ("*.old_final.bam"), file ("*.old_final.bam.bai")
	script:
	"""
	${params.java_path}/java -Xmx16G -jar ${params.abra2_path}/abra2-2.23.jar --in ${params.sequences}/${Sample}*.bam --out ${Sample}.abra.bam --ref ${params.genome} --threads 8 --tmpdir ./ > abra.log
	${params.samtools} sort ${params.sequences}/${Sample}*.bam > ${Sample}.old_final.bam
	${params.samtools} index ${Sample}.old_final.bam > ${Sample}.old_final.bam.bai
	${params.samtools} sort ${Sample}.abra.bam > ${Sample}.final.bam
	${params.samtools} index ${Sample}.final.bam > ${Sample}.final.bam.bai
	"""
}

process hsmetrics_run{
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_hsmetrics.txt'
	input:
		tuple val(Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*_hsmetrics.txt")
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} CollectHsMetrics I= ${finalBam} O= ${Sample}_hsmetrics.txt BAIT_INTERVALS= ${params.bedfile}.interval_list TARGET_INTERVALS= ${params.bedfile}.interval_list R= ${params.genome} VALIDATION_STRINGENCY=LENIENT
	${params.hsmetrics_all} $PWD/Final_Output/hsmetrics.tsv ${Sample} ${Sample}_hsmetrics.txt
	"""
}

process minimap_getitd {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_getitd'
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		path "*_getitd"
	script:
	"""
	${params.samtools} sort ${finalBam} -o ${Sample}.sorted.bam
	${params.samtools} index ${Sample}.sorted.bam
	${params.samtools} view ${Sample}.sorted.bam -b -h chr13 > ${Sample}.chr13.bam
	${params.bedtools} bamtofastq -i ${Sample}.chr13.bam -fq ${Sample}_chr13.fastq
	python ${params.get_itd_path}/getitd.py -reference ${params.get_itd_path}/anno/amplicon.txt -anno ${params.get_itd_path}/anno/amplicon_kayser.tsv -forward_adapter AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -reverse_adapter CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT -nkern 8 ${Sample} ${Sample}_chr13.fastq
	"""
}

process coverage_mosdepth {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*'
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*")
	script:
	"""
	${params.mosdepth_script} ${finalBam} ${Sample}_cov ${params.bedfile}.bed
	${params.extract_COV50_script_path} ${Sample}_cov.regions.bed ${Sample}_median50
	sleep 2s
	"""
}

process freebayes {
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.freebayes.vcf")
	script:
	"""
	${params.freebayes_path} -f ${params.genome} -b ${finalBam} > ${Sample}.freebayes.vcf
	"""
}

process haplotypecaller {
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.haplotypecaller.vcf")
	script:
	"""
	${params.java_path}/java -Xmx10G -jar ${params.GATK38_path} -T HaplotypeCaller -R ${params.genome} -I ${finalBam} -o ${Sample}.haplotypecaller.vcf --dbsnp ${params.site2} 
	"""
}

process strelka {
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.strelka.vcf")
	script:
	"""
	${params.strelka_path}/configureStrelkaGermlineWorkflow.py --bam ${finalBam} --referenceFasta ${params.genome} --targeted --exome --runDir ./
	./runWorkflow.py -m local -j 20
	gunzip -f ./results/variants/variants.vcf.gz
	mv ./results/variants/variants.vcf ${Sample}.strelka.vcf
	"""
}

process platypus {
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file("*.platypus.vcf")
	script:
	"""
	 python2.7 ${params.platypus_path} callVariants --bamFiles=${finalBam[0]} --refFile=${params.genome} --output=${Sample}.platypus.vcf --nCPU=15 --minFlank=0 --filterDuplicates=0 --maxVariants=6 --minReads=6
	"""
}

process vardict {
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.vardict.vcf")
	script:
	"""
	VarDict -G ${params.genome} -f 0.03 -N ${Sample} -b ${finalBam} -c 1 -S 2 -E 3 -g 4 ${params.bedfile} | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${Sample} -E -f 0.03 > ${Sample}.vardict.vcf
	"""
}

process varscan {
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

process pindel {
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

process format_pindel{
	 publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_final.pindel.csv'
	 input:
	 	tuple val (Sample), file(pindelvep), file (bedfile)
	output:
		tuple val (Sample), file ("*")
	script:
	"""
	python3 ${params.format_pindel_wes_script_path} $PWD/Final_Output/${Sample}/${Sample}_cov.regions.bed ${pindelvep} $PWD/Final_Output/${Sample}/${Sample}_final.pindel.csv
	"""
}

process lofreq {
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

process somaticSeq_run {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_somaticseq.vep_annonvar.txt'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '${Sample}.somaticseq.vcf'
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

process cnvnator {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_Cnvnator_calls.tsv'
	input:
		tuple val(Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val(Sample), file ("*_Cnvnator_calls.tsv")
	script:
	"""
	${params.cnvnator_path} ${finalBam} ./ ${Sample}
	"""
}

process cava {
	input:
		tuple val(Sample), file(somaticseqVepAnnovar), file (somaticVcf)
	output:
		tuple val(Sample), file ("*.cava.csv")
	script:
	"""
	${params.cava_path}/cava -c ${params.cava_path}/config_v2.txt -t 10 -i ${somaticVcf} -o ${Sample}.somaticseq
	python3 ${params.cava_script_path} ${Sample}.somaticseq.txt ${Sample}.cava.csv
	"""
}

process merge_csv {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.xlsx'
	input:
		tuple val (Sample), file (pindelVep), file (somaticseqVep), file (somaticVcf), file (DeepVarOut), file (DeepVarVcf), file (cavaCsv)
	output:
		val Sample
	script:
	"""
	python3 ${params.pharma_marker_script} $PWD/Final_Output/${Sample}/${Sample}_somaticseq.vep_annonvar.txt ${params.pharma_input_xlxs} Pharma.tsv
	python3 ${params.merge_csvs_script} ${Sample} ${PWD}/Final_Output/${Sample}/${Sample}.xlsx  ${cavaCsv} $PWD/Final_Output/${Sample}/${Sample}_cov.mosdepth.summary.txt $PWD/Final_Output/${Sample}/${Sample}_cov.regions.bed $PWD/Final_Output/${Sample}/${Sample}_median50 ${pindelVep} ${somaticseqVep} Pharma.tsv ${DeepVarOut}
	sleep 1s
	"""
}

process deepvariant {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.vcf'
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

process ifcnv_run {
	conda '/home/miniconda3/envs/new_base/'
	input:
		val Sample
	output:
		//tuple val (Sample), file ("*.html"), file ("*.tsv")
		val (Sample)
	script:
	"""
	${params.links} $PWD/Final_Output/ ${params.input}
	mkdir ifCNV
	${params.ifcnv} ./ ${params.bedfile}.bed ifCNV

	# Making ifCNV's output directory for each sample
	for i in `cat ${params.input}`
	do
		if [ ! -d $PWD/Final_Output/\${i}/ifCNV ]; then
			mkdir $PWD/Final_Output/\${i}/ifCNV
		fi		
	done

	# Copying output of ifCNV to respective samples
	if [ -f ifCNV/ifCNV.tsv ]; then
		for i in `awk 'NR>1{print \$3}' ifCNV/ifCNV.tsv | awk 'BEGIN{FS="."}{print \$1}' | sort | uniq`	
		do
			cp ifCNV/\${i}.*.html $PWD/Final_Output/\${i}/ifCNV/	
		done
	fi
	"""
}

process GermlineCNVAnalysis {
	conda '/home/diagnostics/.conda/envs/gatk'
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	
	output:
		tuple val (Sample), file("*-genotyped-intervals.vcf.gz"), file("*-genotyped-segments.vcf.gz"), file ("*-denoised-copy-ratios.tsv")
	
	script:
	"""
	# CollectReadCounts
	${params.gatk_path} CollectReadCounts -L ${params.preprocessed_intervals} -R ${params.genome} -imr OVERLAPPING_ONLY -I ${finalBam} --format TSV -O ${Sample}.tsv
	
	# DetermineGermlineContigPloidy_case
	${params.gatk_path} DetermineGermlineContigPloidy --model ${params.ploidy_model} -I ${Sample}.tsv -O . --output-prefix ${Sample}-ploidy --verbosity DEBUG
	
	# GermlineCNVCaller_case
	${params.gatk_path} GermlineCNVCaller --run-mode CASE -I ${Sample}.tsv --contig-ploidy-calls ${Sample}-ploidy-calls --model ${params.cohort_model} --output ${Sample}_gCNV --output-prefix ${Sample} --verbosity DEBUG
	
	# PostprocessGermlineCNVCalls
	${params.gatk_path} PostprocessGermlineCNVCalls --model-shard-path ${params.cohort_model} --calls-shard-path ${Sample}_gCNV/${Sample}-calls --allosomal-contig chrX --allosomal-contig chrY --contig-ploidy-calls ${Sample}-ploidy-calls --sample-index 0 --output-genotyped-intervals ${Sample}-genotyped-intervals.vcf.gz --output-genotyped-segments ${Sample}-genotyped-segments.vcf.gz --sequence-dictionary /home/reference_genomes/hg19_broad/hg19_all.dict --output-denoised-copy-ratios ${Sample}-denoised-copy-ratios.tsv
	"""
}

process deepvariant_gCNV {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.vcf'
	input:
		tuple val (Sample), file (DeepVarOut), file (DeepVarVcf), file(genotypedIntervals), file (genotypedSegments), file(copyRatios)

	output:
		tuple val (Sample), file("${Sample}_deepvar_gCNV_merged.vcf")

	script:
	"""
	bgzip ${DeepVarVcf} 
	bcftools index ${DeepVarVcf}.gz
	bcftools index ${genotypedSegments}
	bcftools concat -a ${genotypedSegments} ${DeepVarVcf}.gz -o ${Sample}_deepvar_gCNV_merged.vcf	
	"""
}

workflow WES {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.map{ it }
		.set { samples_ch }

	main:
		fastqc(samples_ch)
		trimming_trimmomatic(samples_ch) | pair_assembly_pear | mapping_reads | sam_conversion | mark_duplicates
		RealignerTargetCreator(mark_duplicates.out)
		IndelRealigner(RealignerTargetCreator.out.join(mark_duplicates.out)) | BaseRecalibrator
		PrintReads(IndelRealigner.out.join(BaseRecalibrator.out)) | generatefinalbam
		//minimap_getitd(generatefinalbam.out)
		//coverage_mosdepth(generatefinalbam.out)
		//hsmetrics_run(generatefinalbam.out)
		//freebayes(generatefinalbam.out)
		//haplotypecaller(generatefinalbam.out)
		//deepvariant(generatefinalbam.out)
		//strelka(generatefinalbam.out)
		//platypus(generatefinalbam.out)
		//varscan(generatefinalbam.out)
		//lofreq(generatefinalbam.out)
		//pindel(generatefinalbam.out)
		//somaticSeq_run(lofreq.out.join(varscan.out.join(platypus.out.join(strelka.out.join(haplotypecaller.out.join(freebayes.out.join(deepvariant.out.join(generatefinalbam.out))))))))
		//cnvnator(generatefinalbam.out)
		//cava(somaticSeq_run.out)
		//merge_csv(pindel.out.join(somaticSeq_run.out.join(cava.out)))
}

workflow WES_BAMIN {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.map{ it }
		.set { samples_ch }
	main:
		generatefinalbamin(samples_ch)
		minimap_getitd(generatefinalbamin.out)
		coverage_mosdepth(generatefinalbamin.out)
		hsmetrics_run(generatefinalbamin.out)
		freebayes(generatefinalbamin.out)
		haplotypecaller(generatefinalbamin.out)	
		deepvariant(generatefinalbamin.out)
		strelka(generatefinalbamin.out)
		platypus(generatefinalbamin.out)
		varscan(generatefinalbamin.out)
		lofreq(generatefinalbamin.out)
		pindel(generatefinalbamin.out)
		ifcnv_run(generatefinalbamin.out.collect())
		somaticSeq_run(lofreq.out.join(varscan.out.join(platypus.out.join(strelka.out.join(haplotypecaller.out.join(freebayes.out.join(generatefinalbamin.out)))))))
		GermlineCNVAnalysis(generatefinalbamin.out)		
		deepvariant_gCNV(deepvariant.out.join(GermlineCNVAnalysis.out))
		cava(somaticSeq_run.out)
		merge_csv(pindel.out.join(somaticSeq_run.out.join(deepvariant.out.join(cava.out))))
}

workflow WES_NOPAIR {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.map{ it }
		.set { samples_ch }

	main:
		fastqc(samples_ch)
		trimming_trimmomatic(samples_ch) | mapping_both_reads | sam_conversion | mark_duplicates
		RealignerTargetCreator(mark_duplicates.out)
		IndelRealigner(RealignerTargetCreator.out.join(mark_duplicates.out)) | BaseRecalibrator
		PrintReads(IndelRealigner.out.join(BaseRecalibrator.out)) | generatefinalbam
		minimap_getitd(generatefinalbam.out)
		coverage_mosdepth(generatefinalbam.out)
		hsmetrics_run(generatefinalbam.out)
		freebayes(generatefinalbam.out)
		haplotypecaller(generatefinalbam.out)
		deepvariant(generatefinalbam.out)
		strelka(generatefinalbam.out)
		platypus(generatefinalbam.out)
		varscan(generatefinalbam.out)
		lofreq(generatefinalbam.out)
		pindel(generatefinalbam.out)
		somaticSeq_run(lofreq.out.join(varscan.out.join(platypus.out.join(strelka.out.join(haplotypecaller.out.join(freebayes.out.join(deepvariant.out.join(generatefinalbam.out))))))))
		cnvnator(generatefinalbam.out)
		cava(somaticSeq_run.out)
		merge_csv(pindel.out.join(somaticSeq_run.out.join(cava.out)))
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
}
