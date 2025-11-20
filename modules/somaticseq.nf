#!/usr/bin/env nextflow

process SOMATICSEQ {
	label 'process_low'
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*_somaticseq.vep_annonvar.txt'
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '${Sample}.somaticseq.vcf'
	input:
		tuple val (Sample), file (lofreqVcf), file (varscanVcf), file (platypusVcf), file (strelkaVcf), file (haplotypecallerVcf), file (freebayesVcf), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
		path (genome_dir)
		path (genome_fasta)
	output:
		tuple val (Sample), file ("*_somaticseq.vep_annonvar.txt"), file("${Sample}.somaticseq.vcf")
	script:
	"""
	vcf_sorter.sh ${freebayesVcf} ${Sample}.freebayes.sorted.vcf
	vcf_sorter.sh ${platypusVcf} ${Sample}.platypus.sorted.vcf
	vcf_sorter.sh ${haplotypecallerVcf} ${Sample}.haplotypecaller.sorted.vcf

	python3 ${params.splitvcf_path} -infile ${Sample}.platypus.sorted.vcf -snv ${Sample}_platypus_cnvs.vcf -indel ${Sample}_platypus_indels.vcf
	python3 ${params.splitvcf_path} -infile ${Sample}.freebayes.sorted.vcf -snv ${Sample}_freebayes_cnvs.vcf -indel ${Sample}_freebayes_indels.vcf
	python3 ${params.splitvcf_path} -infile ${Sample}.haplotypecaller.sorted.vcf -snv ${Sample}_haplotypecaller_cnvs.vcf -indel ${Sample}_haplotypecaller_indels.vcf

	vcf_sorter.sh ${Sample}_platypus_cnvs.vcf ${Sample}_platypus_cnvs_sort.vcf
	vcf_sorter.sh ${Sample}_platypus_indels.vcf ${Sample}_platypus_indels_sort.vcf
	vcf_sorter.sh ${Sample}_freebayes_cnvs.vcf ${Sample}_freebayes_cnvs_sort.vcf
	vcf_sorter.sh ${Sample}_freebayes_indels.vcf ${Sample}_freebayes_indels_sort.vcf
	vcf_sorter.sh ${Sample}_haplotypecaller_cnvs.vcf ${Sample}_haplotypecaller_cnvs_sort.vcf
	vcf_sorter.sh ${Sample}_haplotypecaller_indels.vcf ${Sample}_haplotypecaller_indels_sort.vcf

	somaticseq_parallel.py --output-directory ${Sample}.somaticseq --genome-reference ${genome_dir}/${genome_fasta} --inclusion-region ${params.bedfile}.bed --threads ${task.cpus} --pass-threshold 0 --lowqual-threshold 0 --algorithm xgboost -minMQ 0 -minBQ 0 -mincaller 0 --dbsnp-vcf /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file ${finalBam} --varscan-vcf ${varscanVcf} --lofreq-vcf ${lofreqVcf} --strelka-vcf ${strelkaVcf} --sample-name ${Sample} --arbitrary-snvs ${Sample}_freebayes_cnvs_sort.vcf ${Sample}_platypus_cnvs_sort.vcf ${Sample}_haplotypecaller_cnvs_sort.vcf --arbitrary-indels ${Sample}_freebayes_indels_sort.vcf ${Sample}_platypus_indels_sort.vcf ${Sample}_haplotypecaller_indels_sort.vcf

	vcf_sorter.sh ${Sample}.somaticseq/Consensus.sSNV.vcf ${Sample}.somaticseq/somaticseq_snv.vcf
	bgzip -c ${Sample}.somaticseq/somaticseq_snv.vcf > ${Sample}.somaticseq/somaticseq_snv.vcf.gz
	${params.bcftools_path} index -t ${Sample}.somaticseq/somaticseq_snv.vcf.gz

	vcf_sorter.sh ${Sample}.somaticseq/Consensus.sINDEL.vcf ${Sample}.somaticseq/somaticseq_indel.vcf
	bgzip -c ${Sample}.somaticseq/somaticseq_indel.vcf > ${Sample}.somaticseq/somaticseq_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.somaticseq/somaticseq_indel.vcf.gz

	${params.bcftools_path} concat -a ${Sample}.somaticseq/somaticseq_snv.vcf.gz ${Sample}.somaticseq/somaticseq_indel.vcf.gz -o ${Sample}.somaticseq.vcf

	sed -i 's/##INFO=<ID=VLK012,Number=6,Type=Integer,Description="Calling decision of the 6 algorithms: VarScan2, LoFreq, Strelka, SnvCaller_0, SnvCaller_1, SnvCaller_2">/##INFO=<ID=VLSFPH,Number=6,Type=String,Description="Calling decision of the 6 algorithms:  VarScan2, LoFreq, Strelka, Freebayes, Platypus, Haplotypecaller">/g' ${Sample}.somaticseq.vcf
	sed -i 's/VLK012/VLSFPH/g' ${Sample}.somaticseq.vcf
	
	# to extract vaf,af,alt and ref count
	extract_somatic.py ${Sample}.somaticseq.vcf ${Sample}.extractedSomaticseq.txt
	#adding vep
	vep_analysis.sh ${Sample}.somaticseq.vcf ${Sample}
	extract_velheader.py ${Sample}_vep_delheaders.txt ${Sample}.extractedvepdelheaders.txt
	
	# for merging extracted somaticsseq and vepheaders
	mergeSomaticVep.py ${Sample}.extractedSomaticseq.txt ${Sample}.extractedvepdelheaders.txt ${Sample}_somaticseq.vep.txt
	sed -i 's/SYMBOL/Gene/g' ${Sample}_somaticseq.vep.txt
	
	#Annotating ${Sample}.somaticseq.vcf using annovar
	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}.somaticseq.vcf --outfile ${Sample}.somaticseq.avinput -withfreq --includeinfo -allsample
	
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.somaticseq.avinput --out ${Sample}.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread ${task.cpus} ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt
	
	#extracting columns  Func.refGene,Gene.refGene,ExonicFunc.refGene,PopFreqMax,InterVar_automated from somaticseq.hg19_multianno.csv and adding them to somaticseq.vep.txt
	extract_annovar.py ${Sample}.somaticseq.hg19_multianno.csv ${Sample}_somaticseq.vep.txt ${Sample}_somaticseq.vep_annonvar.txt
	"""
} 

process SOMATICSEQ_LYMPHOMA {
	label 'process_low'
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*.somaticseq.csv'
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '${Sample}.somaticseq.vcf'
	input:
		tuple val (Sample), file (lofreqVcf), file (varscanVcf), file (platypusVcf), file (strelkaVcf), file (haplotypecallerVcf), file (freebayesVcf), file(DeepSomaticVcf), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
		path (genome_dir)
		path (genome_fasta)
	output:
		tuple val (Sample), file ("*.somaticseq.csv"), file("${Sample}.somaticseq.vcf")
	script:
	"""
	vcf_sorter.sh ${freebayesVcf} ${Sample}.freebayes.sorted.vcf
	vcf_sorter.sh ${platypusVcf} ${Sample}.platypus.sorted.vcf
	vcf_sorter.sh ${haplotypecallerVcf} ${Sample}.haplotypecaller.sorted.vcf
	vcf_sorter.sh ${DeepSomaticVcf} ${Sample}.deepsomatic.sorted.vcf

	python3 ${params.splitvcf_path} -infile ${Sample}.platypus.sorted.vcf -snv ${Sample}_platypus_cnvs.vcf -indel ${Sample}_platypus_indels.vcf
	python3 ${params.splitvcf_path} -infile ${Sample}.freebayes.sorted.vcf -snv ${Sample}_freebayes_cnvs.vcf -indel ${Sample}_freebayes_indels.vcf
	python3 ${params.splitvcf_path} -infile ${Sample}.haplotypecaller.sorted.vcf -snv ${Sample}_haplotypecaller_cnvs.vcf -indel ${Sample}_haplotypecaller_indels.vcf
	python3 ${params.splitvcf_path} -infile ${Sample}.deepsomatic.sorted.vcf -snv ${Sample}_deepsomatic_snvs.vcf -indel ${Sample}_deepsomatic_indels.vcf

	vcf_sorter.sh ${Sample}_platypus_cnvs.vcf ${Sample}_platypus_cnvs_sort.vcf
	vcf_sorter.sh ${Sample}_platypus_indels.vcf ${Sample}_platypus_indels_sort.vcf
	vcf_sorter.sh ${Sample}_freebayes_cnvs.vcf ${Sample}_freebayes_cnvs_sort.vcf
	vcf_sorter.sh ${Sample}_freebayes_indels.vcf ${Sample}_freebayes_indels_sort.vcf
	vcf_sorter.sh ${Sample}_haplotypecaller_cnvs.vcf ${Sample}_haplotypecaller_cnvs_sort.vcf
	vcf_sorter.sh ${Sample}_haplotypecaller_indels.vcf ${Sample}_haplotypecaller_indels_sort.vcf
	vcf_sorter.sh ${Sample}_deepsomatic_snvs.vcf ${Sample}_deepsomatic_snvs_sort.vcf
	vcf_sorter.sh ${Sample}_deepsomatic_indels.vcf ${Sample}_deepsomatic_indels_sort.vcf	

	somaticseq_parallel.py --output-directory ${Sample}.somaticseq --genome-reference ${genome_dir}/${genome_fasta} --inclusion-region ${params.bedfile}.bed --threads ${task.cpus} --pass-threshold 0 --lowqual-threshold 0 --algorithm xgboost -minMQ 0 -minBQ 0 -mincaller 0 --dbsnp-vcf /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file ${finalBam} --varscan-vcf ${varscanVcf} --lofreq-vcf ${lofreqVcf} --strelka-vcf ${strelkaVcf} --sample-name ${Sample} --arbitrary-snvs ${Sample}_freebayes_cnvs_sort.vcf ${Sample}_platypus_cnvs_sort.vcf ${Sample}_haplotypecaller_cnvs_sort.vcf ${Sample}_deepsomatic_snvs_sort.vcf --arbitrary-indels ${Sample}_freebayes_indels_sort.vcf ${Sample}_platypus_indels_sort.vcf ${Sample}_haplotypecaller_indels_sort.vcf ${Sample}_deepsomatic_indels_sort.vcf

	vcf_sorter.sh ${Sample}.somaticseq/Consensus.sSNV.vcf ${Sample}.somaticseq/somaticseq_snv.vcf
	bgzip -c ${Sample}.somaticseq/somaticseq_snv.vcf > ${Sample}.somaticseq/somaticseq_snv.vcf.gz
	${params.bcftools_path} index -t ${Sample}.somaticseq/somaticseq_snv.vcf.gz

	vcf_sorter.sh ${Sample}.somaticseq/Consensus.sINDEL.vcf ${Sample}.somaticseq/somaticseq_indel.vcf
	bgzip -c ${Sample}.somaticseq/somaticseq_indel.vcf > ${Sample}.somaticseq/somaticseq_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.somaticseq/somaticseq_indel.vcf.gz

	${params.bcftools_path} concat -a ${Sample}.somaticseq/somaticseq_snv.vcf.gz ${Sample}.somaticseq/somaticseq_indel.vcf.gz -o ${Sample}.somaticseq.vcf

	sed -i 's/##INFO=<ID=VLK0123,Number=7,Type=Integer,Description="Calling decision of the 7 algorithms: VarScan2, LoFreq, Strelka, SnvCaller_0, SnvCaller_1, SnvCaller_2, SnvCaller_2">/##INFO=<ID=VLSFPHD,Number=7,Type=String,Description="Calling decision of the 7 algorithms:  VarScan2, LoFreq, Strelka, Freebayes, Platypus, Haplotypecaller, Deepsomatic">/g' ${Sample}.somaticseq.vcf
	sed -i 's/VLK0123/VLSFPHD/g' ${Sample}.somaticseq.vcf
	
	
	#Annotating ${Sample}.somaticseq.vcf using annovar
	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}.somaticseq.vcf --outfile ${Sample}.somaticseq.avinput -withfreq --includeinfo -allsample
	
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.somaticseq.avinput --out ${Sample}.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread ${task.cpus} ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt
	
	somaticseqoutput-format.py ${Sample}.somaticseq.hg19_multianno.csv ${Sample}.somaticseq.csv

	"""
} 