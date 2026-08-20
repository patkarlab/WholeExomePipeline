process SOMATICSEQ {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val (Sample), file(bam), file(bamBai), file (lofreqVcf), file (varscanVcf), file (platypusVcf), file (strelkaVcf), file (haplotypecallerVcf), file (freebayesVcf)
		path (GenFile)
		path (GenDir)
		path (bedfile)
		path (known_SNPs)
		path (known_SNPs_index)
	output:
		tuple val (Sample), file("${Sample}_somaticseq_snv.vcf"), file("${Sample}_somaticseq_indel.vcf")
	script:
	"""
	vcf_sorter.sh ${freebayesVcf} ${Sample}.freebayes.sorted.vcf
	vcf_sorter.sh ${platypusVcf} ${Sample}.platypus.sorted.vcf
	vcf_sorter.sh ${haplotypecallerVcf} ${Sample}.haplotypecaller.sorted.vcf

	split_vcf.py -infile ${Sample}.platypus.sorted.vcf -snv ${Sample}_platypus_snvs.vcf -indel ${Sample}_platypus_indels.vcf -genome ${GenFile}
	split_vcf.py -infile ${Sample}.freebayes.sorted.vcf -snv ${Sample}_freebayes_snvs.vcf -indel ${Sample}_freebayes_indels.vcf -genome ${GenFile}
	split_vcf.py -infile ${Sample}.haplotypecaller.sorted.vcf -snv ${Sample}_haplotypecaller_snvs.vcf -indel ${Sample}_haplotypecaller_indels.vcf -genome ${GenFile}

	vcf_sorter.sh ${Sample}_platypus_snvs.vcf ${Sample}_platypus_snvs_sort.vcf
	vcf_sorter.sh ${Sample}_platypus_indels.vcf ${Sample}_platypus_indels_sort.vcf
	vcf_sorter.sh ${Sample}_freebayes_snvs.vcf ${Sample}_freebayes_snvs_sort.vcf
	vcf_sorter.sh ${Sample}_freebayes_indels.vcf ${Sample}_freebayes_indels_sort.vcf
	vcf_sorter.sh ${Sample}_haplotypecaller_snvs.vcf ${Sample}_haplotypecaller_snvs_sort.vcf
	vcf_sorter.sh ${Sample}_haplotypecaller_indels.vcf ${Sample}_haplotypecaller_indels_sort.vcf

	somaticseq_parallel.py \
	--output-directory ./${Sample}.somaticseq \
	--genome-reference ${GenFile} \
	--inclusion-region ${bedfile} \
	--threads ${task.cpus} \
	--pass-threshold 0 \
	--lowqual-threshold 0 \
	--algorithm xgboost  \
	-minMQ 0 \
	-minBQ 0 \
	-mincaller 0 \
	--dbsnp-vcf  ${known_SNPs} \
	single \
	--sample-name ${Sample} \
	--bam-file ${bam} \
	--varscan-vcf ${varscanVcf} \
	--lofreq-vcf ${lofreqVcf} \
	--strelka-vcf ${strelkaVcf}  \
	--arbitrary-snvs ${Sample}_freebayes_snvs_sort.vcf ${Sample}_platypus_snvs_sort.vcf ${Sample}_haplotypecaller_snvs_sort.vcf \
	--arbitrary-indels ${Sample}_freebayes_indels_sort.vcf ${Sample}_platypus_indels_sort.vcf ${Sample}_haplotypecaller_indels_sort.vcf
	
	vcf_sorter.sh ./${Sample}.somaticseq/Consensus.sSNV.vcf ./${Sample}_somaticseq_snv.vcf
	vcf_sorter.sh ./${Sample}.somaticseq/Consensus.sINDEL.vcf ./${Sample}_somaticseq_indel.vcf
	
	"""
	stub:
	"""
	touch ${Sample}_somaticseq_snv.vcf ${Sample}_somaticseq_indel.vcf
	"""
}

process SOMATICSEQ_LYMPHOMA {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val (Sample), file(bam), file(bamBai), file (lofreqVcf), file (varscanVcf), file (platypusVcf), file (strelkaVcf), file (haplotypecallerVcf), file (freebayesVcf), file(DeepSomaticVcf)
		path (GenFile)
		path (GenDir)
		path (bedfile)
                path (known_SNPs)
                path (known_SNPs_index)
	output:
		tuple val (Sample), file("${Sample}_somaticseq_snv.vcf"), file("${Sample}_somaticseq_indel.vcf")
	script:
	"""
	vcf_sorter.sh ${freebayesVcf} ${Sample}.freebayes.sorted.vcf
	vcf_sorter.sh ${platypusVcf} ${Sample}.platypus.sorted.vcf
	vcf_sorter.sh ${haplotypecallerVcf} ${Sample}.haplotypecaller.sorted.vcf
	vcf_sorter.sh ${DeepSomaticVcf} ${Sample}.deepsomatic.sorted.vcf

	split_vcf.py -infile ${Sample}.platypus.sorted.vcf -snv ${Sample}_platypus_snvs.vcf -indel ${Sample}_platypus_indels.vcf -genome ${GenFile}
	split_vcf.py -infile ${Sample}.freebayes.sorted.vcf -snv ${Sample}_freebayes_snvs.vcf -indel ${Sample}_freebayes_indels.vcf -genome ${GenFile}
	split_vcf.py -infile ${Sample}.haplotypecaller.sorted.vcf -snv ${Sample}_haplotypecaller_snvs.vcf -indel ${Sample}_haplotypecaller_indels.vcf -genome ${GenFile}
	split_vcf.py -infile ${Sample}.deepsomatic.sorted.vcf -snv ${Sample}_deepsomatic_snvs.vcf -indel ${Sample}_deepsomatic_indels.vcf -genome ${GenFile}

	vcf_sorter.sh ${Sample}_platypus_snvs.vcf ${Sample}_platypus_snvs_sort.vcf
	vcf_sorter.sh ${Sample}_platypus_indels.vcf ${Sample}_platypus_indels_sort.vcf
	vcf_sorter.sh ${Sample}_freebayes_snvs.vcf ${Sample}_freebayes_snvs_sort.vcf
	vcf_sorter.sh ${Sample}_freebayes_indels.vcf ${Sample}_freebayes_indels_sort.vcf
	vcf_sorter.sh ${Sample}_haplotypecaller_snvs.vcf ${Sample}_haplotypecaller_snvs_sort.vcf
	vcf_sorter.sh ${Sample}_haplotypecaller_indels.vcf ${Sample}_haplotypecaller_indels_sort.vcf
	vcf_sorter.sh ${Sample}_deepsomatic_snvs.vcf ${Sample}_deepsomatic_snvs_sort.vcf
	vcf_sorter.sh ${Sample}_deepsomatic_indels.vcf ${Sample}_deepsomatic_indels_sort.vcf	
	vcf_sorter.sh ${varscanVcf} ${Sample}_varscan_sort.vcf

	somaticseq_parallel.py \
	--output-directory ./${Sample}.somaticseq \
	--genome-reference ${GenFile} \
	--inclusion-region ${bedfile} \
	--threads ${task.cpus} \
	--pass-threshold 0 \
	--lowqual-threshold 0 \
	--algorithm xgboost  \
	-minMQ 0 \
	-minBQ 0 \
	-mincaller 0 \
	--dbsnp-vcf  ${known_SNPs} \
	single \
	--sample-name ${Sample} \
	--bam-file ${bam} \
	--varscan-vcf ${Sample}_varscan_sort.vcf \
	--lofreq-vcf ${lofreqVcf} \
	--strelka-vcf ${strelkaVcf}  \
	--arbitrary-snvs ${Sample}_freebayes_snvs_sort.vcf ${Sample}_platypus_snvs_sort.vcf ${Sample}_haplotypecaller_snvs_sort.vcf ${Sample}_deepsomatic_snvs_sort.vcf \
	--arbitrary-indels ${Sample}_freebayes_indels_sort.vcf ${Sample}_platypus_indels_sort.vcf ${Sample}_haplotypecaller_indels_sort.vcf ${Sample}_deepsomatic_indels_sort.vcf
	
	vcf_sorter.sh ./${Sample}.somaticseq/Consensus.sSNV.vcf ./${Sample}_somaticseq_snv.vcf
	vcf_sorter.sh ./${Sample}.somaticseq/Consensus.sINDEL.vcf ./${Sample}_somaticseq_indel.vcf

	"""
	stub:
	"""
	touch ${Sample}_somaticseq_snv.vcf ${Sample}_somaticseq_indel.vcf
	"""
}
