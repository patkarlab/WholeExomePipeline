process FORMAT_SV_CALLS {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/translocatn", mode: 'copy', pattern: '*.tsv'	
	input:
		tuple val (Sample), file (svaba_somatic_sv_vcf), file (svaba_germline_sv_vcf), file (svaba_somatic_indel_vcf), file (svaba_germline_indel_vcf), file (lumpy_vcf), file (gridss_vcf), file (gridss_somatic_vcf_bgz), file (rm_gridss_somatic_vcf), file (delly_vcf), file (delly_somatic_vcf), file (manta_candidateSV_vcf), file (manta_diploidSV_vcf), file (manta_somaticSV_vcf)
	output:
		tuple val (Sample), file ("${Sample}_svaba_lumpy.tsv"), file ("${Sample}.svaba.somatic.sv.tsv"), file ("${Sample}.svaba.germline.sv.tsv"), file ("${Sample}.svaba.somatic.indel.tsv"), file ("${Sample}.svaba.germline.indel.tsv"), file("${Sample}.gridss.tsv"), file("${Sample}.gridss.somatic.tsv"), file("RM_${Sample}.gridss.somatic.tsv"), file ("${Sample}_delly.tsv"), file ("${Sample}_delly_somatic.tsv"), file ("${Sample}_manta_somaticSV.tsv")
	script:
	"""
	svaba_lumpy_common.py ${svaba_somatic_sv_vcf} ${lumpy_vcf} ${Sample}_common
	svaba_file=\$(basename ${svaba_somatic_sv_vcf} .vcf)
	lumpy_file=\$(basename ${lumpy_vcf} .vcf)
	cat \${svaba_file}.tsv \${lumpy_file}.tsv ${Sample}_common > ${Sample}_svaba_lumpy.tsv

	format_csv.py ${svaba_somatic_sv_vcf} ${Sample}.svaba.somatic.sv.tsv
	format_csv.py ${svaba_germline_sv_vcf} ${Sample}.svaba.germline.sv.tsv
	format_csv.py ${svaba_somatic_indel_vcf} ${Sample}.svaba.somatic.indel.tsv
	format_csv.py ${svaba_germline_indel_vcf} ${Sample}.svaba.germline.indel.tsv
	format_csv.py ${gridss_vcf} ${Sample}.gridss.tsv
	bgzip -d ${gridss_somatic_vcf_bgz}
	format_csv.py ${Sample}.gridss.somatic.vcf ${Sample}.gridss.somatic.tsv
	format_csv.py ${rm_gridss_somatic_vcf} RM_${Sample}.gridss.somatic.tsv
	format_csv.py ${delly_vcf} ${Sample}_delly.tsv
	format_csv.py ${delly_somatic_vcf} ${Sample}_delly_somatic.tsv
	format_manta_csv.py ${manta_somaticSV_vcf} ${Sample}_manta_somaticSV.tsv
	"""
	stub:
	"""
	touch \
		${Sample}_svaba_lumpy.tsv \
		${Sample}.svaba.somatic.sv.tsv \
		${Sample}.svaba.germline.sv.tsv \
		${Sample}.svaba.somatic.indel.tsv \
		${Sample}.svaba.germline.indel.tsv \
		${Sample}.gridss.tsv \
		${Sample}.gridss.somatic.tsv \
		RM_${Sample}.gridss.somatic.tsv \
		${Sample}_delly.tsv \
		${Sample}_delly_somatic.tsv \
		${Sample}_manta_somaticSV.tsv
	"""
}
