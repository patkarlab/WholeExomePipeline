

include { TRIM; MAPBAM; MARK_DUPS; BQSR; APPLY_BQSR; ALIGNMENT_METRICS; INSERT_SIZE_METRICS } from '../modules/fastq_to_bam/main.nf'


adapters = file("${params.adapters}", checkIfExists: true )
genome_loc = file("${params.genome}", checkIfExists: true)
index_files = file("${params.genome_dir}/${params.ind_files}.*")
known_SNPs = file("${params.site2}", checkIfExists: true)
known_SNPs_index = file("${params.site2_idx}", checkIfExists: true)
known_INDELS = file("${params.site1}", checkIfExists: true)
known_INDELS_index = file("${params.site1_idx}", checkIfExists: true)

workflow FASTQTOBAM {
	take:
		samples_ch
	main:
	TRIM(samples_ch, adapters)
	MAPBAM(TRIM.out, genome_loc, index_files)
	MARK_DUPS(MAPBAM.out)
	BQSR(MARK_DUPS.out, genome_loc, index_files, known_SNPs, known_SNPs_index, known_INDELS, known_INDELS_index)
	APPLY_BQSR(MARK_DUPS.out.join(BQSR.out), genome_loc, index_files)
	ALIGNMENT_METRICS(APPLY_BQSR.out, genome_loc, index_files)
	INSERT_SIZE_METRICS(APPLY_BQSR.out)

	emit:
		final_bams_ch = APPLY_BQSR.out
		trimmed_fastq = TRIM.out
}
