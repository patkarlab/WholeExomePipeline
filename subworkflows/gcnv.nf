
include { COLLECT_READ_COUNTS } from '../modules/cnv_call/gcnv/collect_read_counts/main.nf'
include { DETERMINE_PLOIDY }    from '../modules/cnv_call/gcnv/determine_ploidy/main.nf'
include { CALL_GCNV } from '../modules/cnv_call/gcnv/call_gcnv/main.nf'
include { POSTPROCESS_GCNV }    from '../modules/cnv_call/gcnv/postprocess_gcnv/main.nf'

workflow GCNV {
	take:
		bam_ch
		genome_loc
		index_files
		intervals
		ploidy_model
		cohort_model
	main:
	COLLECT_READ_COUNTS(bam_ch, genome_loc, index_files, intervals)
	DETERMINE_PLOIDY(COLLECT_READ_COUNTS.out, ploidy_model)
	CALL_GCNV(DETERMINE_PLOIDY.out, cohort_model)
	POSTPROCESS_GCNV(CALL_GCNV.out, cohort_model )

	emit:
		gcnv_results = POSTPROCESS_GCNV.out
}