#!/usr/bin/env nextflow
nextflow.enable.dsl=2

bedfile = file("${params.bedfile}", checkIfExists: true )
bedfile_zipped = file("${params.bedfile}.*")
genome_loc = file("${params.genome}", checkIfExists: true)
index_files = file("${params.genome_dir}/${params.ind_files}.*")
known_SNPs = file("${params.site2}", checkIfExists: true)
known_SNPs_index = file("${params.site2_idx}", checkIfExists: true)
deepvariant = params.deepvariant
somaticseq = params.somaticseq
dbsnp_somatic = file("${params.dbsnp_somatic}", checkIfExists: true)
cava_config = file("${params.cava_config}", checkIfExists: true)
ensembl_db = file("${params.ensembl_db}", checkIfExists: true)
ensembl_db_index = file("${params.ensembl_db_index}", checkIfExists: true)
known_SNPs_zip = file("${params.site2_zip}", checkIfExists: true)
known_SNPs_tbi = file("${params.site2_zip_idx}", checkIfExists: true)
pharma_input = file("${params.pharma_input}", checkIfExists: true)
preprocessed_intervals = file("${params.preprocessed_intervals}", checkIfExists: true)
ploidy_model = file("${params.ploidy_model}", checkIfExists: true)
cohort_model = file("${params.cohort_model}", checkIfExists: true)

// Adapter Trimming, alignment and GATK BQSR (based on https://github.com/GavinHaLab/fastq_to_bam_paired_snakemake)
include { FASTQTOBAM  } from '../subworkflows/fastq_to_bam.nf'
include { GCNV } from '../subworkflows/gcnv.nf'

// Adapter Trimming, alignment and GATK BQSR (based on https://github.com/GavinHaLab/fastq_to_bam_paired_snakemake)
include { ABRA_BAM } from '../modules/abra/realign/main.nf'
include { ABRA_SORT } from '../modules/samtools/abra_sort/main.nf'
include { COVERAGE_MOSDEPTH } from '../modules/mosdepth/coverage/main.nf'
include { EXTRACT_COV50 } from '../modules/python/extract_mosdepth/main.nf'
include { HSMETRICS ; HSMETRICS_COLLECT } from '../modules/gatk/hsmetrics/main.nf'

// Variant calling
include { FREEBAYES } from '../modules/variant_call/freebayes/main.nf'
include { HAPLOTYPECALLER } from '../modules/variant_call/haplotypecaller/main.nf'
include { STRELKA } from '../modules/variant_call/strelka/main.nf'
include { PLATYPUS } from '../modules/variant_call/platypus/main.nf'
include { VARSCAN } from '../modules/variant_call/varscan/main.nf'
include { LOFREQ } from '../modules/variant_call/lofreq/main.nf'
include { DEEPVARIANT } from '../modules/variant_call/deepvariant/main.nf'
include { PINDEL_FLT3 } from '../modules/variant_call/pindel/pindel_flt3/main.nf'

// Annotation
include { VEP as VEP_DEEPVARIANT ; VEP_ANALYSIS as VEP_PINDEL ; VEP_ANALYSIS as VEP_SOMATICSEQ } from '../modules/vep/annotate/main.nf'
include { ANNOVAR as ANNOVAR_DEEPVARIANT ; ANNOVAR as ANNOVAR_SOMATICSEQ } from '../modules/annovar/annotate/main.nf'
include { CAVA } from '../modules/cava/annotate/main.nf'

// Variant integration 
include { FORMAT_DEEPVARIANT } from '../modules/python/format_deepvariant/main.nf'
include { FORMAT_PINDEL } from '../modules/python/format_pindel/main.nf'
include { FORMAT_SOMATICSEQ } from '../modules/python/format_somaticseq/main.nf'
include { FORMAT_CAVA } from '../modules/python/format_cava/main.nf'
include { SOMATICSEQ } from '../modules/variant_integration/somaticseq/main.nf'
include { SOMATICSEQ_CONCAT } from '../modules/bcftools/concat/main.nf'

// CNV calling
include { IFCNV } from '../modules/cnv_call/ifcnv/main.nf'

// Format output
include { MERGE_DEEPVAR_GCNV } from '../modules/python/merge_deepvar_gcnv/main.nf'
include { MERGE_CSV_WES } from '../modules/python/merge_csv/main.nf'


workflow WES {
	Channel.fromPath(params.input)
		.splitCsv(header:false)
		.map { row ->
			def sample = row[0].trim()
			def r1 = file("${params.sequences}/${sample}_S*_R1_*.fastq.gz", checkIfExists: false)
			def r2 = file("${params.sequences}/${sample}_S*_R2_*.fastq.gz", checkIfExists: false)

			if (!r1 && !r2) {
				r1 = file("${params.sequences}/${sample}*_R1.fastq.gz", checkIfExists: false)
				r2 = file("${params.sequences}/${sample}*_R2.fastq.gz", checkIfExists: false)
			}
			tuple(sample, r1, r2)
		}
		.set { bam_ch }
	main:
		wes_bam_ch = FASTQTOBAM(bam_ch)
		ABRA_BAM(wes_bam_ch.final_bams_ch, bedfile, genome_loc, index_files )
		ABRA_SORT(wes_bam_ch.final_bams_ch.join(ABRA_BAM.out))
		COVERAGE_MOSDEPTH(ABRA_SORT.out.final_bam, bedfile)
		EXTRACT_COV50(COVERAGE_MOSDEPTH.out)
		HSMETRICS(ABRA_SORT.out.final_bam, bedfile, genome_loc, index_files)
		all_hsmetrics = HSMETRICS.out.collect()
		HSMETRICS_COLLECT(all_hsmetrics)
		FREEBAYES(ABRA_SORT.out.final_bam, genome_loc, index_files, bedfile)
		HAPLOTYPECALLER(ABRA_SORT.out.final_bam, genome_loc, index_files, bedfile, known_SNPs, known_SNPs_index)
		STRELKA(ABRA_SORT.out.final_bam, genome_loc, index_files, bedfile, bedfile_zipped)
		PLATYPUS(ABRA_SORT.out.final_bam, genome_loc, index_files, bedfile)
		VARSCAN(ABRA_SORT.out.final_bam, genome_loc, index_files, bedfile)
		DEEPVARIANT(ABRA_SORT.out.final_bam, genome_loc, index_files, bedfile)
		LOFREQ(ABRA_SORT.out.old_bam, genome_loc, index_files, bedfile)
		VEP_DEEPVARIANT(DEEPVARIANT.out)
		ANNOVAR_DEEPVARIANT(DEEPVARIANT.out, deepvariant)
		FORMAT_DEEPVARIANT(DEEPVARIANT.out.join(VEP_DEEPVARIANT.out.join(ANNOVAR_DEEPVARIANT.out)))
		PINDEL_FLT3(ABRA_SORT.out.final_bam, genome_loc, index_files)
		VEP_PINDEL(PINDEL_FLT3.out)
		FORMAT_PINDEL(PINDEL_FLT3.out.join(VEP_PINDEL.out))
		SOMATICSEQ(ABRA_SORT.out.final_bam.join(LOFREQ.out.join(VARSCAN.out.join(PLATYPUS.out.join(STRELKA.out.join(HAPLOTYPECALLER.out.join(FREEBAYES.out)))))), genome_loc, index_files, bedfile, dbsnp_somatic)
		SOMATICSEQ_CONCAT(SOMATICSEQ.out)
		VEP_SOMATICSEQ(SOMATICSEQ_CONCAT.out)
		ANNOVAR_SOMATICSEQ(SOMATICSEQ_CONCAT.out, somaticseq)
		FORMAT_SOMATICSEQ(SOMATICSEQ_CONCAT.out.join(VEP_SOMATICSEQ.out.join(ANNOVAR_SOMATICSEQ.out)))
		IFCNV(ABRA_SORT.out.final_bam, bedfile)
		CAVA(SOMATICSEQ_CONCAT.out,cava_config, genome_loc, index_files, bedfile_zipped, known_SNPs_zip, known_SNPs_tbi, ensembl_db, ensembl_db_index)
		FORMAT_CAVA(CAVA.out)
		GCNV(ABRA_SORT.out.final_bam, genome_loc, index_files, preprocessed_intervals, ploidy_model, cohort_model)
		MERGE_DEEPVAR_GCNV(DEEPVARIANT.out.join(GCNV.out))
		MERGE_CSV_WES(FORMAT_PINDEL.out.join(FORMAT_SOMATICSEQ.out.join(FORMAT_DEEPVARIANT.out.join(CAVA.out.join(COVERAGE_MOSDEPTH.out.join(EXTRACT_COV50.out))))), pharma_input)
}

