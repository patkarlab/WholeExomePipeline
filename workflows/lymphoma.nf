#!/usr/bin/env nextflow
nextflow.enable.dsl=2


bedfile = file("${params.bedfile}", checkIfExists: true )
bedfile_zipped = file("${params.bedfile}.*")
genome_loc = file("${params.genome}", checkIfExists: true)
index_files = file("${params.genome_dir}/${params.ind_files}.*")
known_SNPs = file("${params.site2}", checkIfExists: true)
known_SNPs_index = file("${params.site2_idx}", checkIfExists: true)
coverview_config = file("${params.coverview_config}", checkIfExists: true )
control_bqsr_bam = file("${params.control_bqsr_bam}", checkIfExists: true )
control_bqsr_bamBai = file("${params.control_bqsr_bam}.bai", checkIfExists: true )
control_final_bam = file("${params.control_final_bam}", checkIfExists: true )
control_final_bamBai = file("${params.control_final_bam}.bai", checkIfExists: true )
dbsnp_somatic = file("${params.dbsnp_somatic}", checkIfExists: true)
somaticseq = params.somaticseq
gridss_exclude_list = file("${params.gridss_exclude_list}", checkIfExists: true)
cnvkit_reference = file("${params.cnvkit_reference}", checkIfExists: true)
cnvkit_reference_delIG = file("${params.cnvkit_reference_delIG}", checkIfExists: true)
gene_scatter_list = file("${params.gene_scatter_list}", checkIfExists: true)
cava_config = file("${params.cava_config}", checkIfExists: true)
ensembl_db = file("${params.ensembl_db}", checkIfExists: true)
ensembl_db_index = file("${params.ensembl_db_index}", checkIfExists: true)
known_SNPs_zip = file("${params.site2_zip}", checkIfExists: true)
known_SNPs_tbi = file("${params.site2_zip_idx}", checkIfExists: true)

// Adapter Trimming, alignment and GATK BQSR (based on https://github.com/GavinHaLab/fastq_to_bam_paired_snakemake)
include { FASTQTOBAM  } from '../subworkflows/fastq_to_bam.nf'
include { ABRA_BAM } from '../modules/abra/realign/main.nf'
include { ABRA_SORT } from '../modules/samtools/abra_sort/main.nf'

// HSmetrics calculation
include { HSMETRICS ; HSMETRICS_COLLECT } from '../modules/gatk/hsmetrics/main.nf'

// COVERAGE calculation
include { COVERAGE_COVERVIEW } from '../modules/coverview/coverage/main.nf'
include { FORMAT_COVERVIEW } from '../modules/python/format_coverview/main.nf'
include { COVERAGE_BEDTOOLS } from '../modules/bedtools/coverage/main.nf'

// Variant calling
include { FREEBAYES } from '../modules/variant_call/freebayes/main.nf'
include { HAPLOTYPECALLER } from '../modules/variant_call/haplotypecaller/main.nf'
include { STRELKA } from '../modules/variant_call/strelka/main.nf'
include { PLATYPUS } from '../modules/variant_call/platypus/main.nf'
include { VARSCAN } from '../modules/variant_call/varscan/main.nf'
include { LOFREQ } from '../modules/variant_call/lofreq/main.nf'
include { DEEPSOMATIC } from '../modules/variant_call/deepsomatic/main.nf'
include { PINDEL_FLT3 } from '../modules/variant_call/pindel/pindel_flt3/main.nf'

// // Variant integration 
include { SOMATICSEQ_LYMPHOMA } from '../modules/variant_integration/somaticseq/main.nf'
include { SOMATICSEQ_LYMPHOMA_CONCAT } from '../modules/bcftools/concat/main.nf'

// Annotation
include { VEP_ANALYSIS as VEP_PINDEL ; VEP_ANALYSIS as VEP_SOMATICSEQ } from '../modules/vep/annotate/main.nf'
include { ANNOVAR as ANNOVAR_SOMATICSEQ ; ANNOVAR_LYMPHOMA } from '../modules/annovar/annotate/main.nf'

// // CNV calling
include { IFCNV } from '../modules/cnv_call/ifcnv/main.nf'
include { CNVKIT } from '../modules/cnv_call/cnvkit/main.nf'
include { PLOT_CNVKIT } from '../modules/python/plot_cnvkit/main.nf'

// // Translocation
include { SVABA } from '../modules/sv_call/svaba/main.nf'
include { LUMPY_PREPROCESS } from '../modules/samtools/lumpy_preprocess/main.nf'
include { LUMPY } from '../modules/sv_call/lumpy/main.nf'
include { GRIDSS } from '../modules/sv_call/gridss/main.nf'
include { DELLY } from '../modules/sv_call/delly/main.nf'
include { MANTA } from '../modules/sv_call/manta/main.nf'
include { FORMAT_SV_CALLS } from '../modules/python/format_sv_calls/main.nf'

// // Format output
include { FORMAT_PINDEL } from '../modules/python/format_pindel/main.nf'
include { FORMAT_ANNOVAR } from '../modules/python/format_annovar/main.nf'
include { CAVA } from '../modules/cava/annotate/main.nf'
include { FORMAT_CAVA } from '../modules/python/format_cava/main.nf'
include { MERGE_CSV_LYMPHOMA } from '../modules/python/merge_csv/main.nf'

workflow LYMPHOMA {
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
		lymphoma_bam_ch = FASTQTOBAM(bam_ch)
		ABRA_BAM(lymphoma_bam_ch.final_bams_ch, bedfile, genome_loc, index_files )
		ABRA_SORT(lymphoma_bam_ch.final_bams_ch.join(ABRA_BAM.out))
		COVERAGE_COVERVIEW(ABRA_SORT.out.final_bam, bedfile, coverview_config)
		FORMAT_COVERVIEW(COVERAGE_COVERVIEW.out)
		COVERAGE_BEDTOOLS(ABRA_SORT.out.final_bam, bedfile)	
		HSMETRICS(ABRA_SORT.out.final_bam, bedfile, genome_loc, index_files)
		all_hsmetrics = HSMETRICS.out.collect()
		HSMETRICS_COLLECT(all_hsmetrics)
		FREEBAYES(ABRA_SORT.out.final_bam, genome_loc, index_files, bedfile)
		HAPLOTYPECALLER(ABRA_SORT.out.final_bam, genome_loc, index_files, bedfile, known_SNPs, known_SNPs_index)
		STRELKA(ABRA_SORT.out.final_bam, genome_loc, index_files, bedfile, bedfile_zipped)
		PLATYPUS(ABRA_SORT.out.final_bam, genome_loc, index_files, bedfile)
		VARSCAN(ABRA_SORT.out.final_bam, genome_loc, index_files, bedfile)
		LOFREQ(ABRA_SORT.out.old_bam, genome_loc, index_files, bedfile)
		DEEPSOMATIC(ABRA_SORT.out.final_bam, control_final_bam, control_final_bamBai, genome_loc, index_files, bedfile)
		PINDEL_FLT3(ABRA_SORT.out.final_bam, genome_loc, index_files)
		VEP_PINDEL(PINDEL_FLT3.out)
		FORMAT_PINDEL(PINDEL_FLT3.out.join(VEP_PINDEL.out))
		SOMATICSEQ_LYMPHOMA(ABRA_SORT.out.final_bam.join(LOFREQ.out.join(VARSCAN.out.join(PLATYPUS.out.join(STRELKA.out.join(HAPLOTYPECALLER.out.join(FREEBAYES.out.join(DEEPSOMATIC.out))))))), genome_loc, index_files, bedfile, dbsnp_somatic)																						
		SOMATICSEQ_LYMPHOMA_CONCAT(SOMATICSEQ_LYMPHOMA.out)
		ANNOVAR_LYMPHOMA(SOMATICSEQ_LYMPHOMA_CONCAT.out, somaticseq)
		FORMAT_ANNOVAR(ANNOVAR_LYMPHOMA.out)
		IFCNV(ABRA_SORT.out.final_bam, bedfile)
		CNVKIT(ABRA_SORT.out.final_bam, cnvkit_reference, cnvkit_reference_delIG)
		PLOT_CNVKIT(CNVKIT.out.cnvkit_files, gene_scatter_list)
		CAVA(SOMATICSEQ_LYMPHOMA_CONCAT.out, cava_config, genome_loc, index_files, bedfile, bedfile_zipped, known_SNPs_zip, known_SNPs_tbi, ensembl_db, ensembl_db_index)
		FORMAT_CAVA(CAVA.out)
		SVABA(ABRA_SORT.out.old_bam, genome_loc, index_files, known_SNPs, known_SNPs_index, control_bqsr_bam, control_bqsr_bamBai)
		LUMPY_PREPROCESS(ABRA_SORT.out.old_bam)
		LUMPY(ABRA_SORT.out.old_bam.join(LUMPY_PREPROCESS.out))
		GRIDSS(ABRA_SORT.out.old_bam, genome_loc, index_files, control_bqsr_bam, control_bqsr_bamBai, gridss_exclude_list)
		DELLY(ABRA_SORT.out.old_bam, genome_loc, index_files, control_bqsr_bam, control_bqsr_bamBai)
		MANTA(ABRA_SORT.out.old_bam, genome_loc, index_files, control_bqsr_bam, control_bqsr_bamBai)
		FORMAT_SV_CALLS(SVABA.out.join(LUMPY.out.join(GRIDSS.out.join(DELLY.out.join(MANTA.out)))))
		MERGE_CSV_LYMPHOMA(FORMAT_PINDEL.out.join(FORMAT_ANNOVAR.out.join(FORMAT_CAVA.out.join(COVERAGE_BEDTOOLS.out))))
}

