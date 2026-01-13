#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=
Sample list: ${params.input}
BED file: ${params.bedfile}.bed
Sequences in:${params.sequences}
"""

gene_scatter = file("${params.gene_scatter_list}/exonwise_lymphoma_list.txt", checkIfExists: true )
cnvkit_reference = file("${params.cnvkitRef}", checkIfExists: true )
cnvkit_reference_delIG = file("${params.cnvkitRefDelIG}", checkIfExists: true )
normal_bamfile = file("${params.normal_bam}", checkIfExists: true )
normal_bamBaifile = file("${params.normal_bam}.bai", checkIfExists: true )
genome_loc = file("${params.genome}", checkIfExists: true )
genome_dir = file("${genome_loc.parent}", checkIfExists: true)
genome_fasta = file("${genome_loc.name}")
dbsnp = file("${params.site2}", checkIfExists: true)
bedfile = file("${params.bedfile}.bed", checkIfExists: true)
bed_regions = file("${params.bedfile}_regions.txt", checkIfExists: true)

// Adapter Trimming, alignment and GATK BQSR (based on https://github.com/GavinHaLab/fastq_to_bam_paired_snakemake)
include { FASTQTOBAM; ABRA_BAM } from '../modules/processes.nf'

// HSmetrics calculation
include { HSMETRICS } from '../modules/hsmetrics.nf'

// COVERAGE calculation
include { COVERAGE_BEDTOOLS; COVERAGE_MOSDEPTH; COVERAGE_COVERVIEW} from '../modules/coverage.nf'

// Variant calling
include { FREEBAYES; HAPLOTYPECALLER; STRELKA; PLATYPUS; VARSCAN; DEEPVARIANT; LOFREQ;  PINDEL } from '../modules/variant_call.nf'

// Variant integration 
include { SOMATICSEQ } from '../modules/somaticseq.nf'

// CNV calling
include { IFCNV; GCNV } from '../modules/cnv_call.nf'

// Format output
include { DEEPVARIANT_GCNV; CAVA; MERGE_CSV_WES } from '../modules/format_output.nf'


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
		final_bams_ch = FASTQTOBAM(bam_ch)
		ABRA_BAM(final_bams_ch)
		COVERAGE_MOSDEPTH(ABRA_BAM.out)
		HSMETRICS(ABRA_BAM.out)
		FREEBAYES(ABRA_BAM.out, genome_dir, genome_fasta, bedfile)
		HAPLOTYPECALLER(ABRA_BAM.out, genome_dir, genome_fasta, bedfile, dbsnp)
		STRELKA(ABRA_BAM.out, genome_dir, genome_fasta)
		PLATYPUS(ABRA_BAM.out, genome_dir, genome_fasta, bed_regions)
		VARSCAN(ABRA_BAM.out, genome_dir, genome_fasta, bedfile)
		DEEPVARIANT(ABRA_BAM.out)
		LOFREQ(ABRA_BAM.out, genome_dir, genome_fasta, bedfile)
		PINDEL(ABRA_BAM.out, genome_dir, genome_fasta)
		SOMATICSEQ(LOFREQ.out.join(VARSCAN.out.join(PLATYPUS.out.join(STRELKA.out.join(HAPLOTYPECALLER.out.join(FREEBAYES.out.join(ABRA_BAM.out)))))), genome_dir, genome_fasta)
		IFCNV(ABRA_BAM.out.collect())
		GCNV(ABRA_BAM.out)		
		DEEPVARIANT_GCNV(DEEPVARIANT.out.join(GCNV.out))
		CAVA(SOMATICSEQ.out)
		MERGE_CSV_WES(PINDEL.out.join(SOMATICSEQ.out.join(DEEPVARIANT.out.join(CAVA.out.join(COVERAGE_MOSDEPTH.out.mosdepth_out)))))
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the ${params.output} directory \n" : "Oops .. something went wrong" )
}
