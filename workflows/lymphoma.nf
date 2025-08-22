#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=
Sample list: ${params.input}
BED file: ${params.bedfile}.bed
Sequences in:${params.sequences}
"""

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

// Translocation
include { SVABA; LUMPY; TRANSLOCATION } from '../modules/translocation.nf'

// Format output
include { DEEPVARIANT_GCNV; CAVA; MERGE_CSV} from '../modules/format_output.nf'

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
		final_bams_ch = FASTQTOBAM(bam_ch)
		ABRA_BAM(final_bams_ch)
		COVERAGE_COVERVIEW(ABRA_BAM.out)
		COVERAGE_BEDTOOLS(ABRA_BAM.out)	
		HSMETRICS(ABRA_BAM.out)
		FREEBAYES(ABRA_BAM.out)
		HAPLOTYPECALLER(ABRA_BAM.out)
		STRELKA(ABRA_BAM.out)
		PLATYPUS(ABRA_BAM.out)
		VARSCAN(ABRA_BAM.out)
		DEEPVARIANT(ABRA_BAM.out)
		LOFREQ(ABRA_BAM.out)
		PINDEL(ABRA_BAM.out)
		SOMATICSEQ(LOFREQ.out.join(VARSCAN.out.join(PLATYPUS.out.join(STRELKA.out.join(HAPLOTYPECALLER.out.join(FREEBAYES.out.join(ABRA_BAM.out)))))))																						
		IFCNV(ABRA_BAM.out.collect())
		GCNV(ABRA_BAM.out)		
		DEEPVARIANT_GCNV(DEEPVARIANT.out.join(GCNV.out))
		CAVA(SOMATICSEQ.out)
		SVABA(ABRA_BAM.out)
		LUMPY(ABRA_BAM.out)
		TRANSLOCATION(SVABA.out.join(LUMPY.out))
		MERGE_CSV(PINDEL.out.join(SOMATICSEQ.out.join(DEEPVARIANT.out.join(CAVA.out.join(COVERAGE_BEDTOOLS.out)))))
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the ${params.output} directory \n" : "Oops .. something went wrong" )
}
