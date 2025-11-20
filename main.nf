#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { WES } from './workflows/whole_exome.nf'
include { LYMPHOMA } from './workflows/lymphoma.nf'

//
// WORKFLOW: Run main fastq to bam analysis pipeline
//

workflow WHOLE_EXOME {
	WES ()
}

workflow LYMPHOMA_PANEL {
	LYMPHOMA ()
}
