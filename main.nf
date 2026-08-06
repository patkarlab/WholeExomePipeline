#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=
Sample list: ${params.input}
BED file: ${params.bedfile}
Sequences in:${params.sequences}
"""

include { WES } from './workflows/whole_exome.nf'
include { LYMPHOMA } from './workflows/lymphoma.nf'


workflow WHOLE_EXOME {
	WES ()
}

workflow LYMPHOMA_PANEL {
	LYMPHOMA ()
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the ${params.outdir} directory \n" : "Oops .. something went wrong" )
	log.info ( "Completed at: ${workflow.complete}")
	log.info ( "Total time taken: ${workflow.duration}")
}
