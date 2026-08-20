#! /usr/bin/bash

nextflow -c nextflow.config run main.nf -entry LYMPHOMA_PANEL \
	--bedfile /home/patkarlab-clinical/pipelines/WholeExomePipeline/assets/probes_ok_ACTREC_Lymphoma2p0_TE-92182112_hg38_Rev1_260420063459_sortd.bed \
	-profile docker -resume
