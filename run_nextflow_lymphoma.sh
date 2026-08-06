#!/usr/bin/env bash 

modified_samplesheet=$1

#source activate new_base

nextflow -c nextflow.config run main.nf -entry LYMPHOMA_PANEL \
--input ${modified_samplesheet} -profile docker \
--bedfile /home/diagnostics/pipelines/WholeExomePipeline/assets/TMC_Lymphoma_2.0_master_baseline_hg19_1.04mb_sortd.bed \
-resume -bg

#conda deactivate
