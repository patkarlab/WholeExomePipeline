#!/usr/bin/env bash 



#WHOLE-EXOME

#source activate new_base

#nextflow -c /home/diagnostics/pipelines/WholeExomePipeline/nextflow.config run main.nf -entry WHOLE_EXOME \
#--bedfile /home/diagnostics/pipelines/WholeExomePipeline/bedfiles/Exome_hg19_sortd \
#--output /home/diagnostics/pipelines/WholeExomePipeline/Final_Output/ \
#-resume -bg

#conda deactivate


#Validation - Duplex

source activate new_base

nextflow -c /home/diagnostics/pipelines/WholeExomePipeline/nextflow.config run main.nf -entry WHOLE_EXOME \
--bedfile /home/pipelines/NextSeq_mutation_detector_leukemia/bedfiles/AML_MRD_DUPLEX_probes_sortd \
--output /home/diagnostics/pipelines/WholeExomePipeline/Final_Output/ \
-resume -bg

conda deactivate


#LYMPHOMA

#source activate new_base

#nextflow -c /home/diagnostics/pipelines/WholeExomePipeline/nextflow.config run main.nf -entry LYMPHOMA_PANEL \
#--bedfile /home/diagnostics/pipelines/WholeExomePipeline/bedfiles/lymphoma_probes_21102022_sortd \
#--output /home/diagnostics/pipelines/WholeExomePipeline/Final_Output/ \
#-resume -bg

#conda deactivate

