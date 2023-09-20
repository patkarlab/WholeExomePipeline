#!/usr/bin/env bash 

source activate new_base
./wes.nf -entry WES \
--bedfile /home/pipelines/mutation_detector_nextflow/bedfile/Leukemia_Panel_Myeloid_2023_Feb_hg37_sortd \
-resume -bg
conda deactivate
