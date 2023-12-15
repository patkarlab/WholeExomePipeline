#!/usr/bin/env bash 

source activate new_base
./wes.nf -entry WES \
--bedfile /home/diagnostics/pipelines/WholeExomePipeline/bedfiles/Exome_hg19_sortd \
-resume -bg
conda deactivate
