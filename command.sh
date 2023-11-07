#!/usr/bin/env bash 

source activate new_base
./wes.nf -entry WES \
--bedfile /home/diagnostics/pipelines/WholeExomePipeline/bedfiles/Exome \
-resume -bg
conda deactivate
