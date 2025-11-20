#!/usr/bin/bash

bam_file=$1
bed_file=$2
Sample=$3

/home/programs/CoverView-1.4.3/coverview -i ${bam_file} -b ${bed_file}.bed -c /home/programs/CoverView-1.4.3/config/config.txt -o ${Sample}.coverview
python3 "/home/pipelines/mutation_detector_nextflow/scripts/coverview.py" ${Sample}.coverview_regions.txt Final_Output_Lymphoma_new/${Sample}/${Sample}.coverview_regions.csv