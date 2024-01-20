# WholeExomePipeline

## Description
Whole exome pipeline to identify germline variants.The pipeline trims adapters, creates paired end assembly, maps the reads to the human genome, uses GATK Best Practices to create BAM files. A variety of variant callers are used to call mutations followed by a machine learning based algorithm in SomaticSeq to create a consensus based VCF.Mosdepth tio calculate the coverage in each exome capture region.Detection FLT3-ITD from Illumina sequencing data using Get-ITD. The variants are then annotated with VEP.
## Pipeline Summary

1. Adaptor Trimming (fastq-mcf)

2. Merge paired-end reads (PEAR)

3. Alignment (bwa mem)

4. SAMtools conversion

5. Generating Final BAM files based on GATK Best Practices
a. RealignerTargetCreator
b. BaseRecalibrator
c. PrintReads
d. SAMtools sort and index on aligned recaliberated BAM files

6. Get-ITD.
7. Mosdepth 

8. Variant Calling using
a. Freebayes
b. Platypus
c. Haplotypecaller
d. Varscan
e. Lofreq
f. Strelka
g. Pindel

9. SomaticSeq with inputs from (Lofreq,Strelka,Varscan,Haplotypecaller,Platypus,Freebayes)
10. Annotation (using VEP)
11. CAVA
12. Merge_csv with inputs from (Pindel,SomaticSeq,CAVA)


## Executing the nextfow pipeline
It can be run using 
```
./run_nextflow.sh > script.log
```

The output files will be present in the respective sample directories inside the Final_Output folder
