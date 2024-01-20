# WholeExomePipeline

## Description
Whole exome pipeline to identify germline variants.The pipeline trims adapters, creates paired end assembly, maps the reads to the human genome, uses GATK Best Practices to create BAM files. A variety of variant callers are used to call mutations followed by a machine learning based algorithm in SomaticSeq to create a consensus based VCF.Mosdepth tio calculate the coverage in each exome capture region.Detection FLT3-ITD from Illumina sequencing data using Get-ITD. The variants are then annotated with VEP.
## Pipeline Summary
1. Adaptor Trimming (fastq-mcf)
1. Merge paired-end reads (PEAR)
1. Alignment (bwa mem)
1. SAMtools conversion
1. Generating Final BAM files based on GATK Best Practices
	- RealignerTargetCreator
	- BaseRecalibrator
	- PrintReads
	- SAMtools sort and index on aligned recaliberated BAM files

1. Mosdepth 
1. Variant Calling using
	- Freebayes
	- Platypus
	- Haplotypecaller
	- Varscan
	- Lofreq
	- Strelka
	- Pindel

1. SomaticSeq with inputs from (Lofreq,Strelka,Varscan,Haplotypecaller,Platypus,Freebayes)
1. Annotation (using VEP)
1. CAVA
1. Merge_csv with inputs from (Pindel,SomaticSeq,CAVA)

## Executing the nextfow pipeline
It can be run using 
```
./run_nextflow.sh > script.log
```

The output files will be present in the respective sample directories inside the Final_Output folder
