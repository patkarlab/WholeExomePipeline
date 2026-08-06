process PLATYPUS {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val (Sample), file(bam), file(bamBai)
		path (GenFile)
		path (GenDir)        
        path (bedfile)
	output:
		tuple val(Sample), file ("${Sample}.platypus.vcf")
	script:
	"""
    awk 'BEGIN{FS="\t";OFS=""}{print \$1,":",\$2,"-",\$3}' ${bedfile} > ${bedfile}_sortd_regions.txt
	python2.7 /code/Platypus/bin/Platypus.py callVariants --bamFiles=${bam} --refFile=${GenFile} --output=${Sample}.platypus.vcf --nCPU=${task.cpus} --minFlank=0 --filterDuplicates=0 --maxVariants=6 --minReads=6 --regions=${bedfile}_sortd_regions.txt
	"""
	stub:
	"""
	touch ${Sample}.platypus.vcf
	"""
}
