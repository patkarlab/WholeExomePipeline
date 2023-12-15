#!/usr/bin/env python3

import pandas as pd
import sys
import re
import csv


txt_file = sys.argv[1]
input_xlsx_file = sys.argv[2] #The first file is .xlsx containing the pharmacogenic markers to be mapped
output = sys.argv[3]

map = dict()                    # Dictionary for markers
map_counts = dict()             # Dictionary for the marker occurence
map_comments = dict()   # Dictionary for comments
xlsx_file_data = pd.read_excel(input_xlsx_file, engine='openpyxl', usecols = [0,1,2,3,4,5])
for i in range(len(xlsx_file_data)):
    chromosome = xlsx_file_data.iloc[i, 0]
    position = xlsx_file_data.iloc[i, 1]
    ref = xlsx_file_data.iloc[i, 2]
    alt = xlsx_file_data.iloc[i, 3]
    rsid = xlsx_file_data.iloc[i, 4]
    comment = str (xlsx_file_data.iloc[i, 5]).replace(",", "&")
    comment = comment.lower()

    if re.search('[a-zA-Z]', str(chromosome)):
        chromosome = re.sub ("chr","", chromosome, flags = re.IGNORECASE)
        chromosome = re.sub ("X","23", chromosome, flags = re.IGNORECASE)
        chromosome = re.sub ("Y","y", chromosome, flags = re.IGNORECASE)

    id = str(chromosome) + ':' +''.join(str(position)) + ':' +''.join(str(ref)) + ':' +''.join(str(alt)) 
    map[id] = rsid
    map_counts[id] = 0
    map_comments[id] = comment


outfile = open (output,'w')


print("CHROM","Start","End","REF","ALT","Gene","ID","ref_count","alt_count","VAF","AF","MAX_AF_POPS","ALLELE_NUM","Variant_callers","VARIANT_CLASS","CLIN_SIG","Feature","Feature_type","Consequence","HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","ENSP","Func.refGene","Gene.refGene","ExonicFunc.refGene","PopFreqMax","InterVar_automated","COMMENTS", file=outfile, sep="\t" )


with open(txt_file, 'r') as txt:
    txt_handle = csv.reader(txt, delimiter='\t')
    header = next(txt_handle)
    for str_lines in txt_handle:
        vcf_chr = str_lines[0]
        vcf_start = str_lines[1]
        vcf_end = str_lines[2]
        vcf_ref = str_lines[3]
        vcf_alt = str_lines[4]
        vcf_Gene = str_lines[5]
        vcf_ID= str_lines[6]
        vcf_REF_COUNT = str_lines[7]
        vcf_ALT_COUNT = str_lines[8]
        vcf_VAF = str_lines[9]
        vcf_af = str_lines[10]
        vcf_PopFreqMax = str_lines[11]
        vcf_allele_num= str_lines[12]
        vcf_variant_callers = str_lines[13]
        vcf_VARIANT_CLASS =str_lines[14]
        vcf_CLIN_SIG = str_lines[15]
        vcf_Feature = str_lines[16]
        vcf_Feature_type =str_lines[17]
        vcf_Consequence = str_lines[18]
        vcf_HGVSc = str_lines[19]
        vcf_HGVSp = str_lines[20]
        vcf_cDNA_position = str_lines[21]
        vcf_CDS_position = str_lines[22]
        vcf_Protein_position = str_lines[23]
        vcf_Amino_acids = str_lines[24]
        vcf_Codons =str_lines[25]
        vcf_ENSP = str_lines[26]
        vcf_Func_refGene = str_lines[27]
        vcf_Gene_refGene = str_lines[28]
        vcf_ExonicFunc_refGene = str_lines[29]
        vcf_PopFreqMax =str_lines[30]
        vcf_InterVar_automated= str_lines[31]


        if re.search('[a-zA-Z]', str(vcf_chr)):
            vcf_chr = re.sub ("chr","", vcf_chr, flags = re.IGNORECASE)
            vcf_chr = re.sub ("X","23", vcf_chr, flags = re.IGNORECASE)
            vcf_chr = re.sub ("Y","y", vcf_chr, flags = re.IGNORECASE)

        vcf_id = str(vcf_chr) + ':' +''.join(vcf_start) +  ':' +''.join(vcf_ref) + ':' +''.join(vcf_alt)
        #print(vcf_id)
        
        if vcf_id in map:
            map_counts[vcf_id] = map_counts[vcf_id] + 1
            print (vcf_chr, vcf_start,vcf_end, vcf_ref, vcf_alt, vcf_Gene, vcf_ID,vcf_REF_COUNT, vcf_ALT_COUNT, vcf_VAF, vcf_af, vcf_PopFreqMax, vcf_allele_num,vcf_variant_callers, vcf_VARIANT_CLASS, vcf_CLIN_SIG, vcf_Feature, vcf_Feature_type, vcf_Consequence, vcf_HGVSc, vcf_HGVSp,vcf_cDNA_position, vcf_CDS_position, vcf_Protein_position, vcf_Amino_acids, vcf_Codons, vcf_ENSP, vcf_Func_refGene, vcf_Gene_refGene, vcf_ExonicFunc_refGene, vcf_PopFreqMax ,vcf_InterVar_automated,map_comments[vcf_id].strip("\n"), file=outfile, sep="\t")
            #print(vcf_chr, vcf_start,vcf_variant_callers, no_of_callers, vcf_Gene, vcf_REF_COUNT, vcf_ALT_COUNT, vcf_ref, vcf_alt, vcf_VAF,map[vcf_id], map_comments[vcf_id])


outfile.close()

