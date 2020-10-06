#######################################
####### Snakefile for scRNA-seq #######
#######################################

################
#### IMPORTS ###
################

import os
import pandas as pd
import snakemake

########################
### GLOBAL CONSTANTS ###
########################

#Practice data
prac_samples = list(["ERR523111"])

########################
### INPUT-ONLY RULES ###
########################
rule align:
    input:
        expand("/pellmanlab/stam_niko/data/bam/ERR523111/.{samples}_mockfile.txt", samples=prac_samples)

rule mrk_dups:
    input:
        expand("/pellmanlab/stam_niko/data/bam/ERR523111/mark_dups/{samples}.Aligned.deduped.out.bam", samples=prac_samples)

rule re_sort:
    input:
        expand("/pellmanlab/stam_niko/data/bam/ERR523111/mark_dups/{samples}.Aligned.sorted.deduped.out.bam", samples=prac_samples)

rule split_r:
    input:
        expand("/pellmanlab/stam_niko/data/bam/ERR523111/split_reads/{samples}.Aligned.sorted.deduped.split_r.out.bam", samples=prac_samples)

rule recal:
    input:
        expand("/pellmanlab/stam_niko/data/bam/ERR523111/split_reads/{samples}.Aligned.sorted.deduped.split_r.recal.out.bam", samples=prac_samples)

#rule genotype:
#    input:
#        expand()

#######################
### SNAKEMAKE RULES ###
#######################

#premapping and aligning step for STAR
rule generate_genome_indexes:
    input:
    output:
        ref_dir = "/pellmanlab/stam_niko/STAR/genome/"
    params:
        runMode = "genomeGenerate",
        overhang = config["readlength"] - 1
    threads: config["MAX_THREADS"]
    shell:
        "STAR --genomeFastaFiles {input.in_fasta} --sjdbGTFfile {input.in_gtf} "
        "--runThreadN {threads} --runMode {params.runMode} --sjdbOverhang {params.overhang} "
        "--genomeDir {output.ref_dir} "

#mapping, aligning, and sorting step
rule STAR_alignment:
    input:
        ref_dir = "/pellmanlab/stam_niko/STAR/genome/"
    output:
        mock = "/pellmanlab/stam_niko/data/bam/ERR523111/.{samples}_mockfile.txt"
    params:
        runMode = "alignReads",
    shell:
        "STAR --genomeDir {input.ref_dir} --readFilesIn {input.R1} {input.R2} "

#add read group information and calc tags
rule AddOrReplaceRG:
    input:
        mock = "/pellmanlab/stam_niko/data/bam/ERR523111/.{samples}_mockfile.txt"
    output:
    params:
        rginfo = "--SORT_ORDER queryname --RGID {samples} --RGLB {samples} --RGPL illumina --RGPU {samples} --RGSM {samples} " #Needs to be revisited on real data
    shell:
        "picard AddOrReplaceReadGroups -I {params.bam_in} -O {output.bam_out} {params.rginfo} "

#remove duplicates step        
rule mark_duplicates:
    input:
    output:
    params:
    shell:
        "picard MarkDuplicates --INPUT {input.bam_in} {params.options} --OUTPUT {output.bam_out} --METRICS_FILE {output.metrics} "
       
#sort bam step
rule sort_sam:
    input:
        bam_in = "/pellmanlab/stam_niko/data/bam/ERR523111/mark_dups/{samples}.Aligned.deduped.out.bam"
    output:
    params:
        options = "--CREATE_INDEX true"
    shell:
        "picard SortSam --INPUT {input.bam_in} {params.options} --OUTPUT {output.bam_out} --SORT_ORDER {params.SO} " 
        
#calc tags
rule calctags:
    input:
    output:
    params:
        ref = config["reference_unzip"]
    shell:
        "picard SetNmMdAndUqTags -I {input.bam_in} -O {output.bam_out} --REFERENCE_SEQUENCE {params.ref} "

# Split Reads with N in Cigar
rule splitncigarreads:
    input:
    output:
    params:
        "--create-output-bam-index true"
        #f"--tmp-dir {config['tmp_dir']} "
    shell:
        "gatk SplitNCigarReads --input {input.bam_in} --reference {input.ref} {params} --output {output.bam_out} "
        
# Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context)
rule base_recal:
    input:
    output:
        recal_out = "/pellmanlab/stam_niko/data/bam/ERR523111/split_reads/{samples}.recal_table.txt",  #The output recalibration table file to create  Required.
        bam_out  = "/pellmanlab/stam_niko/data/bam/ERR523111/split_reads/{samples}.Aligned.sorted.deduped.split_r.recal.out.bam",
    params:
        recal = "--read-filter PairedReadFilter",
        applybqsr = "--create-output-bam-index True"
    shell:
        "gatk BaseRecalibrator --input {input.bam_in} --reference {input.ref} --known-sites {input.KS} {params.recal} --output {output.recal_out} "
        "&& gatk ApplyBQSR --input {input.bam_in} --bqsr-recal-file {output.recal_out} --reference {input.ref} {params.applybqsr} --output {output.bam_out} "
