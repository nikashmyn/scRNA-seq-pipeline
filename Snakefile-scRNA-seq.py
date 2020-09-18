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
        in_fasta = config["reference"],
        in_gtf = config["reference_gtf"]
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
        R1 = "/pellmanlab/stam_niko/data/ERR523111/{samples}_1.fastq.gz",
        R2 = "/pellmanlab/stam_niko/data/ERR523111/{samples}_2.fastq.gz",
        ref_dir = "/pellmanlab/stam_niko/STAR/genome/"
    output:
        mock = "/pellmanlab/stam_niko/data/bam/ERR523111/.{samples}_mockfile.txt"
    params:
        names = "/pellmanlab/stam_niko/data/bam/ERR523111/{samples}.",
        runMode = "alignReads",
        outstyle = "BAM Unsorted", #confusion on whether or not must be sorted
        readcmd = "zcat"
    shell:
        "STAR --genomeDir {input.ref_dir} --readFilesIn {input.R1} {input.R2} "
        "--runMode {params.runMode} --outSAMtype {params.outstyle} --readFilesCommand {params.readcmd} "
        "--outFileNamePrefix {params.names} "
        "&& touch {output.mock} "

#add read group information and calc tags
rule AddOrReplaceRG:
    input:
        mock = "/pellmanlab/stam_niko/data/bam/ERR523111/.{samples}_mockfile.txt"
    output:
        bam_out = "/pellmanlab/stam_niko/data/bam/ERR523111/{samples}.Aligned.RG.out.bam",
    params:
        bam_in = "/pellmanlab/stam_niko/data/bam/ERR523111/{samples}.Aligned.out.bam",
        rginfo = "--SORT_ORDER queryname --RGID {samples} --RGLB {samples} --RGPL illumina --RGPU {samples} --RGSM {samples} " #Needs to be revisited on real data
    shell:
        "picard AddOrReplaceReadGroups -I {params.bam_in} -O {output.bam_out} {params.rginfo} "

#remove duplicates step        
rule mark_duplicates:
    input:
        bam_in = "/pellmanlab/stam_niko/data/bam/ERR523111/{samples}.Aligned.RG.out.bam"
    output:
        metrics = "/pellmanlab/stam_niko/data/bam/ERR523111/mark_dups/{samples}.metrics.txt", #File to write duplication metrics to  Required.
        bam_out = "/pellmanlab/stam_niko/data/bam/ERR523111/mark_dups/{samples}.Aligned.deduped.out.bam" #The output file to write marked records to  Required
    params:
        options = f"--CREATE_INDEX true --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER queryname --TMP_DIR {config['tmp_dir']} " #--ASSUME_SORT_ORDER queryname  #not coord sorted not queryname sorted
    shell:
        "picard MarkDuplicates --INPUT {input.bam_in} {params.options} --OUTPUT {output.bam_out} --METRICS_FILE {output.metrics} "
       

#sort bam step
rule sort_sam:
    input:
        bam_in = "/pellmanlab/stam_niko/data/bam/ERR523111/mark_dups/{samples}.Aligned.deduped.out.bam" #Input BAM or SAM file to sort.  Required
    output:
        bam_out =  "/pellmanlab/stam_niko/data/bam/ERR523111/mark_dups/{samples}.Aligned.sorted.deduped.out.bam" #Sorted BAM or SAM output file.  Required.
    params:
        SO = "coordinate", #Sorts primarily according to the SEQ and POS fields of the record.
        options = "--CREATE_INDEX true"
    shell:
        "picard SortSam --INPUT {input.bam_in} {params.options} --OUTPUT {output.bam_out} --SORT_ORDER {params.SO} " 
        
#calc tags
rule calctags:
    input:
        bam_in = "/pellmanlab/stam_niko/data/bam/ERR523111/{samples}.Aligned.sorted.deduped.out.bam"
    output:
        bam_out = "/pellmanlab/stam_niko/data/bam/ERR523111/{samples}.Aligned.sorted.deduped.tagged.out.bam"
    params:
        ref = config["reference_unzip"]
    shell:
        "picard SetNmMdAndUqTags -I {input.bam_in} -O {output.bam_out} --REFERENCE_SEQUENCE {params.ref} "

# Split Reads with N in Cigar
rule splitncigarreads:
    input:
        bam_in =  "/pellmanlab/stam_niko/data/bam/ERR523111/mark_dups/{samples}.Aligned.sorted.deduped.tagged.out.bam", #BAM/SAM/CRAM file containing reads
        ref = config["reference_unzip"] #Reference sequence file and index Required.
    output:
        bam_out = "/pellmanlab/stam_niko/data/bam/ERR523111/split_reads/{samples}.Aligned.sorted.deduped.split_r.out.bam" #Write output to this BAM filename  Required. 
    params:
        "--create-output-bam-index true"
        #f"--tmp-dir {config['tmp_dir']} "
    shell:
        "gatk SplitNCigarReads --input {input.bam_in} --reference {input.ref} {params} --output {output.bam_out} "
        
# Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context)
rule base_recal:
    input:
        bam_in = "/pellmanlab/stam_niko/data/bam/ERR523111/split_reads/{samples}.Aligned.sorted.deduped.split_r.out.bam", #This argument must be specified at least once. Required.
        ref = config["reference_unzip"], #Reference sequence file  Required.
        KS = config["ref_var"] #One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis.  This argument must be specified at least once. Required.
    output:
        bam_out = "/pellmanlab/stam_niko/data/bam/ERR523111/split_reads/{samples}.Aligned.sorted.deduped.split_r.recal.out.bam"  #The output recalibration table file to create  Required.
    params:
        "--create-output-bam-index True",
        "--read-filter PairedReadFilter"
    shell:
        "gatk BaseRecalibrator --input {input.bam_in} --reference {input.ref} --known-sites {input.KS} {params} --output {output.bam_out} "

