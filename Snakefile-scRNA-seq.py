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
        in_fasta = config["reference_unzip"],
        in_gtf = config["reference_gtf_unzip"]
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
        R1 = "/pellmanlab/stam_niko/data/ERR523111/{samples}_1.fastq", #.gz use readcmd
        R2 = "/pellmanlab/stam_niko/data/ERR523111/{samples}_2.fastq", #.gz
        ref_dir = "/pellmanlab/stam_niko/STAR/genome/"
    output:
        mock = "/pellmanlab/stam_niko/data/bam/ERR523111/.{samples}_mockfile.txt"
    params:
        names = "/pellmanlab/stam_niko/data/bam/ERR523111/raw_aligned_bams/{samples}.",
        runMode = "alignReads",
        outstyle = "BAM Unsorted", 
        unmapped = "Within KeepPairs",
        twopass = "Basic",
        tags = "NH HI AS NM MD"
        #readcmd = "zcat"
    shell:
        "STAR --genomeDir {input.ref_dir} --readFilesIn {input.R1} {input.R2} "
        "--runMode {params.runMode} --twopassMode {params.twopass} --outFileNamePrefix {params.names} " #--readFilesCommand {params.readcmd} "
        "--outSAMunmapped {params.unmapped} --outSAMtype {params.outstyle} --outSAMattributes {params.tags} " # For RSEM compatability #like --alignbytype endtoend 
        "&& touch {output.mock} " #For snakemake step tracing

#add read group information and calc tags
rule AddOrReplaceRG:
    input:
        mock = "/pellmanlab/stam_niko/data/bam/ERR523111/.{samples}_mockfile.txt"
    output:
        bam_out = "/pellmanlab/stam_niko/data/bam/ERR523111/raw_aligned_bams/{samples}.Aligned.RG.out.bam",
    params:
        bam_in = "/pellmanlab/stam_niko/data/bam/ERR523111/raw_aligned_bams/{samples}.Aligned.out.bam",
        rginfo = "--SORT_ORDER queryname --RGID {samples} --RGLB {samples} --RGPL illumina --RGPU {samples} --RGSM {samples} " #Needs to be revisited on real data
    shell:
        "picard AddOrReplaceReadGroups -I {params.bam_in} -O {output.bam_out} {params.rginfo} "

#remove duplicates step        
rule mark_duplicates:
    input:
        bam_in = "/pellmanlab/stam_niko/data/bam/ERR523111/raw_aligned_bams/{samples}.Aligned.RG.out.bam"
    output:
        metrics = "/pellmanlab/stam_niko/data/bam/ERR523111/mark_dups/{samples}.metrics.txt", #File to write duplication metrics.
        bam_out = "/pellmanlab/stam_niko/data/bam/ERR523111/mark_dups/{samples}.Aligned.deduped.out.bam"
    params:
        options = f"--CREATE_INDEX true --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER queryname --TMP_DIR {config['tmp_dir']} " 
    shell:
        "picard MarkDuplicates --INPUT {input.bam_in} {params.options} --OUTPUT {output.bam_out} --METRICS_FILE {output.metrics} "
       
#sort bam step
rule sort_sam:
    input:
        bam_in = "/pellmanlab/stam_niko/data/bam/ERR523111/mark_dups/{samples}.Aligned.deduped.out.bam"
    output:
        bam_out =  "/pellmanlab/stam_niko/data/bam/ERR523111/mark_dups/{samples}.Aligned.sorted.deduped.out.bam" 
    params:
        SO = "coordinate", #Downstream tools need coord sorted
        options = "--CREATE_INDEX true"
    shell:
        "picard SortSam --INPUT {input.bam_in} {params.options} --OUTPUT {output.bam_out} --SORT_ORDER {params.SO} " 
        
#calc tags
rule calctags:
    input:
        bam_in = "/pellmanlab/stam_niko/data/bam/ERR523111/mark_dups/{samples}.Aligned.sorted.deduped.out.bam"
    output:
        bam_out = "/pellmanlab/stam_niko/data/bam/ERR523111/mark_dups/{samples}.Aligned.sorted.deduped.tagged.out.bam"
    params:
        ref = config["reference_unzip"]
    shell:
        "picard SetNmMdAndUqTags -I {input.bam_in} -O {output.bam_out} --REFERENCE_SEQUENCE {params.ref} "

# Split Reads with N in Cigar
rule splitncigarreads:
    input:
        bam_in =  "/pellmanlab/stam_niko/data/bam/ERR523111/mark_dups/{samples}.Aligned.sorted.deduped.tagged.out.bam", 
        ref = config["ref_unzip_wdict"] 
    output:
        bam_out = "/pellmanlab/stam_niko/data/bam/ERR523111/split_reads/{samples}.Aligned.sorted.deduped.split_r.out.bam"
    params:
        "--create-output-bam-index true"
        #f"--tmp-dir {config['tmp_dir']} "
    shell:
        "gatk SplitNCigarReads --input {input.bam_in} --reference {input.ref} {params} --output {output.bam_out} "
        
# Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context)
rule base_recal:
    input:
        bam_in = "/pellmanlab/stam_niko/data/bam/ERR523111/split_reads/{samples}.Aligned.sorted.deduped.split_r.out.bam", 
        ref = config["ref_unzip_wdict"], 
        KS = config["ref_var"] #Known polymorphic sites used to exclude regions around known polymorphisms from analysis.
    output:
<<<<<<< HEAD
        recal_out = "/pellmanlab/stam_niko/data/bam/ERR523111/split_reads/{samples}.recal_table.txt",  #The output recalibration table file to create  Required.
        bam_out  = "/pellmanlab/stam_niko/data/bam/ERR523111/split_reads/{samples}.Aligned.sorted.deduped.split_r.recal.out.bam",
=======
        recal_out = "/pellmanlab/stam_niko/data/bam/ERR523111/split_reads/{samples}.recal_table.txt",  #The output recalibration table file to create.
        bam_out  = "/pellmanlab/stam_niko/data/bam/ERR523111/split_reads/{samples}.Aligned.sorted.deduped.split_r.recal.out.bam", #The BAM to build from out table ^
>>>>>>> b5356e03f9f06ad8f29bbf03423f9004bcd54b67
    params:
        recal = "--read-filter PairedReadFilter",
        applybqsr = "--create-output-bam-index True"
    shell:
        "gatk BaseRecalibrator --input {input.bam_in} --reference {input.ref} --known-sites {input.KS} {params.recal} --output {output.recal_out} "
        "&& gatk ApplyBQSR --input {input.bam_in} --bqsr-recal-file {output.recal_out} --reference {input.ref} {params.applybqsr} --output {output.bam_out} "
<<<<<<< HEAD
=======

#BAM input to RSEM differential expression 
### TO DO ### Add BAM input options
#rule rsem_calc_expr:
#    input:
#        R1 = f"/pellmanlab/stam_niko/data/ERR523111/{samples}_1.fastq",
#        R2 = f"/pellmanlab/stam_niko/data/ERR523111/{samples}_2.fastq"
#    params:
#        ref_name = "EnsembleGenome98",
#        sample_name = f"{samples}",
#        options = "--star --calc-pme --calc-ci --output-genome-bam ",
#        num_thrd = "--num-threads 4",
#        seed = "--seed 12 "
#    shell:
#        "rsem-calculate-expression {params.options} {params.num_thrd} {params.seed} "
#        "--paired-end {input.R1} {input.R2} {params.ref_name} {params.sample_name} "
>>>>>>> b5356e03f9f06ad8f29bbf03423f9004bcd54b67
