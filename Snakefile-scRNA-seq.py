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
        expand("/pellmanlab/stam_niko/data/processed_bam/ERR523111/.{samples}_mockfile.txt", samples=prac_samples)

rule mrk_dups:
    input:
        expand("/pellmanlab/stam_niko/data/processed_bam/ERR523111/mark_dups/{samples}.Aligned.toTranscriptome.deduped.out.bam", samples=prac_samples)

rule re_sort:
    input:
        expand("/pellmanlab/stam_niko/data/processed_bam/ERR523111/mark_dups/{samples}.Aligned.toTranscriptome.sorted.deduped.out.bam", samples=prac_samples)

rule rsem:
    input:
        expand("/pellmanlab/stam_niko/data/processed_bam/ERR523111/RSEM/.{samples}_mockfile.rsem.txt", samples=prac_samples)

rule calc_expr:
    input:
        expand("/pellmanlab/stam_niko/data/processed_bam/ERR523111/RSEM/.{samples}_mockfile.rsem_calc.txt", samples=prac_samples)

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
        mock = "/pellmanlab/stam_niko/data/processed_bam/ERR523111/.{samples}_mockfile.txt"
    params:
        names = "/pellmanlab/stam_niko/data/processed_bam/ERR523111/{samples}.",
        runMode = "alignReads",
        outstyle = "BAM Unsorted",
        unmapped = "Within KeepPairs",
        twopass = "Basic",
        tags = "NH HI AS NM MD",
        quant_mode = "TranscriptomeSAM"
        #readcmd = "zcat"
    shell:
        "STAR --genomeDir {input.ref_dir} --readFilesIn {input.R1} {input.R2} "
        "--runMode {params.runMode} --twopassMode {params.twopass} --outFileNamePrefix {params.names} --quantMode {params.quant_mode} " #--readFilesCommand {params.readcmd} "
        "--outSAMunmapped {params.unmapped} --outSAMtype {params.outstyle} --outSAMattributes {params.tags} " # For RSEM compatability
        "&& touch {output.mock} " #For snakemake step tracing

#add read group information
rule AddOrReplaceRG:
    input:
        mock = "/pellmanlab/stam_niko/data/processed_bam/ERR523111/.{samples}_mockfile.txt"
    output:
        bam_out = "/pellmanlab/stam_niko/data/processed_bam/ERR523111/{samples}.Aligned.toTranscriptome.RG.out.bam",
    params:
        bam_in = "/pellmanlab/stam_niko/data/processed_bam/ERR523111/{samples}.Aligned.toTranscriptome.out.bam",
        rginfo = "--SORT_ORDER queryname --RGID {samples} --RGLB {samples} --RGPL illumina --RGPU {samples} --RGSM {samples} " #Needs to be revisited on real data
    shell:
        "picard AddOrReplaceReadGroups -I {params.bam_in} -O {output.bam_out} {params.rginfo} "

#remove duplicates step        
rule mark_duplicates:
    input:
        bam_in = "/pellmanlab/stam_niko/data/processed_bam/ERR523111/{samples}.Aligned.toTranscriptome.RG.out.bam"
    output:
        metrics = "/pellmanlab/stam_niko/data/processed_bam/ERR523111/mark_dups/{samples}.toTranscriptome.metrics.txt", 
        bam_out = "/pellmanlab/stam_niko/data/processed_bam/ERR523111/mark_dups/{samples}.Aligned.toTranscriptome.deduped.out.bam" 
    params:
        options = f"--CREATE_INDEX true --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER queryname --TMP_DIR {config['tmp_dir']} " 
    shell:
        "picard MarkDuplicates --INPUT {input.bam_in} {params.options} --OUTPUT {output.bam_out} --METRICS_FILE {output.metrics} "
       
#sort bam step
rule sort_sam:
    input:
        bam_in = "/pellmanlab/stam_niko/data/processed_bam/ERR523111/mark_dups/{samples}.Aligned.toTranscriptome.deduped.out.bam" #Input BAM or SAM file to sort.  Required
    output:
        bam_out =  "/pellmanlab/stam_niko/data/processed_bam/ERR523111/mark_dups/{samples}.Aligned.toTranscriptome.sorted.deduped.out.bam" #Sorted BAM or SAM output file.  Required.
    params:
        SO = "coordinate", #Sorts primarily according to the SEQ and POS fields of the record.
        options = "--CREATE_INDEX true"
    shell:
        "picard SortSam --INPUT {input.bam_in} {params.options} --OUTPUT {output.bam_out} --SORT_ORDER {params.SO} " 
       
#Prep bam for RSEM post-processing
rule prep_rsem:
    input:
        bam_in = "/pellmanlab/stam_niko/data/processed_bam/ERR523111/mark_dups/{samples}.Aligned.toTranscriptome.sorted.deduped.out.bam",
        fasta = config["reference_unzip"],
        gtf = config["reference_gtf_unzip"]
    output:
        mock = "/pellmanlab/stam_niko/data/processed_bam/ERR523111/RSEM/.{samples}_mockfile.rsem.txt",
    params:
        bam_out = "/pellmanlab/stam_niko/data/processed_bam/ERR523111/RSEM/{samples}.Aligned.toTranscriptome.sorted.deduped.rsem.out",
        transcriptomename = "/pellmanlab/stam_niko/data/processed_bam/ERR523111/RSEM/rsem_ref/Ensembl98",
        transcriptomedir = "/pellmanlab/stam_niko/data/processed_bam/ERR523111/RSEM/rsem_ref/"
    threads: 4
    shell:
        "convert-sam-for-rsem {input.bam_in} {params.bam_out} " #correct bam format for rsem
        "&& mkdir -p {params.transcriptomedir} "
        "&& rsem-prepare-reference -gtf {input.gtf} -p {threads} {input.fasta} {params.transcriptomename} " #needed rsem ref genome
        "&& touch {output.mock} " #For snakemake step tracing

#Run RSEM calculations
rule rsem_calc_expr:
    input:
        mock = "/pellmanlab/stam_niko/data/processed_bam/ERR523111/RSEM/.{samples}_mockfile.rsem.txt"
    output:
        mock = "/pellmanlab/stam_niko/data/processed_bam/ERR523111/RSEM/.{samples}_mockfile.rsem_calc.txt"
    params:
        bam_in = "/pellmanlab/stam_niko/data/processed_bam/ERR523111/RSEM/{samples}.Aligned.toTranscriptome.sorted.deduped.rsem.out.bam",
        ref_name = "/pellmanlab/stam_niko/data/processed_bam/ERR523111/RSEM/rsem_ref/Ensembl98",
        sample_name = "/pellmanlab/stam_niko/data/processed_bam/ERR523111/RSEM/outputs/{samples}",
        sample_dir = "/pellmanlab/stam_niko/data/processed_bam/ERR523111/RSEM/outputs/",
        options = "--paired-end --calc-pme",# --calc-ci --output-genome-bam ", calc-ci giving error maybe bug in version
        num_thrd = "--num-threads 8"
    shell:
        "mkdir -p {params.sample_dir} "
        "&& rsem-calculate-expression {params.options} --alignments {params.num_thrd} "
        "{params.bam_in} {params.ref_name} {params.sample_name} "
        "&& touch {output.mock} "

