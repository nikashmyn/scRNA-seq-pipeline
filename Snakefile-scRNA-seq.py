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
prac_samples_R1 = list(["ERR523111_1"]) #idea: just use names and append the _1 and _2 in the calls when scaling
prac_samples_R2 = list(["ERR523111_2"])
########################
### INPUT-ONLY RULES ###
########################
rule align:
    input:
        expand("/pellmanlab/stam_niko/data/bam/ERR523111/.{samples}_mockfile.txt", samples=prac_samples)

rule mrk_dups:
    input:
        expand("/pellmanlab/stam_niko/data/bam/ERR523111/mark_dups/{samples}.Aligned.sortedByCoord.deduped.out.bam", samples=prac_samples)

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
        mock = temp("/pellmanlab/stam_niko/data/bam/ERR523111/.{samples}_mockfile.txt")
    params:
        names = "/pellmanlab/stam_niko/data/bam/ERR523111/{samples}.",
        runMode = "alignReads",
        outstyle = "BAM SortedByCoordinate",
        readcmd = "zcat"
    shell:
        "STAR --genomeDir {input.ref_dir} --readFilesIn {input.R1} {input.R2} "
        "--runMode {params.runMode} --outSAMtype {params.outstyle} --readFilesCommand {params.readcmd} "
        "--outFileNamePrefix {params.names} "
        "&& touch {output.mock} "

#remove duplicates step        
rule mark_duplicates:
    input:
        mock = "/pellmanlab/stam_niko/data/bam/ERR523111/.{samples}_mockfile.txt",
    output:
        metrics = "/pellmanlab/stam_niko/data/bam/ERR523111/mark_dups/{samples}.metrics.txt", #File to write duplication metrics to  Required.
        bam_out = "/pellmanlab/stam_niko/data/bam/ERR523111/mark_dups/{samples}.Aligned.sortedByCoord.deduped.out.bam" #The output file to write marked records to  Required
    params:
        bam_in = "/pellmanlab/stam_niko/data/bam/ERR523111/{samples}.Aligned.sortedByCoord.out.bam", #In params because output of last step was a mock # Must be coordinate sorted.
        options = "--CREATE_INDEX True --REMOVE_DUPLICATES True --TMP_DIR config['tmp_dir'] "
    shell:
        "picard MarkDuplicates --INPUT {params.bam_in} {params.options} --OUTPUT {output.bam_out} --METRICS_FILE {output.metrics} "
       

#sort bam step
#rule sort_sam:
#    input:
#        bam_in = #Input BAM or SAM file to sort.  Required
#    output:
#        SO = "coordinate", #Sorts primarily according to the SEQ and POS fields of the record.
#        bam_out = #Sorted BAM or SAM output file.  Required.
#    params:
#        "--CREATE_INDEX True"
#    shell:
#        "picard SortSam --INPUT {input.bam_in} {params} --OUTPUT {output.bam_out} --SORT_ORDER {output.SO} " 
        

# Split Reads with N in Cigar
#rule splitncigarreads:
#    input:
#        bam_in = ,#BAM/SAM/CRAM file containing reads  This argument must be specified at least once. Required.   
#        ref = #Reference sequence file  Required.                         
#    output:
#        bam_out = #Write output to this BAM filename  Required. 
#    params:
#        "--create-output-bam-index True",
#        "--tmp-dir  !!!GATKPath ",
#        "--minimum-mapping-quality 20"
#    shell:
#        "gatk SplitNCigarReads --input {input.bam_in} --reference {input.ref} {params} --output {output.bam_out} "
#        
# Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context)
#rule base_recal:
#    input:
#        bam_in = ,#BAM/SAM/CRAM file containing reads  This argument must be specified at least once. Required.
#        ref = ,#Reference sequence file  Required.
#        KS = #One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis.  This argument must be specified at least once. Required.
#    output:
#        bam_out = #The output recalibration table file to create  Required.
#    params:
#        "--create-output-bam-index True",
#        "--read-filter PairedReadFilter",
#        "--tmp-dir !! Temp directory to use.  Default value: null"
#    shell:
#        "gatk BaseRecalibrator --input {input.bam_in} --reference {input.ref} --known-sites {input.KS} {params} --output {output.bam_out} "

