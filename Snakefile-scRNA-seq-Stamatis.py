###########################################################
####### Snakefile for scRNA-seq (Stamatis Specific) #######
###########################################################

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
real_samples = list(["SN218_Run1065_Lane1_190717_Nextera_1A10_L_166350_BCStamatis_Nextera384_190718_P1_A10.demult.bam.qsort.bam"])

########################
### INPUT-ONLY RULES ###
########################
rule index_genome:
    input:
        starindex = "/pellmanlab/stam_niko/refgenomes/STAR/Gencode.v25/", #SAindex",
        rsemindex = "/pellmanlab/stam_niko/refgenomes/RSEM/Gencode.v25/" #genecode.v25"

rule align:
    input:
        expand("/pellmanlab/stam_niko/data/processed_bam/SIS1025f/STAR/.{samples}_mockfile.txt", samples=real_samples)

rule calc_gene_expr:
    input:
        expand("/pellmanlab/stam_niko/data/processed_bam/SIS1025f/RSEM/.{samples}_mockfile.rsem_calc.txt", samples=real_samples)

rule calc_metrics:
    input:
        expand("/pellmanlab/stam_niko/data/processed_bam/SIS1025f/Metrics/{samples}.MultipleMetrics.alignment_summary_metrics", samples=real_samples)
        
#######################
### GENOME INDEXING ###
#######################

#premapping and aligning step for STAR
rule generate_genome_indexes:
    input:
        in_fasta = config["reference_unzip"], #CHANGE TO V25
        in_gtf = config["reference_gtf_unzip"] #CHANGE TO V25
    output:
        ref_dir = "/pellmanlab/stam_niko/refgenomes/STAR/Gencode.v25/gencode.v25" #CHANGE TO V25
    params:
        runMode = "genomeGenerate",
        overhang = config["readlength"] - 1
    threads: config["MAX_THREADS"]
    shell:
        "STAR --genomeFastaFiles {input.in_fasta} --sjdbGTFfile {input.in_gtf} "
        "--runThreadN {threads} --runMode {params.runMode} --sjdbOverhang {params.overhang} "
        "--genomeDir {output.ref_dir} "
        
#Prep bam for RSEM post-processing
rule prep_rsem:
    input:
        fasta = config["reference_unzip"],
        gtf = config["reference_gtf_unzip"]
    output:
        transcriptomedir = "/pellmanlab/stam_niko/refgenomes/RSEM/Gencode.v25/"
    params:
        transcriptomename = "/pellmanlab/stam_niko/refgenomes/RSEM/Gencode.v25/genecode.v25",
    threads: config["MAX_THREADS"]
    shell:
        "rsem-prepare-reference -gtf {input.gtf} -p {threads} --star {input.fasta} {params.transcriptomename} "
        
#####################
### PREPROCESSING ###
#####################

#mapping, aligning, tagging and sorting step
rule STAR_alignment:
    input:
        R1 = "/singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025f/fastqs/{samples}.R1.fastq.gz", #Chance to simplify name here
        R2 = "/singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025f/fastqs/{samples}.R2.fastq.gz",
        ref_dir =  "/pellmanlab/stam_niko/refgenomes/STAR/Gencode.v25/gencode.v25",
        gtf = "/pellmanlab/stam_niko/refgenomes/Gencode/v25/gencode.v25.primary_assembly.annotation.gtf"
    output: #get rid of mockfile by putting one of the output files here and keeping names where it is.
        mock = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/STAR/.{samples}_mockfile.txt" 
    params:
        names = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/STAR/{samples}."
    threads: config["MAX_THREADS"]
    shell:
        """
        STAR \
        --runMode alignReads \
        --runThreadN {threads} \
        --genomeDir {input.ref_dir} \
        --sjdbGTFfile {input.gtf} \
        --twopassMode Basic \
        --alignSJDBoverhangMin 2 \
        --outFilterMismatchNoverLmax 0.1 \
        --outFilterType BySJout \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMatchNminOverLread 0.33 \
        --limitSjdbInsertNsj 1200000 \
        --readFilesIn {input.R1} {input.R2} \
        --readFilesCommand zcat \
        --outFileNamePrefix {params.names} \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs None \
        --alignSoftClipAtReferenceEnds No \
        --quantMode TranscriptomeSAM GeneCounts \
        --outSAMtype BAM SortedByCoordinate \
        --limitBAMsortRAM 171798691840 \
        --outSAMunmapped Within \
        --genomeLoad NoSharedMemory \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --outSAMattributes NH HI AS nM NM MD \
        --outSAMattrRGline ID:rg1 SM:sm1 \
        && touch {output.mock}
        """
        
#MARK DUPLICATES IS CONTRAINDICATED BECAUSE DUPLICATION RATE WAS <5%.
       
########################
### GENE EXPRESSSION ###
########################       
       
#Prep bam for RSEM post-processing
rule rsem_prep_bam:
    input:
        bam_in = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/STAR/{samples}.Aligned.toTranscriptome.out.bam"
    output:
        bam_out = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/RSEM/{samples}.Aligned.toTranscriptome.rsem.out"
    threads: config["MAX_THREADS"]
    shell:
        "convert-sam-for-rsem {input.bam_in} {params.bam_out} " #correct bam format for rsem


#Run RSEM calculations
rule rsem_calc_expr:
    input:
        bam_in = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/RSEM/{samples}.Aligned.toTranscriptome.rsem.out"
    output:
        mock = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/RSEM/.{samples}_mockfile.rsem_calc.txt"
    params:
        ref_name = "/pellmanlab/stam_niko/refgenomes/RSEM/Gencode.v25/genecode.v25",
        sample_name = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/RSEM/outputs/{samples}",
        sample_dir = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/RSEM/outputs/",
        options = "--paired-end --no-bam-output --calc-pme", # --calc-ci ", calc-ci giving error maybe bug in version
        num_thrd = "--num-threads 8"
    shell:
        "mkdir -p {params.sample_dir} "
        "&& rsem-calculate-expression {params.options} --alignments {params.num_thrd} "
        "{input.bam_in} {params.ref_name} {params.sample_name} "
        "&& touch {output.mock} "
 
######################
### METRIC CALLING ###
###################### 

#Collect multiple metrics
rule collect_mult_metrics:
    input:
        mock = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/STAR/.{samples}_mockfile.txt"
    output:
        metrics = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/Metrics/{samples}.MultipleMetrics.alignment_summary_metrics"
    params:
        ref_fasta = config["reference_unzip"],
        bam_in = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/STAR/{samples}.Aligned.sortedByCoord.out.bam",
        metrics = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/Metrics/{samples}.MultipleMetrics"
    shell:
        "picard CollectMultipleMetrics I={params.bam_in} O={params.metrics} R={params.ref_fasta} "
        
####################################
### VARIANT CALLING | GENOTYPING ###
####################################

rule splitncigarreads:
    input:
        bam_in =  "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/STAR/{samples}.Aligned.toTranscriptome.out.bam", 
        ref = config["reference_unzip"] 
    output:
        bam_out = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/VarCall_BAMs/{samples}.Aligned.toTranscriptome.split_r.out.bam" 
    params:
        "--create-output-bam-index true"
    shell:
        "gatk SplitNCigarReads --input {input.bam_in} --reference {input.ref} {params} --output {output.bam_out} "

rule haplotype_variant_calling:
    input:
        bam_in = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/VarCall_BAMs/{samples}.Aligned.toTranscriptome.split_r.out.bam", 
        alleles = config["alleles"],
        fasta = config["reference_unzip"]
    output:
        vcf = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/Variants/{samples}.vcf"
    params:
        rfs = "--num_cpu_threads_per_data_thread {threads} --output_mode EMIT_ALL_SITES --standard_min_confidence_threshold_for_calling 0 -rf DuplicateRead -rf FailsVendorQualityCheck -rf NotPrimaryAlignment -rf BadMate -rf MappingQualityUnavailable -rf UnmappedRead -rf BadCigar",
        options = "--genotyping_mode GENOTYPE_GIVEN_ALLELES --analysis_type HaplotypeCaller --min_mapping_quality_score 30"
    threads: config["MAX_THREADS"]
    shell:
        "gatk --reference_sequence {input.fasta} {params.rfs} {params.options} {input.bam_in} --alleles {input.alleles} "

