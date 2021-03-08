###########################################################
####### Snakefile for scRNA-seq (Stamatis Specific) #######
###########################################################

#TODO: Fix config['MAX_THREADS'] issue where it is setting it to 1 rather than 16
#TODO: make the paths and locations modular

################
#### IMPORTS ###
################

import os
import pandas as pd
import snakemake

########################
### GLOBAL CONSTANTS ###
########################

#GATK_FILTERS = ("-RF MappingQualityReadFilter --minimum-mapping-quality 30 "
#                "-RF OverclippedReadFilter --filter-too-short 50")

#Threads from config
DEFAULT_THREADS = int(config["DEFAULT_THREADS"])
MAX_THREADS = int(config["MAX_THREADS"])

#Chrom Intervals
chrom_intervals = [f'-L chr{i}' for i in range(1, 23)]
alt_intervals = ['-L chrX', '-L chrY', '-L chrM']
chrom_intervals.extend(alt_intervals)
chrom_intervals = set(chrom_intervals)

#Chroms
chroms = [f'chr{i}' for i in range(1, 23)]
alt_chroms = ['chrX', 'chrY', 'chrM']
chroms.extend(alt_chroms)
chroms = set(chroms)

#Lanes
lane = ['Lane1', 'Lane2'] 

#Temp
tmp = "SN218_Run1065_Lane1_190717_Nextera_1E1_L_166282_BCStamatis_Nextera384_190718_P1_E1"

#Temp f experiment
#SIS1025f_samples = pd.read_table("samples/SIS1025f_samples.txt", header=0 )
#all_samples = set(SIS1025f_samples['samples'])

#Whole f experiment
SIS1025f_L1_samples = pd.read_table("samples/SIS1025f_L1_samples.txt", header=0 )
f_L1_samples = set(SIS1025f_L1_samples['samples'])

SIS1025f_L2_samples = pd.read_table("samples/SIS1025f_L2_samples.txt", header=0 )
f_L2_samples = set(SIS1025f_L2_samples['samples'])

#Whole e experiment
SIS1025e_samples = pd.read_table("samples/SIS1025e_samples.txt", header=0 )
e_samples = set(SIS1025e_samples['samples'])

#Whole d experiment
SIS1025d_samples = pd.read_table("samples/SIS1025d_samples.txt", header=0 )
d_samples = set(SIS1025d_samples['samples'])

#Whole b experiment
SIS1025b_samples = pd.read_table("samples/SIS1025b_samples.txt", header=0 )
b_samples = set(SIS1025b_samples['samples'])

#Whole a experiment
SIS1025a_samples = pd.read_table("samples/SIS1025a_samples.txt", header=0 )
a_samples = set(SIS1025a_samples['samples'])

experiments = [ "SIS1025a", "SIS1025b", "SIS1025d", "SIS1025e", "SIS1025f_Lane1", "SIS1025f_Lane2" ]
samples_set = [a_samples, b_samples, d_samples, e_samples, f_L1_samples, f_L2_samples]

#All samples
all_samples = []
all_samples.extend(a_samples)
all_samples.extend(b_samples)
all_samples.extend(d_samples)
all_samples.extend(e_samples)
all_samples.extend(f_L1_samples)
all_samples.extend(f_L2_samples)

########################
### INPUT-ONLY RULES ###
########################
#rule index_genome:
#    input:
#       : starindex = "/pellmanlab/stam_niko/refgenomes/STAR/Gencode.v25/", #SAindex",
#        rsemindex = "/pellmanlab/stam_niko/refgenomes/RSEM/Gencode.v25/" #genecode.v25"

#Add f'...{experimentf}...' to all these string paths
#rule all:
#    input:
#        expand("/pellmanlab/stam_niko/data/processed_bam/SIS1025f/STAR/.{samples}_mockfile.txt", samples=all_samples)

rule align:
    input:
        expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/STAR/.{samples}_mockfile.txt", experiment='SIS1025f_Lane1', samples=f_L1_samples),
        expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/STAR/.{samples}_mockfile.txt", experiment='SIS1025f_Lane2', samples=f_L2_samples),
        expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/STAR/.{samples}_mockfile.txt", experiment='SIS1025e', samples=e_samples),
        #expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/STAR/.{samples}_mockfile.txt", experiment='SIS1025d', samples=d_samples),
        #expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/STAR/.{samples}_mockfile.txt", experiment='SIS1025b', samples=b_samples),
        #expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/STAR/.{samples}_mockfile.txt", experiment='SIS1025a', samples=a_samples)

rule calc_gene_expr:
    input:        
        expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/RSEM/output/.{samples}_mockfile.rsem_calc.txt", experiment='SIS1025f_Lane1', samples=f_L1_samples),
        expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/RSEM/output/.{samples}_mockfile.rsem_calc.txt", experiment='SIS1025f_Lane2', samples=f_L2_samples),
        expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/RSEM/output/.{samples}_mockfile.rsem_calc.txt", experiment='SIS1025e', samples=e_samples),
        #expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/RSEM/output/.{samples}_mockfile.rsem_calc.txt", experiment='SIS1025d', samples=d_samples),
        #expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/RSEM/output/.{samples}_mockfile.rsem_calc.txt", experiment='SIS1025b', samples=b_samples),
        #expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/RSEM/output/.{samples}_mockfile.rsem_calc.txt", experiment='SIS1025a', samples=a_samples)

rule calc_metrics:
    input:
        expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/Metrics/{samples}.MultipleMetrics.alignment_summary_metrics", experiment='SIS1025f_Lane1', samples=f_L1_samples),
        #expand("/pellmanlab/stam_niko/data/processed_bam/SIS1025f/Metrics/{samples}.MultipleMetrics.alignment_summary_metrics", samples=all_samples)

rule change_rgtags:
    input:
        expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/VarCall_BAMs/{samples}.Aligned.sortedByCoord.split_r.RG.out.bam.bai", experiment='SIS1025f_Lane1', samples=tmp),
        #expand("/pellmanlab/stam_niko/data/processed_bam/SIS1025f/VarCall_BAMs/{samples}.Aligned.sortedByCoord.split_r.RG.out.bam.bai", samples=all_samples)

#rule calc_variants:
#    input:
#        expand("/pellmanlab/stam_niko/data/processed_bam/SIS1025f/Variants/{samples}.g.vcf.gz", samples=all_samples)
	#expand("/pellmanlab/stam_niko/data/processed_bam/SIS1025f/Variants/SIS1025f.{chrom}.vcf.gz", chrom=chrom_intervals)

#rule calc_gvcf:
#    input:
#        expand("/pellmanlab/stam_niko/data/processed_bam/SIS1025f/Variants/SIS1025f_{chrom}_allsamples.g.vcf.gz", chrom=chroms)
        
#rule calc_vcf:
#    input:
#        expand("/pellmanlab/stam_niko/data/processed_bam/SIS1025f/Variants/SIS1025f_{chrom}_allsamples.vcf.gz", chrom=chroms)

#rule gen_vcf_jobs:
#    input:
#        expand("/pellmanlab/stam_niko/data/processed_bam/SIS1025f/Variants/mockfile.gengenotypecmds.{lanes}.txt", lanes=lane)


#rule run_vcf_jobs:
#    input:
#        expand("/pellmanlab/stam_niko/data/processed_bam/SIS1025f/Variants/.mockfile.rungenotypecmds.{lanes}.txt", lanes=lane)

rule track_data:
    input:
        expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/Analysis/.mockfile.{experiment}.data_tracking.txt", experiment='SIS1025e', samples=e_samples)

rule run_check_data:
    input:
        [expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/Analysis/.mockfile.{experiment}.data_tracking.txt", experiment=experiments[i]) for i in range(len(experiments))]

rule run_data_aggregation:
    input:
        expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/Analysis/.mockfile.{experiment}.data_aggregation.txt", experiment=experiments)

rule run_data_aggregation_macro:
    input:
        "/pellmanlab/stam_niko/data/processed_bam/aggregated_results/.mockfile.data_aggregation_macro.txt"

rule run_ML:
    input:
        "/pellmanlab/stam_niko/data/processed_bam/aggregated_results/.mockfile.mlscript.txt"

rule run_data_visualization:
    input:
        "/pellmanlab/stam_niko/data/processed_bam/visual_results/.mockfile.visuals.txt"


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
        transcriptomename = "/pellmanlab/stam_niko/refgenomes/RSEM/Gencode.v25/genecode.v25.", #remove end .
    threads: config["MAX_THREADS"]
    shell:
        "rsem-prepare-reference -gtf {input.gtf} -p {threads} --star {input.fasta} {params.transcriptomename} "
        
#####################
### PREPROCESSING ###
#####################

#SN218_Run1065_Lane1_190717_Nextera_1A10_L_166350_BCStamatis_Nextera384_190718_P1_A10.demult.bam.qsort.bam.R1.fastq.gz
#mapping, aligning, tagging and sorting step
rule STAR_alignment:
    input:
        R1 = "/pellmanlab/stam_niko/data/all_fastqs/{samples}.demult.bam.qsort.bam.R1.fastq.gz", #Chance to simplify name here
        R2 = "/pellmanlab/stam_niko/data/all_fastqs/{samples}.demult.bam.qsort.bam.R2.fastq.gz",
        ref_dir =  "/pellmanlab/stam_niko/refgenomes/STAR/Gencode.v25/gencode.v25",
        gtf = "/pellmanlab/stam_niko/refgenomes/Gencode/v25/gencode.v25.primary_assembly.annotation.gtf"
    output: #get rid of mockfile by putting one of the output files here and keeping names where it is.
        mock = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/STAR/.{samples}_mockfile.txt" 
    params:
        names = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/STAR/{samples}.",
        sample = "{samples}"
#    threads: 4 #config["MAX_THREADS"]
    shell:
        """
        STAR \
        --runMode alignReads \
        --runThreadN 1 \
        --genomeDir {input.ref_dir} \
        --sjdbGTFfile {input.gtf} \
        --twopassMode Basic \
        --alignSJDBoverhangMin 2 \
        --outFilterMismatchNoverLmax 0.1 \
        --outFilterType BySJout \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMatchNminOverLread 0.33 \
        --limitSjdbInsertNsj 1200000 \
        --limitBAMsortRAM 171798691840 \
        --readFilesIn {input.R1} {input.R2} \
        --readFilesCommand zcat \
        --outFileNamePrefix {params.names} \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs None \
        --alignSoftClipAtReferenceEnds No \
        --quantMode TranscriptomeSAM GeneCounts \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --genomeLoad NoSharedMemory \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --outSAMattributes NH HI AS nM NM MD \
        --outSAMattrRGline ID:{params.sample} SM:{params.sample} \
        && touch {output.mock} 
        """ 
        #IDEA: instead of touching mockfile, && index and use the index for the snakemake file tracing
        #--limitBAMsortRAM 171,798,691,840       
 
#MARK DUPLICATES IS CONTRAINDICATED BECAUSE DUPLICATION RATE WAS <5%.
       
########################
### GENE EXPRESSSION ###
########################       
       
#Prep bam for RSEM post-processing
rule rsem_prep_bam:
    input:
        mock_in = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/STAR/.{samples}_mockfile.txt"
    output:
        mock = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/RSEM/{samples}.Aligned.toTranscriptome.rsem.mock.txt"
    params:
        bam_in = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/STAR/{samples}.Aligned.toTranscriptome.out.bam",
        bam_out_name = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/RSEM/{samples}.Aligned.toTranscriptome.rsem.out"
    #threads: 4
    shell:
        "convert-sam-for-rsem {params.bam_in} {params.bam_out_name} " #correct bam format for rsem
        "&& touch {output.mock} "

#Run RSEM calculations
rule rsem_calc_expr:
    input:
        mock_in = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/RSEM/{samples}.Aligned.toTranscriptome.rsem.mock.txt"
    output:
        mock = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/RSEM/output/.{samples}_mockfile.rsem_calc.txt"
    params:
        bam_in = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/RSEM/{samples}.Aligned.toTranscriptome.rsem.out.bam",
        ref_name = "/pellmanlab/stam_niko/refgenomes/RSEM/Gencode.v25/genecode.v25.", #remove end .
        sample_name = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/RSEM/output/{samples}",
        options = "--paired-end --no-bam-output --estimate-rspd ", #--calc-pme", #works for 99.9% of samples but some are too small for pme # --calc-ci ", calc-ci giving error maybe bug in version
#        tmp = "--temporary-folder /pellmanlab/stam_niko/data/tmp/"
    #threads: 4
    shell:
        "rsem-calculate-expression {params.options} --alignments " # --num-threads {threads} "
        "{params.bam_in} {params.ref_name} {params.sample_name} "
        "&& touch {output.mock} "
 
######################
### METRIC CALLING ###
###################### 

#Collect multiple metrics
rule collect_mult_metrics:
    input:
        mock = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/STAR/.{samples}_mockfile.txt"
    output:
        metrics = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/Metrics/{samples}.MultipleMetrics.alignment_summary_metrics"
    params:
        ref_fasta = config["reference_unzip"],
        bam_in = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/STAR/{samples}.Aligned.sortedByCoord.out.bam",
        metrics = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/Metrics/{samples}.MultipleMetrics"
    shell:
        "picard CollectMultipleMetrics I={params.bam_in} O={params.metrics} R={params.ref_fasta} "
        
####################################
### VARIANT CALLING | GENOTYPING ###
####################################

rule splitncigarreads:
    input:
        mock_in = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/STAR/.{samples}_mockfile.txt",
        ref = config["ref_unzip_wdict"] 
    output:
        bam_out = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/VarCall_BAMs/{samples}.Aligned.sortedByCoord.split_r.out.bam" 
    params:
        bam_in =  "/pellmanlab/stam_niko/data/processed_bam/{experiment}/STAR/{samples}.Aligned.sortedByCoord.out.bam",
        options = "--create-output-bam-index true", # --read-filter ReassignOneMappingQuality", # -RMQF 255 -RMQT 60",
        intervals = chrom_intervals
    shell:
        "samtools index {params.bam_in} "
        "&& gatk SplitNCigarReads {params.options} {params.intervals} --input {params.bam_in} --reference {input.ref} --output {output.bam_out} "

#add read group information and calc tags
rule AddOrReplaceRG:
    input:
        bam = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/VarCall_BAMs/{samples}.Aligned.sortedByCoord.split_r.out.bam"
    output:
        bam_out = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/VarCall_BAMs/{samples}.Aligned.sortedByCoord.split_r.RG.out.bam"
    params:
        SM = "{samples}".split('.')[0],
    shell:
        "picard AddOrReplaceReadGroups -I {input.bam} -O {output.bam_out} "
        "--RGID {params.SM} --RGLB {params.SM} --RGPL illumina --RGSM {params.SM} --RGPU MiOrHiSeq"

rule index_RG_bams:
    input:
        bam = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/VarCall_BAMs/{samples}.Aligned.sortedByCoord.split_r.RG.out.bam"
    output:
        bam_out = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/VarCall_BAMs/{samples}.Aligned.sortedByCoord.split_r.RG.out.bam.bai"
    shell:
        "samtools index {input.bam} "

#rule haplotype_variant_calling:
#    input:
#        bam_in = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/VarCall_BAMs/{samples}.Aligned.sortedByCoord.split_r.RG.out.bam", 
#        alleles = config["alleles"],
#        fasta = config["ref_unzip_wdict"]
#    output:
#        gvcf = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/Variants/{samples}.g.vcf.gz"
#    params:
	#I need to address the -rf removals used in etai's code
#        options = "-ERC GVCF --standard-min-confidence-threshold-for-calling 0 --minimum-mapping-quality 30"
#    threads: config["MAX_THREADS"]
#    shell:
#        "gatk HaplotypeCaller {params.options} --native-pair-hmm-threads 20 --reference {input.fasta} --alleles {input.alleles} -I {input.bam_in} -O {output.gvcf}"

#rule GVCF_combine_bychr:
#    input:
#        fa = config["ref_unzip_wdict"]
#    output:
#        gvcf = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/Variants/SIS1025f_{chrom}_allsamples.g.vcf.gz"
#    params:
#        gvcfs = [f'-V /pellmanlab/stam_niko/data/processed_bam/SIS1025f/Variants/{i}.g.vcf.gz' for i in all_samples],
#        interval = "-L {chrom}"
#    shell:
#        "gatk CombineGVCFs -R {input.fa} {params.interval} {params.gvcfs} -O {output.gvcf} "

#rule get_vcf_bychr:
#    input:
#        gvcf = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/Variants/SIS1025f_{chrom}_allsamples.g.vcf.gz",
#        fa = config["ref_unzip_wdict"]
#    output:
#        vcf = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/Variants/SIS1025f_{chrom}_allsamples.vcf.gz"
#    params:
#        interval = "-L {chrom}"
#    shell:
#        "gatk GenotypeGVCFs {params.interval} -R {input.fa} -V {input.gvcf} -O {output.vcf} "

#rule count_reads_allelic:
#    input:
#        ref = config["ref_unzip_wdict"],
#        hets = config["alleles"],
#        bam = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/VarCall_BAMs/{samples}.Aligned.sortedByCoord.split_r.out.bam"
#    output:
#        vcf = "/pellmanlab/stam_niko/data/processed_bam/SIS1025f/Variants/{samples}.AD.tsv"
#    shell:
#        "gatk ASEReadCounter -I {input.bam} -V {input.hets} "
#        "-R {input.ref} -O {output.vcf}"

rule etai_genotype_generate_call:
    input:
        [expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/VarCall_BAMs/{samples}.Aligned.sortedByCoord.split_r.RG.out.bam.bai", experiment=experiments[i], samples=samples_set[i]) for i in range(len(experiments))]
    output:
        mockfile = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/Variants/.mockfile.{experiment}.gengenotypecmds.txt"
    params:
        ls = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/VarCall_BAMs/*RG.out.bam",
        bamlist = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/VarCall_BAMs/{experiment}.bamfiles.list",
        script = "/pellmanlab/nikos/Stam_Etai_Scripts/scripts/RPE-1_GRCh38_Genotype_etai.sh",
        out_name = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/Variants/{experiment}"
    shell:
        "ls {params.ls} > {params.bamlist} "
        "&& bash {params.script} {params.bamlist} {params.out_name} 1 "
        "&& touch {output.mockfile}"

rule etai_genotype_run:
    input:
        "/pellmanlab/stam_niko/data/processed_bam/{experiment}/Variants/.mockfile.{experiment}.gengenotypecmds.txt",
    output:
        "/pellmanlab/stam_niko/data/processed_bam/{experiment}/Variants/.mockfile.{experiment}.rungenotypecmds.txt"
    params:
        "/pellmanlab/stam_niko/data/processed_bam/{experiment}/Variants/{experiment}_RPE_hets_GT.UG_jobs.list"
    threads: 25 # didnt put {threads} below because I want the jobs to always be 25 for this part no matter what I put in cmd line
    shell:
        "parallel --jobs {threads} < {params} "
        "&& touch {output}"

####################################
### DATA AGGREGATION | ANALYSIS  ###
####################################

rule check_data:
    input:
        [expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/RSEM/output/.{samples}_mockfile.rsem_calc.txt", experiment=experiments[i], samples=samples_set[i]) for i in range(len(experiments))],
        [expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/Metrics/{samples}.MultipleMetrics.alignment_summary_metrics", experiment=experiments[i], samples=samples_set[i]) for i in range(len(experiments))],
        [expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/Variants/.mockfile.{experiment}.rungenotypecmds.txt", experiment=experiments[i]) for i in range(len(experiments))],
    output:
        mock_out = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/Analysis/.mockfile.{experiment}.data_tracking.txt"
    shell:
        "touch {output.mock_out}"

rule aggregate_data:
    input:
        mock_in = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/Analysis/.mockfile.{experiment}.data_tracking.txt",
    output:
        mock_out = "/pellmanlab/stam_niko/data/processed_bam/{experiment}/Analysis/.mockfile.{experiment}.data_aggregation.txt",
    params:
       script_dir = "/pellmanlab/nikos/Stam_Etai_Scripts",
       expr = "{experiment}",
       wk_dir = "/pellmanlab/stam_niko/data/processed_bam",
       script = "/pellmanlab/nikos/Stam_Etai_Scripts/scripts/Analysis_Etai_Snakemake.R",
       datadir = "/pellmanlab/nikos/Stam_Etai_Data"
    threads: 1 #13 I dont think the multiple cores works
    shell:
        "Rscript {params.script} {params.script_dir} {params.wk_dir} {params.datadir} {params.expr} {threads} "
        "&& touch {output.mock_out} "


rule aggregate_data_macro_and_models:
    input:
        mock_in = [expand("/pellmanlab/stam_niko/data/processed_bam/{experiment}/Analysis/.mockfile.{experiment}.data_aggregation.txt", experiment=experiments[i]) for i in range(len(experiments))],
    output:
        mock_out = "/pellmanlab/stam_niko/data/processed_bam/aggregated_results/.mockfile.data_aggregation_macro.txt", 
    params:
       script_dir = "/pellmanlab/nikos/Stam_Etai_Scripts",
       all_experiments = [experiments[i] for i in range(len(experiments))],
       wk_dir = "/pellmanlab/stam_niko/data/processed_bam",
       script = "/pellmanlab/nikos/Stam_Etai_Scripts/scripts/Data_Aggregation.R",
       datadir = "/pellmanlab/nikos/Stam_Etai_Data"
    shell:
        "Rscript {params.script} {params.script_dir} {params.wk_dir} {params.datadir} {params.all_experiments} "
        "&& touch {output.mock_out} "


rule machine_learning_model:
    input:
        mock_in = "/pellmanlab/stam_niko/data/processed_bam/aggregated_results/.mockfile.data_aggregation_macro.txt",
    output:
        mock_out = "/pellmanlab/stam_niko/data/processed_bam/aggregated_results/.mockfile.mlscript.txt",
    params:
       script_dir = "/pellmanlab/nikos/Stam_Etai_Scripts",
       wk_dir = "/pellmanlab/stam_niko/data/processed_bam",
       script = "/pellmanlab/nikos/Stam_Etai_Scripts/ML/Run_ML.R",
       datadir = "/pellmanlab/nikos/Stam_Etai_Data"
    shell:
        "Rscript {params.script} {params.script_dir} {params.wk_dir} {params.datadir} "
        "&& touch {output.mock_out} "

###########################
### DATA VISUALIZATION  ###
###########################

rule visualization_script:
    input:
        mock_in = "/pellmanlab/stam_niko/data/processed_bam/aggregated_results/.mockfile.mlscript.txt",
    output:
        mock_out = "/pellmanlab/stam_niko/data/processed_bam/visual_results/.mockfile.visuals.txt",
    params:
       script_dir = "/pellmanlab/nikos/Stam_Etai_Scripts",
       wk_dir = "/pellmanlab/stam_niko/data/processed_bam",
       script = "/pellmanlab/nikos/Stam_Etai_Scripts/plots/visual_generator_Nikos.R",
       datadir = "/pellmanlab/nikos/Stam_Etai_Data"
    shell:
        "Rscript {params.script} {params.script_dir} {params.wk_dir} {params.datadir}"
        "&& touch {output.mock_out} "

