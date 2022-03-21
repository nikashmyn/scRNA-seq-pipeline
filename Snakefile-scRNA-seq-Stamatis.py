###########################################################
####### Snakefile for scRNA-seq (Stamatis Specific) #######
###########################################################

os.chdir(config["OUTDIR"]) #change the working directory to the outdir from config
skdir = config["SNAKEDIR"] #pull locations of the git clone from config 
OUTDIR = f'{config["OUTDIR"]}/data' 
datadir = config["DATADIR"]

################
#### IMPORTS ###
################

import os
import pandas as pd
import snakemake

########################
### GLOBAL CONSTANTS ###
########################

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

#Whole a experiment
SIS1025a_samples = pd.read_table(config["a_samples"], header=0 )
a_samples = set(SIS1025a_samples['samples'])

#Whole b experiment
SIS1025b_samples = pd.read_table(config["b_samples"], header=0 )
b_samples = set(SIS1025b_samples['samples'])

#Whole c experiment
SIS1025c_samples = pd.read_table(config["c_samples"], header=0 )
c_samples = set(SIS1025c_samples['samples'])

#Whole d experiment
SIS1025d_samples = pd.read_table(config["d_samples"], header=0 )
d_samples = set(SIS1025d_samples['samples'])

#Whole e experiment
SIS1025e_samples = pd.read_table(config["e_samples"], header=0 )
e_samples = set(SIS1025e_samples['samples'])

#Whole f experiment
SIS1025f_L1_samples = pd.read_table(config["f1_samples"], header=0 )
f_L1_samples = set(SIS1025f_L1_samples['samples'])

SIS1025f_L2_samples = pd.read_table(config["f2_samples"], header=0 )
f_L2_samples = set(SIS1025f_L2_samples['samples'])

#Whole g experiment
SIS1025g_L1_samples = pd.read_table(config["g1_samples"], header=0 )
g_L1_samples = set(SIS1025g_L1_samples['samples'])

SIS1025g_L2_samples = pd.read_table(config["g2_samples"], header=0 )
g_L2_samples = set(SIS1025g_L2_samples['samples'])

#misc early experiments
SIS1025_misc_samples = pd.read_table(config["misc_samples"], header=0)
misc_samples = set(SIS1025_misc_samples['samples'])

#targeted MN experiment
SIS1025_targ_samples = pd.read_table(config["targ_samples"], header=0)
targ_samples = set(SIS1025_targ_samples['samples'])

#list of experiments sets
experiments = [ "SIS1025a", "SIS1025b", "SIS1025c", "SIS1025d", "SIS1025e", "SIS1025f_Lane1", "SIS1025f_Lane2", "SIS1025g_Lane1", "SIS1025g_Lane2", "SIS1025misc", "SIS1025targ" ]
samples_set = [a_samples, b_samples, c_samples, d_samples, e_samples, f_L1_samples, f_L2_samples, g_L1_samples, g_L2_samples, misc_samples, targ_samples]

#All samples
all_samples = []
all_samples.extend(a_samples)
all_samples.extend(b_samples)
all_samples.extend(c_samples)
all_samples.extend(d_samples)
all_samples.extend(e_samples)
all_samples.extend(f_L1_samples)
all_samples.extend(f_L2_samples)
all_samples.extend(g_L1_samples)
all_samples.extend(g_L2_samples)
all_samples.extend(misc_samples)
all_samples.extend(targ_samples)

########################
### INPUT-ONLY RULES ###
########################
#rule run_pipeline:
#    input:
#        expand("{path}/visual_results/.mockfile.visuals.txt", path=OUTDIR)

rule align:
    input:
        [expand("{path}/{experiment}/STAR/.{samples}_mockfile.txt", path=OUTDIR, experiment=experiments[i], samples=samples_set[i]) for i in range(len(experiments))]

rule calc_gene_expr:
    input:
        [expand("{path}/{experiment}/RSEM/output/.{samples}_mockfile.rsem_calc.txt", path=OUTDIR, experiment=experiments[i], samples=samples_set[i]) for i in range(len(experiments))]

rule calc_metrics:
    input:
        [expand("{path}/{experiment}/Metrics/{samples}.MultipleMetrics.alignment_summary_metrics", path=OUTDIR, experiment=experiments[i], samples=samples_set[i]) for i in range(len(experiments))]

rule change_rgtags:
    input:
        [expand("{path}/{experiment}/VarCall_BAMs/{samples}.Aligned.sortedByCoord.split_r.RG.out.bam.bai", path=OUTDIR, experiment=experiments[i], samples=samples_set[i]) for i in range(len(experiments))]

rule run_ASE:
    input:
        [expand("{path}/{experiment}/ASE/{samples}.ASE.table", path=OUTDIR, experiment=experiments[i], samples=samples_set[i]) for i in range(len(experiments))]

rule run_variant_calling:
    input:
        [expand("{path}/{experiment}/Variants/ASE/{samples}.vcf.gz", path=OUTDIR, experiment=experiments[i], samples=samples_set[i]) for i in range(len(experiments))]

rule gen_call_variant:
    input:
        [expand("{path}/{experiment}/Variants/.mockfile.{experiment}.gengenotypecmds.txt", path=OUTDIR, experiment=experiments[i]) for i in range(len(experiments))]

rule run_call_variant:
    input:
        expand("{path}/{experiment}/Variants/.mockfile.{experiment}.rungenotypecmds.{chrs}.txt", path=OUTDIR, experiment=experiments, chrs=chroms)

rule run_check_data:
    input:
        [expand("{path}/{experiment}/Analysis/.mockfile.{experiment}.data_tracking.txt", path=OUTDIR, experiment=experiments[i]) for i in range(len(experiments))]

rule run_data_aggregation:
    input:
        expand("{path}/{experiment}/Analysis/.mockfile.{experiment}.data_aggregation.txt", path=OUTDIR, experiment=experiments)

rule run_data_aggregation_macro:
    input:
        expand("{path}/aggregated_results/.mockfile.data_aggregation_macro.txt", path=OUTDIR)

#######################
### GENOME INDEXING ###
#######################

#premapping and aligning step for STAR
rule generate_genome_indexes:
    input:
        in_fasta = config["reference_unzip"], 
        in_gtf = config["reference_gtf_unzip"] 
    output:
        ref_dir = f"{datadir}/refgenomes/STAR/Gencode.v25/gencode.v25"
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
        transcriptomedir = f"{datadir}/refgenomes/RSEM/Gencode.v25/"
    params:
        transcriptomename = f"{datadir}/refgenomes/RSEM/Gencode.v25/genecode.v25.", 
    threads: config["MAX_THREADS"]
    shell:
        "rsem-prepare-reference -gtf {input.gtf} -p {threads} --star {input.fasta} {params.transcriptomename} "
        
#####################
### PREPROCESSING ###
#####################

#mapping, aligning, tagging and sorting step
rule STAR_alignment:
    input:
        R1 = "%s/all_fastqs/%s.R1.fastq.gz" % (datadir, "{samples}"),         
        R2 = "%s/all_fastqs/%s.R2.fastq.gz" % (datadir, "{samples}"),
        ref_dir =  f"{datadir}/refgenomes/STAR/Gencode.v25/gencode.v25",
        gtf = f'{datadir}/refgenomes/Gencode/v25/gencode.v25.primary_assembly.annotation.gtf',
    output:
        mock = "{path}/{experiment}/STAR/.{samples}_mockfile.txt"
    params:
        names = "{path}/{experiment}/STAR/{samples}.",
        sample = "{samples}"
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
 
#MARK DUPLICATES IS CONTRAINDICATED BECAUSE DUPLICATION RATE WAS <5%.
       
########################
### GENE EXPRESSSION ###
########################       
       
#Prep bam for RSEM post-processing
rule rsem_prep_bam:
    input:
        mock_in = "{path}/{experiment}/STAR/.{samples}_mockfile.txt"
    output:
        mock = "{path}/{experiment}/RSEM/{samples}.Aligned.toTranscriptome.rsem.mock.txt"
    params:
        bam_in = "{path}/{experiment}/STAR/{samples}.Aligned.toTranscriptome.out.bam",
        bam_out_name = "{path}/{experiment}/RSEM/{samples}.Aligned.toTranscriptome.rsem.out"
    shell:
        "convert-sam-for-rsem {params.bam_in} {params.bam_out_name} " 
        "&& touch {output.mock} "

#Run RSEM calculations
rule rsem_calc_expr:
    input:
        mock_in = "{path}/{experiment}/RSEM/{samples}.Aligned.toTranscriptome.rsem.mock.txt"
    output:
        mock = "{path}/{experiment}/RSEM/output/.{samples}_mockfile.rsem_calc.txt"
    params:
        bam_in = "{path}/{experiment}/RSEM/{samples}.Aligned.toTranscriptome.rsem.out.bam",
        ref_name = f"{datadir}/refgenomes/RSEM/Gencode.v25/genecode.v25.", 
        sample_name = "{path}/{experiment}/RSEM/output/{samples}",
        options = "--paired-end --no-bam-output --estimate-rspd ", 
    shell:
        "rsem-calculate-expression {params.options} --alignments "
        "{params.bam_in} {params.ref_name} {params.sample_name} "
        "&& touch {output.mock} "
 
######################
### METRIC CALLING ###
###################### 

#Collect multiple metrics for each bam file
rule collect_mult_metrics:
    input:
        mock = "{path}/{experiment}/STAR/.{samples}_mockfile.txt"
    output:
        metrics = "{path}/{experiment}/Metrics/{samples}.MultipleMetrics.alignment_summary_metrics"
    params:
        ref_fasta = config["reference_unzip"],
        bam_in = "{path}/{experiment}/STAR/{samples}.Aligned.sortedByCoord.out.bam",
        metrics = "{path}/{experiment}/Metrics/{samples}.MultipleMetrics"
    shell:
        "picard CollectMultipleMetrics I={params.bam_in} O={params.metrics} R={params.ref_fasta} "
        
####################################
### VARIANT CALLING | GENOTYPING ###
####################################

#reads that are "split", subread segments aligned to multiple locations, are aligned to one location with the rest of the read parts being hardclipped. 
rule splitncigarreads:
    input:
        mock_in = "{path}/{experiment}/STAR/.{samples}_mockfile.txt",
        ref = config["ref_unzip_wdict"] 
    output:
        bam_out = "{path}/{experiment}/VarCall_BAMs/{samples}.Aligned.sortedByCoord.split_r.out.bam" 
    params:
        bam_in =  "{path}/{experiment}/STAR/{samples}.Aligned.sortedByCoord.out.bam",
        options = "--create-output-bam-index true",
        intervals = chrom_intervals
    shell:
        "samtools index {params.bam_in} "
        "&& gatk SplitNCigarReads {params.options} {params.intervals} --input {params.bam_in} --reference {input.ref} --output {output.bam_out} "

#add read group information and calc tags as annotations. These annotations are used in the variant calling later. 
rule AddOrReplaceRG:
    input:
        bam = "{path}/{experiment}/VarCall_BAMs/{samples}.Aligned.sortedByCoord.split_r.out.bam"
    output:
        bam_out = "{path}/{experiment}/VarCall_BAMs/{samples}.Aligned.sortedByCoord.split_r.RG.out.bam"
    params:
        SM = "{samples}".split('.')[0],
    shell:
        "picard AddOrReplaceReadGroups -I {input.bam} -O {output.bam_out} "
        "--RGID {params.SM} --RGLB {params.SM} --RGPL illumina --RGSM {params.SM} --RGPU MiOrHiSeq"

#After being annotated with read groups etc the new file needs to be reindexed.
rule index_RG_bams:
    input:
        bam = "{path}/{experiment}/VarCall_BAMs/{samples}.Aligned.sortedByCoord.split_r.RG.out.bam"
    output:
        bam_out = "{path}/{experiment}/VarCall_BAMs/{samples}.Aligned.sortedByCoord.split_r.RG.out.bam.bai"
    shell:
        "samtools index {input.bam} "


rule allele_specific_expression:
    input:
        bam_index = "{path}/{experiment}/VarCall_BAMs/{samples}.Aligned.sortedByCoord.split_r.RG.out.bam.bai",
        bam = "{path}/{experiment}/VarCall_BAMs/{samples}.Aligned.sortedByCoord.split_r.RG.out.bam",
        fasta = config["ref_unzip_wdict"],
        ref_SNPs = config['ref_het_SNPs'],
    output:
        table = "{path}/{experiment}/ASE/{samples}.ASE.table"
    shell:
        "gatk ASEReadCounter -R {input.fasta} -I {input.bam} -V {input.ref_SNPs} -O {output.table} "

#Annotated indexed bams are now ready to be variant called with unified genotyper (UG). This script writes files with the UG calls. UG is fastest because it doesn't construct allles like heplotype caller.
rule etai_genotype_generate_call:
    input:
        [expand("{path}/{experiment}/VarCall_BAMs/{samples}.Aligned.sortedByCoord.split_r.RG.out.bam.bai", path=OUTDIR, experiment=experiments[i], samples=samples_set[i]) for i in range(len(experiments))]
    output:
        mockfile = "{path}/{experiment}/Variants/.mockfile.{experiment}.gengenotypecmds.txt"
    params:
        ls = "{path}/{experiment}/VarCall_BAMs/*RG.out.bam",
        hets = config["alleles"],
        ref = config["ref_unzip_wdict"],
        GATK = config["GATK_path"],
        Picard = config["Picard_path"],
        bamlist = "{path}/{experiment}/VarCall_BAMs/{experiment}.bamfiles.list",
        script = f"{skdir}/all_scripts/scripts/RPE-1_GRCh38_Genotype_nikos.sh", #*_etai.sh
        out_name = "{path}/{experiment}/Variants/{experiment}"
    shell:
        "ls {params.ls} > {params.bamlist} "
        "&& bash {params.script} {params.bamlist} {params.out_name} 1 {params.ref} {params.hets} {params.GATK} {params.Picard} "
        "&& touch {output.mockfile}"

#Now the the files with UG calls are ready we can parallelize the runs. These UG calls create one vcf file (for each chr) from all the bams in the experiment.
rule etai_genotype_run:
    input:
        "{path}/{experiment}/Variants/.mockfile.{experiment}.gengenotypecmds.txt",
    output:
        "{path}/{experiment}/Variants/.mockfile.{experiment}.rungenotypecmds.{chrs}.txt"
    params:
        "{path}/{experiment}/Variants/{experiment}_RPE_hets_GT.UG_jobs.{chrs}.sh"
    log:
        "{path}/{experiment}/Variants/{experiment}.rungenotypecmds.{chrs}.log"
    shell:
        "sh {params} 2> {log} " # "parallel --jobs {threads} < {params} "
        "&& touch {output}"

####################################
### DATA AGGREGATION | ANALYSIS  ###
####################################

#Before moving on to the analysis we check the required files have been created. 
rule check_data:
    input:
        [expand("{path}/{experiment}/RSEM/output/.{samples}_mockfile.rsem_calc.txt", path=OUTDIR, experiment=experiments[i], samples=samples_set[i]) for i in range(len(experiments))],
        [expand("{path}/{experiment}/Metrics/{samples}.MultipleMetrics.alignment_summary_metrics", path=OUTDIR, experiment=experiments[i], samples=samples_set[i]) for i in range(len(experiments))],
        expand("{path}/{experiment}/Variants/.mockfile.{experiment}.rungenotypecmds.{chrs}.txt", path=OUTDIR, experiment=experiments, chrs=chroms),
    output:
        mock_out = "{path}/{experiment}/Analysis/.mockfile.{experiment}.data_tracking.txt"
    shell:
        "touch {output.mock_out}"

#We now use this script to aggregate the data into R matrices for each experiment.
rule aggregate_data:
    input:
        mock_in = [expand("{path}/{experiment}/Analysis/.mockfile.{experiment}.data_tracking.txt", path=OUTDIR, experiment=experiments[i]) for i in range(len(experiments))]
    output:
        mock_out = "{path}/{experiment}/Analysis/.mockfile.{experiment}.data_aggregation.txt",
    params:
        script_dir = f"{skdir}/all_scripts",
        expr = "{experiment}",
        wk_dir = f"{OUTDIR}",
        script = f"{skdir}/all_scripts/scripts/Analysis_Etai_Snakemake.R",
        datadir = f"{datadir}"
    threads: 1
    shell:
        "Rscript {params.script} {params.script_dir} {params.wk_dir} {params.datadir} {params.expr} {threads} "
        "&& touch {output.mock_out} "

rule aggregate_data_macro_and_models:
    input:
        mock_in = [expand("{path}/{experiment}/Analysis/.mockfile.{experiment}.data_aggregation.txt", path=OUTDIR, experiment=experiments[i]) for i in range(len(experiments))],
    output:
        mock_out = f"{OUTDIR}/aggregated_results/.mockfile.data_aggregation_macro.txt", 
    params:
        script_dir = f"{skdir}/all_scripts",
        all_experiments = [experiments[i] for i in range(len(experiments))],
        wk_dir = f"{OUTDIR}",
        script = f"{skdir}/all_scripts/scripts/Data_Aggregation.R",
        datadir = f"{datadir}"
    shell:
        "Rscript {params.script} {params.script_dir} {params.wk_dir} {params.datadir} {params.all_experiments} "
        "&& touch {output.mock_out} "

