
########################
### Passed arguments ###
########################

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/rerun_6_9_2021/data", "/pellmanlab/nikos/Stam_Etai_Data", "SIS1025a",  "SIS1025b", "SIS1025c", "SIS1025d", "SIS1025e", "SIS1025f_Lane1", "SIS1025f_Lane2", "SIS1025g_Lane1", "SIS1025g_Lane2", "SIS1025misc", "SIS1025targ")
print(args)
scriptsdir <- args[1]
dirpath <- args[2]
datadir <- args[3]
mk_agg_dir <- sprintf("mkdir -p %s/ML_data", dirpath)
system(mk_agg_dir)

#set seed before loading in TF and Keras
configs <- readRDS(sprintf("%s/param_config_list.rds", datadir))
configs$seed <- 205
set.seed(configs$seed)

########################
### loading packages ###
########################

require(keras)
require(tensorflow)
require(data.table)
require(readxl)
require(dplyr)
require(zoo)
require(MASS)

###############################
### Input variable and data ###
###############################

#let's list the relevant samples for downstream analysis:
exclude <- c("113015_T12", "061516_T12", "070516_LookRNAseq", "072516_LookRNAseq", 
             "071816_LookRNAseq", "071816_LookRNAseq", "161005_LookRNAseq", 
             "061316_LookRNAseq", "161005_LookRNAseq_control", "Jinyu", "170726_MCF10A_control",
             "040416_ATCC", "161017_LookRNAseq", "062716_LookRNAseq", "080816LookRNAseq", "080816_LookRNAseq",
             "071416_LookRNAseq", "161003_LookRNAseq")

#lets load all of our processed data for the ML model
rsemtpm <- readRDS(file = sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", dirpath))
controlSampleIDs <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))
controlSampleIDs2 <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs2.rds", dirpath))
controlIDs <- readRDS(file = sprintf("%s/aggregated_results/controlIDs.rds", dirpath))
adt <- adt.default <- data.table(readRDS(file = sprintf("%s/aggregated_results/adt.rds", dirpath)))
adt.na <- data.table(readRDS(file = sprintf("%s/aggregated_results/adt.na.rds", dirpath)))
ganno <- data.table(readRDS(file = sprintf("%s/aggregated_results/ganno.rds", dirpath)))
centromeres <- data.table(readRDS(file = sprintf("%s/centromeres.rds", datadir)))
arms <- readRDS(file = sprintf("%s/CN_data/CN_predictions.byarm.rds", dirpath))

#order my data objects 
adt <- adt.default <- adt[order(seqnames, start, end)]
adt <- cbind( adt[,c(1:4)], setcolorder(adt[,-c(1:4)], order(colnames(adt[,-c(1:4)]))) )
adt.na <- adt.na[order(seqnames, start, end)]
adt.na <- cbind( adt.na[,c(1:4)], setcolorder(adt.na[,-c(1:4)], order(colnames(adt.na[,-c(1:4)]))) )
ganno <- ganno[order(seqnames, start, end)]
stopifnot(sum(ganno$id != adt$id) == 0)

#get samples to run based on samples in df
anno <- data.table(read.csv( sprintf("%s/work_in_progress/annotation_list.csv", datadir)))
samples_to_use <- anno[!LookRNAseq.Exp %in% exclude]$WTA.plate
columns <- colnames(adt)[-c(1:4)]
samples_to_use <- c(intersect(columns, samples_to_use))
samples_to_use <- c(intersect(columns, anno$WTA.plate))

#get high quality samples from raw tpm values
all_QC <- readRDS(file = sprintf("%s/aggregated_results/all_QC.rds", dirpath))
high_qc_ids <- as.character(all_QC[which(all_QC$th5 >= quantile(all_QC$th5, c(.10))),]$id)
high_qc_ids <- intersect(high_qc_ids, samples_to_use)

#Variant matrix in the coding and UTR regions
coding <- readRDS(sprintf("%s/aggregated_results/ASE.coding.rds", dirpath))
#Low var counts exclusion. Note: both distributions of var counts are fairly similar.
#high_varqc_idsA <- names(which(colSums(coding$cnts.A[,-c(1:2)]) > quantile(colSums(coding$cnts.A[,-c(1:2)]), c(.10))))
#high_varqc_idsB <- names(which(colSums(coding$cnts.B[,-c(1:2)]) > quantile(colSums(coding$cnts.B[,-c(1:2)]), c(.10))))
#high_varqc_ids <- intersect(high_varqc_idsA, high_varqc_idsB)
#high_qc_ids <- intersect(high_qc_ids, high_varqc_ids)

#columns annotations used in file naming convention for visuals
col_anno <- data.table(read.csv( sprintf("%s/work_in_progress/annotation_list.csv", datadir)))
dim(col_anno)
col_anno <- col_anno[ WTA.plate %in% high_qc_ids]
dim(col_anno)
col_anno[Pairs == "NA", Pairs := NA]

# Temporary: train only on high quality samples 
samples_to_use <- high_qc_ids

###########################
### Running OLR Scripts ###
###########################

#Let's find small chrs 
sort(table(adt$seqnames))
#arms_to_exclude <- c("10q", "Xp", "Xq", "21p", "21q", "18p", "18q", "13q", "22p", "22q", "20p", "20q")
arms_to_exclude <- c("Xp", "Xq")
#arms_to_exclude <- c()

#Run ML scripts
source(sprintf("%s/ML/data_preparation.base_model.R", scriptsdir)) #seeded for reproducibility
source(sprintf("%s/ML/training.base_model2.R", scriptsdir)) #seeded for reproducibility
#source(sprintf("%s/ML/predict.base_model.on_our_data.R", scriptsdir))
source(sprintf("%s/ML/predict.Stamatis_lab_meeting_Jan08_2019.R", scriptsdir))

