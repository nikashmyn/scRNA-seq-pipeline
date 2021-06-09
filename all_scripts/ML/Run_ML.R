
#################################################################################
# Passed arguments:
#################################################################################

#args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/data/processed_bam", "/pellmanlab/nikos/Stam_Etai_Data")
args <- commandArgs(trailingOnly = TRUE)
print(args)
scriptsdir <- args[1]
dirpath <- args[2]
datadir <- args[3]
mk_agg_dir <- sprintf("mkdir -p %s/ML_data", dirpath)
system(mk_agg_dir)

#################################################################################
# loading packages:
#################################################################################

require(keras)
require(tensorflow)
require(data.table)
require(readxl)
require(dplyr)
require(zoo)
require(MASS)

###########################
# Input variable and data #
###########################

#let's list the relevant samples for downstream analysis:
exclude <- c("113015_T12", "061516_T12", "070516_LookRNAseq", "072516_LookRNAseq", 
             "071816_LookRNAseq", "071816_LookRNAseq", "161005_LookRNAseq", 
             "061316_LookRNAseq", "161005_LookRNAseq_control", "Jinyu", "170726_MCF10A_control",
             "040416_ATCC", "161017_LookRNAseq", "062716_LookRNAseq", "080816LookRNAseq", "080816_LookRNAseq",
             "071416_LookRNAseq", "161003_LookRNAseq")

#let's list the chrs we do not want to use for training and testing:
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
anno <- data.table(readRDS(sprintf("%s/work_in_progress/Annotation_list_long_vnikos_210129.rds", datadir) ))
samples_to_use <- anno[!LookRNAseq.Exp %in% exclude]$WTA.plate
columns <- colnames(adt)[-c(1:4)]
samples_to_use <- c(intersect(columns, samples_to_use))

#get high quality samples from raw tpm values
high_qc_ids <- names(which(colSums(rsemtpm>5)>4000))

#columns annotations used in file naming convention for visuals
col_anno <- data.table(readRDS( sprintf("%s/work_in_progress/Annotation_list_long_vnikos_210129.rds", datadir) ))
col_anno <- col_anno[ WTA.plate %in% high_qc_ids]
col_anno[Pairs == "NA", Pairs := NA]


#################################################################################
# Running OLR Scripts:
#################################################################################

#Let's find small chrs 
sort(table(adt$seqnames))
arms_to_exclude <- c("10q", "Xp", "Xq", "21p", "21q", "18p", "18q", "13q", "22p", "22q", "20p", "20q")
#arms_to_exclude <- c("Xp", "Xq")

#Run ML scripts
source(sprintf("%s/ML/data_preparation.base_model.R", scriptsdir))
source(sprintf("%s/ML/training.base_model2.R", scriptsdir))
#source(sprintf("%s/ML/predict.base_model.on_our_data.R", scriptsdir))
source(sprintf("%s/ML/predict.Stamatis_lab_meeting_Jan08_2019.R", scriptsdir))


