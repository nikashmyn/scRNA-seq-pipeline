#args <- c("/pellmanlab/nikos/Stam_Etai_Scripts", "/pellmanlab/stam_niko/data/processed_bam", "/pellmanlab/nikos/Stam_Etai_Data")
args <- commandArgs(trailingOnly = TRUE)
print(args)
scriptsdir <- args[1]
dirpath <- args[2]
datadir <- args[3]
mk_agg_dir <- sprintf("mkdir -p %s/ML_data", dirpath)
system(mk_agg_dir)

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
#rsemtpm <- readRDS(file = sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", dirpath))
#controlSampleIDs <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))
#controlSampleIDs2 <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs2.rds", dirpath))
#controlIDs <- readRDS(file = sprintf("%s/aggregated_results/controlIDs.rds", dirpath))
#adt <- adt.default <- data.table(readRDS(file = sprintf("%s/aggregated_results/adt.rds", dirpath)))
#adt.na <- data.table(readRDS(file = sprintf("%s/aggregated_results/adt.na.rds", dirpath)))
#coding <- readRDS(file = sprintf("%s/aggregated_results/ASE.coding.rds", dirpath))
#rsemtpm <- readRDS(file = sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", dirpath))
#centromeres <- data.table(readRDS(file = sprintf("%s/centromeres.rds", datadir)))
#arms <- readRDS(file = sprintf("%s/CN_data/CN_predictions.byarm.rds", dirpath))

#order my data objects 
#adt <- adt.default <- adt[order(seqnames, start, end)]
#adt <- cbind( adt[,c(1:4)], setcolorder(adt[,-c(1:4)], order(colnames(adt[,-c(1:4)]))) )
#adt.na <- adt.na[order(seqnames, start, end)]
#adt.na <- cbind( adt.na[,c(1:4)], setcolorder(adt.na[,-c(1:4)], order(colnames(adt.na[,-c(1:4)]))) )
#ganno <- ganno[order(seqnames, start, end)]
#stopifnot(sum(ganno$id != adt$id) == 0)

rsemtpm_e <- readRDS("/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2/data/rsemtpm.rds")
adt_e <- adt_e.default <- readRDS("/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2/data/adt.rds")
adt.na_e <- readRDS("/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2/data/adt.na.rds")
ganno_e <- readRDS("/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2/data/ganno.rds")
arms_e <- readRDS("/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2/CN_predictions.byarm.rds")
arms_e <- arms_e[c("raw.A", "raw.B", "MA50", "MA75", "MA100", "SSM.TE")]
coding_e <- readRDS("/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2/data/ASE.coding.rds")

#intersect arms
#columns_to_keep <- intersect(unlist(colnames(arms[["MA50"]])), unlist(colnames(arms_e[["MA50"]])))
#arms$MA50 <- arms$MA50[, ..columns_to_keep]
#arms_e$MA50 <- arms_e$MA50[, ..columns_to_keep]
#columns_to_keep <- intersect(unlist(colnames(arms[["MA75"]])), unlist(colnames(arms_e[["MA75"]])))
#arms$MA75 <- arms$MA75[, ..columns_to_keep]
#arms_e$MA75 <- arms_e$MA75[, ..columns_to_keep]
#columns_to_keep <- intersect(unlist(colnames(arms[["MA100"]])), unlist(colnames(arms_e[["MA100"]])))
#arms$MA100 <- arms$MA100[, ..columns_to_keep]
#arms_e$MA100 <- arms_e$MA100[, ..columns_to_keep]
#columns_to_keep <- intersect(unlist(colnames(arms[["raw.A"]])), unlist(colnames(arms_e[["raw.A"]])))
#arms$raw.A <- arms$raw.A[, ..columns_to_keep]
#arms_e$raw.A <- arms_e$raw.A[, ..columns_to_keep]
#columns_to_keep <- intersect(unlist(colnames(arms[["raw.B"]])), unlist(colnames(arms_e[["raw.B"]])))
#arms$raw.B <- arms$raw.B[, ..columns_to_keep]
#arms_e$raw.B <- arms_e$raw.B[, ..columns_to_keep]
#columns_to_keep <- intersect(unlist(colnames(arms[["SSM.TE"]])), unlist(colnames(arms_e[["SSM.TE"]])))
#arms$SSM.TE <- arms$SSM.TE[, ..columns_to_keep]
#arms_e$SSM.TE <- arms_e$SSM.TE[, ..columns_to_keep]

#intersect adt
#columns_to_keep <- append(colnames(adt)[1:4], columns_to_keep[-c(1)])
#rows_adt <- intersect(adt$id, adt_e$id)
#adt <- adt[which(adt$id %in% rows_adt), ..columns_to_keep]
#adt_e <- adt_e[which(adt_e$id %in% rows_adt), ..columns_to_keep]

#intersect ganno
#rows_to_keep <- intersect(ganno$id, ganno_e$id)
#ganno <- ganno[which(ganno$id %in% rows_to_keep),]
#ganno_e <- ganno_e[which(ganno_e$id %in% rows_to_keep),]
#
#intersect coding$A
#rows_to_keep <- intersect(coding_e$A$id, coding$A$id)
#coding_e$A <- coding_e$A[which(coding_e$A$id %in% rows_to_keep),]
#coding_e$A <- coding_e$A[which(coding_e$A$id %in% rows_to_keep),]

#intersect coding$B
#rows_to_keep <- intersect(coding_e$B$id, coding$B$id)
#coding_e$B <- coding_e$B[which(coding_e$B$id %in% rows_to_keep),]
#coding_e$B <- coding_e$B[which(coding_e$B$id %in% rows_to_keep),]

#Use Etai data sets
rsemtpm <- rsemtpm_e; adt <- adt.default <- adt_e; adt.na <- adt.na_e; ganno <- ganno_e; arms <- arms_e; coding <- coding_e;

#get samples to run based on samples in df
#anno <- data.table(readRDS(sprintf("%s/work_in_progress/Annotation_list_long_vnikos_210129.rds", datadir) ))
#samples_to_use <- anno[!LookRNAseq.Exp %in% exclude]$WTA.plate
#columns <- colnames(adt)[-c(1:4)]
#samples_to_use <- c(intersect(columns, samples_to_use))

anno <- data.table(read_xlsx(path = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2/data/samples_annotation.xlsx"))
samples_to_use <- anno[!LookRNAseq.Exp %in% exclude]$WTA.plate

#get high quality samples from raw tpm values
high_qc_ids <- names(which(colSums(rsemtpm>5)>4000))

#columns annotations used in file naming convention for visuals
#col_anno <- data.table(readRDS(sprintf("%s/work_in_progress/Annotation_list_long_vnikos_210129.rds", datadir) ))
col_anno <- data.table(read_xlsx(path = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2/data/samples_annotation.xlsx"))
col_anno <- col_anno[ WTA.plate %in% high_qc_ids]
col_anno[Pairs == "NA", Pairs := NA]

#Let's find small chrs 
sort(table(adt$seqnames))
arms_to_exclude <- c("10q", "Xp", "Xq", "21p", "21q", "18p", "18q", "13q", "22p", "22q", "20p", "20q")

#dirpath adjustment to save seperately from normal models
#dirpath <- sprintf("%s/test", dirpath)

#Run ML scripts: M %s/ML/... was run before 
source(sprintf("%s/original_ML/data_preparation.base_model.R", scriptsdir))
source(sprintf("%s/original_ML/training.base_model.R", scriptsdir))
#source(sprintf("%s/ML/predict.base_model.on_our_data.R", scriptsdir))
source(sprintf("%s/original_ML/predict.Stamatis_lab_meeting_Jan08_2019.R", scriptsdir))





#### Data Comparisons ###
dirpath <- "/pellmanlab/stam_niko/data/processed_bam"
olr_model_fname <- sprintf("%s/NN/olr_model_v1.AVG%d.rds", dirpath, NUM_OF_FEATURES)
olr_model_nikos <- readRDS(olr_model_fname)


dirpath <- "/pellmanlab/stam_niko/data/processed_bam/test"
olr_model_fname <- sprintf("%s/NN/olr_model_v1.AVG%d.rds", dirpath, NUM_OF_FEATURES)
olr_model_etai2 <- readRDS(olr_model_fname)


#
##Samples to compare
#samples <- c("161130_B6", "170202_B7", "190701_9F", "190628_4B")
#
##arms
#arms.SSM <- unlist(arms[["raw.B"]][,"161130_B6"]) #SSM.TE MA100 raw.B 
#arms1.SSM <- unlist(arms_e[["raw.B"]][,"161130_B6"])
#plot(x = arms.SSM, y = arms1.SSM)
#
##adt
#rows_adt <- intersect(adt$id, adt_e$id)
#adt_short <- adt[which(adt$id %in% rows_adt), c("id", "160627_B1")]
#adt_e_short <- adt_e[which(adt_e$id %in% rows_adt), c("id", "160627_B1")]
#adt_short <- adt_short[order(id),]; adt_e_short <- adt_e_short[order(id),]
#plot( unlist(adt_short), unlist(adt_e_short) )
#
##output models
#olr_model_etai <- readRDS("/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2/NN/olr_model_v1.AVG50.rds")
#olr_model_nikos <- readRDS(olr_model_fname)

