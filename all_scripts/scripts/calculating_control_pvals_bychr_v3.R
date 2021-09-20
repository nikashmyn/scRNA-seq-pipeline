#Calculating p-values script by Nikos 8/18/2021

######################
### Source Scripts ###
######################

require(data.table)
require(readxl)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
#install.packages("factoextra")
library(factoextra) # clustering algorithms & visualization
source(sprintf('%s/scripts/fromRawTPMtoExprsRatio.R', scriptsdir))

################
#### Imports ###
################

#Annotion list with information on every cell and IDs of control samples
anno <- data.table(read.csv( sprintf("%s/work_in_progress/Annotation_list_long_vnikos_1_9_21.csv", datadir)))
controlSampleIDs <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))

#QC Data Frame and high QC IDs based on 5 read cnt threshold for genes
all_QC <- readRDS(file = sprintf("%s/aggregated_results/all_QC.rds", dirpath))
high_qc_ids <- as.character(all_QC[which(all_QC$th5 >= quantile(all_QC$th5, c(.10))),]$id) #Excludes bottom 10%

#Data Frames with TPM and Variant information
#adt <- readRDS(sprintf("%s/aggregated_results/adt.rds", dirpath)) #Normalized TPM matrix
coding <- readRDS(sprintf("%s/aggregated_results/ASE.coding.rds", dirpath)) #Raw Variant matrix in the coding and UTR regions

#import raw gene cnts and gene ranges
geneRanges <- readRDS(sprintf("%s/geneRanges.rds", datadir))
rsemtpm <- readRDS(sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", dirpath))

#new adt object without gene normalization
adt_tmp <- fromRawTPMtoExprsRatio(rsemtpm, geneRanges, controlSampleIDs, max_pos_col = 6, plusOne = 1, 
                              zerosAsNA = F, normBySd = F, doNormalizeByControls = F,
                              maxExp = 1200, quantileOfExp = 0.99, minDetectionLevel = 5, minNumOfSamplesToDetect = 5 )$adt

############################################
### Re-Engineer raw variant count Matrix ###
############################################

#Bin adt 
adt_tmp <- adt_tmp[order(seqnames, start)]

#group TPM by chr 
adt.bychr <- adt_tmp[,-c(1,3,4,5,6)] %>%
  group_by(seqnames) %>%
  summarise_all(mean, na.rm = TRUE)
adt.bychr <- data.table(adt.bychr)
adt.bychr <- adt.bychr[order(seqnames)]

#normalize adt columns
adt_rowmeans <- rowMeans(adt.bychr[,..controlSampleIDs], na.rm = T)
adt_rowmeans_mat_reciprocal <- data.table(replicate(ncol(adt.bychr[,-c(1)]), adt_rowmeans))
adt.bychr <- cbind(adt.bychr[,c(1)] , (adt.bychr[,-c(1)] - adt_rowmeans_mat_reciprocal))

adt.bychr <- cbind(adt.bychr[,c(1)], 2^adt.bychr[,-c(1)])
adt3_colmeans <- colMeans(adt.bychr[,-c(1)], na.rm = T)
adt3_colmeans_mat_reciprocal <- data.table(1/(t(replicate(nrow(adt.bychr), adt3_colmeans))))
adt3_colmeans_mat_reciprocal[adt3_colmeans_mat_reciprocal == Inf] <- 0 
adt.bychr <- cbind(adt.bychr[,c(1)] , (adt.bychr[,-c(1)] * adt3_colmeans_mat_reciprocal))

#Read in data
golden_samples <- readRDS(sprintf("%s/ML_data/golden_set_ids.rds", dirpath))

#Read in grouped sample information
hand_samples <- readRDS(sprintf("%s/aggregated_results/grouped_control_aneuploidies.rds", dirpath))

#change loss and gain examples to hand picked grouped samples
golden_samples$loss <- hand_samples[which(hand_samples$CN == 1),c(1,3)]
golden_samples$gain <- hand_samples[which(hand_samples$CN == 3),c(1,3)]
setnames(golden_samples$loss, old = "ID", new = "sample_id")
setnames(golden_samples$gain, old = "ID", new = "sample_id")

#add chr column to golden samples by adjusting arm strings
golden_samples$normal$chr <- gsub('.{1}$', '', golden_samples$normal$arm)
golden_samples$normal$chr <- paste0("chr", golden_samples$normal$chr)
golden_samples$loss$chr <- paste0("chr", golden_samples$loss$chr)
golden_samples$gain$chr <- paste0("chr", golden_samples$gain$chr)

#Change name convention
adt_bychr <- adt.bychr
setnames(adt_bychr, "seqnames", "chr")
adt_bychr$chr <- as.character(adt_bychr$chr)
adt_bychr <- as.data.frame(adt_bychr, col.names = colnames(adt_bychr))

### get adt values for golden set regions ###

#Get loss samples' tpm by arm and by chr
adt_golden_loss_bychr <- c()
for (i in 1:nrow(golden_samples$loss)) {
  adt_golden_loss_bychr <- append(adt_golden_loss_bychr, as.numeric(as.character(adt_bychr[which(adt_bychr$chr %in% golden_samples$loss$chr[i]), golden_samples$loss$sample_id[i]])) )
}
golden_samples$loss_bychr <- data.table(cbind(golden_samples$loss[,c("sample_id", "chr")], adt_golden_loss_bychr))
colnames(golden_samples$loss_bychr)[3] <- "tpm"
golden_samples$loss_bychr$tpm <- sapply(golden_samples$loss_bychr$tpm, as.numeric)

#Get loss samples summary stats
golden_samples$loss_stats_bychr <- c()
golden_samples$loss_stats_bychr$mean <- mean(golden_samples$loss_bychr$tpm)
golden_samples$loss_stats_bychr$sd <- sd(golden_samples$loss_bychr$tpm)
golden_samples$loss_stats_bychr$n <- length(golden_samples$loss_bychr$tpm)

#get normal samples' tpm by arm and by chr
adt_golden_normal_bychr <- c()
for (i in 1:nrow(golden_samples$normal)) {
  adt_golden_normal_bychr <- append(adt_golden_normal_bychr, as.numeric(as.character(adt_bychr[which(adt_bychr$chr %in% golden_samples$normal$chr[i]), golden_samples$normal$sample_id[i]])))
}
golden_samples$normal_bychr <- data.table(cbind(golden_samples$normal[,c("sample_id", "chr")], adt_golden_normal_bychr))
colnames(golden_samples$normal_bychr)[3] <- "tpm"
golden_samples$normal_bychr$tpm <- sapply(golden_samples$normal_bychr$tpm, as.numeric)

#Get normal samples summary stats
golden_samples$normal_stats_bychr <- c()
golden_samples$normal_stats_bychr$mean <- mean(golden_samples$normal_bychr$tpm)
golden_samples$normal_stats_bychr$sd <- sd(golden_samples$normal_bychr$tpm)
golden_samples$normal_stats_bychr$n <- length(golden_samples$normal_bychr$tpm)

#Get gain samples' tpm by arm and by chr
adt_golden_gain_bychr <- c()
for (i in 1:nrow(golden_samples$gain)) {
  adt_golden_gain_bychr <- append(adt_golden_gain_bychr, as.numeric(as.character(adt_bychr[which(adt_bychr$chr %in% golden_samples$gain$chr[i]), golden_samples$gain$sample_id[i]])))
}
golden_samples$gain_bychr <- data.table(cbind(golden_samples$gain[,c("sample_id", "chr")], adt_golden_gain_bychr))
colnames(golden_samples$gain_bychr)[3] <- "tpm"
golden_samples$gain_bychr$tpm <- sapply(golden_samples$gain_bychr$tpm, as.numeric)

#Get gain samples summary stats
golden_samples$gain_stats_bychr <- c()
golden_samples$gain_stats_bychr$mean <- mean(golden_samples$gain_bychr$tpm)
golden_samples$gain_stats_bychr$sd <- sd(golden_samples$gain_bychr$tpm)
golden_samples$gain_stats_bychr$n <- length(golden_samples$gain_bychr$tpm)

#Get control samples' tpm by arm and by chr
adt_control_normal_bychr <- data.frame()
for (i in 1:length(controlSampleIDs)) {
  adt_control_normal_bychr <- rbind( adt_control_normal_bychr, cbind( controlSampleIDs[i], cbind( data.frame(adt_bychr$chr), data.frame(unname(adt_bychr[,controlSampleIDs[i]])) ) ) )
}
#colnames(adt_control_normal_bychr) <- adt_control_normal_bychr[c(1),]; adt_control_normal_bychr <- adt_control_normal_bychr[-c(1),]
colnames(adt_control_normal_bychr) <- c("sample", "chr", "tpm")
adt_control_normal_bychr$tpm <- sapply(adt_control_normal_bychr$tpm, as.numeric)
golden_samples$control_bychr <- data.table(adt_control_normal_bychr)

#Get control samples summary stats
golden_samples$control_stats_bychr <- c()
golden_samples$control_stats_bychr$mean <- mean(golden_samples$control_bychr$tpm)
golden_samples$control_stats_bychr$sd <- sd(golden_samples$control_bychr$tpm)
golden_samples$control_stats_bychr$n <- length(golden_samples$control_bychr$tpm)

#Get arm specific control sample summary stats
golden_samples$control_stats_bychr$chr_means <- golden_samples$control_bychr[,-c(1)] %>% group_by(chr) %>% summarise_all(mean, na.rm = TRUE)
golden_samples$control_stats_bychr$chr_sds <- golden_samples$control_bychr[,-c(1)] %>% group_by(chr) %>% summarise_all(sd, na.rm = TRUE)
golden_samples$control_stats_bychr$chr_ns <- golden_samples$control_bychr[,-c(1)] %>% group_by(chr) %>% summarise_all(length)

#save adjusted golden samples object
saveRDS(golden_samples, file = sprintf("%s/aggregated_results/golden_set_tpms.rds", dirpath))

#############################################################################
### Get an adt-like matrix of z-scores from the labeled sample statistics ###
#############################################################################

#z-score for loss sample statistics
adt_center_loss_bychr <- cbind( adt_bychr[,c(1)], sweep(adt_bychr[,-c(1)], 1, golden_samples$loss_stats_bychr$mean, '-') )
adt_zs_loss_bychr <- cbind( adt_bychr[,c(1)], sweep(adt_center_loss_bychr[,-c(1)], 1, golden_samples$loss_stats_bychr$sd, '/') )

#for normal sample statistics
adt_center_normal_bychr <- cbind( adt_bychr[,c(1)], sweep(adt_bychr[,-c(1)], 1, golden_samples$normal_stats_bychr$mean, '-') )
adt_zs_normal_bychr <- cbind( adt_bychr[,c(1)], sweep(adt_center_normal_bychr[,-c(1)], 1, golden_samples$normal_stats_bychr$sd, '/') )

#for gain sample statistics
adt_center_gain_bychr <- cbind( adt_bychr[,c(1)], sweep(adt_bychr[,-c(1)], 1, golden_samples$gain_stats_bychr$mean, '-') )
adt_zs_gain_bychr <- cbind( adt_bychr[,c(1)], sweep(adt_center_gain_bychr[,-c(1)], 1, golden_samples$gain_stats_bychr$sd, '/') )

#got control samples statistics
adt_center_control_bychr <- cbind( adt_bychr[,c(1)], sweep(adt_bychr[,-c(1)], 1, golden_samples$control_stats_bychr$mean, '-') )
adt_zs_control_bychr <- cbind( adt_bychr[,c(1)], sweep(adt_center_control_bychr[,-c(1)], 1, golden_samples$control_stats_bychr$sd, '/') )

#z-score for control samples statistics
adt_center_control_bychr_specific <- adt_bychr[,-c(1)] - c(golden_samples$control_stats_bychr$chr_means[,2])
adt_zs_control_bychr_specific <- cbind( adt_bychr[,c(1)], adt_center_control_bychr_specific / c(golden_samples$control_stats_bychr$chr_sds[,2]) )

### Get pval matrix by chr ###
pval_matrix_loss_bychr <- c(); pval_matrix_normal_bychr <- c(); pval_matrix_gain_bychr <- c(); pval_matrix_control_bychr <- c(); pval_matrix_control_bychr_specific <- c();
for (i in 1:nrow(adt_bychr)) {
  row_loss <- c(); row_normal <- c(); row_gain <- c(); row_control <- c(); row_control_specific <- c();
  for (j in 2:ncol(adt_bychr)) {
    row_loss <- append( row_loss, 2*pnorm( -abs(adt_zs_loss_bychr[i,j]) ) ) #Should these be one or two tail?
    row_normal <- append( row_normal, 2*pnorm( -abs(adt_zs_normal_bychr[i,j]) ) ) #pval is multiplied by two cause each z-score could be above or below the class specific population stats
    row_gain <- append( row_gain, 2*pnorm( -abs(adt_zs_gain_bychr[i,j]) ) ) #Should these be one or two tail?
    row_control <- append( row_control, 2*pnorm( -abs(adt_zs_control_bychr[i,j]) ) )
    row_control_specific <- append( row_control_specific, 2*pnorm( -abs(adt_zs_control_bychr_specific[i,j]) ))
  }
  pval_matrix_loss_bychr <- cbind(pval_matrix_loss_bychr, data.table(row_loss))
  pval_matrix_normal_bychr <- cbind(pval_matrix_normal_bychr, data.table(row_normal))
  pval_matrix_gain_bychr <- cbind(pval_matrix_gain_bychr, data.table(row_gain))
  pval_matrix_control_bychr <- cbind(pval_matrix_control_bychr, data.table(row_control))
  pval_matrix_control_bychr_specific <- cbind(pval_matrix_control_bychr_specific, data.table(row_control_specific))
}
pval_matrix_loss_bychr <- data.table(t(pval_matrix_loss_bychr))
pval_matrix_normal_bychr <- data.table(t(pval_matrix_normal_bychr))
pval_matrix_gain_bychr <- data.table(t(pval_matrix_gain_bychr))
pval_matrix_control_bychr <- data.table(t(pval_matrix_control_bychr))
pval_matrix_control_bychr_specific <- data.table(t(pval_matrix_control_bychr_specific))

pval_matrix_loss_bychr <- cbind(adt_bychr$chr, pval_matrix_loss_bychr)
pval_matrix_normal_bychr <- cbind(adt_bychr$chr, pval_matrix_normal_bychr)
pval_matrix_gain_bychr <- cbind(adt_bychr$chr, pval_matrix_gain_bychr)
pval_matrix_control_bychr <- cbind(adt_bychr$chr, pval_matrix_control_bychr)
pval_matrix_control_bychr_specific <- cbind(adt_bychr$chr, pval_matrix_control_bychr_specific)

colnames(pval_matrix_loss_bychr) <- colnames(adt_bychr)
colnames(pval_matrix_normal_bychr) <- colnames(adt_bychr)
colnames(pval_matrix_gain_bychr) <- colnames(adt_bychr)
colnames(pval_matrix_control_bychr) <- colnames(adt_bychr)
colnames(pval_matrix_control_bychr_specific) <- colnames(adt_bychr)

saveRDS(pval_matrix_loss_bychr, file = sprintf("%s/aggregated_results/pval_matrix_loss_bychr.rds", dirpath))
saveRDS(pval_matrix_normal_bychr, file = sprintf("%s/aggregated_results/pval_matrix_normal_bychr.rds", dirpath))
saveRDS(pval_matrix_gain_bychr, file = sprintf("%s/aggregated_results/pval_matrix_gain_bychr.rds", dirpath))
saveRDS(pval_matrix_control_bychr, file = sprintf("%s/aggregated_results/pval_matrix_control_bychr.rds", dirpath))
saveRDS(pval_matrix_control_bychr_specific, file = sprintf("%s/aggregated_results/pval_matrix_control_bychr_specific.rds", dirpath))

##############################
### Get significant chroms ###
##############################

#reorder chroms by number
pval_matrix_control_bychr$chr <- sub(pattern = "chr", replacement = "", x = pval_matrix_control_bychr$chr)
pval_matrix_control_bychr$chr <- as.integer(sub(pattern = "X", replacement = "23", x = pval_matrix_control_bychr$chr))
pval_matrix_control_bychr <- pval_matrix_control_bychr[order(chr)]

pval_matrix_control_bychr_noname <- setkey(pval_matrix_control_bychr)[,-c(1)]

#setkey(pval_matrix_control_bychr_noname)
inds <- which(pval_matrix_control_bychr_noname < .05, arr.ind = TRUE, useNames = TRUE)
#test <- subset(pval_matrix_control_bychr_noname, pval_matrix_control_bychr_noname < .05 )
#test <- apply(pval_matrix_control_bychr_noname, 2, function(r) which(r < .05))
rows <- rownames(pval_matrix_control_bychr_noname)[inds[,1]]
cols <- colnames(pval_matrix_control_bychr_noname)[inds[,2]]
non_diploid_chrs <- cbind(sample = cols, chr = rows) #already sorted
non_diploid_chrs_vals <- c()
for (i in 1:nrow(inds)) {
  ind1 <- inds[i,1]; ind2 <- inds[i,2];
  tmp <- unlist(pval_matrix_control_bychr_noname[ind1,])
  tmp2 <- tmp[ind2]
  non_diploid_chrs_vals <- append(non_diploid_chrs_vals, tmp2)
}
non_diploid_chrs_w_vals <- cbind(non_diploid_chrs, non_diploid_chrs_vals)
rownames(non_diploid_chrs_w_vals) <- NULL
colnames(non_diploid_chrs_w_vals)[3] <- "Pval"
non_diploid_chrs_w_vals[,c("chr")] <- paste0("chr", non_diploid_chrs_w_vals[,c("chr")])

saveRDS(non_diploid_chrs_w_vals, file = sprintf("%s/aggregated_results/non_diploid_chromosomes.rds", dirpath))
