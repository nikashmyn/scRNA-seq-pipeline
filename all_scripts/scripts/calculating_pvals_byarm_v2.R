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
anno <- data.table(read.csv( sprintf("%s/work_in_progress/annotation_list.csv", datadir)))
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

#list of preset parameters for the experiment
configs <- readRDS(sprintf("%s/param_config_list.rds", datadir))

#new adt object without gene normalization
adt_tmp <- fromRawTPMtoExprsRatio(rsemtpm, geneRanges, controlSampleIDs, max_pos_col = 6, plusOne = 1, 
                              zerosAsNA = F, normBySd = F, doNormalizeByControls = F,
                              maxExp = 1200, quantileOfExp = 0.99, minDetectionLevel = configs$minDetectionLevel, minNumOfSamplesToDetect = configs$minNumOfSamplesToDetect )$adt

############################################
### Re-Engineer raw variant count Matrix ###
############################################

#order my data objects 
adt_ord <- adt_tmp[order(seqnames, start, end)]
adt_ord <- data.table(cbind( adt_tmp[,c(1:4)], setcolorder(adt_tmp[,-c(1:4)], order(colnames(adt_tmp[,-c(1:4)]))) ))
ganno_ord <- data.table(ganno[order(seqnames, start, end), c("id", "seqnames", "arm", "start", "end")])
stopifnot(sum(ganno_ord$id != adt$id) == 0)
adt_noanno <- adt_ord[,-c(1:4)]

#setkey
armRanges <- readRDS( sprintf("%s/armRanges.rds", datadir) )[,c("seqnames", "start", "end", "arm")]
setkeyv(adt_ord, c("seqnames", "start", "end"))
setkeyv(armRanges, c("seqnames", "start", "end"))
adt_ord_arm <- foverlaps(adt_ord, armRanges)
setnames(adt_ord_arm, old = c("i.start", "i.end", "start", "end"), new = c("start", "end", "arm_start", "arm_end"))
adt_ord_arm <- cbind(adt_ord_arm[,c("id", "seqnames", "arm", "start", "end", "arm_start", "arm_end", "strand", "width")], adt_ord_arm[,-c("id", "seqnames", "arm", "start", "end", "arm_start", "arm_end", "strand", "width")])

#group TPM by chr 
adt.byarm <- adt_ord_arm[,-c(1,2,4,5,6,7,8,9)] %>%
  group_by(arm) %>%
  summarise_all(mean, na.rm = TRUE)
adt.byarm <- data.table(adt.byarm)

#normalize adt columns
adt_rowmeans <- rowMeans(adt.byarm[,..controlSampleIDs], na.rm = T)
adt_rowmeans_mat_reciprocal <- data.table(replicate(ncol(adt.byarm[,-c(1)]), adt_rowmeans))
adt.byarm <- cbind(adt.byarm[,c(1)] , (adt.byarm[,-c(1)] - adt_rowmeans_mat_reciprocal))

adt.byarm <- cbind(adt.byarm[,c(1)], 2^adt.byarm[,-c(1)])
adt3_colmeans <- colMeans(adt.byarm[,-c(1)], na.rm = T)
adt3_colmeans_mat_reciprocal <- data.table(1/(t(replicate(nrow(adt.byarm), adt3_colmeans))))
adt3_colmeans_mat_reciprocal[adt3_colmeans_mat_reciprocal == Inf] <- 0 
adt.byarm <- cbind(adt.byarm[,c(1)] , (adt.byarm[,-c(1)] * adt3_colmeans_mat_reciprocal))

#Read in data
golden_samples <- readRDS(sprintf("%s/ML_data/golden_set_ids.rds", dirpath))

#Read in grouped sample information
hand_samples <- readRDS(sprintf("%s/aggregated_results/grouped_control_aneuploidies.rds", dirpath))

#
hand_samples_x2 <- rbind(hand_samples, hand_samples)
hand_samples_p <- paste(hand_samples$chr, "p", sep = "")
hand_samples_q <- paste(hand_samples$chr, "q", sep = "")
hand_samples_x2$chr <- append(hand_samples_p, hand_samples_q)
hand_samples <- hand_samples_x2
setnames(hand_samples, "chr", "arm")
hand_samples <- hand_samples[-which(hand_samples$arm %in% c("13p", "14p", "15p", "22p")),] #remove 13p, 14p, 15p, 22p because they are too small to have genes covering them. 

#change loss and gain examples to hand picked grouped samples
golden_samples$loss <- hand_samples[which(hand_samples$CN == 1),c(1,3)]
golden_samples$gain <- hand_samples[which(hand_samples$CN == 3),c(1,3)]
setnames(golden_samples$loss, old = "ID", new = "sample_id")
setnames(golden_samples$gain, old = "ID", new = "sample_id")

#add chr column to golden samples by adjusting arm strings
#golden_samples$normal$chr <- gsub('.{1}$', '', golden_samples$normal$arm)
#golden_samples$normal$chr <- paste0("chr", golden_samples$normal$chr)
#golden_samples$loss$chr <- paste0("chr", golden_samples$loss$chr)
#golden_samples$gain$chr <- paste0("chr", golden_samples$gain$chr)

#Change name convention
adt_byarm <- adt.byarm
#setnames(adt_byarm, "seqnames", "chr")
adt_byarm$arm <- as.character(adt_byarm$arm)
adt_byarm <- as.data.frame(adt_byarm, col.names = colnames(adt_bychr))

### get adt values for golden set regions ###

#Get loss samples' tpm by arm and by chr
adt_golden_loss_byarm <- c()
for (i in 1:nrow(golden_samples$loss)) {
  adt_golden_loss_byarm <- append(adt_golden_loss_byarm, as.numeric(as.character(adt_byarm[which(adt_byarm$arm %in% golden_samples$loss$arm[i]), golden_samples$loss$sample_id[i]])))
}
golden_samples$loss_byarm <- data.table(cbind(golden_samples$loss[,c("sample_id", "arm")], adt_golden_loss_byarm))
colnames(golden_samples$loss_byarm)[3] <- "tpm"
golden_samples$loss_byarm$tpm <- sapply(golden_samples$loss_byarm$tpm, as.numeric)

#Get loss samples summary stats
golden_samples$loss_stats_byarm <- c()
golden_samples$loss_stats_byarm$mean <- mean(golden_samples$loss_byarm$tpm)
golden_samples$loss_stats_byarm$sd <- sd(golden_samples$loss_byarm$tpm)
golden_samples$loss_stats_byarm$n <- length(golden_samples$loss_byarm$tpm)

#get normal samples' tpm by arm and by chr
adt_golden_normal_byarm <- c()
for (i in 1:nrow(golden_samples$normal)) {
  adt_golden_normal_byarm <- append(adt_golden_normal_byarm, as.numeric(as.character(adt_byarm[which(adt_byarm$arm %in% golden_samples$normal$arm[i]), golden_samples$normal$sample_id[i]])))
}
golden_samples$normal_byarm <- data.table(cbind(golden_samples$normal[,c("sample_id", "arm")], adt_golden_normal_byarm))
colnames(golden_samples$normal_byarm)[3] <- "tpm"
golden_samples$normal_byarm$tpm <- sapply(golden_samples$normal_byarm$tpm, as.numeric)

#Get normal samples summary stats
golden_samples$normal_stats_byarm <- c()
golden_samples$normal_stats_byarm$mean <- mean(golden_samples$normal_byarm$tpm)
golden_samples$normal_stats_byarm$sd <- sd(golden_samples$normal_byarm$tpm)
golden_samples$normal_stats_byarm$n <- length(golden_samples$normal_byarm$tpm)

#Get gain samples' tpm by arm and by chr
adt_golden_gain_byarm <- c()
for (i in 1:nrow(golden_samples$gain)) {
  adt_golden_gain_byarm <- append(adt_golden_gain_byarm, as.numeric(as.character(adt_byarm[which(adt_byarm$arm %in% golden_samples$gain$arm[i]), golden_samples$gain$sample_id[i]])))
}
golden_samples$gain_byarm <- data.table(cbind(golden_samples$gain[,c("sample_id", "arm")], adt_golden_gain_byarm))
colnames(golden_samples$gain_byarm)[3] <- "tpm"
golden_samples$gain_byarm$tpm <- sapply(golden_samples$gain_byarm$tpm, as.numeric)

#Get gain samples summary stats
golden_samples$gain_stats_byarm <- c()
golden_samples$gain_stats_byarm$mean <- mean(golden_samples$gain_byarm$tpm)
golden_samples$gain_stats_byarm$sd <- sd(golden_samples$gain_byarm$tpm)
golden_samples$gain_stats_byarm$n <- length(golden_samples$gain_byarm$tpm)

#Get control samples' tpm by arm and by chr
adt_control_normal_byarm <- data.table() #"sample", "arm", "tpm")
for (i in 1:length(controlSampleIDs)) {
  adt_control_normal_byarm <- rbind( adt_control_normal_byarm, cbind( rep(controlSampleIDs[i], length(adt_byarm$arm)), cbind( data.frame(adt_byarm$arm), data.frame(unname(adt_byarm[,controlSampleIDs[i]])) ) ) )
}
colnames(adt_control_normal_byarm) <- c("sample", "arm", "tpm")
#colnames(adt_control_normal_byarm) <- adt_control_normal_byarm[c(1),];
#adt_control_normal_byarm <- adt_control_normal_byarm[-c(1),]
adt_control_normal_byarm$tpm <- sapply(adt_control_normal_byarm$tpm, as.numeric)
golden_samples$control_byarm <- data.table(adt_control_normal_byarm)

#Get control samples summary stats
golden_samples$control_stats_byarm <- c()
golden_samples$control_stats_byarm$mean <- mean(golden_samples$control_byarm$tpm)
golden_samples$control_stats_byarm$sd <- sd(golden_samples$control_byarm$tpm)
golden_samples$control_stats_byarm$n <- length(golden_samples$control_byarm$tpm)

#Get arm specific control sample summary stats
golden_samples$control_stats_byarm$arm_means <- golden_samples$control_byarm[,-c(1)] %>% group_by(arm) %>% summarise_all(mean, na.rm = TRUE)
golden_samples$control_stats_byarm$arm_sds <- golden_samples$control_byarm[,-c(1)] %>% group_by(arm) %>% summarise_all(sd, na.rm = TRUE)
golden_samples$control_stats_byarm$arm_ns <- golden_samples$control_byarm[,-c(1)] %>% group_by(arm) %>% summarise_all(length)

#save adjusted golden samples object
golden_tpms <- golden_samples
saveRDS(golden_tpms, file = sprintf("%s/aggregated_results/golden_set_tpms_byarm.rds", dirpath))

#############################################################################
### Get an adt-like matrix of z-scores from the labeled sample statistics ###
#############################################################################

#z-score for loss sample statistics
adt_center_loss_byarm <- cbind( adt_byarm[,c(1)], sweep(adt_byarm[,-c(1)], 1, golden_samples$loss_stats_byarm$mean, '-') )
adt_zs_loss_byarm <- cbind( adt_byarm[,c(1)], sweep(adt_center_loss_byarm[,-c(1)], 1, golden_samples$loss_stats_byarm$sd, '/') )

#z-score for normal sample statistics
adt_center_normal_byarm <- cbind( adt_byarm[,c(1)], sweep(adt_byarm[,-c(1)], 1, golden_samples$normal_stats_byarm$mean, '-') )
adt_zs_normal_byarm <- cbind( adt_byarm[,c(1)], sweep(adt_center_normal_byarm[,-c(1)], 1, golden_samples$normal_stats_byarm$sd, '/') )

#z-score for gain sample statistics
adt_center_gain_byarm <- cbind( adt_byarm[,c(1)], sweep(adt_byarm[,-c(1)], 1, golden_samples$gain_stats_byarm$mean, '-') )
adt_zs_gain_byarm <- cbind( adt_byarm[,c(1)], sweep(adt_center_gain_byarm[,-c(1)], 1, golden_samples$gain_stats_byarm$sd, '/') )

#z-score for control samples statistics
adt_center_control_byarm <- cbind( adt_byarm[,c(1)], sweep(adt_byarm[,-c(1)], 1, golden_samples$control_stats_byarm$mean, '-') )
adt_zs_control_byarm <- cbind( adt_byarm[,c(1)], sweep(adt_center_control_byarm[,-c(1)], 1, golden_samples$control_stats_byarm$sd, '/') )

#z-score for control samples statistics
adt_center_control_byarm_specific <- adt_byarm[,-c(1)] - c(golden_samples$control_stats_byarm$arm_means[,2])
adt_zs_control_byarm_specific <- cbind( adt_byarm[,c(1)], adt_center_control_byarm_specific / c(golden_samples$control_stats_byarm$arm_sds[,2]) )

### Get pval matrix by arm ###
pval_matrix_loss_byarm <- c(); pval_matrix_normal_byarm <- c(); pval_matrix_gain_byarm <- c(); pval_matrix_control_byarm <- c(); pval_matrix_control_byarm_specific <- c();
for (i in 1:nrow(adt_byarm)) {
  row_loss <- c(); row_normal <- c(); row_gain <- c(); row_control <- c(); row_control_specific <- c();
  for (j in 2:ncol(adt_byarm)) {
    row_loss <- append( row_loss, 2*pnorm( -abs(adt_zs_loss_byarm[i,j]) ) ) #Should these be one or two tail?
    row_normal <- append( row_normal, 2*pnorm( -abs(adt_zs_normal_byarm[i,j]) ) ) #pval is multiplied by two cause each z-score could be above or below the class specific population stats
    row_gain <- append( row_gain, 2*pnorm( -abs(adt_zs_gain_byarm[i,j]) ) ) #Should these be one or two tail?
    row_control <- append( row_control, 2*pnorm( -abs(adt_zs_control_byarm[i,j]) ) )
    row_control_specific <- append( row_control_specific, 2*pnorm( -abs(adt_zs_control_byarm_specific[i,j]) ))
  }
  pval_matrix_loss_byarm <- cbind(pval_matrix_loss_byarm, data.table(row_loss))
  pval_matrix_normal_byarm <- cbind(pval_matrix_normal_byarm, data.table(row_normal))
  pval_matrix_gain_byarm <- cbind(pval_matrix_gain_byarm, data.table(row_gain))
  pval_matrix_control_byarm <- cbind(pval_matrix_control_byarm, data.table(row_control))
  pval_matrix_control_byarm_specific <- cbind(pval_matrix_control_byarm_specific, data.table(row_control_specific))
}
pval_matrix_loss_byarm <- data.table(t(pval_matrix_loss_byarm))
pval_matrix_normal_byarm <- data.table(t(pval_matrix_normal_byarm))
pval_matrix_gain_byarm <- data.table(t(pval_matrix_gain_byarm))
pval_matrix_control_byarm <- data.table(t(pval_matrix_control_byarm))
pval_matrix_control_byarm_specific <- data.table(t(pval_matrix_control_byarm_specific))

pval_matrix_loss_byarm <- cbind(adt_byarm$arm, pval_matrix_loss_byarm)
pval_matrix_normal_byarm <- cbind(adt_byarm$arm, pval_matrix_normal_byarm)
pval_matrix_gain_byarm <- cbind(adt_byarm$arm, pval_matrix_gain_byarm)
pval_matrix_control_byarm <- cbind(adt_byarm$arm, pval_matrix_control_byarm)
pval_matrix_control_byarm_specific <- cbind(adt_byarm$arm, pval_matrix_control_byarm_specific)

colnames(pval_matrix_loss_byarm) <- colnames(adt_byarm)
colnames(pval_matrix_normal_byarm) <- colnames(adt_byarm)
colnames(pval_matrix_gain_byarm) <- colnames(adt_byarm)
colnames(pval_matrix_control_byarm) <- colnames(adt_byarm)
colnames(pval_matrix_control_byarm_specific) <- colnames(adt_byarm)

saveRDS(pval_matrix_loss_byarm, file = sprintf("%s/aggregated_results/pval_matrix_loss_byarm.rds", dirpath))
saveRDS(pval_matrix_normal_byarm, file = sprintf("%s/aggregated_results/pval_matrix_normal_byarm.rds", dirpath))
saveRDS(pval_matrix_gain_byarm, file = sprintf("%s/aggregated_results/pval_matrix_gain_byarm.rds", dirpath))
saveRDS(pval_matrix_control_byarm, file = sprintf("%s/aggregated_results/pval_matrix_control_byarm.rds", dirpath))
saveRDS(pval_matrix_control_byarm_specific, file = sprintf("%s/aggregated_results/pval_matrix_control_byarm_specific.rds", dirpath))
