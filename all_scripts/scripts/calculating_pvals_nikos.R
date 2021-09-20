#Calculating p-values script by Nikos 4/2/2021

#Set path prefixs
#args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/data/processed_bam", "/pellmanlab/nikos/Stam_Etai_Data")
#args <- commandArgs(trailingOnly = TRUE)
#print(args)
#scriptsdir <- args[1]
#dirpath <- args[2]
#datadir <- args[3]

#Read in data
#rsemtpm <- readRDS(file = sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", dirpath))
#controlSampleIDs <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))
#controlSampleIDs2 <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs2.rds", dirpath))
#controlIDs <- readRDS(file = sprintf("%s/aggregated_results/controlIDs.rds", dirpath))
#adt <- adt.default <- data.table(readRDS(file = sprintf("%s/aggregated_results/adt.rds", dirpath)))
#adt.na <- data.table(readRDS(file = sprintf("%s/aggregated_results/adt.na.rds", dirpath)))
#ganno <- data.table(readRDS(file = sprintf("%s/aggregated_results/ganno.rds", dirpath)))
#centromeres <- data.table(readRDS(file = sprintf("%s/centromeres.rds", datadir)))
#arms <- readRDS(file = sprintf("%s/CN_data/CN_predictions.byarm.rds", dirpath))
golden_samples <- readRDS(sprintf("%s/ML_data/golden_set_ids.rds", dirpath))

#add chr column to golden samples by adjusting arm strings
golden_samples$loss$chr <- gsub('.{1}$', '', golden_samples$loss$arm)
golden_samples$loss$chr <- paste0("chr", golden_samples$loss$chr)
golden_samples$normal$chr <- gsub('.{1}$', '', golden_samples$normal$arm)
golden_samples$normal$chr <- paste0("chr", golden_samples$normal$chr)
golden_samples$gain$chr <- gsub('.{1}$', '', golden_samples$gain$arm)
golden_samples$gain$chr <- paste0("chr", golden_samples$gain$chr)

#order my data objects 
adt <- adt.default <- adt[order(seqnames, start, end)]
adt <- cbind( adt[,c(1:4)], setcolorder(adt[,-c(1:4)], order(colnames(adt[,-c(1:4)]))) )
adt.na <- adt.na[order(seqnames, start, end)]
adt.na <- cbind( adt.na[,c(1:4)], setcolorder(adt.na[,-c(1:4)], order(colnames(adt.na[,-c(1:4)]))) )
ganno <- ganno[order(seqnames, start, end)]
stopifnot(sum(ganno$id != adt$id) == 0)

#column center adt | This did not really help.
#adt_noanno <- adt[,-c(1:4)]
#adt_noanno <- scale(adt_noanno, center=T, scale=F)

#exponentiate and divide by col mean
adt_noanno <- adt[,-c(1:4)]
adt_noanno <- 2^adt_noanno

means <- colMeans(adt_noanno)
adt_noanno <- sweep(adt_noanno, 2, means, '/')

#add arm annotation to adt
adt <- cbind(ganno, adt_noanno)
adt_arm <- cbind(adt$arm, adt[,-c(1:9)])
colnames(adt_arm)[1] <- "arm"
adt_chr <- cbind(adt$seqnames, adt[,-c(1:9)])
colnames(adt_chr)[1] <- "chr"

### group adt object by arm & chr ###
adt_byarm <- adt_arm %>% 
               group_by(arm) %>%
                 summarise_all(mean, na.rm = TRUE)

adt_bychr <- adt_chr %>% 
               group_by(chr) %>%
                 summarise_all(mean, na.rm = TRUE)

### get adt values for golden set regions ###

#Get loss samples' tpm by arm and by chr
adt_golden_loss_byarm <- c()
adt_golden_loss_bychr <- c()
for (i in 1:nrow(golden_samples$loss)) {
  adt_golden_loss_byarm <- append(adt_golden_loss_byarm, as.numeric(as.character(adt_byarm[which(adt_byarm$arm %in% golden_samples$loss$arm[i]), golden_samples$loss$sample_id[i]])))
  adt_golden_loss_bychr <- append(adt_golden_loss_bychr, as.numeric(as.character(adt_bychr[which(adt_bychr$chr %in% golden_samples$loss$chr[i]), golden_samples$loss$sample_id[i]])))
}
golden_samples$loss_byarm <- data.table(cbind(golden_samples$loss[,c("sample_id", "arm")], adt_golden_loss_byarm))
golden_samples$loss_bychr <- data.table(cbind(golden_samples$loss[,c("sample_id", "chr")], adt_golden_loss_bychr))
colnames(golden_samples$loss_byarm)[3] <- "tpm"
colnames(golden_samples$loss_bychr)[3] <- "tpm"
golden_samples$loss_byarm$tpm <- sapply(golden_samples$loss_byarm$tpm, as.numeric)
golden_samples$loss_bychr$tpm <- sapply(golden_samples$loss_bychr$tpm, as.numeric)

#Get loss samples summary stats
golden_samples$loss_stats_byarm <- c()
golden_samples$loss_stats_bychr <- c()
golden_samples$loss_stats_byarm$mean <- mean(golden_samples$loss_byarm$tpm)
golden_samples$loss_stats_bychr$mean <- mean(golden_samples$loss_bychr$tpm)
golden_samples$loss_stats_byarm$sd <- sd(golden_samples$loss_byarm$tpm)
golden_samples$loss_stats_bychr$sd <- sd(golden_samples$loss_bychr$tpm)
golden_samples$loss_stats_byarm$n <- length(golden_samples$loss_byarm$tpm)
golden_samples$loss_stats_bychr$n <- length(golden_samples$loss_bychr$tpm)

#get normal samples' tpm by arm and by chr
adt_golden_normal_byarm <- c()
adt_golden_normal_bychr <- c()
for (i in 1:nrow(golden_samples$normal)) {
  adt_golden_normal_byarm <- append(adt_golden_normal_byarm, as.numeric(as.character(adt_byarm[which(adt_byarm$arm %in% golden_samples$normal$arm[i]), golden_samples$normal$sample_id[i]])))
  adt_golden_normal_bychr <- append(adt_golden_normal_bychr, as.numeric(as.character(adt_bychr[which(adt_bychr$chr %in% golden_samples$normal$chr[i]), golden_samples$normal$sample_id[i]])))
}
golden_samples$normal_byarm <- data.table(cbind(golden_samples$normal[,c("sample_id", "arm")], adt_golden_normal_byarm))
golden_samples$normal_bychr <- data.table(cbind(golden_samples$normal[,c("sample_id", "chr")], adt_golden_normal_bychr))
colnames(golden_samples$normal_byarm)[3] <- "tpm"
colnames(golden_samples$normal_bychr)[3] <- "tpm"
golden_samples$normal_byarm$tpm <- sapply(golden_samples$normal_byarm$tpm, as.numeric)
golden_samples$normal_bychr$tpm <- sapply(golden_samples$normal_bychr$tpm, as.numeric)

#Get normal samples summary stats
golden_samples$normal_stats_byarm <- c()
golden_samples$normal_stats_bychr <- c()
golden_samples$normal_stats_byarm$mean <- mean(golden_samples$normal_byarm$tpm)
golden_samples$normal_stats_bychr$mean <- mean(golden_samples$normal_bychr$tpm)
golden_samples$normal_stats_byarm$sd <- sd(golden_samples$normal_byarm$tpm)
golden_samples$normal_stats_bychr$sd <- sd(golden_samples$normal_bychr$tpm)
golden_samples$normal_stats_byarm$n <- length(golden_samples$normal_byarm$tpm)
golden_samples$normal_stats_bychr$n <- length(golden_samples$normal_bychr$tpm)

#Get gain samples' tpm by arm and by chr
adt_golden_gain_byarm <- c()
adt_golden_gain_bychr <- c()
for (i in 1:nrow(golden_samples$gain)) {
  adt_golden_gain_byarm <- append(adt_golden_gain_byarm, as.numeric(as.character(adt_byarm[which(adt_byarm$arm %in% golden_samples$gain$arm[i]), golden_samples$gain$sample_id[i]])))
  adt_golden_gain_bychr <- append(adt_golden_gain_bychr, as.numeric(as.character(adt_bychr[which(adt_bychr$chr %in% golden_samples$gain$chr[i]), golden_samples$gain$sample_id[i]])))
}
golden_samples$gain_byarm <- data.table(cbind(golden_samples$gain[,c("sample_id", "arm")], adt_golden_gain_byarm))
golden_samples$gain_bychr <- data.table(cbind(golden_samples$gain[,c("sample_id", "chr")], adt_golden_gain_bychr))
colnames(golden_samples$gain_byarm)[3] <- "tpm"
colnames(golden_samples$gain_bychr)[3] <- "tpm"
golden_samples$gain_byarm$tpm <- sapply(golden_samples$gain_byarm$tpm, as.numeric)
golden_samples$gain_bychr$tpm <- sapply(golden_samples$gain_bychr$tpm, as.numeric)

#Get gain samples summary stats
golden_samples$gain_stats_byarm <- c()
golden_samples$gain_stats_bychr <- c()
golden_samples$gain_stats_byarm$mean <- mean(golden_samples$gain_byarm$tpm)
golden_samples$gain_stats_bychr$mean <- mean(golden_samples$gain_bychr$tpm)
golden_samples$gain_stats_byarm$sd <- sd(golden_samples$gain_byarm$tpm)
golden_samples$gain_stats_bychr$sd <- sd(golden_samples$gain_bychr$tpm)
golden_samples$gain_stats_byarm$n <- length(golden_samples$gain_byarm$tpm)
golden_samples$gain_stats_bychr$n <- length(golden_samples$gain_bychr$tpm)

#Get control samples' tpm by arm and by chr
adt_control_normal_byarm <- c("sample", "arm", "tpm")
adt_control_normal_bychr <- c("sample", "chr", "tpm")
for (i in 1:length(controlSampleIDs)) {
  adt_control_normal_byarm <- rbind( adt_control_normal_byarm, cbind( controlSampleIDs[i], cbind( data.frame(adt_byarm$arm), data.frame(unname(adt_byarm[,controlSampleIDs[i]])) ) ) )
  adt_control_normal_bychr <- rbind( adt_control_normal_bychr, cbind( controlSampleIDs[i], cbind( data.frame(adt_bychr$chr), data.frame(unname(adt_bychr[,controlSampleIDs[i]])) ) ) )
}
colnames(adt_control_normal_byarm) <- adt_control_normal_byarm[c(1),]; adt_control_normal_byarm <- adt_control_normal_byarm[-c(1),]
colnames(adt_control_normal_bychr) <- adt_control_normal_bychr[c(1),]; adt_control_normal_bychr <- adt_control_normal_bychr[-c(1),]
adt_control_normal_byarm$tpm <- sapply(adt_control_normal_byarm$tpm, as.numeric)
adt_control_normal_bychr$tpm <- sapply(adt_control_normal_bychr$tpm, as.numeric)
golden_samples$control_byarm <- data.table(adt_control_normal_byarm)
golden_samples$control_bychr <- data.table(adt_control_normal_bychr)

#Get control samples summary stats
golden_samples$control_stats_byarm <- c()
golden_samples$control_stats_bychr <- c()
golden_samples$control_stats_byarm$mean <- mean(golden_samples$control_byarm$tpm)
golden_samples$control_stats_byarm$sd <- sd(golden_samples$control_byarm$tpm)
golden_samples$control_stats_byarm$n <- length(golden_samples$control_byarm$tpm)
golden_samples$control_stats_bychr$mean <- mean(golden_samples$control_bychr$tpm)
golden_samples$control_stats_bychr$sd <- sd(golden_samples$control_bychr$tpm)
golden_samples$control_stats_bychr$n <- length(golden_samples$control_bychr$tpm)

#save adjusted golden samples object
saveRDS(golden_samples, file = sprintf("%s/aggregated_results/golden_set_tpms.rds", dirpath))

#############################################################################
### Get an adt-like matrix of z-scores from the labeled sample statistics ###
#############################################################################

#z-score for loss sample statistics
adt_center_loss_byarm <- cbind( adt_byarm[,c(1)], sweep(adt_byarm[,-c(1)], 1, golden_samples$loss_stats_byarm$mean, '-') )
adt_zs_loss_byarm <- cbind( adt_byarm[,c(1)], sweep(adt_center_loss_byarm[,-c(1)], 1, golden_samples$loss_stats_byarm$sd, '/') )

adt_center_loss_bychr <- cbind( adt_bychr[,c(1)], sweep(adt_bychr[,-c(1)], 1, golden_samples$loss_stats_bychr$mean, '-') )
adt_zs_loss_bychr <- cbind( adt_bychr[,c(1)], sweep(adt_center_loss_bychr[,-c(1)], 1, golden_samples$loss_stats_bychr$sd, '/') )

#for normal sample statistics
adt_center_normal_byarm <- cbind( adt_byarm[,c(1)], sweep(adt_byarm[,-c(1)], 1, golden_samples$normal_stats_byarm$mean, '-') )
adt_zs_normal_byarm <- cbind( adt_byarm[,c(1)], sweep(adt_center_normal_byarm[,-c(1)], 1, golden_samples$normal_stats_byarm$sd, '/') )

adt_center_normal_bychr <- cbind( adt_bychr[,c(1)], sweep(adt_bychr[,-c(1)], 1, golden_samples$normal_stats_bychr$mean, '-') )
adt_zs_normal_bychr <- cbind( adt_bychr[,c(1)], sweep(adt_center_normal_bychr[,-c(1)], 1, golden_samples$normal_stats_bychr$sd, '/') )

#for gain sample statistics
adt_center_gain_byarm <- cbind( adt_byarm[,c(1)], sweep(adt_byarm[,-c(1)], 1, golden_samples$gain_stats_byarm$mean, '-') )
adt_zs_gain_byarm <- cbind( adt_byarm[,c(1)], sweep(adt_center_gain_byarm[,-c(1)], 1, golden_samples$gain_stats_byarm$sd, '/') )

adt_center_gain_bychr <- cbind( adt_bychr[,c(1)], sweep(adt_bychr[,-c(1)], 1, golden_samples$gain_stats_bychr$mean, '-') )
adt_zs_gain_bychr <- cbind( adt_bychr[,c(1)], sweep(adt_center_gain_bychr[,-c(1)], 1, golden_samples$gain_stats_bychr$sd, '/') )

#got control samples statistics
adt_center_control_byarm <- cbind( adt_byarm[,c(1)], sweep(adt_byarm[,-c(1)], 1, golden_samples$control_stats_byarm$mean, '-') )
adt_zs_control_byarm <- cbind( adt_byarm[,c(1)], sweep(adt_center_control_byarm[,-c(1)], 1, golden_samples$control_stats_byarm$sd, '/') )

adt_center_control_bychr <- cbind( adt_bychr[,c(1)], sweep(adt_bychr[,-c(1)], 1, golden_samples$control_stats_bychr$mean, '-') )
adt_zs_control_bychr <- cbind( adt_bychr[,c(1)], sweep(adt_center_control_bychr[,-c(1)], 1, golden_samples$control_stats_bychr$sd, '/') )

### Get pval matrix by arm ###
pval_matrix_loss_byarm <- c(); pval_matrix_normal_byarm <- c(); pval_matrix_gain_byarm <- c(); pval_matrix_control_byarm <- c();
for (i in 1:nrow(adt_byarm)) {
  row_loss <- c(); row_normal <- c(); row_gain <- c(); row_control <- c();
  for (j in 2:ncol(adt_byarm)) {
    row_loss <- append( row_loss, 2*pnorm( -abs(adt_zs_loss_byarm[i,j]) ) ) #Should these be one or two tail?
    row_normal <- append( row_normal, 2*pnorm( -abs(adt_zs_normal_byarm[i,j]) ) ) #pval is multiplied by two cause each z-score could be above or below the class specific population stats
    row_gain <- append( row_gain, 2*pnorm( -abs(adt_zs_gain_byarm[i,j]) ) ) #Should these be one or two tail?
    row_control <- append( row_control, 2*pnorm( -abs(adt_zs_control_byarm[i,j]) ) )
  }
  pval_matrix_loss_byarm <- cbind(pval_matrix_loss_byarm, data.table(row_loss))
  pval_matrix_normal_byarm <- cbind(pval_matrix_normal_byarm, data.table(row_normal))
  pval_matrix_gain_byarm <- cbind(pval_matrix_gain_byarm, data.table(row_gain))
  pval_matrix_control_byarm <- cbind(pval_matrix_control_byarm, data.table(row_control))
}
pval_matrix_loss_byarm <- data.table(t(pval_matrix_loss_byarm))
pval_matrix_normal_byarm <- data.table(t(pval_matrix_normal_byarm))
pval_matrix_gain_byarm <- data.table(t(pval_matrix_gain_byarm))
pval_matrix_control_byarm <- data.table(t(pval_matrix_control_byarm))
pval_matrix_loss_byarm <- cbind(adt_byarm$arm, pval_matrix_loss_byarm)
pval_matrix_normal_byarm <- cbind(adt_byarm$arm, pval_matrix_normal_byarm)
pval_matrix_gain_byarm <- cbind(adt_byarm$arm, pval_matrix_gain_byarm)
pval_matrix_control_byarm <- cbind(adt_byarm$arm, pval_matrix_control_byarm)
colnames(pval_matrix_loss_byarm) <- colnames(adt_byarm)
colnames(pval_matrix_normal_byarm) <- colnames(adt_byarm)
colnames(pval_matrix_gain_byarm) <- colnames(adt_byarm)
colnames(pval_matrix_control_byarm) <- colnames(adt_byarm)

saveRDS(pval_matrix_loss_byarm, file = sprintf("%s/aggregated_results/pval_matrix_loss_byarm.rds", dirpath))
saveRDS(pval_matrix_normal_byarm, file = sprintf("%s/aggregated_results/pval_matrix_normal_byarm.rds", dirpath))
saveRDS(pval_matrix_gain_byarm, file = sprintf("%s/aggregated_results/pval_matrix_gain_byarm.rds", dirpath))
saveRDS(pval_matrix_control_byarm, file = sprintf("%s/aggregated_results/pval_matrix_control_byarm.rds", dirpath))

### Get pval matrix by chr ###
pval_matrix_loss_bychr <- c(); pval_matrix_normal_bychr <- c(); pval_matrix_gain_bychr <- c()
for (i in 1:nrow(adt_bychr)) {
  row_loss <- c(); row_normal <- c(); row_gain <- c()
  for (j in 2:ncol(adt_bychr)) {
    row_loss <- append( row_loss, 2*pnorm( -abs(adt_zs_loss_bychr[i,j]) ) ) #Should these be one or two tail?
    row_normal <- append( row_normal, 2*pnorm( -abs(adt_zs_normal_bychr[i,j]) ) ) #pval is multiplied by two cause each z-score could be above or below the class specific population stats
    row_gain <- append( row_gain, 2*pnorm( -abs(adt_zs_gain_bychr[i,j]) ) ) #Should these be one or two tail?
  }
  pval_matrix_loss_bychr <- cbind(pval_matrix_loss_bychr, data.table(row_loss))
  pval_matrix_normal_bychr <- cbind(pval_matrix_normal_bychr, data.table(row_normal))
  pval_matrix_gain_bychr <- cbind(pval_matrix_gain_bychr, data.table(row_gain))
}
pval_matrix_loss_bychr <- data.table(t(pval_matrix_loss_bychr))
pval_matrix_normal_bychr <- data.table(t(pval_matrix_normal_bychr))
pval_matrix_gain_bychr <- data.table(t(pval_matrix_gain_bychr))
pval_matrix_loss_bychr <- cbind(adt_bychr$chr, pval_matrix_loss_bychr)
pval_matrix_normal_bychr <- cbind(adt_bychr$chr, pval_matrix_normal_bychr)
pval_matrix_gain_bychr <- cbind(adt_bychr$chr, pval_matrix_gain_bychr)
colnames(pval_matrix_loss_bychr) <- colnames(adt_bychr)
colnames(pval_matrix_normal_bychr) <- colnames(adt_bychr)
colnames(pval_matrix_gain_bychr) <- colnames(adt_bychr)

saveRDS(pval_matrix_loss_bychr, file = sprintf("%s/aggregated_results/pval_matrix_loss_bychr.rds", dirpath))
saveRDS(pval_matrix_normal_bychr, file = sprintf("%s/aggregated_results/pval_matrix_normal_bychr.rds", dirpath))
saveRDS(pval_matrix_gain_bychr, file = sprintf("%s/aggregated_results/pval_matrix_gain_bychr.rds", dirpath))

#We can use the pval_martix_normal to show arms that are clearly not a part of the normal CN distribution. 
#H_0 = tpm for an arm belongs to the normal CN dist.
#H_a = tpm for an arm does not belong to the normal CN dist
#low p-vals show arms where we should reject the null hypothesis in favor of the alternative. 
#Or a low p-val shows an arm most likely does not belong to the normal CN dist. 
#aka an arm with a low pval is likely an abnormal CN.

##############################
### Get significant chroms ###
##############################

pval_matrix_normal_bychr_noname <- pval_matrix_normal_bychr[,-c(1)]
setkey(pval_matrix_normal_bychr_noname)
inds <- which(pval_matrix_normal_bychr_noname < .05, arr.ind = TRUE, useNames = TRUE)
rows <- rownames(pval_matrix_normal_bychr_noname)[inds[,1]]
cols <- colnames(pval_matrix_normal_bychr_noname)[inds[,2]]
non_diploid_chrs <- cbind(sample = cols, chr = rows) #already sorted
non_diploid_chrs_vals <- c()
for (i in 1:nrow(inds)) {
  ind1 <- inds[i,1]; ind2 <- inds[i,2];
  tmp <- unlist(pval_matrix_normal_bychr_noname[ind1,])
  tmp2 <- tmp[ind2]
  non_diploid_chrs_vals <- append(non_diploid_chrs_vals, tmp2)
}
non_diploid_chrs_w_vals <- cbind(non_diploid_chrs, non_diploid_chrs_vals)
rownames(non_diploid_chrs_w_vals) <- NULL
colnames(non_diploid_chrs_w_vals)[3] <- "Pval"

saveRDS(non_diploid_chrs_w_vals, file = sprintf("%s/aggregated_results/non_diploid_chromosomes.rds", dirpath))

###############
### VISUALS ###
###############

#destDir <- sprintf("%s/visual_results/pval_plots", dirpath)
#system(sprintf("mkdir -p %s", destDir))
#
##main data to be plotted byarm
#boxplots.vals <- c()
#boxplots.vals$loss <- golden_samples$loss_byarm$tpm; boxplots.vals$normal <- golden_samples$normal_byarm$tpm; boxplots.vals$gain <- golden_samples$gain_byarm$tpm; boxplots.vals$control <- golden_samples$control_byarm$tpm;
#dist_colors <- c("deepskyblue4", "darkorchid4", "darkred", "darkgreen")
#point_colors <- c("green", "yellow")
#
##Plot figures at arm level
#for ( i in 2:ncol(pval_matrix_normal_byarm) ) { #2:3) { #
#  myid <- colnames(pval_matrix_normal_byarm)[i]
#  pdfFile <- sprintf("%s/%s.pval_boxplot_byarm.pdf", destDir, myid)
#  pdf(file = pdfFile, width = 10, height = 14)
#  for ( j in 1:nrow(pval_matrix_normal_byarm) ) {
#    boxplot( boxplots.vals, xlab="CN State", ylab="Normalized Ratio of Arm Average TPM",
#             pch = 21, bg = dist_colors, axes = F, frame.plot = FALSE, col = dist_colors)
#    grid(col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
#    mtext(side = 1, text = names(boxplots.vals), at = c(1,2,3,4),
#          col = "grey20", line = 1, cex = 0.9)
#    mtext(side = 2, text = c(0.5, 1.0, 1.5, 2.0), at = c(0.5, 1.0, 1.5, 2.0),
#          col = "grey20", line = 1, cex = 0.9)
#    abline(a = as.numeric(as.character(adt_byarm[j,i])), b=0, col="yellow")
#    title( sprintf("%s | %s", colnames(pval_matrix_normal_byarm[j,..i]), pval_matrix_normal_byarm[j,1]) )
#    text(x=1, y = 2.15, labels = sprintf("P-Value = %s", signif(pval_matrix_loss_byarm[j,..i], digits = 2)), col="black", cex=.8)
#    text(x=2, y = 2.15, labels = sprintf("P-Value = %s", signif(pval_matrix_normal_byarm[j,..i], digits = 2)), col="black", cex=.8 )
#    text(x=3, y = .15, labels = sprintf("P-Value = %s", signif(pval_matrix_gain_byarm[j,..i], digits = 2)), col="black", cex=.8 )
#    text(x=4, y = .15, labels = sprintf("P-Value = %s", signif(pval_matrix_control_byarm[j,..i], digits = 2)), col="black", cex=.8 )
#  }
#  dev.off()
#}
#
##main data to be plotted bychr
#boxplots.vals <- c()
#boxplots.vals$loss <- golden_samples$loss_bychr$tpm; boxplots.vals$normal <- golden_samples$normal_bychr$tpm; boxplots.vals$gain <- golden_samples$gain_bychr$tpm; boxplots.vals$control <- golden_samples$control_bychr$tpm;
#plot_colors <- c("deepskyblue4", "darkorchid4", "darkred", "darkgreen")
#
##Plot figures at chr level
#for ( i in 2:ncol(pval_matrix_normal_bychr) ) { #2:3) {
#  myid <- colnames(pval_matrix_normal_bychr)[i]
#  pdfFile <- sprintf("%s/%s.pval_boxplot_bychr.pdf", destDir, myid)
#  pdf(file = pdfFile, width = 10, height = 14)
#  for ( j in 1:nrow(pval_matrix_normal_bychr) ) {
#    boxplot( boxplots.vals, xlab="CN State", ylab="Normalized Ratio of Chr Average TPM",
#             pch = 21, bg = dist_colors, axes = F, frame.plot = FALSE, col = dist_colors)
#    grid(col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
#    mtext(side = 1, text = names(boxplots.vals), at = c(1,2,3),
#          col = "grey20", line = 1, cex = 0.9)
#    mtext(side = 2, text = c(0.5, 1.0, 1.5, 2.0), at = c(0.5, 1.0, 1.5, 2.0),
#          col = "grey20", line = 1, cex = 0.9)
#    abline(a = as.numeric(as.character(adt_bychr[j,i])), b=0, col="green")
#    title( sprintf("%s | %s", colnames(pval_matrix_normal_bychr[j,..i]), pval_matrix_normal_bychr[j,1]) )
#    text(x=1, y = 2.1, labels = sprintf("P-Value = %s", signif(pval_matrix_loss_bychr[j,..i], digits = 2)), col="black", cex=.8)
#    text(x=2, y = 2.1, labels = sprintf("P-Value = %s", signif(pval_matrix_normal_bychr[j,..i], digits = 2)), col="black", cex=.8 )
#    text(x=3, y = .25, labels = sprintf("P-Value = %s", signif(pval_matrix_gain_bychr[j,..i], digits = 2)), col="black", cex=.8 )
#  }
#  dev.off()
#}