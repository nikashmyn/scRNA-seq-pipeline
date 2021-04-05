#Calculating p-values script by Nikos 4/2/2021

#Set path prefixs
#args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/data/processed_bam", "/pellmanlab/nikos/Stam_Etai_Data")
args <- commandArgs(trailingOnly = TRUE)
print(args)
scriptsdir <- args[1]
dirpath <- args[2]
datadir <- args[3]

#Read in data
rsemtpm <- readRDS(file = sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", dirpath))
controlSampleIDs <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))
controlSampleIDs2 <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs2.rds", dirpath))
controlIDs <- readRDS(file = sprintf("%s/aggregated_results/controlIDs.rds", dirpath))
adt <- adt.default <- data.table(readRDS(file = sprintf("%s/aggregated_results/adt.rds", dirpath)))
adt.na <- data.table(readRDS(file = sprintf("%s/aggregated_results/adt.na.rds", dirpath)))
ganno <- data.table(readRDS(file = sprintf("%s/aggregated_results/ganno.rds", dirpath)))
centromeres <- data.table(readRDS(file = sprintf("%s/centromeres.rds", datadir)))
arms <- readRDS(file = sprintf("%s/CN_data/CN_predictions.byarm.rds", dirpath))
golden_samples <- readRDS(sprintf("%s/ML_data/golden_set_ids.rds", dirpath))

#order my data objects 
adt <- adt.default <- adt[order(seqnames, start, end)]
adt <- cbind( adt[,c(1:4)], setcolorder(adt[,-c(1:4)], order(colnames(adt[,-c(1:4)]))) )
adt.na <- adt.na[order(seqnames, start, end)]
adt.na <- cbind( adt.na[,c(1:4)], setcolorder(adt.na[,-c(1:4)], order(colnames(adt.na[,-c(1:4)]))) )
ganno <- ganno[order(seqnames, start, end)]
stopifnot(sum(ganno$id != adt$id) == 0)

#add arm annotation to adt
adt <- cbind(ganno, adt[,-c(1:4)])
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

#Get loss samples' tpm
adt_golden_loss <- c()
for (i in 1:nrow(golden_samples$loss)) {
  adt_golden_loss <- append(adt_golden_loss, as.numeric(as.character(adt_byarm[which(adt_byarm$arm %in% golden_samples$loss$arm[i]), golden_samples$loss$sample_id[i]])))
} 
golden_samples$loss$tpm <- data.table(adt_golden_loss)
loss_stats <- c()
loss_stats$mean <- mean(golden_samples$loss$tpm)
loss_stats$sd <- sd(golden_samples$loss$tpm)
loss_stats$n <- length(golden_samples$loss$tpm)

#get normal samples' tpm
adt_golden_normal <- c()
for (i in 1:nrow(golden_samples$normal)) {
  adt_golden_normal <- append(adt_golden_normal, as.numeric(as.character(adt_byarm[which(adt_byarm$arm %in% golden_samples$normal$arm[i]), golden_samples$normal$sample_id[i]])))
}
golden_samples$normal$tpm <- data.table(adt_golden_normal)
normal_stats <- c()
normal_stats$mean <- mean(golden_samples$normal$tpm)
normal_stats$sd <- sd(golden_samples$normal$tpm)
normal_stats$n <- length(golden_samples$normal$tpm)

#Get gain samples' tpm
adt_golden_gain <- c()
for (i in 1:nrow(golden_samples$gain)) {
  adt_golden_gain <- append(adt_golden_gain, as.numeric(as.character(adt_byarm[which(adt_byarm$arm %in% golden_samples$gain$arm[i]), golden_samples$gain$sample_id[i]])))
}
golden_samples$gain$tpm <- data.table(adt_golden_gain)
gain_stats <- c()
gain_stats$mean <- mean(golden_samples$gain$tpm)
gain_stats$sd <- sd(golden_samples$gain$tpm)
gain_stats$n <- length(golden_samples$gain$tpm)

pval_matrix_loss <- c(); pval_matrix_normal <- c(); pval_matrix_gain <- c()
for (i in 1:nrow(adt_byarm)) {
  row_loss <- c(); row_normal <- c(); row_gain <- c()
  for (j in 2:ncol(adt_byarm)) {
    row_loss <- append(row_loss, 1-pnorm(as.numeric(as.character(adt_byarm[i,j])), mean=loss_stats$mean, sd=loss_stats$sd) )
    row_normal <- append(row_normal, 1-pnorm(as.numeric(as.character(adt_byarm[i,j])), mean=normal_stats$mean, sd=normal_stats$sd) )
    row_gain <- append(row_gain, 1-pnorm(as.numeric(as.character(adt_byarm[i,j])), mean=gain_stats$mean, sd=gain_stats$sd) )
  }
  pval_matrix_loss <- cbind(pval_matrix_loss, data.table(row_loss))
  pval_matrix_normal <- cbind(pval_matrix_normal, data.table(row_normal))
  pval_matrix_gain <- cbind(pval_matrix_gain, data.table(row_gain))
}
pval_matrix_loss <- data.table(t(pval_matrix_loss))
pval_matrix_normal <- data.table(t(pval_matrix_normal))
pval_matrix_gain <- data.table(t(pval_matrix_gain))
pval_matrix_loss <- cbind(adt_byarm$arm, pval_matrix_loss)
pval_matrix_normal <- cbind(adt_byarm$arm, pval_matrix_normal)
pval_matrix_gain <- cbind(adt_byarm$arm, pval_matrix_gain)
colnames(pval_matrix_loss) <- colnames(adt_byarm)
colnames(pval_matrix_normal) <- colnames(adt_byarm)
colnames(pval_matrix_gain) <- colnames(adt_byarm)

#The code framework is done. I need to think about the stats I am calculating now. 
#Could this work as a one-sample t-test? I think no because the distribution i want to use is the population not the sample.
#I could do a two sample t-test to compare two means but I dont have the standard deviation of the sample... wait wait
#I could use the distr for a chr before I bin and compare to the distrs of the unbinned golden data?
t.test(golden_samples$normal$tpm, mu=as.numeric(as.character(adt_byarm[2,2])))

#rownames(pval_matrix) <- rownames(adt_byarm)


