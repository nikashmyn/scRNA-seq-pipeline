
######################
### Source Scripts ###
######################

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

#SIS1025g Experiment
SISg_IDs <- readRDS(file = sprintf("%s/work_in_progress/SIS1025g_sample_IDs.rds", datadir))
SISg_IDs <- SISg_IDs[SISg_IDs %in% colnames(rsemtpm)]

#QC for all experiments
all_QC <- readRDS(file = sprintf("%s/aggregated_results/all_QC.rds", wkpath))

controlSampleIDs_all <- anno[ (mainGroup == "A_level_control" | mainGroup == "B_level_control" | event == "control")]$WTA.plate
controlSampleIDs_all <- controlSampleIDs_all[controlSampleIDs_all %in% colnames(rsemtpm)]
controlSampleIDs_all <- controlSampleIDs_all[ which( controlSampleIDs_all %in% all_QC$id[ which(all_QC$th5 > 6000)] ) ]

controlSampleIDs_old <- anno[ (mainGroup == "A_level_control" | mainGroup == "B_level_control" ) ]$WTA.plate
controlSampleIDs_old <- controlSampleIDs_old[controlSampleIDs_old %in% colnames(rsemtpm)]
controlSampleIDs_old <- controlSampleIDs_old[ which( controlSampleIDs_old %in% all_QC$id[ which(all_QC$th5 > 6000)] ) ]

#Data Frames with TPM and Variant information
#adt <- readRDS(sprintf("%s/aggregated_results/adt.rds", dirpath)) #Normalized TPM matrix
coding <- readRDS(sprintf("%s/aggregated_results/ASE.coding.rds", dirpath)) #Raw Variant matrix in the coding and UTR regions

#import raw gene cnts and gene ranges
geneRanges <- readRDS(sprintf("%s/geneRanges.rds", datadir))
rsemtpm <- readRDS(sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", dirpath))
rsemtpm_controls <- rsemtpm[,controlSampleIDs_all]

#cols <- sample(ncol(rsemtpm_controls), 200)
rows <- sample(nrow(rsemtpm_controls), 10000)
#
sampled_rsemtpm <- rsemtpm_controls[rows,] #cols

m_sampled_rsemtpm <- data.table(cbind( rep(rownames(sampled_rsemtpm), ncol(sampled_rsemtpm)), melt(sampled_rsemtpm)))
setnames(m_sampled_rsemtpm, c("gene", "id", "tpm"))
m_sampled_rsemtpm_no0 <- m_sampled_rsemtpm[m_sampled_rsemtpm$tpm > 0,]
m_rows <- sample(nrow(m_sampled_rsemtpm_no0), 10000)

#data.table and merge
m_sampledx2_rsemtpm <- data.table(m_sampled_rsemtpm_no0[m_rows,])
sampled_QC <- data.table(all_QC[which(all_QC$id %in% m_sampledx2_rsemtpm$id),])
setkey(m_sampledx2_rsemtpm, "id")
setkey(sampled_QC, "id")
QC_TPM <- merge(sampled_QC, m_sampledx2_rsemtpm)

#Make the column for fraction below a value
dist <- sort(QC_TPM$tpm)
ecdf_fun <- function(x,perc) ecdf(x)(perc)
QC_TPM$TPM_quantile <- ecdf_fun(dist, QC_TPM$tpm)

#Get rid of the top 5% of the values
#upp_lim_QC_TPM <- QC_TPM[which(QC_TPM$TPM_quantile < .95),]

scatter_tpm <- ggplot(data = QC_TPM, mapping = aes(x = tpm, y = TPM_quantile)) +
  geom_point(mapping = aes(x = tpm, y = TPM_quantile)) +
  #xlim(c(0,100)) + 
  labs(x = "TPM", y = "TPM_quantile", title = sprintf("TPM vs Quantile of TPM"))

#Now calculate mean and variance of control samples by gene
means <- rowMeans(rsemtpm_controls)
vars <- rowVars(as.matrix(rsemtpm_controls))
genes <- rownames(rsemtpm_controls)
Control_cell_stats <- data.table(id = genes, means = means, vars = vars)
Control_cell_stats_red <- Control_cell_stats[which(Control_cell_stats$means < quantile(Control_cell_stats$means, .99)),]
Control_cell_stats_red <- Control_cell_stats_red[which(Control_cell_stats_red$vars < quantile(Control_cell_stats_red$vars, .99)),]

#plot control cell stats
scatter_tpm2 <- ggplot(data = Control_cell_stats_red, mapping = aes(x = means, y = vars)) +
  geom_point(mapping = aes(x = means, y = vars)) +
  #xlim(c(0,1000)) + 
  labs(x = "Mean TPM", y = "Variance", title = sprintf("Mean TPM vs TPM Variance by gene in control cells"))


#rank genes by tpm for a control cell
control_TPM_cell <- data.table(id = rownames(rsemtpm_controls), TPM = rsemtpm_controls[,"210721_4F"])
control_TPM_cell <- copy(QC_TPM)
setnames(control_TPM_cell, "tpm", "TPM")
setkey(control_TPM_cell, TPM)
control_TPM_cell_n0 <- control_TPM_cell[which(control_TPM_cell$TPM != 0),]
control_TPM_cell_n0$log_tpm <- log2(control_TPM_cell_n0$TPM + 1)
frac <- c()
for (i in 1:length(control_TPM_cell_n0$TPM)) {
  frac <- append(frac, (sum(control_TPM_cell_n0$TPM[1:i])/sum(control_TPM_cell_n0$TPM)) ) 
}
control_TPM_cell_n0$frac <- frac

#plot this
scatter_tpm3 <- ggplot(data = control_TPM_cell_n0, mapping = aes(x = TPM, y = frac)) +
  geom_point(mapping = aes(x = TPM, y = frac)) +
  #xlim(c(0,1000)) + 
  labs(x = "TPM", y = "sum(TPM) / total(TPM)", title = sprintf("sum of TPM by total TPM as a function of TPM"))


#plot this
scatter_tpm4 <- ggplot(data = control_TPM_cell_n0, mapping = aes(x = log_tpm, y = frac)) +
  geom_point(mapping = aes(x = log_tpm, y = frac)) +
  #xlim(c(0,1000)) + 
  labs(x = "log2(TPM)", y = "sum(TPM) / total(TPM)", title = sprintf("sum of TPM by total TPM as a function of log2(TPM)"))


#Display plots together
grid.arrange(scatter_tpm, scatter_tpm2, scatter_tpm3, scatter_tpm4, nrow=c(2), ncol=c(2))





#######################################################




#Based on the new the graphic lets redo these graphic with with a reduced set of genes
QC_TPM_2 <- QC_TPM[which(QC_TPM$tpm > 2^5 & QC_TPM$tpm < 2^10),]

#Now calculate mean and variance of control samples by gene
rsemtpm_controls_red <- rsemtpm_controls[which(rownames(rsemtpm_controls) %in% QC_TPM_2$gene),]
means <- rowMeans(rsemtpm_controls_red)
vars <- rowVars(as.matrix(rsemtpm_controls_red))
genes <- rownames(rsemtpm_controls_red)
Control_cell_stats <- data.table(id = genes, means = means, vars = vars)
Control_cell_stats_red <- Control_cell_stats[which(Control_cell_stats$means < quantile(Control_cell_stats$means, .99)),]
Control_cell_stats_red <- Control_cell_stats_red[which(Control_cell_stats_red$vars < quantile(Control_cell_stats_red$vars, .99)),]

#plot control cell stats
scatter_tpm5 <- ggplot(data = Control_cell_stats_red, mapping = aes(x = means, y = vars)) +
  geom_point(mapping = aes(x = means, y = vars)) +
  xlim(c(0,250)) + ylim(c(0,6000)) +
  labs(x = "Mean TPM", y = "Variance", title = sprintf("Mean TPM vs TPM Variance by gene in control cells"))

grid.arrange(scatter_tpm, scatter_tpm2, scatter_tpm3, scatter_tpm4, scatter_tpm5, nrow=c(3), ncol=c(2))



#######################################################################



#average of every gene in the control cells
single_cell <- sample(colnames(rsemtpm), 1)
single_cell %in% controlSampleIDs

single_cell <- "210503_10B"

single_cell_tpm <- log2(rsemtpm[,single_cell] + 1)
rsemtpm_controls_means <- log2(rowMeans(rsemtpm_controls) + 1)

df <- data.table(sc_tpm = single_cell_tpm, cnt_tpm = rsemtpm_controls_means)

scatter_tpm6 <- ggplot(data = df, mapping = aes(x = sc_tpm, y = cnt_tpm)) +
  geom_point(mapping = aes(x = sc_tpm, y = cnt_tpm)) +
  #xlim(c(0,10000)) + ylim(c(0,10000)) +
  labs(x = "Mean TPM", y = "Variance", title = sprintf("Mean TPM vs TPM Variance by gene in control cells"))
scatter_tpm6



Control_cell_stats_red$means <- log2(Control_cell_stats_red$means+1)
Control_cell_stats_red$vars <- log2(Control_cell_stats_red$vars+1)

scatter_tpm5 <- ggplot(data = Control_cell_stats_red, mapping = aes(x = means, y = vars)) +
  geom_point(mapping = aes(x = means, y = vars)) +
  #xlim(c(0,250)) + ylim(c(0,6000)) +
  labs(x = "Mean TPM", y = "Variance", title = sprintf("Mean TPM vs TPM Variance by gene in control cells"))


