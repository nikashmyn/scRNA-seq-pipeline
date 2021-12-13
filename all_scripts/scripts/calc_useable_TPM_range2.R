####### Global Variables ########

#This version of the script is the best version.

#args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/rerun_6_9_2021/data", "/pellmanlab/nikos/Stam_Etai_Data", "SIS1025a",  "SIS1025b", "SIS1025c", "SIS1025d", "SIS1025e", "SIS1025f_Lane1", "SIS1025f_Lane2", "SIS1025g_Lane1", "SIS1025g_Lane2", "SIS1025misc", "SIS1025targ")
args <- commandArgs(trailingOnly = TRUE)
print(args)
scriptsdir <- args[1]
wkpath <- dirpath <- args[2]
datadir <- args[3]
experiments <- args[4:length(args)]
mk_agg_dir <- sprintf("mkdir -p %s/aggregated_results", wkpath)
system(mk_agg_dir)

######################
### Source Scripts ###
######################

library(tidyverse)  # data manipulation
source(sprintf('%s/scripts/fromRawTPMtoExprsRatio.R', scriptsdir))

################
#### Imports ###
################

#Annotion list with information on every cell and IDs of control samples
anno <- data.table(read.csv( sprintf("%s/work_in_progress/Annotation_list_long_vnikos_1_9_21.csv", datadir)))
controlSampleIDs <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))

rsemtpm <- readRDS(sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", dirpath))

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
rsemtpm_controls <- rsemtpm[,controlSampleIDs_all]

rsemtpm_control_means <- data.table(gene = rownames(rsemtpm), tpm = rowMeans(rsemtpm_controls))
rsemtpm_control_means <- rsemtpm_control_means[which(rsemtpm_control_means$tpm > 0),]

#Make the column for fraction below a value
dist <- sort(rsemtpm_control_means$tpm)
ecdf_fun <- function(x,perc) ecdf(x)(perc)
rsemtpm_control_means$TPM_quantile <- ecdf_fun(dist, rsemtpm_control_means$tpm)

#Get rid of the top 5% of the values
#upp_lim_QC_TPM <- QC_TPM[which(QC_TPM$TPM_quantile < .95),]

scatter_tpm <- ggplot(data = rsemtpm_control_means, mapping = aes(x = tpm, y = TPM_quantile)) +
  geom_point(mapping = aes(x = tpm, y = TPM_quantile)) +
  #xlim(c(0,100)) + 
  labs(x = "TPM", y = "Quantile of TPM", title = sprintf("TPM vs Quantile of TPM"))


#Now calculate mean and variance of control samples by gene
means <- log2(rowMeans(rsemtpm_controls) + 1)
vars <- log2(rowVars(as.matrix(rsemtpm_controls)) + 1)
genes <- rownames(rsemtpm_controls)
Control_cell_stats <- data.table(id = genes, means = means, vars = vars)
Control_cell_stats <- Control_cell_stats[which(Control_cell_stats$means > 0),]
#Control_cell_stats_red <- Control_cell_stats[which(Control_cell_stats$means < quantile(Control_cell_stats$means, .99)),]
#Control_cell_stats_red <- Control_cell_stats_red[which(Control_cell_stats_red$vars < quantile(Control_cell_stats_red$vars, .99)),]

#plot control cell stats
scatter_tpm2 <- ggplot(data = Control_cell_stats, mapping = aes(x = means, y = vars)) +
  geom_point(mapping = aes(x = means, y = vars)) +
  #xlim(c(0,800)) + ylim(c(0,150000)) +
  labs(x = "log2(<TPM>)", y = "log2(S^2)", title = sprintf("Mean vs Variance of control TPM values by gene in log scale"))


#rank genes by tpm for a control cell
#control_TPM_cell <- data.table(id = rownames(rsemtpm_controls), TPM = rsemtpm_controls[,"210721_4F"])
control_TPM <- copy(rsemtpm_control_means)
setnames(control_TPM, "tpm", "TPM")
setkey(control_TPM, TPM)
control_TPM$TPM <- as.numeric(control_TPM$TPM)
control_TPM$log_tpm <- log2(control_TPM$TPM + 1)
frac <- c()
for (i in 1:length(control_TPM$TPM)) {
  frac <- append(frac, (sum(control_TPM$TPM[1:i])/sum(control_TPM$TPM)) ) 
}
control_TPM$frac <- frac

#plot this
scatter_tpm3 <- ggplot(data = control_TPM, mapping = aes(x = TPM, y = frac)) +
  geom_point(mapping = aes(x = TPM, y = frac)) +
  #xlim(c(0,1000)) + 
  labs(x = "<TPM>", y = "sum <TPM>_i (1-i) / sum <TPM>_i (1-n)", title = sprintf("sum to TPM_i by sum of all TPM as a function of TPM"))


#plot this
scatter_tpm4 <- ggplot(data = control_TPM, mapping = aes(x = log_tpm, y = frac)) +
  geom_point(mapping = aes(x = log_tpm, y = frac)) +
  #xlim(c(0,1000)) + 
  labs(x = "log2(<TPM>)", y = "sum <TPM>_i (1-i) / sum <TPM>_i (1-n)", title = sprintf("sum to TPM_i by sum of all TPM as a function of log2(TPM)"))


#Display plots together
#grid.arrange(scatter_tpm, scatter_tpm2, scatter_tpm3, scatter_tpm4, nrow=c(2), ncol=c(2))





#######################################################




#Based on the new the graphic lets redo these graphic with with a reduced set of genes
rsemtpm_control_means_2 <- rsemtpm_control_means[which(rsemtpm_control_means$tpm > 5 & rsemtpm_control_means$tpm < 2^10),]

#Now calculate mean and variance of control samples by gene
rsemtpm_controls_red <- rsemtpm_controls[which(rownames(rsemtpm_controls) %in% rsemtpm_control_means_2$gene),]
means <- log2(rowMeans(rsemtpm_controls_red) + 1)
vars <- log2(rowVars(as.matrix(rsemtpm_controls_red)) + 1)
genes <- rownames(rsemtpm_controls_red)
Control_cell_stats2 <- data.table(id = genes, means = means, vars = vars)
Control_cell_stats2 <- Control_cell_stats2[which(Control_cell_stats2$means > 0),]
#Control_cell_stats_red <- Control_cell_stats[which(Control_cell_stats$means < quantile(Control_cell_stats$means, .99)),]
#Control_cell_stats_red <- Control_cell_stats_red[which(Control_cell_stats_red$vars < quantile(Control_cell_stats_red$vars, .99)),]

#plot control cell stats
scatter_tpm5 <- ggplot(data = Control_cell_stats2, mapping = aes(x = means, y = vars)) +
  geom_point(mapping = aes(x = means, y = vars)) +
  #xlim(c(0,800)) + ylim(c(0,150000)) +
  labs(x = "log2(<TPM>)", y = "log2(S^2)", title = sprintf("Mean vs Variance limited by mean TPM values (5-1024)"))

#grid.arrange(scatter_tpm, scatter_tpm2, scatter_tpm3, scatter_tpm4, scatter_tpm5, nrow=c(3), ncol=c(2))



#Now select genes mean variance regression
x <- Control_cell_stats$means
y <- Control_cell_stats$vars
df <- data.table(x = x, y = y)

#plot(x, y)
res <- lm(y ~ x) # init = "lts"
#abline(res)

resi <- residuals(res)
stan_resi <- rstandard(res)

#standardized residuals |r| > 2 should be excluded
df_red <- df[which(abs(stan_resi) < 2),]

#plot control cell stats
scatter_tpm6 <- ggplot(data = df_red, mapping = aes(x = x, y = y)) +
  geom_point(mapping = aes(x = x, y = y)) +
  geom_abline(intercept = res$coefficients[1], slope = res$coefficients[2]) +
  #xlim(c(0,800)) + ylim(c(0,150000)) +
  labs(x = "log2(<TPM>)", y = "log2(S^2)", title = sprintf("Mean vs Variance limited by linear fit standardized residuals (|sd(r)| < 2"))

grid.arrange(scatter_tpm, scatter_tpm2, scatter_tpm3, scatter_tpm4, scatter_tpm5, scatter_tpm6, nrow=c(3), ncol=c(2))

#######################################################################



