
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
anno <- data.table(read.csv( sprintf("%s/work_in_progress/annotation_list.csv", datadir)))
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

controlSampleIDs_SISg <- anno[ ( event == "control") ]$WTA.plate
controlSampleIDs_SISg <- intersect(controlSampleIDs_SISg, SISg_IDs)
controlSampleIDs_SISg <- controlSampleIDs_SISg[controlSampleIDs_SISg %in% colnames(rsemtpm)]
controlSampleIDs_SISg <- controlSampleIDs_SISg[ which( controlSampleIDs_SISg %in% all_QC$id[ which(all_QC$th5 > 6000)] ) ]

#Data Frames with TPM and Variant information
#adt <- readRDS(sprintf("%s/aggregated_results/adt.rds", dirpath)) #Normalized TPM matrix
coding <- readRDS(sprintf("%s/aggregated_results/ASE.coding.rds", dirpath)) #Raw Variant matrix in the coding and UTR regions

#import raw gene cnts and gene ranges
geneRanges <- readRDS(sprintf("%s/geneRanges.rds", datadir))
rsemtpm <- readRDS(sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", dirpath))

#list of preset parameters for the experiment
configs <- readRDS(sprintf("%s/param_config_list.rds", datadir))
#configs$minDetectionLevel <- 25

#new adt object without gene normalization
adt3 <- fromRawTPMtoExprsRatio(rsemtpm, geneRanges, controlSampleIDs, max_pos_col = 6, plusOne = 1, 
                              zerosAsNA = F, normBySd = F, doNormalizeByControls = F,
                              maxExp = 1200, quantileOfExp = 0.99, minDetectionLevel = configs$minDetectionLevel, minNumOfSamplesToDetect = configs$minNumOfSamplesToDetect )$adt

### Exclude noisy chromosomes a priori ###
#exclude_chrs <- c("chr13", "chr18", "chr21", "chr22", "chrX")
exclude_chrs <- c()
if (length(exclude_chrs) > 0) { adt <- adt[-which(adt$seqnames %in% exclude_chrs),] }

############################################
### Re-Engineer raw variant count Matrix ###
############################################

#reorder adt
adt <- adt[order(seqnames, start)]
ordered_cols <- colnames(adt[,-c(1:6)])[order(colnames(adt[,-c(1:6)]))]
adt <- cbind(adt[,c(1:6)], adt[,..ordered_cols])

#group TPM by chr 
adt.bychr <- adt[,-c(1,3,4,5,6)] %>%
  group_by(seqnames) %>%
  summarise_all(mean, na.rm = TRUE)
adt.bychr <- data.table(adt.bychr)
adt.bychr$seqnames <- sub(pattern = "chr", replacement = "", x = adt.bychr$seqnames)
adt.bychr$seqnames <- as.integer(sub(pattern = "X", replacement = "23", x = adt.bychr$seqnames))
adt.bychr <- adt.bychr[order(seqnames)]

#normalize adt columns
adt_rowmeans <- rowMeans(adt.bychr[,..controlSampleIDs_all], na.rm = T)
adt_rowmeans_mat_reciprocal <- data.table(replicate(ncol(adt.bychr[,-c(1)]), adt_rowmeans))
adt.bychr_all <- cbind(adt.bychr[,c(1)] , (adt.bychr[,-c(1)] - adt_rowmeans_mat_reciprocal))
adt.bychr_all <- cbind(adt.bychr_all[,c(1)], 2^adt.bychr_all[,-c(1)])

#normalize adt columns to old control samples
adt_rowmeans <- rowMeans(adt.bychr[,..controlSampleIDs_old], na.rm = T)
adt_rowmeans_mat_reciprocal <- data.table(replicate(ncol(adt.bychr[,-c(1)]), adt_rowmeans))
adt.bychr_old <- cbind(adt.bychr[,c(1)] , (adt.bychr[,-c(1)] - adt_rowmeans_mat_reciprocal))
adt.bychr_old <- cbind(adt.bychr_old[,c(1)], 2^adt.bychr_old[,-c(1)])

#normalize adt columns to new control samples
adt_rowmeans <- rowMeans(adt.bychr[,..controlSampleIDs_SISg], na.rm = T)
adt_rowmeans_mat_reciprocal <- data.table(replicate(ncol(adt.bychr[,-c(1)]), adt_rowmeans))
adt.bychr_new <- cbind(adt.bychr[,c(1)] , (adt.bychr[,-c(1)] - adt_rowmeans_mat_reciprocal))
adt.bychr_new <- cbind(adt.bychr_new[,c(1)], 2^adt.bychr_new[,-c(1)])

#bind the th5 values with the chr1 TPM ratio values
TPM <- data.table(id = colnames(adt.bychr_old), TPM = unlist(adt.bychr_old[c(2),]) )
th5 <- data.table(all_QC[,c("id", "th5")])
setkey(TPM, "id")
setkey(th5, "id")
dt <- merge(TPM, th5)
colors <- ifelse( test = (dt$id %in% SISg_IDs), "SIS1025g", "other")
dt <- cbind(dt, colors)

#Manually input fragmented experiment IDs
fragmented_exp_ids <- c("210428_8A", "210428_8B", "210428_8C", "210428_8E", "210428_8F", "210428_8G", "210503_10D", "210503_10E", "210503_10F", "210503_10G", "210503_9A", "210503_9B", "210503_9C", "210503_9D", "210503_9E", "210503_9F", "210503_9G", "210503_9H", "210503_10A", "210503_10B", "210503_10C", "210510_2A", "210510_2B", "210510_2C", "210510_2D")

#reduce data table by the manual IDs
dt$colors[which(dt$id %in% fragmented_exp_ids)] <- "Targeted Experiment"

scatter_th5 <- ggplot(data = dt, mapping = aes(x = th5, y = TPM, color = colors)) +
  geom_point(mapping = aes(x = th5, y = TPM, color = colors)) +
  #xlim(c(0,2)) + ylim(c(0,2)) +
  labs(x = "Genes with at least 5 Reads", y = "TPM Ratio", title = sprintf("Expression Scatter Plot"))

scatter_ids <- ggplot(data = dt, mapping = aes(x = id, y = TPM, color = colors)) +
  geom_point(mapping = aes(x = id, y = TPM, color = colors)) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  geom_hline(yintercept=c(.5,1,1.5), linetype='dotted', col = 'grey')+
  #xlim(c(0,2)) + ylim(c(0,2)) + 
  labs(x = "Cells (Chronological)", y = "TPM Ratio", title = sprintf("Expression Scatter Plot"))

grid.arrange(scatter_th5, scatter_ids, nrow=c(2))

### After seeing the log relationship between TPM and th5 we need to model it ###

#first we should make the relationship linear by logging the TPM values. 
dt$log_TPM <- log2(dt$TPM)

#now we can linearly fit the data 
linear_fit <- lm(dt$log_TPM ~ dt$th5)
m <- coef(linear_fit)[2]
b <- coef(linear_fit)[1]

#we can now normalize the log_TPM accounting for the relationship with th5
# log(TPM) - log(e^(m*th5)) / b = 1 => log(TPM) - (m*th5)/ b = 1
dt$log_TPM_shift <- ((dt$log_TPM) - (m*dt$th5)) / b  

#Now we can undo the log
dt$TPM_shift <- 2^(dt$log_TPM_shift)

#Now we can plot the same relationship and see if it worked
scatter_th5 <- ggplot(data = dt, mapping = aes(x = th5, y = dt$log_TPM_shift, color = colors)) +
  geom_point(mapping = aes(x = th5, y = dt$log_TPM_shift, color = colors)) +
  #xlim(c(0,2)) + ylim(c(0,2)) +
  labs(x = "Genes with at least 5 Reads", y = "TPM Ratio", title = sprintf("Expression Scatter Plot"))

scatter_ids <- ggplot(data = dt, mapping = aes(x = id, y = dt$log_TPM_shift, color = colors)) +
  geom_point(mapping = aes(x = id, y = dt$log_TPM_shift, color = colors)) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  geom_hline(yintercept=c(.5,1,1.5), linetype='dotted', col = 'grey')+
  #xlim(c(0,2)) + ylim(c(0,2)) + 
  labs(x = "Cells (Chronological)", y = "TPM Ratio", title = sprintf("Expression Scatter Plot"))

grid.arrange(scatter_th5, scatter_ids, nrow=c(2))
