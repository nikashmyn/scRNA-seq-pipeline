###################################################
### 2D Contour Map (TPM Ratio vs Variant Ratio) ###
###################################################

########################
### Passed Arguments ###
########################

args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/rerun_6_9_2021/data", "/pellmanlab/nikos/Stam_Etai_Data")
#args <- commandArgs(trailingOnly = TRUE)
print(args)
scriptsdir <- args[1]
dirpath <- args[2]
datadir <- args[3]
mk_agg_dir <- sprintf("mkdir -p %s/aggregated_results", dirpath)
system(mk_agg_dir)

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
anno <- data.table(read.csv(sprintf("%s/work_in_progress/Annotation_list_long_vnikos_11_6_21.csv", datadir) ))
controlSampleIDs <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))

#control samples with aneuploidies
#c("161228_A2", "")

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
adt <- fromRawTPMtoExprsRatio(rsemtpm, geneRanges, controlSampleIDs, max_pos_col = 6, plusOne = 1, 
                              zerosAsNA = F, normBySd = F, doNormalizeByControls = F,
                              maxExp = 1200, quantileOfExp = 0.99, minDetectionLevel = 5, minNumOfSamplesToDetect = 5 )$adt

#Var cnts per 50 genomically sequential genes for allele A and B
var_mat_A <- coding$cnts.A
var_mat_B <- coding$cnts.B
var_mat <- cbind(var_mat_A[,c(1:2)], var_mat_A[,-c(1:2)] + var_mat_B[,-c(1:2)])

### Exclude noisy chromosomes a priori ###
exclude_chrs <- c("chr13", "chr18", "chr21", "chr22", "chrX")
adt <- adt[-which(adt$seqnames %in% exclude_chrs),]
var_mat <- var_mat[-which(var_mat$seqnames %in% exclude_chrs),]
var_mat$seqnames <- sub(pattern = "chr", replacement = "", x = var_mat$seqnames)
var_mat$seqnames <- as.integer(sub(pattern = "X", replacement = "23", x = var_mat$seqnames))
var_mat <- var_mat[order(seqnames, bin)]

############################################
### Re-Engineer raw variant count Matrix ###
############################################

#Bin adt 
adt <- adt[order(seqnames, start)]

#group TPM by chr 
adt.bychr <- adt[,-c(1,3,4,5,6)] %>%
  group_by(seqnames) %>%
  summarise_all(mean, na.rm = TRUE)
adt.bychr <- data.table(adt.bychr)
adt.bychr$seqnames <- sub(pattern = "chr", replacement = "", x = adt.bychr$seqnames)
adt.bychr$seqnames <- as.integer(sub(pattern = "X", replacement = "23", x = adt.bychr$seqnames))
adt.bychr <- adt.bychr[order(seqnames)]

#group var cnts by chr
var_mat.bychr <- var_mat[,-c(2)] %>%
  group_by(seqnames) %>%
  summarise_all(mean, na.rm = TRUE)
#var_mat.bychr$seqnames <- as.numeric(c(1:23))
var_mat.bychr <- data.table(var_mat.bychr)
var_mat.bychr <- var_mat.bychr[order(seqnames)]

#normalize adt columns
adt_rowmeans <- rowMeans(adt.bychr[,..controlSampleIDs], na.rm = T)
adt_rowsds <- rowSds(as.matrix(adt.bychr[,..controlSampleIDs]), na.rm = T)
adt_rowmeans_mat_reciprocal <- data.table(replicate(ncol(adt.bychr[,-c(1)]), adt_rowmeans))
adt_rowsds_mat_reciprocal <- data.table(replicate(ncol(adt.bychr[,-c(1)]), adt_rowsds))
adt.bychr <- cbind(adt.bychr[,c(1)] , (adt.bychr[,-c(1)] - adt_rowmeans_mat_reciprocal))
adt.bychr <- cbind(adt.bychr[,c(1)] , (adt.bychr[,-c(1)] / adt_rowsds_mat_reciprocal))

adt3_colmeans <- colMeans(adt.bychr[,-c(1)], na.rm = T)
adt3_colsds <- colSds(as.matrix(adt.bychr[,-c(1)]), na.rm = T)
adt3_colmeans_mat_reciprocal <- data.table((t(replicate(nrow(adt.bychr), adt3_colmeans))))
adt3_colsds_mat_reciprocal <- data.table((t(replicate(nrow(adt.bychr), adt3_colsds))))
adt.bychr <- cbind(adt.bychr[,c(1)] , (adt.bychr[,-c(1)] - adt3_colmeans_mat_reciprocal))
adt.bychr <- cbind(adt.bychr[,c(1)] , (adt.bychr[,-c(1)] / adt3_colsds_mat_reciprocal))

#Exponentiate to unlog
adt.bychr <- cbind(adt.bychr[,c(1)], 2^adt.bychr[,-c(1)])

#Normalize variant dataframe
dt2go <- copy(var_mat.bychr[, -c(1)])

#log before mean:
dt2go <- log2(dt2go + 1)

var_rowmeans <- rowMeans(dt2go[,..controlSampleIDs], na.rm = T)
var_rowsds <- rowSds(as.matrix(dt2go[,..controlSampleIDs]), na.rm = T)
var_rowmeans_mat_reciprocal <- data.table((replicate(ncol(dt2go), var_rowmeans)))
var_rowsds_mat_reciprocal <- data.table((replicate(ncol(dt2go), var_rowsds)))
var_mat.bychr <- cbind(var_mat.bychr[,c(1)], (dt2go - var_rowmeans_mat_reciprocal))
var_mat.bychr <- cbind(var_mat.bychr[,c(1)], (var_mat.bychr[,-c(1)] / var_rowsds_mat_reciprocal))

var_colmeans <- colMeans(var_mat.bychr[,-c(1)], na.rm = T)
var_colsds <- colSds(as.matrix(var_mat.bychr[,-c(1)]), na.rm = T)
var_colmeans_mat_reciprocal <- data.table((t(replicate(nrow(var_mat.bychr), var_colmeans))))
var_colsds_mat_reciprocal <- data.table((t(replicate(nrow(var_mat.bychr), var_colsds))))
var_mat.bychr <- cbind(var_mat.bychr[,c(1)] , (var_mat.bychr[,-c(1)] - var_colmeans_mat_reciprocal))
var_mat.bychr <- cbind(var_mat.bychr[,c(1)] , (var_mat.bychr[,-c(1)] / var_colsds_mat_reciprocal))

#Exponentiate to unlog
var_mat.bychr <- cbind(var_mat.bychr[,c(1)], 2^var_mat.bychr[,-c(1)])

#check cols are the same
colIDs <- intersect(colnames(var_mat.bychr), colnames(adt.bychr))

#Check that binning is the same
adt.bychr <- adt.bychr[,..colIDs]; var_mat.bychr <- var_mat.bychr[,..colIDs];
stopifnot(dim(var_mat.bychr)[2] == dim(adt.bychr)[2])

###################
### Contour Map ###
###################

#get families from anno list
setkey(anno, WTA.plate)
families <- sort(table(anno[Pairs != "NA" & Pairs != "mother"]$Pairs), decreasing = T)

#Set seed for reproducability
set.seed(123)

#Prepare Visual Objects
myfamily = names(families)[24] #43, 24, 142, 30, 25
message(myfamily)
myids <- anno[Pairs %in% myfamily]$WTA.plate
visual.data <- data.table(TPM = as.numeric(flatten(adt.bychr[,..myids])))
visual.data$VAR <- as.numeric(flatten(var_mat.bychr[,..myids]))
visual.data$chr <- as.factor(rep(unique(adt.bychr$seqnames), length = length(visual.data$VAR)))
visual.data$cell <- rep(myids, each=23-length(exclude_chrs))

#Second Draft Plots
base.visual <- ggplot(data = visual.data, mapping = aes(x = TPM, y = VAR, color = chr, shape = cell)) +
  geom_point(mapping = aes(x = TPM, y = VAR, color = chr, shape = cell)) +
  xlim(c(0,2)) + ylim(c(0,2)) +
  labs(x = "TPM Ratio", y = "SNP Count Ratio", title = sprintf("Expression Contour Plot | %s ", myfamily))
#geom_abline(aes(intercept = 0, slope = 2/3, alpha = .25)) +  # linetype = "dotted", size = .1
#geom_abline(aes(intercept = 0, slope = 3/2, alpha = .25))    # linetype = "dotted", size = .1

####################
### Cluster Plot ###
####################
#shorten data name
df <- visual.data[, c("TPM", "VAR")]

#Visualize optimal number of clusters
wss <- fviz_nbclust(df, kmeans, method = "wss")
silhouette <- fviz_nbclust(df, kmeans, method = "silhouette")

#Extract optimal number of means
optimal_means_data <- fviz_nbclust(df, kmeans, method = "silhouette")$data
optimal_means <-as.numeric(optimal_means_data$clusters[which.max(optimal_means_data$y)])

#Visualize cluster map using optimal number of means
cluster.data <- kmeans(df, centers = optimal_means, nstart = 25)
cluster.visual <- fviz_cluster(cluster.data, data = df, geom = "point")

####################
### Final Visual ###
####################

grid.arrange(arrangeGrob(base.visual), arrangeGrob(cluster.visual, silhouette, ncol=1, nrow=2), widths=c(2,1))
