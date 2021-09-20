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
var_mat <- cbind( var_mat_A[,c(1:2)], var_mat_A[,-c(1:2)] + var_mat_B[,-c(1:2)] )
vardiff_mat <- cbind( var_mat_A[,c(1:2)], abs(var_mat_A[,-c(1:2)] - var_mat_B[,-c(1:2)]) )

#allelediff_mat <- cbind(var_mat_A[,c(1:2)], var_mat_A[,-c(1:2)] - var_mat_B[,-c(1:2)])

### Exclude noisy chromosomes a priori ###
exclude_chrs <- c("chr13", "chr18", "chr21", "chr22", "chrX")
adt <- adt[-which(adt$seqnames %in% exclude_chrs),]

var_mat <- var_mat[-which(var_mat$seqnames %in% exclude_chrs),]
var_mat$seqnames <- sub(pattern = "chr", replacement = "", x = var_mat$seqnames)
var_mat$seqnames <- as.integer(sub(pattern = "X", replacement = "23", x = var_mat$seqnames))
var_mat <- var_mat[order(seqnames, bin)]

vardiff_mat <- vardiff_mat[-which(vardiff_mat$seqnames %in% exclude_chrs),]
vardiff_mat$seqnames <- sub(pattern = "chr", replacement = "", x = vardiff_mat$seqnames)
vardiff_mat$seqnames <- as.integer(sub(pattern = "X", replacement = "23", x = vardiff_mat$seqnames))
vardiff_mat <- vardiff_mat[order(seqnames, bin)]

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
var_mat.bychr <- data.table(var_mat.bychr)
var_mat.bychr <- var_mat.bychr[order(seqnames)]

#group vardiff by chr
vardiff_mat.bychr <- vardiff_mat[,-c(2)] %>%
  group_by(seqnames) %>%
  summarise_all(mean, na.rm = TRUE)
vardiff_mat.bychr <- data.table(vardiff_mat.bychr)
vardiff_mat.bychr <- vardiff_mat.bychr[order(seqnames)]

#normalize adt columns
adt_rowmeans <- rowMeans(adt.bychr[,..controlSampleIDs], na.rm = T)
adt_rowmeans_mat_reciprocal <- data.table(replicate(ncol(adt.bychr[,-c(1)]), adt_rowmeans))
adt.bychr <- cbind(adt.bychr[,c(1)] , (adt.bychr[,-c(1)] - adt_rowmeans_mat_reciprocal))

adt.bychr <- cbind(adt.bychr[,c(1)], 2^adt.bychr[,-c(1)])
adt3_colmeans <- colMeans(adt.bychr[,-c(1)], na.rm = T)
adt3_colmeans_mat_reciprocal <- data.table(1/(t(replicate(nrow(adt.bychr), adt3_colmeans))))
adt3_colmeans_mat_reciprocal[adt3_colmeans_mat_reciprocal == Inf] <- 0 
adt.bychr <- cbind(adt.bychr[,c(1)] , (adt.bychr[,-c(1)] * adt3_colmeans_mat_reciprocal))

#Normalize variant dataframe
dt2go <- copy(var_mat.bychr[, -c(1)])
#log before mean:
dt2go <- log2(dt2go + 1)
#normalize rows
var_rowmeans <- rowMeans(dt2go[,..controlSampleIDs], na.rm = T)
var_rowmeans_mat_reciprocal <- data.table((replicate(ncol(dt2go), var_rowmeans)))
var_mat.bychr <- cbind(var_mat.bychr[,c(1)], (dt2go - var_rowmeans_mat_reciprocal))
#normalize columns
var_mat.bychr <- cbind(var_mat.bychr[,c(1)], 2^var_mat.bychr[,-c(1)])
var_colmeans <- colMeans(var_mat.bychr[,-c(1)], na.rm = T)
var_colmeans_mat_reciprocal <- data.table(1/(t(replicate(nrow(var_mat.bychr), var_colmeans))))
var_colmeans_mat_reciprocal[var_colmeans_mat_reciprocal == Inf] <- 0 
var_mat.bychr <- cbind(var_mat.bychr[,c(1)] , (var_mat.bychr[,-c(1)] * var_colmeans_mat_reciprocal))

#Normalize vardiff dataframe
dt2go <- copy(vardiff_mat.bychr[, -c(1)])
#log before mean:
dt2go <- log2(dt2go + 1)
#normalize rows
vardiff_rowmeans <- rowMeans(dt2go[,..controlSampleIDs], na.rm = T)
vardiff_rowmeans_mat_reciprocal <- data.table((replicate(ncol(dt2go), vardiff_rowmeans)))
vardiff_mat.bychr <- cbind(vardiff_mat.bychr[,c(1)], (dt2go - vardiff_rowmeans_mat_reciprocal))
#normalize columns
vardiff_mat.bychr <- cbind(vardiff_mat.bychr[,c(1)], 2^vardiff_mat.bychr[,-c(1)])
vardiff_colmeans <- colMeans(vardiff_mat.bychr[,-c(1)], na.rm = T)
vardiff_colmeans_mat_reciprocal <- data.table(1/(t(replicate(nrow(vardiff_mat.bychr), vardiff_colmeans))))
vardiff_colmeans_mat_reciprocal[vardiff_colmeans_mat_reciprocal == Inf] <- 0 
vardiff_mat.bychr <- cbind(vardiff_mat.bychr[,c(1)] , (vardiff_mat.bychr[,-c(1)] * vardiff_colmeans_mat_reciprocal))

#check cols are the same
colIDs <- intersect(colnames(var_mat.bychr), colnames(adt.bychr))

#Check that binning is the same
adt.bychr <- adt.bychr[,..colIDs]; var_mat.bychr <- var_mat.bychr[,..colIDs]; vardiff_mat.bychr <- vardiff_mat.bychr[,..colIDs];
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
visual.data$DIFF <- as.numeric(flatten(vardiff_mat.bychr[,..myids]))
visual.data$chr <- as.factor(rep(unique(adt.bychr$seqnames), length = length(visual.data$VAR)))
visual.data$cell <- rep(myids, each=23-length(exclude_chrs))

####################
### Cluster Plot ###
####################
#shorten data name
df <- visual.data[, c("TPM", "VAR", "DIFF")]

#Visualize optimal number of clusters
wss <- fviz_nbclust(df, kmeans, method = "wss")
silhouette <- fviz_nbclust(df, kmeans, method = "silhouette")

#Extract optimal number of means
optimal_means_data <- fviz_nbclust(df, kmeans, method = "silhouette")$data
optimal_means <-as.numeric(optimal_means_data$clusters[which.max(optimal_means_data$y)])

#Visualize cluster map using optimal number of means
#cluster.data <- kmeans(df, centers = optimal_means, nstart = 25)
cluster.data <- kmeans(df, centers = 3, nstart = 25)

####################
### Final Visual ###
####################

#TPM/VAR plane plot
VAR.visual <- ggplot(data = visual.data, mapping = aes(x = TPM, y = VAR, color = chr, shape = cell)) +
  geom_point(mapping = aes(x = TPM, y = VAR, color = chr, shape = cell)) +
  xlim(c(0,2)) + ylim(c(0,2)) +
  labs(x = "TPM Ratio", y = "SNP Count Ratio", title = sprintf("Expression Contour Plot | %s ", myfamily))

#TPM/Diff plane plot
DIFF.visual <- ggplot(data = visual.data, mapping = aes(x = TPM, y = DIFF, color = chr, shape = cell)) +
  geom_point(mapping = aes(x = TPM, y = DIFF, color = chr, shape = cell)) +
  #xlim(c(0,2)) + ylim(c(0,2)) +
  labs(x = "TPM Ratio", y = "Allelic Fraction Ratio", title = sprintf("Expression Contour Plot | %s ", myfamily))

#Cluster visual
cluster.visual <- fviz_cluster(cluster.data, data = df, geom = "point")

grid.arrange(arrangeGrob(VAR.visual), arrangeGrob(cluster.visual, DIFF.visual, ncol=1, nrow=2), widths=c(2,1))






#########################
### Single Chr Visual ###
#########################

#Prepare Visual Objects
chrom = 12
chr.data <- data.table()
ids_to_vis <- which(colnames(adt.bychr) %in% high_qc_ids)
chr.data$TPM <- as.numeric(adt.bychr[c(chrom), ..ids_to_vis])
chr.data$VAR <- as.numeric(var_mat.bychr[c(chrom), ..ids_to_vis])
chr.data$chr <- as.factor(rep(unique(adt.bychr$seqnames[c(chrom)]), length = length(chr.data$VAR)))
chr.data <- data.table(chr.data)

#Second Draft Plots
chr.visual <- ggplot(data = chr.data, mapping = aes(x = TPM, y = VAR, color = chr)) +
  geom_point(mapping = aes(x = TPM, y = VAR, color = chr)) +
  #xlim(c(0,2)) + ylim(c(0,2)) +
  labs(x = "TPM Ratio", y = "SNP Count Ratio", title = sprintf("Expression Contour Plot | chr%s ", chrom))

####################
### Cluster Plot ###
####################
#shorten data name
df <- chr.data[, c("TPM", "VAR")]
if (length(which(df$TPM > 2)) > 0) {df <- df[-which(df$TPM > 2),]}

#Visualize optimal number of clusters
wss <- fviz_nbclust(df, kmeans, method = "wss")
silhouette <- fviz_nbclust(df, kmeans, method = "silhouette")

#Extract optimal number of means
optimal_means_data <- fviz_nbclust(df, kmeans, method = "silhouette")$data
optimal_means <-as.numeric(optimal_means_data$clusters[which.max(optimal_means_data$y)])

#Visualize cluster map using optimal number of means
cluster.data <- kmeans(df, centers = optimal_means, nstart = 25)
cluster.visual <- fviz_cluster(cluster.data, data = df, geom = "point")

grid.arrange(arrangeGrob(chr.visual), arrangeGrob(cluster.visual, silhouette, ncol=1, nrow=2), widths=c(2,1))
