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
adt <- readRDS(sprintf("%s/aggregated_results/adt.rds", dirpath)) #Normalized TPM matrix
coding <- readRDS(sprintf("%s/aggregated_results/ASE.coding.rds", dirpath)) #Raw Variant matrix in the coding and UTR regions

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

#creates a column in adt with the bin number. putting genes in groups of 50.
bins <- lapply(unique(adt$seqnames), function(chr) unlist(lapply(1:ceiling(nrow(adt[seqnames==chr])/50), rep, 50 ))[1:nrow(adt[seqnames==chr])])
#bins <- lapply(unique(adt$seqnames), function(chr) unlist(rep(c(1:floor(nrow(adt[seqnames==chr])/50)), times=c(rep(50, floor(nrow(adt[seqnames==chr])/50)-1), 50+(nrow(adt[seqnames==chr]) - floor(nrow(adt[seqnames==chr])/50)*50)) ))[1:nrow(adt[seqnames==chr])])
names(bins) <- unique(adt$seqnames)
bins2 <- rbindlist(lapply(names(bins), function(chr) data.table(seqnames = chr, bin = bins[[chr]])))
stopifnot(sum(bins2$seqnames != adt$seqnames)==0)
adt2 <- cbind(bin = bins2$bin, adt)

#Aggregate values by bin
adt.bin <- adt2[,-c(2,4,5)] %>% 
  group_by(seqnames, add=T) %>%
  group_by(bin, add=T) %>%
  summarise_all(mean, na.rm = TRUE)
adt.bin <- data.table(adt.bin)
adt.bin$seqnames <- sub(pattern = "chr", replacement = "", x = adt.bin$seqnames)
adt.bin$seqnames <- as.integer(sub(pattern = "X", replacement = "23", x = adt.bin$seqnames))
adt.bin <- adt.bin[order(seqnames, bin)]

#check cols are the same
colIDs <- intersect(colnames(var_mat), colnames(adt.bin))
adt.bin <- adt.bin[,..colIDs]; var_mat <- var_mat[,..colIDs];

#check rows are the same
if (dim(adt.bin) != dim(var_mat)) {
  adt.bin <- adt.bin[-which(adt.bin$bin != var_mat$bin)[1],]
}

#Check that binning is the same
stopifnot(adt.bin$bin == var_mat$bin)

#group TPM by chr 
adt.bychr <- adt.bin[,-c(2)] %>%
  group_by(seqnames) %>%
  summarise_all(mean, na.rm = TRUE)
adt.bychr <- data.table(adt.bychr)
adt.bychr <- adt.bychr[order(seqnames)]

#group var cnts by chr
var_mat.bychr <- var_mat[,-c(2)] %>%
  group_by(seqnames) %>%
  summarise_all(mean, na.rm = TRUE)
#var_mat.bychr$seqnames <- as.numeric(c(1:23))
var_mat.bychr <- data.table(var_mat.bychr)
var_mat.bychr <- var_mat.bychr[order(seqnames)]

#normalize adt columns
adt.bychr <- cbind(adt.bychr[,c(1)], 2^adt.bychr[,-c(1)])
adt3_colmeans <- colMeans(adt.bychr[,-c(1)], na.rm = T)
adt3_colmeans_mat_reciprocal <- data.table(1/(t(replicate(nrow(adt.bychr), adt3_colmeans))))
adt3_colmeans_mat_reciprocal[adt3_colmeans_mat_reciprocal == Inf] <- 0 
adt.bychr <- cbind(adt.bychr[,c(1)] , (adt.bychr[,-c(1)] * adt3_colmeans_mat_reciprocal))
#saveRDS(adt.bin, file = sprintf("%s/aggregated_results/adt.bin.rds", dirpath))

#Normalize variant dataframe
var_rowmeans <- rowMeans(var_mat.bychr[,..controlSampleIDs], na.rm = T)
var_rowmeans_mat_reciprocal <- data.table(1/(replicate(ncol(var_mat.bychr[,-c(1)]), var_rowmeans)))
var_rowmeans_mat_reciprocal[var_rowmeans_mat_reciprocal == Inf] <- 0 
var_mat.bychr <- cbind(var_mat.bychr[,c(1)] , (var_mat.bychr[,-c(1)] * var_rowmeans_mat_reciprocal))

var_colmeans <- colMeans(var_mat.bychr[,-c(1)], na.rm = T)
var_colmeans_mat_reciprocal <- data.table(1/(t(replicate(nrow(var_mat.bychr), var_colmeans))))
var_colmeans_mat_reciprocal[var_colmeans_mat_reciprocal == Inf] <- 0 
var_mat.bychr <- cbind(var_mat.bychr[,c(1)] , (var_mat.bychr[,-c(1)] * var_colmeans_mat_reciprocal))

###################
### Contour Map ###
###################

#get families from anno list
setkey(anno, WTA.plate)
families <- sort(table(anno[Pairs != "NA" & Pairs != "mother"]$Pairs), decreasing = T)

#Set seed for reproducability
set.seed(123)

#Prepare Visual Objects
myfamily = names(families)[30] #43, 24, 142, 30
message(myfamily)
myids <- anno[Pairs %in% myfamily]$WTA.plate
visual.data <- data.table(TPM = as.numeric(flatten(adt.bychr[,..myids])))
visual.data$VAR <- as.numeric(flatten(var_mat.bychr[,..myids]))
visual.data$chr <- as.factor(rep(unique(adt.bin$seqnames), length = length(visual.data$VAR)))
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
