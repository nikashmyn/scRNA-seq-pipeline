
#This is the best version as of now. 
#I changed the variant binning so that it calculates the mean number of hits over the SNPs in a gene. 
#It then aggregates the means for every gene into a value by chr. 

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
ganno <- data.table(readRDS(file = sprintf("%s/aggregated_results/ganno.rds", dirpath))) #called arms in visualization script

#Data Frames with TPM and Variant information
#adt <- readRDS(sprintf("%s/aggregated_results/adt.rds", dirpath)) #Normalized TPM matrix
coding <- readRDS(sprintf("%s/aggregated_results/ASE.coding.rds", dirpath)) #Raw Variant matrix in the coding and UTR regions

#import raw gene cnts and gene ranges
geneRanges <- readRDS(sprintf("%s/geneRanges.rds", datadir))
rsemtpm <- readRDS(sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", dirpath))

#list of preset parameters for the experiment
configs <- readRDS(sprintf("%s/param_config_list.rds", datadir))

#new adt object without gene normalization
adt <- fromRawTPMtoExprsRatio(rsemtpm, geneRanges, controlSampleIDs, max_pos_col = 6, plusOne = 1, 
                              zerosAsNA = F, normBySd = F, doNormalizeByControls = F,
                              maxExp = 1200, quantileOfExp = 0.99, minDetectionLevel = configs$minDetectionLevel, minNumOfSamplesToDetect = configs$minNumOfSamplesToDetect )$adt

#I need to open coding[bins.A] and take the mean instead of the sum by gene 
#After taking the mean by gene sum by bin as usual and feed as varA_mat
varA_mat.bySNP <- data.table(copy(coding$A))
setnames(varA_mat.bySNP, old = "id", new = "snp_id")
varB_mat.bySNP <- data.table(copy(coding$B))
setnames(varB_mat.bySNP, old = "id", new = "snp_id")

A <- foverlaps(x = ganno[, c("id", "seqnames", "start", "end", "arm")], y = varA_mat.bySNP)
setcolorder(A, neworder = c("id", "seqnames", "arm", "i.start", "i.end", "snp_id", "start", "end"))
setnames(A, old = c("i.start", "i.end", "start", "end"), new = c("start", "end", "snp_start", "snp_end"))
B <- foverlaps(x = ganno[, c("id", "seqnames", "start", "end", "arm")], y = varB_mat.bySNP)
setcolorder(B, neworder = c("id", "seqnames", "arm", "i.start", "i.end", "snp_id", "start", "end"))
setnames(B, old = c("i.start", "i.end", "start", "end"), new = c("start", "end", "snp_start", "snp_end"))

#group vars by gene 
varA_mat.bygene <- A[,-c(2:8)] %>%
  group_by(id) %>%
  summarise_all(mean, na.rm = TRUE)
varA_mat <- merge(ganno, varA_mat.bygene, by = "id")

#group vars by gene 
varB_mat.bygene <- B[,-c(2:8)] %>%
  group_by(id) %>%
  summarise_all(mean, na.rm = TRUE)
varB_mat <- merge(ganno, varB_mat.bygene, by = "id")

#Var cnts per 50 genomically sequential genes for allele A and B
var_mat_A <- copy(varA_mat); var_mat_B <- copy(varB_mat);
var_mat <- cbind( var_mat_A[,c(1:9)], var_mat_A[,-c(1:9)] + var_mat_B[,-c(1:9)] )

### Exclude noisy chromosomes a priori ###
#exclude_chrs <- c("chr13", "chr18", "chr21", "chr22", "chrX")
exclude_chrs <- c()
if (length(exclude_chrs) > 0) { adt <- adt[-which(adt$seqnames %in% exclude_chrs),] }

if (length(exclude_chrs) > 0) { var_mat <- var_mat[-which(var_mat$seqnames %in% exclude_chrs),] }
var_mat$seqnames <- sub(pattern = "chr", replacement = "", x = var_mat$seqnames)
var_mat$seqnames <- as.integer(sub(pattern = "X", replacement = "23", x = var_mat$seqnames))
var_mat <- var_mat[order(seqnames)]

if (length(exclude_chrs) > 0) { varA_mat <- varA_mat[-which(varA_mat$seqnames %in% exclude_chrs),] }
varA_mat$seqnames <- sub(pattern = "chr", replacement = "", x = varA_mat$seqnames)
varA_mat$seqnames <- as.integer(sub(pattern = "X", replacement = "23", x = varA_mat$seqnames))
varA_mat <- varA_mat[order(seqnames)]

if (length(exclude_chrs) > 0) { varB_mat <- varB_mat[-which(varB_mat$seqnames %in% exclude_chrs),] }
varB_mat$seqnames <- sub(pattern = "chr", replacement = "", x = varB_mat$seqnames)
varB_mat$seqnames <- as.integer(sub(pattern = "X", replacement = "23", x = varB_mat$seqnames))
varB_mat <- varB_mat[order(seqnames)]

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
var_mat.bychr <- var_mat[,-c(1,3:9)] %>%
  group_by(seqnames) %>%
  summarise_all(sum, na.rm = TRUE)
var_mat.bychr <- data.table(var_mat.bychr)
var_mat.bychr <- var_mat.bychr[order(seqnames)]

#group varA by chr
varA_mat.bychr <- varA_mat[,-c(1,3:9)] %>%
  group_by(seqnames) %>%
  summarise_all(sum, na.rm = TRUE)
varA_mat.bychr <- data.table(varA_mat.bychr)
varA_mat.bychr <- varA_mat.bychr[order(seqnames)]

#group varB by chr
varB_mat.bychr <- varB_mat[,-c(1,3:9)] %>%
  group_by(seqnames) %>%
  summarise_all(sum, na.rm = TRUE)
varB_mat.bychr <- data.table(varB_mat.bychr)
varB_mat.bychr <- varB_mat.bychr[order(seqnames)]

#normalize adt columns
adt_rowmeans <- rowMeans(adt.bychr[,..controlSampleIDs], na.rm = T)
adt_rowmeans_mat_reciprocal <- data.table(replicate(ncol(adt.bychr[,-c(1)]), adt_rowmeans))
adt.bychr <- cbind(adt.bychr[,c(1)] , (adt.bychr[,-c(1)] - adt_rowmeans_mat_reciprocal))

#normalize adt columns to old control samples
#adt_rowmeans <- rowMeans(adt.bychr[,..controlSampleIDs2], na.rm = T)
#adt_rowmeans_mat_reciprocal <- data.table(replicate(ncol(adt.bychr[,-c(1)]), adt_rowmeans))
#adt.bychr <- cbind(adt.bychr[,c(1)] , (adt.bychr[,-c(1)] - adt_rowmeans_mat_reciprocal))

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

#Normalize varA SNP counts to get rid of cell wide effects. 
dtgo <- copy(varA_mat.bychr[, -c(1)] + .0001)
#log before mean:
dt2go <- log2(dtgo + 1)
#normalize rows
varA_rowmeans <- rowMeans(dt2go[,..controlSampleIDs], na.rm = T)
varA_rowmeans_mat_reciprocal <- data.table((replicate(ncol(dt2go), varA_rowmeans)))
dt3go <- cbind(varA_mat.bychr[,c(1)], (dt2go - varA_rowmeans_mat_reciprocal))
#normalize columns
dt4go <- cbind(varA_mat.bychr[,c(1)], 2^dt3go[,-c(1)]) #add.0001 so that when both alleles are zero you dont get NA.  
varA_colmeans <- colMeans(dt4go[,-c(1)], na.rm = T)
dtgo_biased <- data.table(1/t(replicate(nrow(dt4go), varA_colmeans)))
dtgo_biased[dtgo_biased == Inf] <- 0 
dt5go_A <- cbind(varA_mat.bychr[,c(1)] , (dtgo * dtgo_biased))

#Normalize varB SNP counts to get rid of cell wide effects. 
dtgo <- copy(varB_mat.bychr[, -c(1)] + .0001)
#log before mean:
dt2go <- log2(dtgo + 1)
#normalize rows
varB_rowmeans <- rowMeans(dt2go[,..controlSampleIDs], na.rm = T)
varB_rowmeans_mat_reciprocal <- data.table((replicate(ncol(dt2go), varB_rowmeans)))
dt3go <- cbind(varB_mat.bychr[,c(1)], (dt2go - varB_rowmeans_mat_reciprocal))
#normalize columns
dt4go <- cbind(varB_mat.bychr[,c(1)], 2^dt3go[,-c(1)])
varB_colmeans <- colMeans(dt4go[,-c(1)], na.rm = T)
dtgo_biased <- data.table(1/t(replicate(nrow(dt4go), varB_colmeans)))
dtgo_biased[dtgo_biased == Inf] <- 0 
dt5go_B <- cbind(varB_mat.bychr[,c(1)] , (dtgo * dtgo_biased))

#Normalize varA SNP ratio dataframe
dt2go <- copy(varA_mat.bychr[, -c(1)])
#log before mean:
dt2go <- log2(dt2go + 1)
#normalize rows
varA_rowmeans <- rowMeans(dt2go[,..controlSampleIDs], na.rm = T)
varA_rowmeans_mat_reciprocal <- data.table((replicate(ncol(dt2go), varA_rowmeans)))
varA_mat.bychr <- cbind(varA_mat.bychr[,c(1)], (dt2go - varA_rowmeans_mat_reciprocal))
#normalize columns
varA_mat.bychr <- cbind(varA_mat.bychr[,c(1)], 2^varA_mat.bychr[,-c(1)])
varA_colmeans <- colMeans(varA_mat.bychr[,-c(1)], na.rm = T)
varA_colmeans_mat_reciprocal <- data.table(1/(t(replicate(nrow(varA_mat.bychr), varA_colmeans))))
varA_colmeans_mat_reciprocal[varA_colmeans_mat_reciprocal == Inf] <- 0 
varA_mat.bychr <- cbind(varA_mat.bychr[,c(1)] , (varA_mat.bychr[,-c(1)] * varA_colmeans_mat_reciprocal))

#Normalize varB SNP ratio dataframe 
dt2go <- copy(varB_mat.bychr[, -c(1)])
#log before mean:
dt2go <- log2(dt2go + 1)
#normalize rows
varB_rowmeans <- rowMeans(dt2go[,..controlSampleIDs], na.rm = T)
varB_rowmeans_mat_reciprocal <- data.table((replicate(ncol(dt2go), varB_rowmeans)))
varB_mat.bychr <- cbind(varB_mat.bychr[,c(1)], (dt2go - varB_rowmeans_mat_reciprocal))
#normalize columns
varB_mat.bychr <- cbind(varB_mat.bychr[,c(1)], 2^varB_mat.bychr[,-c(1)])
varB_colmeans <- colMeans(varB_mat.bychr[,-c(1)], na.rm = T)
varB_colmeans_mat_reciprocal <- data.table(1/(t(replicate(nrow(varB_mat.bychr), varB_colmeans))))
varB_colmeans_mat_reciprocal[varB_colmeans_mat_reciprocal == Inf] <- 0 
varB_mat.bychr <- cbind(varA_mat.bychr[,c(1)] , (varB_mat.bychr[,-c(1)] * varB_colmeans_mat_reciprocal))

#check cols are the same
colIDs <- intersect(colnames(var_mat.bychr), colnames(adt.bychr))

#Check that binning is the same
adt.bychr <- adt.bychr[,..colIDs]; var_mat.bychr <- var_mat.bychr[,..colIDs]; varA_mat.bychr <- varA_mat.bychr[,..colIDs]; varB_mat.bychr <- varB_mat.bychr[,..colIDs];
stopifnot(dim(var_mat.bychr)[2] == dim(adt.bychr)[2])

###################
### Add Alleles ###
###################

#vector difference between alleles cnts per chr
allele_diff_mat.bychr <- cbind(varA_mat.bychr[,c(1)], varA_mat.bychr[,-c(1)] - varB_mat.bychr[,-c(1)])
#absolute difference between allele cnts per chr
abs_allele_diff_mat.bychr <- cbind(varA_mat.bychr[,c(1)], abs(varA_mat.bychr[,-c(1)] - varB_mat.bychr[,-c(1)]))
#total SNP counts over both alleles taken after separate allele normalization. ~same as var_mat.bychr. 
var_mat_sep.bychr <- cbind(varA_mat.bychr[,c(1)], (varA_mat.bychr[,-c(1)] + varB_mat.bychr[,-c(1)])/2 )
#Allelic Fraction by dividing each allele by the sum of the two alleles.
allele_A <- data.table(cbind(dt5go_A[,c(1)], (dt5go_A[,-c(1)])/(dt5go_A[,-c(1)] + dt5go_B[,-c(1)]) ))
allele_B <- data.table(cbind(dt5go_B[,c(1)], (dt5go_B[,-c(1)])/(dt5go_A[,-c(1)] + dt5go_B[,-c(1)]) ))
allele_frac_mat.bychr <- list(A = allele_A, B = allele_B)

####################################
### Save normalized data objects ###
####################################

#TPM object normalized with chr bins
saveRDS(adt.bychr, file=sprintf("%s/aggregated_results/adt.bychr.rds", dirpath))
#vector difference between alleles cnts per chr
saveRDS(allele_diff_mat.bychr, file=sprintf("%s/aggregated_results/allele_diff_mat.bychr.rds", dirpath))
#absolute difference between allele cnts per chr
saveRDS(abs_allele_diff_mat.bychr, file=sprintf("%s/aggregated_results/abs_allele_diff_mat.bychr.rds", dirpath))
#total SNP counts over both alleles normalized together after summing cnts. 
saveRDS(var_mat.bychr, file=sprintf("%s/aggregated_results/var_mat.bychr.rds", dirpath))
#total SNP counts over both alleles taken after separate allele normalization. ~same as var_mat.bychr. 
saveRDS(var_mat_sep.bychr, file=sprintf("%s/aggregated_results/var_mat_sep.bychr.rds", dirpath))
#Allele fractions normalized on the level.
saveRDS(allele_frac_mat.bychr, file=sprintf("%s/aggregated_results/allele_frac_mat.bychr.rds", dirpath))



