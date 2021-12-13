#normalize_tpm_and_vars_bychr3.R

#This version treats the variants as CZ suggests

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

#list of preset parameters for the experiment
configs <- readRDS(sprintf("%s/param_config_list.rds", datadir))

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
#varA_mat <- copy(var_mat_A); varB_mat <- copy(var_mat_B);


#### Exclude noisy chromosomes a priori ###
#exclude_chrs <- configs[["chr_to_excl"]]
#
#if (length(exclude_chrs) > 0) { varA_mat <- varA_mat[-which(varA_mat$seqnames %in% exclude_chrs),] }
varA_mat$seqnames <- sub(pattern = "chr", replacement = "", x = varA_mat$seqnames)
varA_mat$seqnames <- as.integer(sub(pattern = "X", replacement = "23", x = varA_mat$seqnames))
varA_mat <- varA_mat[order(seqnames)]

#if (length(exclude_chrs) > 0) { varB_mat <- varB_mat[-which(varB_mat$seqnames %in% exclude_chrs),] }
varB_mat$seqnames <- sub(pattern = "chr", replacement = "", x = varB_mat$seqnames)
varB_mat$seqnames <- as.integer(sub(pattern = "X", replacement = "23", x = varB_mat$seqnames))
varB_mat <- varB_mat[order(seqnames)]

############################################
### Re-Engineer raw variant count Matrix ###
############################################

#group varA by chr
varA_mat.bychr <- varA_mat[,-c(1,3:9)] %>%
  group_by(seqnames) %>%
  summarise_all(mean, na.rm = TRUE)
varA_mat.bychr <- data.table(varA_mat.bychr)
varA_mat.bychr <- varA_mat.bychr[order(seqnames)]

#group varB by chr
varB_mat.bychr <- varB_mat[,-c(1,3:9)] %>%
  group_by(seqnames) %>%
  summarise_all(mean, na.rm = TRUE)
varB_mat.bychr <- data.table(varB_mat.bychr)
varB_mat.bychr <- varB_mat.bychr[order(seqnames)]

#Normalize varA SNP counts to get rid of cell wide effects. 
#We need a specific treatment for when both alleles are zero  
dt2go_A <- copy(varA_mat.bychr[, -c(1)])

#Normalize varB SNP counts to get rid of cell wide effects. 
dt2go_B <- copy(varB_mat.bychr[, -c(1)])

#Get the allelic fractions
dt3go_A <- dt2go_A / (dt2go_A + dt2go_B) 
dt3go_B <- dt2go_B / (dt2go_A + dt2go_B) 

#add back chr information
allele_frac_A.bychr <- cbind(varA_mat.bychr[, c(1)], dt3go_A)
allele_frac_B.bychr <- cbind(varB_mat.bychr[, c(1)], dt3go_B)

#combine the information
allele_frac.bychr <- list(A = allele_frac_A.bychr, B = allele_frac_B.bychr)


####################################
### Save normalized data objects ###
####################################

#Allele fractions normalized on the level.
saveRDS(allele_frac.bychr, file=sprintf("%s/aggregated_results/allele_frac.bychr.rds", dirpath))




