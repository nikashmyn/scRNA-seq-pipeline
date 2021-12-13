#normalize_tpm_and_vars_bychr4.R

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
varA_mat.bychr <- A[,-c(1,3:8)] %>%
  group_by(seqnames) %>%
  summarise_all(sum, na.rm = TRUE)
varA_mat.bychr <- data.table(varA_mat.bychr)
varA_mat.bychr <- varA_mat.bychr[order(seqnames)]

#group vars by gene 
varB_mat.bychr <- B[,-c(1,3:8)] %>%
  group_by(seqnames) %>%
  summarise_all(sum, na.rm = TRUE)
varB_mat.bychr <- data.table(varB_mat.bychr)
varB_mat.bychr <- varB_mat.bychr[order(seqnames)]

#### Exclude noisy chromosomes a priori ###
#exclude_chrs <- configs[["chr_to_excl"]]
#
#if (length(exclude_chrs) > 0) { varA_mat <- varA_mat[-which(varA_mat$seqnames %in% exclude_chrs),] }
varA_mat.bychr$seqnames <- sub(pattern = "chr", replacement = "", x = varA_mat.bychr$seqnames)
varA_mat.bychr$seqnames <- as.integer(sub(pattern = "X", replacement = "23", x = varA_mat.bychr$seqnames))
varA_mat.bychr <- varA_mat.bychr[order(seqnames)]

#if (length(exclude_chrs) > 0) { varB_mat <- varB_mat[-which(varB_mat$seqnames %in% exclude_chrs),] }
varB_mat.bychr$seqnames <- sub(pattern = "chr", replacement = "", x = varB_mat.bychr$seqnames)
varB_mat.bychr$seqnames <- as.integer(sub(pattern = "X", replacement = "23", x = varB_mat.bychr$seqnames))
varB_mat.bychr <- varB_mat.bychr[order(seqnames)]

############################################
### Re-Engineer raw variant count Matrix ###
############################################

#Normalize varA SNP counts to get rid of cell wide effects. 
#We need a specific treatment for when both alleles are zero  
dt2go_A <- copy(varA_mat.bychr[, -c(1)])

#Normalize varB SNP counts to get rid of cell wide effects. 
dt2go_B <- copy(varB_mat.bychr[, -c(1)])

#Get the allelic fractions
dt3go_A <- dt2go_A / (dt2go_A + dt2go_B) 
dt3go_B <- dt2go_B / (dt2go_A + dt2go_B) 

#add back chr information
allele_cnts_A.bychr <- cbind(varA_mat.bychr[, c(1)], dt3go_A)
allele_cnts_B.bychr <- cbind(varB_mat.bychr[, c(1)], dt3go_B)

#combine the information
allele_cnts.bychr <- list(A = allele_cnts_A.bychr, B = allele_cnts_B.bychr)


####################################
### Save normalized data objects ###
####################################

#Allele cnts normalized on the chr level.
saveRDS(allele_cnts.bychr, file=sprintf("%s/aggregated_results/allele_cnts.bychr.rds", dirpath))


############################################
### Part 2 allelic fraction on bin level ###
############################################

#cints with bins
varA_mat <- data.table(copy(coding$bins.A))
varB_mat <- data.table(copy(coding$bins.B))

#add the chr to the bin so that all bins are unique
varA_mat$bin <- paste(varA_mat$seqnames, varA_mat$bin, sep = "_")
varB_mat$bin <- paste(varB_mat$seqnames, varB_mat$bin, sep = "_")

#group vars by bin 
varA_mat.bybin <- varA_mat[,-c(1,3:7)] %>%
  group_by(bin) %>%
  summarise_all(sum, na.rm = TRUE)

#group vars by bin 
varB_mat.bybin <- varB_mat[,-c(1,3:7)] %>%
  group_by(bin) %>%
  summarise_all(sum, na.rm = TRUE)

#turn cnts by bin into allelic fraction
varA_frac.bybin <- cbind(varA_mat.bybin[,1], varA_mat.bybin[,-1] / (varA_mat.bybin[,-1] + varB_mat.bybin[,-1]))
varB_frac.bybin <- cbind(varB_mat.bybin[,1], varB_mat.bybin[,-1] / (varA_mat.bybin[,-1] + varB_mat.bybin[,-1]))

#Turn NaNs into NAs
varA_frac.bybin[varA_frac.bybin == "NaN"] <- NA
varB_frac.bybin[varB_frac.bybin == "NaN"] <- NA

#split the bin column into two columns
varA_frac.bybin <- data.table(cbind(str_split_fixed(varA_frac.bybin$bin, "_", 2), varA_frac.bybin[,-1]))
varB_frac.bybin <- data.table(cbind(str_split_fixed(varB_frac.bybin$bin, "_", 2), varB_frac.bybin[,-1]))

#adjust the column names
colnames(varA_frac.bybin)[1:2] <- c("seqnames", "bin"); varA_frac.bybin$bin <- as.numeric(varA_frac.bybin$bin)
colnames(varB_frac.bybin)[1:2] <- c("seqnames", "bin"); varB_frac.bybin$bin <- as.numeric(varB_frac.bybin$bin)

#combine into one object for storage
allele_frac.bybin <- list(A = varA_frac.bybin, B = varB_frac.bybin)

#Allele fractions normalized on the bin level.
saveRDS(allele_frac.bybin, file=sprintf("%s/aggregated_results/allele_frac.bybin.rds", dirpath))


