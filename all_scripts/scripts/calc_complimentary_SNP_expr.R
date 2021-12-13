
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

#Data Frames with TPM and Variant information
#adt <- readRDS(sprintf("%s/aggregated_results/adt.rds", dirpath)) #Normalized TPM matrix
coding <- readRDS(sprintf("%s/aggregated_results/ASE.coding.rds", dirpath)) #Raw Variant matrix in the coding and UTR regions

#import raw gene cnts and gene ranges
geneRanges <- readRDS(sprintf("%s/geneRanges.rds", datadir))
rsemtpm <- readRDS(sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", dirpath))

#list of preset parameters for the experiment
configs <- readRDS(sprintf("%s/param_config_list.rds", datadir))
#configs$minDetectionLevel <- 50

#Var cnts per 50 genomically sequential genes for allele A and B
varA_mat <- var_mat_A <- coding$A
varB_mat <- var_mat_B <- coding$B

### Exclude noisy chromosomes a priori ###
#exclude_chrs <- c("chr13", "chr18", "chr21", "chr22", "chrX")
exclude_chrs <- c()

if (length(exclude_chrs) > 0) { varA_mat <- varA_mat[-which(varA_mat$seqnames %in% exclude_chrs),] }
varA_mat$seqnames <- sub(pattern = "chr", replacement = "", x = varA_mat$seqnames)
varA_mat$seqnames <- as.integer(sub(pattern = "X", replacement = "23", x = varA_mat$seqnames))
varA_mat <- varA_mat[order(seqnames, start, end)]

if (length(exclude_chrs) > 0) { varB_mat <- varB_mat[-which(varB_mat$seqnames %in% exclude_chrs),] }
varB_mat$seqnames <- sub(pattern = "chr", replacement = "", x = varB_mat$seqnames)
varB_mat$seqnames <- as.integer(sub(pattern = "X", replacement = "23", x = varB_mat$seqnames))
varB_mat <- varB_mat[order(seqnames, start, end)]

######################################################
### Get complimentary expression DF from two cells ###
######################################################

Get_SNP_info_dfs <- function(sample) {
  cell <- list(A = varA_mat[,..sample], B = varB_mat[,..sample])
  return(cell)
}

Get_comp_expr2 <- function(cell_annos, cell1, cell2) {
  #frac is minimum difference in expression fraction needed to be classified as complimentary
  cell1 <- cell1 #prevent NaNs
  cell2 <- cell2 #and ensures bins with no expression in either cell are 50/50
  abs_diff_cells_A <- abs(cell1$A - cell2$A)
  abs_sum_cells_A <- abs(cell1$A + cell2$A)
  abs_diff_cells_B <- abs(cell1$B - cell2$B)
  abs_sum_cells_B <- abs(cell1$B + cell2$B)
  abs_diff_cells_A <- cbind(cell_annos, abs_diff_cells_A)
  abs_diff_cells_B <- cbind(cell_annos, abs_diff_cells_B)
  abs_diff_cells_red_A <- abs_diff_cells_A[which(abs_sum_cells_A >= 1),]
  abs_diff_cells_red_B <- abs_diff_cells_B[which(abs_sum_cells_B >= 1),]
  colnames(abs_diff_cells_red_A) <- c("seqnames", "A")
  colnames(abs_diff_cells_red_B) <- c("seqnames", "B")
  abs_diff_cells_red <- list(A = abs_diff_cells_red_A, B = abs_diff_cells_red_B)
  return(abs_diff_cells_red)
}

get_frac_comp <- function (abs_diff_cells_red, chr) {
  abs_diff_cells_red_chr_A <- unlist(abs_diff_cells_red[["A"]][which(abs_diff_cells_red[["A"]]$seqnames == chr),"A"])
  abs_diff_cells_red_chr_B <- unlist(abs_diff_cells_red[["B"]][which(abs_diff_cells_red[["B"]]$seqnames == chr),"B"])
  result <- list(chr = sprintf("chr%s", chr), Frac_A = mean(abs_diff_cells_red_chr_A), Frac_B = mean(abs_diff_cells_red_chr_B))
  return(result)
}

macro_fun <- function(sample1_name, sample2_name, cell_annos, chr) {
  
  #Get information for each sample
  sample1 <- Get_SNP_info_dfs(sample1_name) 
  sample2 <- Get_SNP_info_dfs(sample2_name) 
  
  #use information to get complimentary expresssion 
  comp_bins_binary <- Get_comp_expr2(cell_annos, sample1, sample2)
  
  #Synthesize results
  frac_comp_bins <- data.frame(append(c(sample1_name, sample2_name), get_frac_comp(comp_bins_binary, chr)))
  colnames(frac_comp_bins) <- c("cell1", "cell2", "chr", "A", "B")
  
  #return results
  return(frac_comp_bins)
  
}


sample1_name <- "210503_10F" #210503_10B, ! 210503_10F, ! 210503_9F, ! 210428_8F, 210104_2B #monosomy 210201_6G 210621_9D #trisomy 210621_2E #control 170202_A1
sample2_name <- "210503_10G" #210503_10C, ! 210503_10G, ! 210503_9G, ! 210428_8G, 210104_2C #disomy 210201_6H 210621_9E #disomy 210621_2F #control 170202_A2


samples_targ_exp <- data.frame(c("210428_8B", "210428_8C"), c("210428_8F", "210428_8G"), c("210503_10D", "210503_10E"), c("210503_10F", "210503_10G"), c("210503_9A", "210503_9B"), c("210503_9C", "210503_9D"), c("210503_9F", "210503_9G"), c("210503_9H", "210503_10A"), c("210503_10B", "210503_10C"))
cell_annos <- varA_mat[,c('seqnames')]

comp_info_df <- data.frame()
for (i in 1:ncol(samples_targ_exp)) {
  for (j in c(1:12,14:17,20)) {
    tmp <- macro_fun(samples_targ_exp[1,i], samples_targ_exp[2,i], cell_annos, j)
    comp_info_df <- rbind(comp_info_df, tmp)
  }
}

A_th <- quantile(comp_info_df$A, .95)
B_th <- quantile(comp_info_df$B, .95)

outlier_complimentary_df <- comp_info_df[union(which(comp_info_df$A > A_th), which(comp_info_df$B > B_th)),]


