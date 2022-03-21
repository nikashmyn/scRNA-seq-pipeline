#########################################
#### Input Desired script Directories ###
#########################################
#
#args <- commandArgs(trailingOnly = TRUE)
##args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/rerun_6_9_2021/data", "/pellmanlab/nikos/Stam_Etai_Data", "SIS1025a",  "SIS1025b", "SIS1025c", "SIS1025d", "SIS1025e", "SIS1025f_Lane1", "SIS1025f_Lane2", "SIS1025g_Lane1", "SIS1025g_Lane2", "SIS1025misc", "SIS1025targ")
#print(args)
#scriptsdir <- args[1]
#wkpath <- dirpath <- args[2]
#datadir <- args[3]
#experiments <- args[4:length(args)]

###############
### Imports ###
###############

#Annotation list with information on every cell and IDs of control samples
anno <- data.table(read.csv( sprintf("%s/work_in_progress/annotation_list.csv", datadir)))
controlSampleIDs <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))

#QC Data Frame and high QC IDs based on 5 read cnt threshold for genes
#all_QC <- readRDS(file = sprintf("%s/aggregated_results/all_QC.rds", dirpath))
#high_qc_ids <- as.character(all_QC[which(all_QC$th5 >= quantile(all_QC$th5, c(.10))),]$id) #Excludes bottom 10%

#import samples that were manually determined to be aneuploid
hand_samples <- readRDS(sprintf("%s/aggregated_results/grouped_control_aneuploidies.rds", dirpath))

#import gene and arm annotations
geneRanges <- readRDS(sprintf("%s/geneRanges_Nikos.rds", datadir))
armRanges <- readRDS(sprintf("%s/armRanges.rds", datadir))

#Read in raw TPM matrix
TPM <- readRDS(file=sprintf("%s/aggregated_results/TPM.bygene.rds", dirpath))

#use function to bin each 
TPM_byarm <- readRDS(file=sprintf("%s/aggregated_results/TPM.inv_var.byarm.rds", dirpath))

#use function to bin each 
TPM_bychr <- readRDS(file=sprintf("%s/aggregated_results/TPM.inv_var.bychr.rds", dirpath))

#################
### Functions ###
#################

#function to get pvals through z-test
pval_ztest <- function(x, pop_distr){
  #x is a vector of points to check against population statistics
  #pop_distr is a reference vector from which to pull statistics
  mu <- mean(pop_distr)
  sd <- sd(pop_distr)
  #calc z-score for all values in x against stats from pop_distr
  z_score <- ((x - mu) / sd)
  #calculate a p-value for each element in the vector of z-scores based on normal distr
  pvals <- 2*pnorm( -abs(z_score) ) #*2 for two tailed test
  #return a p-value vector
  return(pvals)
}

#function to get pvals for a matrix 
run_ztest_for_matrix <- function(mat, pop_distr, anno_cols = 1){
  #run pval_ztest function on non-annotation columns of input matrix
  changeCols <- colnames(mat)[-c(1:anno_cols)]
  pval_mat <- copy(mat)
  pval_mat[,(changeCols):= lapply(.SD, pval_ztest, pop_distr = pop_distr), .SDcols = changeCols] 
  #return pval matrix 
  return(pval_mat)
}

get_vals_from_mat <- function(samples_df, mat) {
  vals <- c()
  for(i in 1:nrow(samples_df)) {
    cell <- samples_df$ID[i]
    chr <- paste0("chr", samples_df$chr[i]); if(chr == "chr23") {chr <- "chrX"};
    row <- which(mat$chr == chr)
    vals <- append(vals, unlist(mat[row, ..cell]))
  }
  samples_vals <- cbind(samples_df, vals)
  return(samples_vals)
}

#Get the factors used to do global normalization 
global_factors <- function(mat, anno_cols = 6) {
  #Copy dataframe and remove anno columns
  cols <- c(1:anno_cols)
  dt2go <- copy(data.table(mat)[,-..cols])
  log2go <- log2(dt2go)
  log2go[log2go == -Inf] <- NA
  
  #<(log(TPM_i)>_genes - <<log(TPM_i)>_ctrl>_genes
  cell_means <- colMeans(log2go, na.rm = T)
  ctrl_means <- mean(colMeans(log2go[,..controlSampleIDs], na.rm = T), na.rm = T)
  cell_TPM_diff <- 2^(cell_means - ctrl_means)
  
  #return normed ratios
  return(cell_TPM_diff)
}

##################################################
### Create ref distributions for each CN state ###
##################################################

#split ref samples by CN state
mono_samples <- hand_samples[which(hand_samples$CN == 1),]
ctrl_samples <- data.table(cbind(ID = rep(controlSampleIDs, each=24), CN = 2, chr = rep(c(1:9, "10a", "10b", 11:23), length(controlSampleIDs))))
tri_samples <- hand_samples[which(hand_samples$CN == 3),]

#run get_vals_from_mat function on for each ref group
mono_TPM <- get_vals_from_mat(samples_df=mono_samples, mat=TPM_bychr)
ctrl_TPM <- get_vals_from_mat(samples_df=ctrl_samples, mat=TPM_bychr)
tri_TPM <- get_vals_from_mat(samples_df=tri_samples, mat=TPM_bychr)
ref_TPM <- rbind(mono_TPM, ctrl_TPM, tri_TPM)

#save ref TPM values
saveRDS(ref_TPM, file = sprintf("%s/aggregated_results/ref_TPM.rds", dirpath))

#get global norm factor as QC
TPM_global_factors <- global_factors(TPM)
saveRDS(TPM_global_factors, file = sprintf("%s/aggregated_results/TPM_global_norm_factor.rds", dirpath))

#Pull norm factors for ref samples
mono_factors <- data.table(cbind(ID = names(TPM_global_factors[which(names(TPM_global_factors) %in% mono_TPM$ID)]), factor = TPM_global_factors[which(names(TPM_global_factors) %in% mono_TPM$ID)]))
ctrl_factors <- data.table(cbind(ID = names(TPM_global_factors[which(names(TPM_global_factors) %in% ctrl_TPM$ID)]), factor = TPM_global_factors[which(names(TPM_global_factors) %in% ctrl_TPM$ID)]))
tri_factors <- data.table(cbind(ID = names(TPM_global_factors[which(names(TPM_global_factors) %in% tri_TPM$ID)]), factor = TPM_global_factors[which(names(TPM_global_factors) %in% tri_TPM$ID)]))

#add norm factors to ref data tables
setkey(mono_TPM, ID); setkey(mono_factors, ID); mono_TPM <- merge(mono_TPM, mono_factors);
setkey(ctrl_TPM, ID); setkey(ctrl_factors, ID); ctrl_TPM <- merge(ctrl_TPM, ctrl_factors);
setkey(tri_TPM, ID); setkey(tri_factors, ID); tri_TPM <- merge(tri_TPM, tri_factors);

#run pval calculations on ref tables
TPM_bychr_mono_pvals <- run_ztest_for_matrix(mat=TPM_bychr, pop_distr = mono_TPM$vals, anno_cols = 1)
TPM_bychr_ctrl_pvals <- run_ztest_for_matrix(mat=TPM_bychr, pop_distr = ctrl_TPM$vals, anno_cols = 1)
TPM_bychr_tri_pvals <- run_ztest_for_matrix(mat=TPM_bychr, pop_distr = tri_TPM$vals, anno_cols = 1)

#save pvals for each state
saveRDS(TPM_bychr_mono_pvals, file = sprintf("%s/aggregated_results/pval_matrix_loss_bychr.rds", dirpath))
saveRDS(TPM_bychr_ctrl_pvals, file = sprintf("%s/aggregated_results/pval_matrix_control_bychr.rds", dirpath))
saveRDS(TPM_bychr_tri_pvals, file = sprintf("%s/aggregated_results/pval_matrix_gain_bychr.rds", dirpath))



