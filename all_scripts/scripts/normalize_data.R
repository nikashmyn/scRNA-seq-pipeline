##########################
### Script Explanation ###
##########################

#-------------------------------------------------------
#This script was written by Nikos Mynhier 

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

#################################
### Read in ASE and TPM files ###
#################################

#Read in data binned by genomic coordinates
TPM_nolim <- readRDS(file=sprintf("%s/aggregated_results/TPM.nolim.rds", dirpath))
TPM <- readRDS(file=sprintf("%s/aggregated_results/TPM.bygene.rds", dirpath))
ASE <- readRDS(file=sprintf("%s/aggregated_results/ASE.bygene.rds", dirpath))
#AS_TPM <- readRDS(file=sprintf("%s/aggregated_results/AS-TPM.bygene.rds", dirpath))

#Control samples
controlSampleIDs <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs.rds", wkpath))

#######################################
### Function to normalize a matrix ####
#######################################

norm_data_table <- function(mat, anno_cols = 6) {
   #Copy dataframe and remove anno columns
   cols <- c(1:anno_cols)
   dt2go <- copy(data.table(mat)[,-..cols])
   
   #Divide by row means to get ratio to normal 
   rowmeans <- rowMeans(dt2go[,..controlSampleIDs], na.rm = T)
   rowmeans_mat_reciprocal <- data.table((replicate(ncol(dt2go), rowmeans)))
   dt3go <- (dt2go / rowmeans_mat_reciprocal)
   dt3go[dt3go == "NaN"] <- NA
   dt3go[dt3go == Inf] <- NA
   
   #add annotations to ratios
   ratios <- cbind(data.table(mat)[,..cols], dt3go)
   
   #return normed ratios
   return(ratios)
}

##################################
### aggregate TPM whole genome ###
##################################

adjust_by_inv_var <- function(mat, anno_cols = 6) {
   
   #read in geneRanges
   geneRanges <- readRDS(sprintf("%s/geneRanges_Nikos.rds", datadir))
   
   #Ensure that chr columns in character and sort mat
   mat <- data.table(mat)
   mat$chr <- as.character(mat$chr) #turn factor chrs to character chrs
   setkey(mat, chr, start, end)
   
   #get the inv variance of each row #TODO: only get var over control cells?
   cols <- c(1:anno_cols)
   gene_inv_variance <- 1/rowVars(as.matrix(mat[,..controlSampleIDs]), na.rm = T)
   gene_inv_variance[gene_inv_variance == Inf] <- 0
   gene_inv_variance[is.na(gene_inv_variance)] <- 0
   
   #merge weights with mat for aggregation
   mat_w <- data.table(cbind(weight = gene_inv_variance, mat))
   
   #aggregate with mean weighted by inv variance
   cols <- c(1:(anno_cols+1))
   sum_cols <- colnames(mat_w[,-..cols])
   mat_w_cell <- data.table(mat_w %>% summarise_at(sum_cols, funs(weighted.mean(., weight, na.rm = T))))
   
   #adjust the input matrix
   cols <- c(1:anno_cols)
   adjustment_mat <- data.table(t((replicate(nrow(mat), unlist(mat_w_cell)))))
   mat_adjusted <- (mat[,-..cols] / adjustment_mat)
   mat_out <- cbind(mat[,..cols], mat_adjusted)
   
   return(mat_out)
}


###############################################################
### "Global" Normalization: Cell-specific effect correction ###
###############################################################

global_normalization <- function(mat, anno_cols = 6) {
   #Copy dataframe and remove anno columns
   cols <- c(1:anno_cols)
   dt2go <- copy(data.table(mat)[,-..cols])
   log2go <- log2(dt2go)
   log2go[log2go == -Inf] <- NA
   
   #<(log(TPM_i)>_genes - <<log(TPM_i)>_ctrl>_genes
   cell_means <- colMeans(log2go, na.rm = T)
   ctrl_means <- mean(colMeans(log2go[,..controlSampleIDs], na.rm = T), na.rm = T)
   cell_TPM_diff <- 2^(cell_means - ctrl_means)
   
   #scale each cell by differences in average TPM and add annotations
   cell_TPM_diff_df <- data.table(t((replicate(nrow(dt2go), cell_TPM_diff))))
   dt3go <- (dt2go / cell_TPM_diff_df)
   mat_out <- cbind(mat[,..cols], dt3go)
   
   #return normed ratios
   return(mat_out)
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


global_adjustment <- function(mat, anno_cols = 6) {
   
   #xi
   cols <- c(1:anno_cols)
   x <- copy(data.table(mat)[,-..cols])
   lnx <- as.matrix(log2(x))
   lnx[lnx == -Inf] <- NA
   
   #Ci
   c <- rowMeans(data.table(mat)[,..controlSampleIDs])
   lnc <- log2(c)
   lnc[lnc == -Inf] <- NA
   
   #linear fit between lnx and lnc
   model_mat <- c()
   for(i in 1:ncol(lnx)){
      model <- lm(lnx[,i] ~ lnc)$coefficients
      model_mat <- rbind(model_mat, model)
   }
   
   #change rownames to cell ids and colnames to slope and intercept
   model_mat <- data.table(cbind(cell = colnames(lnx), model_mat))
   setnames(model_mat, old = c("(Intercept)", "lnc"), new = c("intercept", "slope"))
   
   #change columns to numeric
   changeCols <- c("intercept", "slope")
   model_mat[,(changeCols):= lapply(.SD, as.numeric), .SDcols = changeCols] #converts all columns in "changeCols" to numeric for allele A
   
   out_mat <- c()
   for(i in 1:ncol(lnx)){
      adjusted_cell <- (lnx[,i] - model_mat$intercept[i]) / (model_mat$slope[i])
      out_mat <- cbind(out_mat, adjusted_cell)
   }
   
   colnames(out_mat) <- colnames(lnx)
   out_mat_2 <- 2^out_mat
   out <- data.table(cbind(mat[,..cols], out_mat_2))
   
   #return normed ratios
   return(out)
}

########################################################################
### Split chr10 into two parts by adding additional chr factor level ###
########################################################################

split_chr10_anno <- function(mat){
   
   #change chr column to characters instead of factors
   mat$chr <- as.character(unlist(mat$chr))

   #get genes before and after 61Mb
   mat_10 <- mat[which(mat$chr == "chr10")]
   mat_10a_genes <- mat_10$gene_id[which(mat_10$start < 61000000)]
   mat_10b_genes <- mat_10$gene_id[which(mat_10$start >= 61000000)]

   #change regions before 61Mb to a and after 61Mb to b
   mat$chr[which(mat$gene_id %in% mat_10a_genes)] <- "chr10a"
   mat$chr[which(mat$gene_id %in% mat_10b_genes)] <- "chr10b"
   
   return(mat)
   
}

################################
### Norm all data structures ###
################################

#read in geneRanges to change chr entries
geneRanges <- readRDS(sprintf("%s/geneRanges_Nikos.rds", datadir))

#split chr10 into two chr entries
geneRanges2 <- split_chr10_anno(geneRanges)
saveRDS(geneRanges2, sprintf("%s/geneRanges_Nikos.rds", datadir))

#split chr10 into two chr entries
TPM_nolim2 <- split_chr10_anno(TPM_nolim)
saveRDS(TPM_nolim2, file=sprintf("%s/aggregated_results/TPM.nolim.rds", dirpath))

#split chr10 into two chr entries
ASE2 <- list()
ASE2$A <- split_chr10_anno(ASE$A)
ASE2$B <- split_chr10_anno(ASE$B)
ASE2$AF <- split_chr10_anno(ASE$AF)
ASE2$TC <- split_chr10_anno(ASE$TC)
saveRDS(ASE2, file=sprintf("%s/aggregated_results/ASE.bygene.rds", dirpath))

#split chr10 into two chr entries
TPM2 <- split_chr10_anno(TPM)

#calculate TPM ratios and do global adjustment
TPM_normed <- norm_data_table(TPM2)
TPM_adjusted <- adjust_by_inv_var(TPM_normed)
saveRDS(TPM_adjusted, file=sprintf("%s/aggregated_results/TPM.bygene.normed.rds", dirpath))

## or #
#
##new approach
#TPM_adjusted_2 <- global_adjustment(TPM)
#TPM_normed_2 <- norm_data_table(TPM_adjusted_2)
#saveRDS(TPM_normed_2, file=sprintf("%s/aggregated_results/TPM.bygene.normed.rds", dirpath))


