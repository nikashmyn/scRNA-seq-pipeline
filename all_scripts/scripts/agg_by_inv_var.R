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

TPM_normed <- readRDS(file=sprintf("%s/aggregated_results/TPM.bygene.normed.rds", dirpath))
#ASE_normed <- readRDS(file=sprintf("%s/aggregated_results/ASE.bygene.normed.rds", dirpath))
#AS_TPM_normed <- readRDS(file=sprintf("%s/aggregated_results/AS-TPM.bygene.normed.rds", dirpath))
ASE <- readRDS(file=sprintf("%s/aggregated_results/ASE.bygene.rds", dirpath))
controlSampleIDs <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))

##########################################################################
### Function to aggregate with inverse variance by genomic coordinates ###
##########################################################################

agg_by_tpm_inv_var <- function(genomic_region = 10000000, mat, anno_cols = 6) {
  
  #read in geneRanges
  geneRanges <- readRDS(sprintf("%s/geneRanges_Nikos.rds", datadir))
  
  #Ensure that chr columns in character and sort mat
  mat <- data.table(mat)
  mat$chr <- as.character(mat$chr) #turn factor chrs to character chrs
  setkey(mat, chr, start, end)
  
  #Get the max length of each chromosome
  chr_lengths <- c()
  for (i in unique(geneRanges$chr)) {
    chr_lengths <- append(chr_lengths, max(geneRanges$end[which(geneRanges$chr == i)]))
  }
  chr_seq_all <- data.table(); ind <- 0;
  for (i in 1:length(unique(geneRanges$chr))){
    num_bins <- floor(chr_lengths[i]/genomic_region)
    chr_seq <- seq(0, (num_bins*genomic_region), by = genomic_region)
    chr_seq[length(chr_seq)] <- chr_lengths[i]
    ind <- ind + length(chr_seq) - 1
    chr_seq_mat <- cbind(chr = as.character(unique(geneRanges$chr)[i]), cbind(start = chr_seq[1:(length(chr_seq)-1)], end = chr_seq[2:length(chr_seq)]))
    chr_seq_all <- rbind(chr_seq_all, chr_seq_mat)
  }
  chr_seq_all$start <- as.numeric(chr_seq_all$start); chr_seq_all$end <- as.numeric(chr_seq_all$end); 
  chr_seq_all <- cbind(bin = c(1:ind), chr_seq_all)
  setkey(chr_seq_all, chr, start, end)
  
  #Prep data table for foverlaps
  changeCols <- c("start", "end") #stores (non-anno) columns that need to be converted to numeric values
  mat[,(changeCols):= lapply(.SD, as.numeric), .SDcols = changeCols] #converts all columns in "changeCols" to numeric for allele A
  setkey(mat, chr, start, end)
  
  #foverlaps mat with genomic coordinate bins
  #mat_bin <- foverlaps(mat, chr_seq_all, mult = "first")
  mat_bin <- foverlaps(chr_seq_all, mat, nomatch=NA) #
  
  #Relabel columns and make bin numeric
  #setnames(mat_bin, old=c("start", "end", "i.start", "i.end"), new=c("bin_start", "bin_end", "gene_start", "gene_end"))
  setnames(mat_bin, old=c("start", "end", "i.start", "i.end"), new=c("gene_start", "gene_end", "bin_start", "bin_end"))
  mat_bin <- cbind(mat_bin[,c("chr", "bin", "bin_start", "bin_end", "gene_id", "gene_start", "gene_end", "width", "strand")], mat_bin[,-c("chr", "bin", "bin_start", "bin_end", "gene_id", "gene_start", "gene_end", "width", "strand")])
  
  #get the inv variance of each row #TODO: only get var over control cells?
  cols <- c(1:anno_cols)
  gene_inv_variance <- 1/rowVars(as.matrix(mat_bin[,..controlSampleIDs]), na.rm = T)
  gene_inv_variance[gene_inv_variance == Inf] <- 0
  gene_inv_variance[is.na(gene_inv_variance)] <- 0
  
  ##plot(variance by avg tpm in control)
  #vars <- unlist(rowVars(as.matrix(TPM[,-c(1:6)]), na.rm = T))
  #means <- unlist(rowMeans(as.matrix(TPM[,-c(1:6)]), na.rm = T))
  #norm_vars <- unlist(rowVars(as.matrix(tpm[,-c(1:6)]), na.rm = T))
  #norm_means <- unlist(rowMeans(as.matrix(tpm[,-c(1:6)]), na.rm = T))
  #plot(log2(means), log2(vars))
  #plot(log2(means), log2(norm_vars))
  #plot(log2(norm_means), log2(norm_vars))
  
  #merge weights with mat for aggregation
  mat_w <- data.table(cbind(weight = gene_inv_variance, mat_bin))

  #aggregate with mean weighted by inv variance
  cols <- c(1:(anno_cols+4))
  sum_cols <- colnames(mat_w[,-..cols])
  mat_w_binned <- data.table(mat_w %>% 
    group_by(bin) %>% 
    summarise_at(sum_cols, funs(weighted.mean(., weight, na.rm = T))))
  
  #merge mat with annos by bin
  setkey(chr_seq_all, bin)
  setkey(mat_w_binned, bin)
  mat_binned2 <- merge(chr_seq_all, mat_w_binned)
  mid <- (mat_binned2$start + mat_binned2$end) / 2
  mat_binned3 <- cbind(mat_binned2[,c(1:4)], cbind(mid, mat_binned2[,-c(1:4)]))
  
  #make sure all value columns are numeric
  changeCols <- colnames(mat_binned2)[-c(1:5)]
  mat_binned3[,(changeCols):= lapply(.SD, as.numeric), .SDcols = changeCols] #converts all columns in "changeCols" to numeric for allele A
  
  return(mat_binned3)
}

############################################
### aggregate ASE by genomic coordinates ###
############################################

#use function to bin each 
ASE_bybin <- list()
#ASE_bybin$A <- agg_by_tpm_inv_var(tpm = ASE_normed$A, mat = ASE_normed$A)
#ASE_bybin$B <- agg_by_tpm_inv_var(tpm = ASE_normed$B, mat = ASE_normed$B)
ASE_bybin$AF <- agg_by_tpm_inv_var(mat = ASE$AF)
#ASE_bybin$TC <- agg_by_tpm_inv_var(tpm = ASE_normed$TC, mat = ASE_normed$TC)

saveRDS(ASE_bybin, file=sprintf("%s/aggregated_results/ASE.inv_var.bybin.rds", dirpath))

############################################
### aggregate TPM by genomic coordinates ###
############################################

#use function to bin each 
TPM_bybin <- agg_by_tpm_inv_var(mat = TPM_normed)

saveRDS(TPM_bybin, file=sprintf("%s/aggregated_results/TPM.inv_var.bybin.rds", dirpath))

###############################################
### aggregate AS-TPM by genomic coordinates ###
###############################################

#AS_TPM_bybin <- list()
#AS_TPM_bybin$A <- agg_by_tpm_inv_var(tpm = AS_TPM_normed$A, mat = AS_TPM_normed$A)
#AS_TPM_bybin$B <- agg_by_tpm_inv_var(tpm = AS_TPM_normed$B, mat = AS_TPM_normed$B)
#
#saveRDS(AS_TPM_bybin, file=sprintf("%s/aggregated_results/AS-TPM.inv_var.bygene.rds", dirpath))

############################
#### Calculate AS-TPM v2 ###
############################

AS_TPM_bybin <- list()
AS_TPM_bybin$A <- cbind(ASE_bybin$AF[,c(1:5)], (TPM_bybin[,-c(1:5)] * (2*ASE_bybin$AF[,-c(1:5)]) ))
AS_TPM_bybin$B <- cbind(ASE_bybin$AF[,c(1:5)], (TPM_bybin[,-c(1:5)] * (2*(1-ASE_bybin$AF[,-c(1:5)])) ))

saveRDS(AS_TPM_bybin, file=sprintf("%s/aggregated_results/AS-TPM.inv_var.bybin.rds", dirpath))

###################################################################
### Function to aggregate with inverse variance over chromosome ###
###################################################################

agg_inv_var_byarm_or_chr <- function(mat, size = "chr", anno_cols = 6) {
  
  #Ensure that chr columns in character and sort mat
  mat <- data.table(mat)
  mat$chr <- as.character(mat$chr) #turn factor chrs to character chrs
  setkey(mat, chr, start, end)
  
  #If size is arm read in armRanges and merge anno
  if(size == "arm") {
    
    message("aggregating matrix over chr arms based on variance")
    
    #read in arm annotations
    armRanges <- data.table(readRDS(file = sprintf("%s/armRanges.rds", datadir)))
    setnames(armRanges, old = c("seqnames"), new = c("chr"))
    setkey(armRanges, chr, start, end)
    
    #foverlaps chr and arm annotations
    mat_arm <- foverlaps(mat, armRanges)
    cols <- c("chr", "start", "end", "width", "strand")
    mat_arm <- mat_arm[,-..cols]
    setnames(mat_arm, old = c("i.start", "i.end", "i.width", "i.strand"), new = c("start", "end", "width", "strand"))
    
    #change the order of the columns to match orginal form
    mat <- data.table(cbind(mat_arm[,c("gene_id", "arm", "start", "end", "width", "strand")], mat_arm[,-c("gene_id", "arm", "start", "end", "width", "strand")]))
    
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
    mat_byarm <- data.table(mat_w %>% 
                              group_by(arm) %>% 
                              summarise_at(sum_cols, funs(weighted.mean(., weight, na.rm = T))))
    
    #make sure all value columns are numeric
    changeCols <- colnames(mat_byarm)[-c(1)]
    mat_byarm[,(changeCols):= lapply(.SD, as.numeric), .SDcols = changeCols] #converts all columns in "changeCols" to numeric for allele A
    
    #return agg by chr matrix
    return(mat_byarm)
  }
  
  if(size == "chr") {
    
    message("aggregating matrix over chrs based on variance")
    
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
    mat_bychr <- data.table(mat_w %>% 
                              group_by(chr) %>% 
                              summarise_at(sum_cols, funs(weighted.mean(., weight, na.rm = T))))
    
    #make sure all value columns are numeric
    changeCols <- colnames(mat_bychr)[-c(1)]
    mat_bychr[,(changeCols):= lapply(.SD, as.numeric), .SDcols = changeCols] #converts all columns in "changeCols" to numeric for allele A
    
    #return agg by chr matrix
    return(mat_bychr)
  }

}

############################
### aggregate ASE by chr ###
############################

#use function to bin each 
ASE_bychr <- list()
#ASE_bychr$A <- agg_inv_var_byarm_or_chr(mat = ASE_normed$A, size = "chr")
#ASE_bychr$B <- agg_inv_var_byarm_or_chr(mat = ASE_normed$B, size = "chr")
ASE_bychr$AF <- agg_inv_var_byarm_or_chr(mat = ASE$AF, size = "chr")
#ASE_bychr$TC <- agg_inv_var_byarm_or_chr(mat = ASE_normed$TC, size = "chr")
saveRDS(ASE_bychr, file=sprintf("%s/aggregated_results/ASE.inv_var.bychr.rds", dirpath))

##TODO:BY arm may need to adjust foverlaps by switching which data is x and y so that there is a bin for each arm no matter what
##use function to bin each 
#ASE_byarm <- list()
##ASE_byarm$A <- agg_inv_var_byarm_or_chr(mat = ASE_normed$A, size = "arm")
##ASE_byarm$B <- agg_inv_var_byarm_or_chr(mat = ASE_normed$B, size = "arm")
#ASE_byarm$AF <- agg_inv_var_byarm_or_chr(mat = ASE_normed$AF, size = "arm")
##ASE_byarm$TC <- agg_inv_var_byarm_or_chr(mat = ASE_normed$TC, size = "arm")
#saveRDS(ASE_byarm, file=sprintf("%s/aggregated_results/ASE.inv_var.byarm.rds", dirpath))

############################
### aggregate TPM by chr ###
############################

#use function to bin each 
TPM_bychr <- agg_inv_var_byarm_or_chr(mat = TPM_normed, size = "chr")
saveRDS(TPM_bychr, file=sprintf("%s/aggregated_results/TPM.inv_var.bychr.rds", dirpath))

#factors <- global_factors(TPM_bychr, anno_cols = 1)
#
##use function to bin each 
#TPM_byarm <- agg_inv_var_byarm_or_chr(mat = TPM_normed, size = "arm")
#saveRDS(TPM_byarm, file=sprintf("%s/aggregated_results/TPM.inv_var.byarm.rds", dirpath))

###############################
### aggregate AS-TPM by chr ###
###############################

#AS_TPM_bychr <- list()
#AS_TPM_bychr$A <- agg_inv_var_byarm_or_chr(mat = AS_TPM_normed$A, size = "chr")
#AS_TPM_bychr$B <- agg_inv_var_byarm_or_chr(mat = AS_TPM_normed$B, size = "chr")
#saveRDS(AS_TPM_bychr, file=sprintf("%s/aggregated_results/AS-TPM.inv_var.bychr.rds", dirpath))
#
#AS_TPM_byarm <- list()
#AS_TPM_byarm$A <- agg_inv_var_byarm_or_chr(mat = AS_TPM_normed$A, size = "arm")
#AS_TPM_byarm$B <- agg_inv_var_byarm_or_chr(mat = AS_TPM_normed$B, size = "arm")
#saveRDS(AS_TPM_byarm, file=sprintf("%s/aggregated_results/AS-TPM.inv_var.byarm.rds", dirpath))


############################
#### Calculate AS-TPM v2 ###
############################

AS_TPM_bychr <- list()
AS_TPM_bychr$A <- cbind(ASE_bychr$AF[,c(1:5)], (TPM_bychr[,-c(1:5)] * (2*ASE_bychr$AF[,-c(1:5)]) ))
AS_TPM_bychr$B <- cbind(ASE_bychr$AF[,c(1:5)], (TPM_bychr[,-c(1:5)] * (2*(1-ASE_bychr$AF[,-c(1:5)])) ))

saveRDS(AS_TPM_bychr, file=sprintf("%s/aggregated_results/AS-TPM.inv_var.bychr.rds", dirpath))


