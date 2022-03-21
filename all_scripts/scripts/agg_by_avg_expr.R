##########################
### Script Explanation ###
##########################

#-------------------------------------------------------
#This script was written by Nikos Mynhier 

########################################
### Input Desired script Directories ###
########################################

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/rerun_6_9_2021/data", "/pellmanlab/nikos/Stam_Etai_Data", "SIS1025a",  "SIS1025b", "SIS1025c", "SIS1025d", "SIS1025e", "SIS1025f_Lane1", "SIS1025f_Lane2", "SIS1025g_Lane1", "SIS1025g_Lane2", "SIS1025misc", "SIS1025targ")
print(args)
scriptsdir <- args[1]
wkpath <- dirpath <- args[2]
datadir <- args[3]
experiments <- args[4:length(args)]

#################################
### Read in ASE and TPM files ###
#################################

TPM <- readRDS(file=sprintf("%s/aggregated_results/TPM.bygene.rds", dirpath))
ASE <- readRDS(file=sprintf("%s/aggregated_results/ASE.bygene.rds", dirpath))
AS_TPM <- readRDS(file=sprintf("%s/aggregated_results/AS-TPM.bygene.rds", dirpath))

##########################################################################
### Function to aggregate with inverse variance by genomic coordinates ###
##########################################################################

agg_by_avg_expr <- function(genomic_region = 20000000, tpm, mat, anno_cols = 6) {
  
  #Ensure that chr columns in character and sort mat
  mat <- data.table(mat)
  mat$chr <- as.character(mat$chr) #turn factor chrs to character chrs
  setkey(mat, chr, start, end)
  
  #make tpm data into a data.table
  tpm <- data.table(tpm)
  
  #Get the max length of each chromosome
  chr_lengths <- c()
  for (i in unique(mat$chr)) {
    chr_lengths <- append(chr_lengths, max(mat$end[which(mat$chr == i)]))
  }
  chr_seq_all <- data.table(); ind <- 0;
  for (i in 1:length(unique(mat$chr))){
    num_bins <- ceiling(chr_lengths[i]/genomic_region)
    chr_seq <- seq(0, (num_bins*genomic_region), by = genomic_region)
    ind <- ind + length(chr_seq) - 1
    chr_seq_mat <- cbind(chr = unique(mat$chr)[i], cbind(start = chr_seq[1:(length(chr_seq)-1)], end = chr_seq[2:length(chr_seq)]))
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
  mat_bin <- foverlaps(mat, chr_seq_all, mult = "first")
  
  #Relabel columns and make bin numeric
  setnames(mat_bin, old=c("start", "end", "i.start", "i.end"), new=c("bin_start", "bin_end", "gene_start", "gene_end"))
  
  #get the inv variance of each row
  avg_ctrl_gene_expr <- rowMeans(tpm[,..controlSampleIDs], na.rm = T)
  #avg_ctrl_gene_expr[avg_ctrl_gene_expr <= 1] <- 0
  
  #merge weights with mat for aggregation
  mat_w <- data.table(cbind(weight = avg_ctrl_gene_expr, mat_bin))
  
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
ASE_bybin$A <- agg_by_avg_expr(tpm = TPM, mat = ASE$A)
ASE_bybin$B <- agg_by_avg_expr(tpm = TPM, mat = ASE$B)
ASE_bybin$AF <- agg_by_avg_expr(tpm = TPM, mat = ASE$AF)
ASE_bybin$TC <- agg_by_avg_expr(tpm = TPM, mat = ASE$TC)

saveRDS(ASE_bybin, file=sprintf("%s/aggregated_results/ASE.avg_expr.bybin.rds", dirpath))

############################################
### aggregate TPM by genomic coordinates ###
############################################

#use function to bin each 
TPM_bybin <- agg_by_avg_expr(tpm = TPM, mat = TPM)

saveRDS(TPM_bybin, file=sprintf("%s/aggregated_results/TPM.avg_expr.bybin.rds", dirpath))

###############################################
### aggregate AS-TPM by genomic coordinates ###
###############################################

AS_TPM_bybin <- list()
AS_TPM_bybin$A <- agg_by_avg_expr(tpm = TPM, mat = AS_TPM$A)
AS_TPM_bybin$B <- agg_by_avg_expr(tpm = TPM, mat = AS_TPM$B)

saveRDS(AS_TPM_bybin, file=sprintf("%s/aggregated_results/AS-TPM.avg_expr.bybin.rds", dirpath))

###########################
### Calculate AS-TPM v2 ###
###########################

AS_TPM_bybin_v2 <- list()
AS_TPM_bybin_v2$A <- cbind(ASE_bybin$AF[,c(1:5)], (TPM_bybin[,-c(1:5)] * ASE_bybin$AF[,-c(1:5)]))
AS_TPM_bybin_v2$B <- cbind(ASE_bybin$AF[,c(1:5)], (TPM_bybin[,-c(1:5)] * (1-ASE_bybin$AF[,-c(1:5)])))

saveRDS(AS_TPM_bybin_v2, file=sprintf("%s/aggregated_results/AS-TPM.avg_expr.bybin.v2.rds", dirpath))
