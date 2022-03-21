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

################################
### read in exploration data ###
################################

TPM_list <- readRDS(file=sprintf("%s/aggregated_results/TPM_list.rds", dirpath))

plot(log2(rowMeans(TPM_list$raw_tpm[,..controlSampleIDs])), log2(TPM_list$ratio_ctrl_var))

plot(log2(rowVars(as.matrix(TPM_list$ratio_tpm[,-c(1:6)]))), log2(TPM_list$ratio_ctrl_var))

plot(log2(rowMeans(as.matrix(TPM_list$ratio_tpm[,-c(1:6)]))), log2(rowMeans(as.matrix(TPM_list$ratio_tpm[,..controlSampleIDs]))))

plot(log2(unlist(TPM_list$raw_tpm[,c(100)])), log2(TPM_list$raw_ctrl_means)) #, xlim = c(0,10000), ylim = c(0,10000))

#mean AF on y in ctrl and mean tpm in ctrl in x
plot(TPM_list$raw_ctrl_means, rowMeans(ASE$AF[,..controlSampleIDs]))

#plot mean ratio_tpm vs mean tpm
plot(rowMeans(TPM_list$raw_tpm[,-c(1:6)]), rowMeans(TPM_list$ratio_tpm[,-c(1:6)]), xlim = c(0,2000), ylim = c(0,10))
#Conclusions:
#1) this reiterates what we know about TPM ratio, that the ratio is more stable at higher expression levels
#2) This begs the question what is the variance as a function of tpm ratio

#plot variance as function of tpm ratio
plot(rowMeans(TPM_list$ratio_tpm[,-c(1:6)]), TPM_list$ratio_ctrl_var, xlim = c(0,10), ylim = c(0,100))
#Conclusions:
#1) there is a contingent of points that are high tpm ratio and most likely lowly expressed genes with low variance
#2) this would heavily weight these lowly expressed genes which we know to be inaccurate
#3) There are also genes with low tpm ratio and high variance also likely from lowly expressed genes
#3) Perhaps if we re-implement the inv vars with a greater threshold we would eliminate these genes
#   and that prevent us from using the inv vars method.
#4) Inv vars method is not working because the lowly variant genes hace a wide spread of TPM ratios
#5) the trail up and to the side at the edges of the plots are driven by the lowly expressed genes which have higher and lowly TPM ratios
#6) The assumption that lowly expressed genes are highly variant is apparently not true but how can that be? Biological effect? 

#TODO: Would be interesting to explore why and which lowly expressed genes can have high ratio and low variance

#get rid of lowly expressed genes
TPM_list <- readRDS(file=sprintf("%s/aggregated_results/TPM_list.rds", dirpath))
genes <- which(rowMeans(TPM_list$raw_tpm[,..controlSampleIDs]) >= 35)
TPM_list$raw_tpm <- TPM_list$raw_tpm[genes,]
TPM_list$ratio_tpm <- TPM_list$ratio_tpm[genes,]
TPM_list$raw_ctrl_means <- TPM_list$raw_ctrl_means[genes]
TPM_list$ratio_ctrl_var <- TPM_list$ratio_ctrl_var[genes]


#re-plot mean ratio_tpm vs mean tpm with less lowly expressed genes
plot(rowMeans(TPM_list$raw_tpm[,-c(1:6)]), rowMeans(TPM_list$ratio_tpm[,-c(1:6)]), xlim = c(0,2000), ylim = c(0,10))
#re-plot variance as function of tpm ratio with less lowly expressed genes
plot(rowMeans(TPM_list$ratio_tpm[,-c(1:6)]), TPM_list$ratio_ctrl_var, ylim = c(0,20)) # xlim = c(0,5)

#Conclusions:
#1) As you eliminate the lowly expressed genes you do reduce the two problematic point groups. but not eliminate.
#2) There persist a few genes that have high expression & high TPM ratios with very low variance
#   throughout the experimental cells. Are these perhaps genes responding to the experimental treatment?

#Question: what is the relationship between raw expression and variance rather than ratio expression and variance
#plot variance as a function of raw tpm
plot(rowMeans(TPM_list$raw_tpm[,-c(1:6)]), TPM_list$ratio_ctrl_var)
plot(rowMeans(TPM_list$ratio_tpm[,-c(1:6)]), TPM_list$ratio_ctrl_var)
#Conclusions:
#1) There is a relationship between variance and raw expression but not variance and tpm ratio.
#   That means that its better to do inv variance weighting with tpm ratio.


#Question: Which genes have TPM ratio greater than 2?
test <- TPM_list$raw_tpm[which(rowMeans(TPM_list$ratio_tpm[,-c(1:6)]) > 2),]



#OVERALL CONCLUSIONS
#1) The reason the inv variance weighing wasnt working is because there were high and low tpm ratio
#   points with low variance coming from lowly expressed genes that were distorting the signal.
#2) The main incorrect assumption was that inv variance would be the same as summing the points together
#   because lowly expressed genes were more variant but that is not always true. These genes were consistently
#   expressed at higher and lower ratios than the controls in the treated cells.


#QUESTION
#   Without just cutting off a significant number of genes how do I down weight the genes
#   with low variance but extreme TPM ratios?
#POSSIBLE ANSWERS
#1) Just use a cutoff as I did or just do a simple sum. Essentially weighing by expression. (similar to what is done currently)
#2) Differential Expression? Which genes are deferentially expressed throughout the experimental cells
#   any cell that is too consistently expressed will be problematic to determining changes in CN. 
#NOTE:
#1) I tried to do this EARLY on in the project ~1yr ago. I was able to find a set of genes what was more 
#   differentially expressed but I had to cut this short due to time constraints.
#2) When detecting changes in TPM ratio you cannot weigh genes that are extremely consistent too heavily
#   perhaps they are up or down regulated no matter the CN state. 


#SEPARATE NOTE
#1) The reason why the plots I did that were weighed based on expression were not exactly the same as just a 
#   simple avg was because the weights that were acting on the ratios were determined by the avg expression
#   across all or control cells rather than each cell's individual gene expression. 

###################################################################################################
# AF Exploration #
###################################################################################################

#read in ASE information
ASE <- readRDS(file=sprintf("%s/aggregated_results/ASE.bygene.rds", dirpath))
ASE_bychr <- readRDS(file=sprintf("%s/aggregated_results/ASE.inv_var.bychr.rds", dirpath))

#plot rowmeans vs row variance by gene
VAR_A_bygene <- ASE$AF[!which(ASE$AF$chr %in% c("chr10", "chr12", "chrX")),-c(1:6)]
plot(rowMeans(VAR_A_bygene, na.rm = T), rowVars(as.matrix(VAR_A_bygene), na.rm = T))

#plot rowmeans vs row variance by gene for chr10 specifically
VAR_A_bygene <- ASE$AF[which(ASE$AF$chr %in% c("chr12")),-c(1:6)]
VAR_A_bygene_2 <- ASE$AF[which(ASE$AF$chr %in% c("chr10")),]
VAR_A_bygene_3 <- VAR_A_bygene_2[which(VAR_A_bygene_2$start > 61000000), -c(1:6)]
plot(rowMeans(VAR_A_bygene_3, na.rm = T), rowVars(as.matrix(VAR_A_bygene_3), na.rm = T))

#plot rowmeans vs row variance by chr
VAR_A <- ASE_bychr$AF[,-c(1)]
plot(rowMeans(VAR_A, na.rm = T), rowVars(as.matrix(VAR_A), na.rm = T))

test <- cbind(ASE_bychr$AF[,c(1)], rowMeans(VAR_A, na.rm = T))

######################################################################################
### Exploring the chr16 cell AF and TPM ###
######################################################################################

TPM <- readRDS(file=sprintf("%s/aggregated_results/TPM.bygene.rds", dirpath))
ASE <- readRDS(file=sprintf("%s/aggregated_results/ASE.bygene.rds", dirpath))

TPM_spec <- unlist(TPM[which(TPM$chr == "chr2"),"210621_9D"])
AF_spec <- unlist(ASE$AF[which(ASE$AF$chr == "chrX"),"210621_9D"])

AF_spec_alt <- unlist(ASE$AF[which(ASE$AF$chr == "chr1"),"190628_4B"])

plot(log2(unlist(TPM[which(TPM$chr == "chr11"),"210621_9D"])), log2(unlist(rowMeans(TPM[which(TPM$chr == "chr11"),..controlSampleIDs]))), xlim = c(0,15), ylim = c(0,15))

TPM_normed <- readRDS(file=sprintf("%s/aggregated_results/TPM.bygene.normed.rds", dirpath))
mean(unlist(TPM_normed[which(TPM$chr == "chr11"),"210621_9D"]), na.rm = T)

test <- ASE_all_cells[which(ASE_all_cells$start == 18400940),c(1:3,which(colnames(ASE_all_cells) %in% c("210621_9D", "210621_9D_L2")))]


ASE$AF[which(ASE$AF$gene_id == "ENSG00000134333.13"), "210621_9D"]

View(ASE$AF[which(ASE$AF$chr == "chr11"),"210621_9D"])


######################################################################################
### TPM adjustment ###
######################################################################################

TPM <- readRDS(file=sprintf("%s/aggregated_results/TPM.nolim.rds", dirpath))
TPM <- readRDS(file=sprintf("%s/aggregated_results/TPM.bygene.rds", dirpath))

#Get the factors used to do global normalization 
global_factors <- function(mat, anno_cols = 6) {
  
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
  
  #return normed ratios
  return(model_mat)
}

factors <- global_factors(TPM)
setkey(factors, intercept)
boxplot(factors$intercept, outline = F)
boxplot(factors$slope, outline = F)

num1 <- .99
num2 <- .05
cell1 <- factors$cell[quantile(probs = num1, c(1:length(factors$cell)))]
cell2 <- factors$cell[quantile(probs = num2, c(1:length(factors$cell)))]

#sum of log TPM
cell1_TPM <- unlist(log2(TPM[,..cell1]))
cell1_TPM[cell1_TPM == -Inf] <- NA
mean(cell1_TPM, na.rm = T)

#sum of log TPM ctrl
ctrl_TPM <- log2(rowMeans(TPM[,..controlSampleIDs]))
ctrl_TPM[ctrl_TPM == -Inf] <- NA
mean(ctrl_TPM, na.rm = T)

#difference in scaling factor
cell2_TPM <- unlist(log2(TPM[,..cell2]))
cell2_TPM[cell2_TPM == -Inf] <- NA
mean(cell2_TPM, na.rm = T)

#y - x = b
dif <- mean(cell1_TPM, na.rm = T) - mean(ctrl_TPM, na.rm = T)
cell1_TPM_adjsuted <- cell1_TPM - dif
plot(ctrl_TPM, cell1_TPM_adjsuted, xlim = c(0,20), ylim = c(0,20))
plot(ctrl_TPM, cell1_TPM, xlim = c(0,20), ylim = c(0,20))
abline(a = 0, b = 1)



ctrl_TPM <- log2(rowMeans(TPM[,..controlSampleIDs]))

cell1_TPM <- unlist(log2(TPM[,..cell1]))
cell1_TPM[cell1_TPM == -Inf] <- NA



plot(ctrl_TPM, cell1_TPM, xlim = c(-15,15), ylim = c(-5,15),  main = unlist(factors[which(factors$cell == cell1),])) #quantile(probs = c(num1), factors)
abline(lm(cell1_TPM ~ ctrl_TPM), col = "red")
plot(ctrl_TPM, cell2_TPM, xlim = c(-15,15), ylim = c(-5,15), main = unlist(factors[which(factors$cell == cell2),])) #quantile(probs = c(num2), factors)
abline(lm(cell2_TPM ~ ctrl_TPM), col = "red")

cell1_TPM_adjust <- (cell1_TPM - factors$intercept[which(factors$cell == cell1)]) / (factors$slope[which(factors$cell == cell1)])

plot(ctrl_TPM, cell1_TPM_adjust, xlim = c(-15,15), ylim = c(-5,15), main = lm(cell1_TPM_adjust ~ ctrl_TPM)$coefficients) #quantile(probs = c(num1), factors)
abline(lm(cell1_TPM_adjust ~ ctrl_TPM), col = "red")

plot(ctrl_TPM, unlist(out[,c("210618_4C")]), xlim = c(0,10), ylim = c(0,10))
abline(lm(unlist(out[,c("210618_4C")]) ~ ctrl_TPM), col = "red")

plot(cell1_TPM, cell2_TPM, xlim = c(-15,15), ylim = c(-5,15), main = lm(cell2_TPM ~ cell1_TPM)$coefficients)
abline(lm(cell2_TPM ~ cell1_TPM), col = "red")

###################################################################################################
# AS-TPM Exploration of chr10 and 12 #
###################################################################################################

#read in ASE and TPM info
ASE <- readRDS(file=sprintf("%s/aggregated_results/ASE.bygene.rds", dirpath))
TPM <- readRDS(file=sprintf("%s/aggregated_results/TPM.nolim.rds", dirpath))

#calculate AS_TPM by genes
ASE_genes <- ASE$AF$gene_id
TPM_red <- TPM[which(TPM$gene_id %in% ASE_genes),]
AS_TPM_bygene_A <- cbind(TPM_red[,c(1:6)], ASE$AF[,-c(1:6)]*TPM_red[,-c(1:6)])
AS_TPM_bygene_B <- cbind(TPM_red[,c(1:6)], ((1-ASE$AF[,-c(1:6)])*TPM_red[,-c(1:6)]))
TPM_red_noanno <- TPM_red[which(TPM_red$chr == "chr12"),-c(1:6)]

#get chr specific 
AS_TPM_A <- AS_TPM_bygene_A[which(AS_TPM_bygene_A$chr == "chr12"),-c(1:6)]
AS_TPM_A[abs(AS_TPM_A) == Inf] <- NA
AS_TPM_B <- AS_TPM_bygene_B[which(AS_TPM_bygene_B$chr == "chr12"),-c(1:6)]
AS_TPM_B[abs(AS_TPM_B) == Inf] <- NA

#get colmeans for each cell
AS_TPM_A_means <- colMeans(AS_TPM_A, na.rm = T)
AS_TPM_B_means <- colMeans(AS_TPM_B, na.rm = T)
TPM_red_noanno_means <- colMeans(TPM_red_noanno, na.rm = T)

#plot AS_TPM_A cs AS_TPM_B
plot(log2(AS_TPM_A_means), log2(AS_TPM_B_means))

###################################################################################################

#get chr specific 
TPM_red_noanno <- TPM_red[which(TPM_red$chr == "chr10"),]
AS_TPM_A <- AS_TPM_bygene_A[which(AS_TPM_bygene_A$chr == "chr10"),]
AS_TPM_B <- AS_TPM_bygene_B[which(AS_TPM_bygene_B$chr == "chr10"),]

#Just greater than 61Mb
TPM_red_noanno_2 <- TPM_red_noanno[which(TPM_red_noanno$start > 61000000),-c(1:6)]
AS_TPM_A_2 <- AS_TPM_A[which(AS_TPM_A$start > 61000000),..controlSampleIDs]
AS_TPM_A_2[abs(AS_TPM_A_2) == Inf] <- NA
AS_TPM_B_2 <- AS_TPM_B[which(AS_TPM_B$start > 61000000),..controlSampleIDs]
AS_TPM_B_2[abs(AS_TPM_B_2) == Inf] <- NA

#plot AS_TPM_A cs AS_TPM_B
plot(log2(rowMeans(AS_TPM_A_2, na.rm = T)), log2(rowMeans(AS_TPM_B_2, na.rm = T)))

###################################################################################################

TPM_red <- TPM[which(TPM$gene_id %in% ASE_genes),]
AF <- ASE$AF

TPM_rednoanno <- TPM_red[which(TPM_red$chr == "chr10"),-c(1:6)]
AF_noanno <- AF[which(AF$chr == "chr10"),-c(1:6)]

plot(log2(rowMeans(TPM_rednoanno, na.rm = T)), rowMeans(AF_noanno, na.rm = T))

###################################################################################################

TPM_red <- TPM[which(TPM$gene_id %in% ASE_genes),]
AF <- ASE$AF

TPM_rednoanno <- TPM_red[which(TPM_red$chr == "chr10"),]
AF_noanno <- AF[which(AF$chr == "chr10"),]

TPM_61 <- TPM_rednoanno[which(TPM_rednoanno$start > 61000000),-c(1:6)]
AF_61 <- AF_noanno[which(AF_noanno$start > 61000000),-c(1:6)]

plot(log2(rowMeans(TPM_61, na.rm = T)), rowMeans(AF_61, na.rm = T))

###################################################################################################

cell = "210525_4G"


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



