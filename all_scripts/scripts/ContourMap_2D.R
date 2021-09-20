###################################################
### 2D Contour Map (TPM Ratio vs Variant Ratio) ###
###################################################

### Goals ###
#We want to see clustering of chromosomes based on CN state
#We want to see how exclusions effect clustering inorder to determine best parameters

### Passed Arguments ###

args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/rerun_6_9_2021/data", "/pellmanlab/nikos/Stam_Etai_Data")
#args <- commandArgs(trailingOnly = TRUE)
print(args)
scriptsdir <- args[1]
dirpath <- args[2]
datadir <- args[3]
mk_agg_dir <- sprintf("mkdir -p %s/aggregated_results", dirpath)
system(mk_agg_dir)

### Source Scripts ###

require(data.table)
require(readxl)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
#install.packages("factoextra")
library(factoextra) # clustering algorithms & visualization
#source( sprintf('%s/scripts/runRNAseqPipelineFor_BF9VP.utils.R', scriptsdir) )
#source( sprintf('%s/scripts/generate_baseline_datasets_for_ML_and_stats_analysis.R', scriptsdir) )
#source( sprintf('%s/plots/class_prob_to_cn_level.R', scriptsdir) )
#
#
#### Imports ###
#
##centromeres object is called inside AddFeaturesandPhasing function
#centromeres <- readRDS( sprintf("%s/centromeres.rds", datadir) )
#
##raw TPM matrix
#rsemtpm <- readRDS(sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", wkpath))
#
##haplotype phased variants
#mym <- readRDS(file = sprintf("%s/aggregated_results/haplotypePhasedADs.rds", wkpath))
#
##Add features to snps
#ASE <- readRDS(file = sprintf("%s/aggregated_results/ASE.rds", wkpath))
#
##Allele level data:
#alleles.all <- prepareAlleleData(m = ASE, useFraction = F , controlSampleIDs = controlSampleIDs)
#saveRDS(alleles.all, file = sprintf("%s/aggregated_results/alleles.all.rds", wkpath))
#
#alleles.all.frac <- prepareAlleleData(m = ASE, useFraction = T , controlSampleIDs = controlSampleIDs)
#saveRDS(alleles.all.frac, file = sprintf("%s/aggregated_results/alleles.all.frac.rds", wkpath))

#Annotion list with information on every cell and IDs of control samples
anno <- data.table(read.csv(sprintf("%s/work_in_progress/Annotation_list_long_vnikos_11_6_21.csv", datadir) ))
controlSampleIDs <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))

#control samples with aneuploidies
#c("161228_A2", "")

#QC Data Frame and high QC IDs based on 5 read cnt threshold for genes
all_QC <- readRDS(file = sprintf("%s/aggregated_results/all_QC.rds", dirpath))
high_qc_ids <- as.character(all_QC[which(all_QC$th5 >= quantile(all_QC$th5, c(.10))),]$id) #Excludes bottom 10%

#Standard set of annotations for all genes
#geneRanges <- readRDS(sprintf("%s/geneRanges.rds", datadir))

#list of preset parameters for the experiment
#configs <- readRDS(sprintf("%s/param_config_list.rds", datadir))
#configs$minDetectionLevel <- 5 #5, 50,
#configs$minNumOfSamplesToDetect <- 5 #5, 5, 

##normalize tpm matrix, limit high and low expression
#generate_adt_adt0_adt.na_and_nonzeros_data(dirpath = wkpath, th = 4, scriptsdir = scriptsdir, datadir = datadir, normBySd = F) 
#
##creates coding_snp rds object for ML
#generate_coding_snps_data(dirpath = wkpath)
#
##creates ASE.coding.rds and ASE.noncoding.rds
#get_snps_fraction_bins_zscore(saveMe = T, dirpath = wkpath)
#
##creates the nonzeros.zs.bin50.rds and nonzeros.bin50.rds objects
#nonzeros <- compute_nonZero_bins_zscore(saveMe = T, dirpath = wkpath)
#
##Creates arm level annotations
#generate_ganno2(dirpath = wkpath, datapath = datadir)

#Data Frames with TPM and Variant information
adt <- readRDS(sprintf("%s/aggregated_results/adt.rds", dirpath)) #Normalized TPM matrix
coding <- readRDS(sprintf("%s/aggregated_results/ASE.coding.rds", dirpath)) #Raw Variant matrix in the coding and UTR regions

#Var cnts per 50 genomically sequential genes for allele A and B
var_mat_A <- coding$cnts.A
var_mat_B <- coding$cnts.B
var_mat <- cbind(var_mat_A[,c(1:2)], var_mat_A[,-c(1:2)] + var_mat_B[,-c(1:2)])

### Exclude noisy chromosomes a priori ###
exclude_chrs <- c("chr13", "chr18", "chr21", "chr22", "chrX")
#exclude_chrs <- c("chrX")
adt <- adt[-which(adt$seqnames %in% exclude_chrs),]
var_mat <- var_mat[-which(var_mat$seqnames %in% exclude_chrs),]
var_mat$seqnames <- sub(pattern = "chr", replacement = "", x = var_mat$seqnames)
var_mat$seqnames <- as.integer(sub(pattern = "X", replacement = "23", x = var_mat$seqnames))
var_mat <- var_mat[order(seqnames, bin)]

### Re-Engineer raw variant count Matrix ###
#Step 1: Get Binary Variant Matrix based on thresholds. No pass = 0 ; Pass = 1 [2D Binary Matrix]
#Step 2: Calculate the fraction of cells that pass threshold for each variant [1D Vector of Fractions]
#Step 3: Plot frequency distribution of binned fractions  [Histogram of 1D Vector of Fractions]
#Step 4: Based on distribution of those fractions decide a cutoff fraction to exclude variants seen infrequently in cells
#NOTE: Step 4 could be skipped if that exclusions is determined to be unnecessary. 
#Step 5: Normalize cells by dividing values by the mean of the control cells for each remaining variant. 
#Step 6: Normalize variants by dividing values by the mean/median of variants for each cell
#Bin adt 
adt <- adt[order(seqnames, start)]

#creates a column in adt with the bin number. putting genes in groups of 50.
bins <- lapply(unique(adt$seqnames), function(chr) unlist(lapply(1:ceiling(nrow(adt[seqnames==chr])/50), rep, 50 ))[1:nrow(adt[seqnames==chr])])
#bins <- lapply(unique(adt$seqnames), function(chr) unlist(rep(c(1:floor(nrow(adt[seqnames==chr])/50)), times=c(rep(50, floor(nrow(adt[seqnames==chr])/50)-1), 50+(nrow(adt[seqnames==chr]) - floor(nrow(adt[seqnames==chr])/50)*50)) ))[1:nrow(adt[seqnames==chr])])
names(bins) <- unique(adt$seqnames)
bins2 <- rbindlist(lapply(names(bins), function(chr) data.table(seqnames = chr, bin = bins[[chr]])))
stopifnot(sum(bins2$seqnames != adt$seqnames)==0)
adt2 <- cbind(bin = bins2$bin, adt)
adt3 <- cbind(adt2[,c(1:5)], 2^adt2[,-c(1:5)])

#normalize adt columns
adt3_colmeans <- colMeans(adt3[,-c(1:5)], na.rm = T)
adt3_colmeans_mat_reciprocal <- data.table(1/(t(replicate(nrow(adt3), adt3_colmeans))))
adt3_colmeans_mat_reciprocal[adt3_colmeans_mat_reciprocal == Inf] <- 0 
adt3 <- cbind(adt3[,c(1:5)] , (adt3[,-c(1:5)] * adt3_colmeans_mat_reciprocal))

#Aggregate values by bin
adt.bin <- adt3[,-c(2,4,5)] %>% 
  group_by(seqnames, add=T) %>%
  group_by(bin, add=T) %>%
  summarise_all(mean, na.rm = TRUE)
adt.bin <- data.table(adt.bin)
adt.bin$seqnames <- sub(pattern = "chr", replacement = "", x = adt.bin$seqnames)
adt.bin$seqnames <- as.integer(sub(pattern = "X", replacement = "23", x = adt.bin$seqnames))
adt.bin <- adt.bin[order(seqnames, bin)]
saveRDS(adt.bin, file = sprintf("%s/aggregated_results/adt.bin.rds", dirpath))

#check rows are the same
if (dim(adt.bin)[1] != dim(var_mat)[1]) {
  adt.bin <- adt.bin[-which(adt.bin$bin != var_mat$bin)[1],]
}

#remove errorneous bins 
#adt.bin <- adt.bin[-which(unlist(rowMeans(adt.bin[,..controlSampleIDs]))*2>2.5),]
#var_mat <- var_mat[-which(unlist(rowMeans(adt.bin[,..controlSampleIDs]))*2>2.5),]

#Check that binning is the same
stopifnot(adt.bin$bin == var_mat$bin)

#Normalize variant dataframe
var_rowmeans <- rowMeans(var_mat[,..controlSampleIDs], na.rm = T)
var_rowmeans_mat_reciprocal <- data.table(1/(replicate(ncol(var_mat[,-c(1:2)]), var_rowmeans)))
var_rowmeans_mat_reciprocal[var_rowmeans_mat_reciprocal == Inf] <- 0 
var_mat <- cbind(var_mat[,c(1:2)] , (var_mat[,-c(1:2)] * var_rowmeans_mat_reciprocal))

var_colmeans <- colMeans(var_mat[,-c(1:2)], na.rm = T)
var_colmeans_mat_reciprocal <- data.table(1/(t(replicate(nrow(var_mat), var_colmeans))))
var_colmeans_mat_reciprocal[var_colmeans_mat_reciprocal == Inf] <- 0 
var_mat <- cbind(var_mat[,c(1:2)] , (var_mat[,-c(1:2)] * var_colmeans_mat_reciprocal))

#check cols are the same
colIDs <- intersect(colnames(var_mat), colnames(adt.bin))
adt.bin <- adt.bin[,..colIDs]; var_mat <- var_mat[,..colIDs];

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

##Flatten dataframes
#flat_adt.bin <- flatten(adt.bychr[,..controlSampleIDs])
#flat_var_mat <- flatten(var_mat.bychr[,..controlSampleIDs])
#c(5,8,9,10,11,13,14,15,23) , c(1,12,16:18)

#I Need to revisit the calculation/aggregation of the adt object because it doesnt seem to look accurate.
#This also wont work because I need to combine the alleles before I plot

#Notes: 
#It seems that the few number of aneuploidies is not seriously affecting the means and therefore not likely responsible for high spread. 

#Methods for finding news exclusion ideas:
#go after outliers in the scatter plot and see why they are outside the clusters. 

#TODO: Things to try and increasing clustering of copy number states
#1) Restrict to less variable chroms
#2) Restrict to high qc samples
#3) Cluster on log scale?
#4) New gene exclusion criteria
#   a) Keep only highly expressed genes
#   b) Keep only differentially expressed genes
#5) Use only biallelic germline SNPs for total CN information even though you lose phasing.

#6) implement this clustering 8/11/21 https://uc-r.github.io/kmeans_clustering#prep 

#adt.bychr.control <- adt.bychr[,..controlSampleIDs]
#adt.bychr.control.1_5 <- adt.bychr[c(1,12,16:18),..controlSampleIDs]

#THE LAST THING I DID WAS RUN WITH 5 MIN hits for genes

#Plot 
#plot(2*unlist(flat_adt.bin), 2*unlist(flat_var_mat), xlim = c(0,5), ylim = c(0,5))
#plot(unlist(adt.bychr[,c(552)])*2, unlist(var_mat.bychr[,c(552)])*2, xlim = c(0,5), ylim = c(0,5))
#plot(unlist(adt.bychr[c(4),..high_qc_ids])*2, unlist(var_mat.bychr[c(4),..high_qc_ids])*2, xlim = c(0,5), ylim = c(0,5))
#plot(unlist(rowMeans(adt.bin[,..controlSampleIDs]))*2, unlist(rowMeans(var_mat[,..controlSampleIDs]))*2, xlim = c(0,5), ylim = c(0,5))

setkey(anno, WTA.plate)
families <- sort(table(anno[Pairs != "NA" & Pairs != "mother"]$Pairs), decreasing = T)

#Prepare Visual Objects
myfamily = names(families)[25] #43, 24, 
message(myfamily)
myids <- anno[Pairs %in% myfamily]$WTA.plate
visual.data <- data.table(TPM = as.numeric(flatten(adt.bychr[,..myids])))
visual.data$VAR <- as.numeric(flatten(var_mat.bychr[,..myids]))
visual.data$chr <- as.factor(rep(unique(adt.bin$seqnames), length = length(visual.data$VAR)))
visual.data$cell <- rep(myids, each=23-length(exclude_chrs))

#Second Draft Plots
base_visual <- ggplot(data = visual.data, mapping = aes(x = TPM, y = VAR, color = chr, shape = cell)) +
                      geom_point(mapping = aes(x = TPM, y = VAR, color = chr, shape = cell)) +
                      xlim(c(0,2)) + ylim(c(0,2)) +
                      labs(x = "TPM Ratio", y = "SNP Count Ratio", title = sprintf("Expression Contour Plot | %s ", myfamily)) +
                      geom_abline(aes(intercept = 0, slope = 2/3, alpha = .25)) +  # linetype = "dotted", size = .1
                      geom_abline(aes(intercept = 0, slope = 3/2, alpha = .25))    # linetype = "dotted", size = .1
grid.arrange(base_visual)

#make_path <- sprintf("mkdir -p %s/byfamily/", destDir)
#system(make_path)
#
#foreach(i = c(1:length(names(families)))) %dopar% {
#  #for(i in 3:4) {
#  myfamily = names(families)[i]
#  message(myfamily)
#  myids <- anno[Pairs %in% myfamily]$WTA.plate
#  pdf(file = sprintf("%s/byfamily/%s.ContourMap.pdf", destDir, myfamily), width = 18, height = 12)
#
#  #print(p)
#  dev.off()
#}
#
#
#flat_var_cnts <- flatten(coding$A[,-c(1:4)])
#length(flat_var_cnts[flat_var_cnts >= 2]) #<- 0
#flat_var_cnts[flat_var_cnts >= 2] <- 1
#length(flat_var_cnts[flat_var_cnts == 1])
#sum(as.numeric(flat_var_cnts))


#### Normalize seperately ###
#
##Normalize Alelle A dataframe
#var_A_rowmeans <- rowMeans(var_mat_A[,..controlSampleIDs], na.rm = T)
#var_A_rowmeans_mat_reciprocal <- data.table(1/(replicate(ncol(var_mat_A[,-c(1:2)]), var_A_rowmeans)))
#var_A_rowmeans_mat_reciprocal[var_A_rowmeans_mat_reciprocal == Inf] <- 0 
#var_mat_A <- cbind(var_mat_A[,c(1:2)] , (var_mat_A[,-c(1:2)] * var_A_rowmeans_mat_reciprocal))
#
#var_A_colmeans <- colMeans(var_mat_A[,-c(1:2)], na.rm = T)
#var_A_colmeans_mat_reciprocal <- data.table(1/(t(replicate(nrow(var_mat_A), var_A_colmeans))))
#var_A_colmeans_mat_reciprocal[var_A_colmeans_mat_reciprocal == Inf] <- 0 
#var_mat_A <- cbind(var_mat_A[,c(1:2)] , (var_mat_A[,-c(1:2)] * var_A_colmeans_mat_reciprocal))
#var_mat_A <- var_mat_A[order(seqnames, bin)]
#
##Normalize Alelle B dataframe
#var_B_rowmeans <- rowMeans(var_mat_B[,..controlSampleIDs], na.rm = T)
#var_B_rowmeans_mat_reciprocal <- data.table(1/(replicate(ncol(var_mat_B[,-c(1:2)]), var_B_rowmeans)))
#var_B_rowmeans_mat_reciprocal[var_B_rowmeans_mat_reciprocal == Inf] <- 0 
#var_mat_B <- cbind(var_mat_B[,c(1:2)] , (var_mat_B[,-c(1:2)] * var_B_rowmeans_mat_reciprocal))
#
#var_B_colmeans <- colMeans(var_mat_B[,-c(1:2)], na.rm = T)
#var_B_colmeans_mat_reciprocal <- data.table(1/(t(replicate(nrow(var_mat_B), var_B_colmeans))))
#var_B_colmeans_mat_reciprocal[var_B_colmeans_mat_reciprocal == Inf] <- 0 
#var_mat_B <- cbind(var_mat_B[,c(1:2)] , (var_mat_B[,-c(1:2)] * var_B_colmeans_mat_reciprocal))
#var_mat_B <- var_mat_B[order(seqnames, bin)]
#
##Aggregaye by chr seperately
#var_mat_A.bychr <- var_mat_A[,-c(2)] %>%
#  group_by(seqnames) %>%
#  summarise_all(mean, na.rm = TRUE)
#var_mat_A.bychr <- data.table(var_mat_A.bychr)
#var_mat_A.bychr <- var_mat_A.bychr[order(seqnames)]
#
#var_mat_B.bychr <- var_mat_B[,-c(2)] %>%
#  group_by(seqnames) %>%
#  summarise_all(mean, na.rm = TRUE)
#var_mat_B.bychr <- data.table(var_mat_B.bychr)
#var_mat_B.bychr <- var_mat_B.bychr[order(seqnames)]
