#### Compare Objects Script ####
#This is used to compare Nikos Mynhier's data objects with Etai Jacob's

etaipath <- "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2/data"
nikospath <- "/pellmanlab/stam_niko/data/processed_bam/aggregated_results"

compdir_nikos <- "/pellmanlab/stam_niko/data/comparable_datasets/nikos_dataset"
compdir_etai <- "/pellmanlab/stam_niko/data/comparable_datasets/etai_dataset"

source('/pellmanlab/nikos/Stam_Etai_Scripts/runRNAseqPipelineFor_BF9VP.utils.R')
source('/pellmanlab/nikos/Stam_Etai_Scripts/generate_baseline_datasets_for_ML_and_stats_analysis.R')


###TODO:
#all objects need to be reordered to compare
#adt_nikos needs start and stop columns fixed, this is not always true. look into. 
#ASE needs allelic feature calculations to be fixed

#####################################################
### Random Sample of Genes and samples to compare ###
#####################################################
random_samples <- sample(x = colnames(tpm_nikos_new), size = 50, replace = FALSE)
random_genes <- sample(x = rownames(tpm_nikos_new), size = 50, replace = FALSE)

#################################
### Annotation list Breakdown ###
#################################
anno <- data.table(read.csv("/pellmanlab/nikos/Stam_Etai_Data/work_in_progress/Annotation_list_long_vnikos_210129.csv"))
#171122, 171218 samples are jinyu samples. can be taken out. 
#171205, 180227, 180228, 180301 are LCMlook samples. 
#170323 are Miseq samples
#Miseq samples can be taken out. new and old. 

Miseq_samples <- rbindlist(lapply(1:length(anno$Seq.platform), function(i)  {
  res <- grep(pattern = "Miseq", ignore.case = T, x = anno$Seq.platform[i])
  data.table(id = anno$WTA.plate[i], 
             idx = ifelse(length(res) > 0, res, -1)) }))
Miseq_samples <- Miseq_samples[idx>0]$id

Jinyu_samples <- rbindlist(lapply(1:length(anno$LookRNAseq.Exp), function(i)  {
  res <- grep(pattern = "Jinyu", ignore.case = T, x = anno$LookRNAseq.Exp[i])
  data.table(id = anno$WTA.plate[i], 
             idx = ifelse(length(res) > 0, res, -1)) }))
Jinyu_samples <- Jinyu_samples[idx>0]$id

LCMlook_samples <- rbindlist(lapply(1:length(anno$LookRNAseq.Exp), function(i)  {
  res <- grep(pattern = "LCMlook", ignore.case = T, x = anno$LookRNAseq.Exp[i])
  data.table(id = anno$WTA.plate[i], 
             idx = ifelse(length(res) > 0, res, -1)) }))
LCMlook_samples <- LCMlook_samples[idx>0]$id

Excluded_samples <- union(Miseq_samples, Jinyu_samples)
Excluded_samples2 <- union(Excluded_samples, LCMlook_samples)

###########################################
### single experiment tpm objects: Done ###
###########################################
#These are the same cells in two different lanes. They have similar tpm.
tpm_SISf_L1_nikos <- readRDS("/pellmanlab/stam_niko/data/processed_bam/SIS1025f_Lane1/expr_results/SIS1025f_Lane1_rsemtpm.rds")
tpm_SISf_L2_nikos <- readRDS("/pellmanlab/stam_niko/data/processed_bam/SIS1025f_Lane2/expr_results/SIS1025f_Lane2_rsemtpm.rds")

#checking for a samples
tpm_SISe_nikos <- readRDS("/pellmanlab/stam_niko/data/processed_bam/SIS1025e/expr_results/SIS1025e_rsemtpm.rds")
colnames(tpm_SISe_nikos)

The_sample <- rbindlist(lapply(1:length(colnames(tpm_SISe_nikos)), function(i)  {
  res <- grep(pattern = "181015_Nextera_1F10", ignore.case = T, x = colnames(tpm_SISe_nikos)[i])
  data.table(id = colnames(tpm_SISe_nikos)[i], 
             idx = ifelse(length(res) > 0, res, -1)) }))
The_sample <- The_sample[idx>0]$id

###################
### tpm objects ###
###################
# if you skip to last 50 values you see that the values are fairly similar as with the big tpm objects.
# I need to exclude some old samples from tpm_etai
# 39-153 inclusive are old samples. that is 115 samples. Anno v15 only has 106
tpm_etai <- readRDS(sprintf("%s/rsemtpm.rds", etaipath))
tpm_nikos <- readRDS(sprintf("%s/all_experiments_rsemtpm.rds", nikospath))

Etai_Samples <- colnames(tpm_etai)
Nikos_Samples <- colnames(tpm_nikos)

Missing_Samples <- setdiff(Etai_Samples, Nikos_Samples)
All_Samples_v15 <- anno$WTA.plate

rows <- c()
for (i in 1:length(Missing_Samples)) {
  for (j in 1:length(All_Samples_v15)) {
    if (Missing_Samples[i]==All_Samples_v15[j]) {
      rows <- append(rows, j)
    }
  }
}

anno_missing_samples <- anno[rows, ]
#write.csv(anno_missing_samples, "/pellmanlab/nikos/Stam_Etai_Data/work_in_progress/Annotation_list_excluded_samples.csv")

columns <- c()
for (i in 1:length(Missing_Samples)) {
  for (j in 1:length(colnames(tpm_etai))) {
    if (Missing_Samples[i]==colnames(tpm_etai)[j]) {
      columns <- append(columns, j)
    }
  }
}
tpm_etai_new <- tpm_etai[,-columns]
tpm_etai_new <- tpm_etai_new[,order(colnames(tpm_etai_new))]

columns <- c()
for (i in 1:length(colnames(tpm_etai_new))) {
  for (j in 1:length(colnames(tpm_nikos))) {
    if (colnames(tpm_etai_new)[i]==colnames(tpm_nikos)[j]) {
      columns <- append(columns, j)
    }
  }
}
tpm_nikos_new <- tpm_nikos[,columns]
tpm_nikos_new <- tpm_nikos_new[,order(colnames(tpm_nikos_new))]

rows_missing <- setdiff( rownames(tpm_etai_new), rownames(tpm_nikos_new) )
rows <- c()
for (i in 1:length(rows_missing)) {
  for (j in 1:length(rownames(tpm_etai_new))) {
    if (rows_missing[i]==rownames(tpm_etai_new)[j]) {
      rows <- append(rows, j)
    }
  }
}
tpm_etai_new <- tpm_etai_new[-rows,]

minDetectionLevel = 5 #option1: 5 #original: 5
minNumOfSamplesToDetect = 3 #option1: 10 #original: 3
wi <- which(rowSums(tpm_nikos_new[, controlSampleIDs] >= minDetectionLevel) >= minNumOfSamplesToDetect)
tpm_nikos_new <- tpm_nikos_new[wi,]

#import geneRanges annotations
geneRanges <- readRDS(sprintf("%s/geneRanges.rds", nikospath))
genes <- geneRanges$id

#Intersect and reduce geneRanges and tpm genes
final_genes <- intersect(rownames(tpm_nikos_new), genes)
tpm_nikos_new <- tpm_nikos_new[final_genes,]
tpm_etai_new <- tpm_etai_new[final_genes,]
geneRanges <- geneRanges[geneRanges$id %in% final_genes]

#export final objects
mk_agg_dir <- sprintf("mkdir -p %s/aggregated_results", compdir_nikos)
system(mk_agg_dir)
mk_agg_dir <- sprintf("mkdir -p %s/aggregated_results", compdir_etai)
system(mk_agg_dir)
saveRDS(tpm_nikos_new, sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", compdir_nikos))
saveRDS(tpm_etai_new,  sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", compdir_etai) )

#Final dfs
rsemtpm_nikos <- as.data.frame(cbind(geneRanges, tpm_nikos_new))
rsemtpm_etai <- as.data.frame(cbind(geneRanges, tpm_etai_new))

### compare dfs ###
#abs_tpm_difference <- abs(tpm_etai_new - tpm_nikos_new) 
#avg_bysample_abs_tpm_difference <- colMeans(abs_tpm_difference)
#median_avg_bysample_abs_tpm_difference <- median(avg_bysample_abs_tpm_difference)
#print(median_avg_bysample_abs_tpm_difference)
#boxplot(avg_bysample_abs_tpm_difference)
#
#tpm_difference <- tpm_etai_new - tpm_nikos_new 
#avg_bysample_tpm_difference <- colMeans(tpm_difference)
#median_avg_bysample_tpm_difference <- median(avg_bysample_tpm_difference)
#print(median_avg_bysample_tpm_difference)
#boxplot(avg_bysample_tpm_difference)


#Final tpm dataframe
tpm_difference <- tpm_etai_new - tpm_nikos_new 
rsemtpm_dif <- cbind(geneRanges, tpm_difference)

#reorder df
rsemtpm_dif <- rsemtpm_dif[order(seqnames, start)]

#creates a column in rsemtpm with the bin number. putting genes in groups of 50.
bins <- lapply(unique(rsemtpm_dif$seqnames), function(chr) unlist(lapply(1:ceiling(nrow(rsemtpm_dif[seqnames==chr])/50), rep, 50))[1:nrow(rsemtpm_dif[seqnames==chr])])
names(bins) <- unique(rsemtpm_dif$seqnames)
bins2 <- rbindlist(lapply(names(bins), function(chr) data.table(seqnames = chr, bin = bins[[chr]])))
stopifnot(sum(bins2$seqnames != rsemtpm_dif$seqnames)==0)
rsemtpm_dif2 <- cbind(bin = bins2$bin, rsemtpm_dif)
rsemtpm_dif3 <- rsemtpm_dif2[, lapply(.SD, function(x) mean(x, na.rm = T)), .SDcols = c(colnames(rsemtpm_dif2)), by = c("seqnames", "bin")]
rsemtpm_dif4 <- rsemtpm_dif3[,-c(3:9)]

#Get just chr1
chrom <- "chr1"
rsemtpm_dif4_chr1 <- rsemtpm_dif4[rsemtpm_dif4$seqnames == chrom]
rsemtpm_dif4_chr1 <- data.matrix(t(rsemtpm_dif4_chr1[,-c(1:2)]))
boxplot(rsemtpm_dif4_chr1, outline = F)

#expected size 12238 x 797
expected_size_flat_v <- (12238 * 797)
library(purrr)
tpm_etai_new <- as.data.frame.logical(tpm_etai_new)
tpm_nikos_new <- as.data.frame.logical(tpm_nikos_new)
tpm_etai_new_flat <- flatten(tpm_etai_new)
tpm_nikos_new_flat <- flatten(tpm_nikos_new)

rsemtpm_nikos_chr1 <- as.data.frame(data.table(rsemtpm_nikos)[rsemtpm_nikos$seqnames == chrom])
rsemtpm_etai_chr1 <- as.data.frame(data.table(rsemtpm_etai)[rsemtpm_etai$seqnames==chrom])
rsemtpm_nikos_chr1_flat <- flatten(rsemtpm_nikos_chr1[,-c(1:6)])
rsemtpm_etai_chr1_flat <- flatten(rsemtpm_etai_chr1[,-c(1:6)])

plot.default(x = rsemtpm_nikos_chr1_flat, y = rsemtpm_etai_chr1_flat,
             xlim = c(0, 2000), ylim = c(0, 2000))



#write.csv(tpm_difference, "/pellmanlab/stam_niko/data/processed_bam/misc/Difference_between_tpm_matrices.csv")
#avgdif <- c()
#for ( i in 1:(length(rownames(tpm_difference))) ) {
#  avgdif <- append(avgdif, mean(tpm_difference[i,], na.rm = TRUE))
#}
#
#for (i in 1:length(avgdif)) {
#  if (is.nan(avgdif[i]) == F) {
#    if (avgdif[i] == Inf) {
#      avgdif[i] <- NaN
#    }
#  }
#}
#
#meandif <- mean(avgdif, na.rm = TRUE)


# TODO:
# 1) scatter plot of my points vs etai's points 
# 2) see if I do a z-score adt object if it looks the same as etai's given adt object. could also use scatter plot. 

###################################
### adt object: Mostly Complete ###
###################################
#creates adt rds objects for ML
#The based a small sample of test genes and samples the values are very similar. like the tpm matrix.
adt_etai <- readRDS(sprintf("%s/adt.rds", etaipath))
adt_nikos <- readRDS(sprintf("%s/adt.rds", nikospath))
#adt_nikos_zs <- readRDS(sprintf("%s/adt.zscore.rds", nikospath))
adt_etai <- as.data.frame(adt_etai)
adt_nikos <- as.data.frame(adt_nikos)
#adt_nikos_zs <- as.data.frame(adt_nikos_zs)

controls_e <- adt_etai[,controlSampleIDs]
controls_n <- adt_nikos[rownames(adt_etai),controlSampleIDs]

rownames(adt_etai) <- adt_etai[,1]
rownames(adt_nikos) <- adt_nikos[,1]

rownames(adt_nikos_zs) <- adt_nikos_zs[,1]

adt_etai_sampled <- adt_etai[random_genes, random_samples]
adt_nikos_sampled <- adt_nikos[random_genes, random_samples]
adt_nikos_zs_sampled <- adt_nikos_zs[random_genes, random_samples]

adt_nikos_short <- adt_nikos[rownames(adt_etai) , -c(1:4)]
adt_etai_short <- adt_etai[, -c(1:4)]

shared_columns <- intersect(colnames(adt_etai_short), colnames(adt_nikos_short))

adt_nikos_short <- adt_nikos_short[, shared_columns]
adt_etai_short <- adt_etai_short[, shared_columns]

adt_nikos_short_ordered <- adt_nikos_short[order(rownames(adt_nikos_short)),order(colnames(adt_nikos_short))]
adt_etai_short_ordered <- adt_etai_short[order(rownames(adt_etai_short)),order(colnames(adt_etai_short))]


adt_difference <- (adt_etai_short_ordered - adt_nikos_short_ordered)

##############################
### ASE object: Incomplete ###
##############################
# These variants are notably different but they should be exactly the same. 
ASE_etai <- readRDS(sprintf("%s/ASE.rds", etaipath))
ASE_nikos <- readRDS(sprintf("%s/ASE.rds", nikospath))

ASE_etai_A <- as.data.frame(ASE_etai$A)
ASE_nikos_A <- as.data.frame(ASE_nikos$A)

columns <- intersect(colnames(ASE_nikos_A), colnames(ASE_etai_A))
ASE_etai$A <- as.data.frame(ASE_etai$A)[,columns]
ASE_nikos$A <- as.data.frame(ASE_nikos$A)[,columns]
ASE_etai$B <- as.data.frame(ASE_etai$B)[,columns]
ASE_nikos$B <- as.data.frame(ASE_nikos$B)[,columns]

saveRDS(ASE_etai, sprintf("%s/aggregated_results/ASE.rds", compdir_etai))
saveRDS(ASE_nikos, sprintf("%s/aggregated_results/ASE.rds", compdir_nikos))

#################################
# Generate comparable data sets #
#################################

#controls
controlSampleIDs <- anno[ (mainGroup == "A_level_control" | mainGroup == "B_level_control") & th1 > 6000 ]$WTA.plate
controlSampleIDs <- controlSampleIDs[controlSampleIDs %in% colnames(all_rsemtpm)]
saveRDS(controlSampleIDs, sprintf("%s/aggregated_results/controlSampleIDs.rds", compdir_etai))
saveRDS(controlSampleIDs, sprintf("%s/aggregated_results/controlSampleIDs.rds", compdir_nikos))

controlSampleIDs2 <- anno[ (mainGroup == "A_level_control") & th1 > 6000 ]$WTA.plate
controlSampleIDs2 <- controlSampleIDs2[controlSampleIDs2 %in% colnames(all_rsemtpm)]
saveRDS(controlSampleIDs2, sprintf("%s/aggregated_results/controlSampleIDs2.rds", compdir_etai))
saveRDS(controlSampleIDs2, sprintf("%s/aggregated_results/controlSampleIDs2.rds", compdir_nikos))

#geneRanges
load("/pellmanlab/stam_niko/refgenomes/Gencode/v25/gencode.v25.primary_assembly.ERCC92.annotation.geneRanges.RData")
names(geneRanges)[1] <- "id"
saveRDS(geneRanges, sprintf("%s/aggregated_results/geneRanges.rds", compdir_etai))
saveRDS(geneRanges, sprintf("%s/aggregated_results/geneRanges.rds", compdir_nikos))

#creates adt rds objects for ML
tmppath <- "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2"
configs <- readRDS(sprintf("%s/data/param_config_list.rds", tmppath))

#To me it looks like etai's df is standardized as well
generate_adt_adt0_adt.na_and_nonzeros_data(dirpath = compdir_nikos, scriptsdir = "/pellmanlab/nikos/Stam_Etai_Scripts", th = 4) #normBySd = T would standardize the data
generate_adt_adt0_adt.na_and_nonzeros_data(dirpath = compdir_etai, scriptsdir = "/pellmanlab/nikos/Stam_Etai_Scripts", th = 4)

#creates coding_snp rds object for ML
generate_coding_snps_data(dirpath = compdir_nikos)
generate_coding_snps_data(dirpath = compdir_etai)

#creates ASE.coding.rds and ASE.noncoding.rds
get_snps_fraction_bins_zscore(saveMe = T, dirpath = compdir_nikos)
get_snps_fraction_bins_zscore(saveMe = T, dirpath = compdir_etai)

#creates the nonzeros.zs.bin50.rds and nonzeros.bin50.rds objects
nonzeros <- compute_nonZero_bins_zscore(saveMe = T, dirpath = compdir_nikos)
nonzeros <- compute_nonZero_bins_zscore(saveMe = T, dirpath = compdir_etai)

coding <- readRDS(file = sprintf("%s/ASE.coding.rds", nikospath))

