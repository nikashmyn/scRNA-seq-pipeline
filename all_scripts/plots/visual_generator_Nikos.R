#################################################################################
# Passed arguments:
#################################################################################
#args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/rerun_6_9_2021/data", "/pellmanlab/nikos/Stam_Etai_Data")
args <- commandArgs(trailingOnly = TRUE)
print(args)
scriptsdir <- args[1]
dirpath <- args[2]
datadir <- args[3]
mk_vis_dir <- sprintf("mkdir -p %s/visual_results", dirpath)
system(mk_vis_dir)

#################################################################################
# loading data:
#################################################################################

require(data.table)
require(gtools)
require(readxl)
require(ggplot2)
require(matrixStats)
require(reshape2)
require(gridExtra)
require(dplyr)
source(sprintf("%s/scripts/CNML.R", scriptsdir))
source(sprintf("%s/plots/ClusterMap_byfamily.R", scriptsdir))
source(sprintf("%s/plots/class_prob_to_cn_level.R", scriptsdir))
source(sprintf("%s/plots/plot_cn_prediction_probabilities.R", scriptsdir))
source(sprintf("%s/plots/calc_and_plot_intermediate_pvalue.R", scriptsdir))
source(sprintf("%s/plots/plot_raw_data_and_prediction_boxplots2.R", scriptsdir))
source(sprintf("%s/plots/plot_barplots_of_AllelicAndExpBiasPerSamples.R", scriptsdir))
source(sprintf('%s/plots/plot_cn_and_intermediate_prediction_latest_model_intermediate_bars.R', scriptsdir))
source(sprintf('%s/plots/plot_pvals_and_tpm_distributions_bychr.R', scriptsdir))

#############################Parameters##################################

#Sample families to exclude
exclude <- c("113015_T12", "061516_T12", "070516_LookRNAseq", "072516_LookRNAseq", 
             "071816_LookRNAseq", "071816_LookRNAseq", "161005_LookRNAseq", 
             "061316_LookRNAseq", "161005_LookRNAseq_control", "Jinyu", "170726_MCF10A_control",
             "040416_ATCC", "161017_LookRNAseq", "062716_LookRNAseq", "080816LookRNAseq", "080816_LookRNAseq",
             "071416_LookRNAseq", "161003_LookRNAseq")

#Window size used in models
NUM_OF_FEATURES <- 100

##########################Data for Raw Plots###########################
message("Loading data...")

#log-space centered tpm object
adt <- readRDS(sprintf("%s/aggregated_results/adt.rds", dirpath))
adt.na <- readRDS(sprintf("%s/aggregated_results/adt.na.rds", dirpath))
adt <- adt.default <- adt[order(seqnames, start, end)]
adt <- cbind( adt[,c(1:4)], setcolorder(adt[,-c(1:4)], order(colnames(adt[,-c(1:4)]))) )

#TMP: Create a fake object to feed as ML preds. 
adt_fake <- adt[,-c(1:4)]
adt_fake[,] <- 0
adt_fake <- cbind(adt[,c(1:4)], adt_fake)

#raw tpm object
rsemtpm <- readRDS(sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", dirpath))
#high_qc_ids <- names(which(colSums(rsemtpm>5)>4000)) #original criteria from etai
all_QC <- readRDS(file = sprintf("%s/aggregated_results/all_QC.rds", dirpath))
high_qc_ids <- as.character(all_QC[which(all_QC$th5 >= quantile(all_QC$th5, c(.10))),]$id)

#Variant matrix in the coding and UTR regions
coding <- readRDS(sprintf("%s/aggregated_results/ASE.coding.rds", dirpath))
#Low var counts exclusion. Note: both distributions of var counts are fairly similar.
high_varqc_idsA <- names(which(colSums(coding$cnts.A[,-c(1:2)]) > quantile(colSums(coding$cnts.A[,-c(1:2)]), c(.10))))
high_varqc_idsB <- names(which(colSums(coding$cnts.B[,-c(1:2)]) > quantile(colSums(coding$cnts.B[,-c(1:2)]), c(.10))))
high_varqc_ids <- intersect(high_varqc_idsA, high_varqc_idsB)
high_qc_ids <- intersect(high_qc_ids, high_varqc_ids)

anno <- data.table(read.csv( sprintf("%s/work_in_progress/Annotation_list_long_vnikos_1_9_21.csv", datadir)))
#samples_to_use <- anno[!LookRNAseq.Exp %in% exclude]$WTA.plate
columns <- colnames(adt)[-c(1:4)]
#samples_to_use <- c(intersect(columns, samples_to_use))
samples_to_use <- c(intersect(columns, anno$WTA.plate))
high_qc_ids <- intersect(high_qc_ids, samples_to_use)

#columns annotations used in file naming convention for visuals
col_anno <- data.table(read.csv( sprintf("%s/work_in_progress/Annotation_list_long_vnikos_1_9_21.csv", datadir)))
dim(col_anno)
excluded_by_qc <- col_anno$WTA.plate[-which(col_anno$WTA.plate %in% high_qc_ids)]
col_anno <- col_anno[ WTA.plate %in% high_qc_ids]
dim(col_anno)
col_anno[Pairs == "NA", Pairs := NA]

#fraction of nonzero tpm (binned) object
nonzeros.zs <- readRDS(sprintf("%s/aggregated_results/nonzeros.zs.bin50.rds", dirpath))

#Gene and cell annotaions
geneRanges <- readRDS(sprintf("%s/geneRanges.rds", datadir))
controlSampleIDs <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))
controlSampleIDs2 <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs2.rds", dirpath))
controlIDs <- readRDS(file = sprintf("%s/aggregated_results/reduced_controlIDs.rds", dirpath))
arms <- ganno <- data.table(readRDS(file = sprintf("%s/aggregated_results/ganno.rds", dirpath))) #called arms in visualization script
centromeres <- readRDS( sprintf("%s/centromeres.rds", datadir) )

#Why???
controlSampleIDs <- controlSampleIDs[-which(controlSampleIDs %in% "170425_A3")]
controlSampleIDs2 <- controlSampleIDs2[-which(controlSampleIDs2 %in% "170425_A3")]

##########################Data for OLR Plots###########################
message("Loading more data...")

#Main ML Predictions object
ourpreds <- readRDS(file = sprintf("%s/ML_data/preds.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES))
predicted_samples <- names(ourpreds$preds)
samples_to_use <- c(intersect(predicted_samples, samples_to_use))

#Auxiliary ML Data Objects
interstat.chr <- data.table(readRDS(file = sprintf("%s/ML_data/interstat.chrlevel.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES)))
interstat.arm <- data.table(readRDS(file = sprintf("%s/ML_data/interstat.armlevel.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES)))
fracstat <- data.table(readRDS(file = sprintf("%s/ML_data/fracstat.chrlevel.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES)))

##########################Data for Barcharts###########################
message("Loading more data...")

#SSM, MA, and Allelic CN information
Ms.chr <- readRDS(file = sprintf("%s/CN_data/CN_predictions.bychr.rds", dirpath))

#total expression by chr binned OLR
OLR_preds_wchr <- cbind(ourpreds[["cns"]][,c(2)], ourpreds[["cns"]][,-c(1:4)])
OLR_preds_bychr <- OLR_preds_wchr %>% 
                     group_by(seqnames) %>%
                       summarise_all(mean, na.rm = TRUE)
chr.TE <- setnames(OLR_preds_bychr, old = "seqnames", new = "bin_id")
chr.TE <- as.data.frame(chr.TE, stringsAsFactors = FALSE)
chr.TE$bin_id <- as.character(chr.TE$bin_id)

##Allele expression (option B (raw aggs by gene)):
#tmp.chr.Af <- copy(Ms.chr$raw.A)
#setnames(tmp.chr.Af, old = "seqnames", new = "bin_id")
#tmp.chr.Bf <- copy(Ms.chr$raw.B)
#setnames(tmp.chr.Bf, old = "seqnames", new = "bin_id")
#tmp.chr.Af <- cbind(tmp.chr.Af[,1], tmp.chr.Af[,-1]/(tmp.chr.Af[, -1] + tmp.chr.Bf[, -1]))
#tmp.chr.Bf <- cbind(tmp.chr.Bf[,1], tmp.chr.Bf[,-1]/(tmp.chr.Af[, -1] + tmp.chr.Bf[, -1]))
#IDs <- intersect(colnames(chr.TE)[-1], colnames(tmp.chr.Af[-1]))
#samples_to_use <- c(intersect(samples_to_use, IDs))

#Allele expression (option B (raw aggs by gene)):
allele_frac_mat.bychr <- readRDS(file=sprintf("%s/aggregated_results/allele_frac_mat.bychr.rds", dirpath))
chr.Af <- copy(allele_frac_mat.bychr$A)
setnames(chr.Af, old = "seqnames", new = "bin_id")
chr.Af$bin_id <- as.character(paste0("chr", chr.Af$bin_id))
chr.Af$bin_id[23] <- "chrX" #Change chr23 to chrX
chr.Bf <- copy(allele_frac_mat.bychr$B)
setnames(chr.Bf, old = "seqnames", new = "bin_id")
chr.Bf$bin_id <- as.character(paste0("chr", chr.Bf$bin_id))
chr.Bf$bin_id[23] <- "chrX" #Change chr23 to chrX
IDs <- intersect(colnames(chr.TE)[-1], colnames(chr.Af[-1]))
samples_to_use <- c(intersect(samples_to_use, IDs))

#Change the loss category to reflect CN = 0 as well as CN = 1
rawA <- copy(Ms.chr$raw.A)
setnames(rawA, old = "seqnames", new = "bin_id")
rawB <- copy(Ms.chr$raw.B)
setnames(rawB, old = "seqnames", new = "bin_id")
indices <- data.frame(stringsAsFactors = FALSE)
for (i in 1:ncol(rawA[,-c(1)])) {
  for (j in 1:nrow(rawA[,-c(1)])) {
    if (!is.nan(unlist(rawA[,-c(1)][j,..i]))) {
      if (rawA[,-c(1)][j,..i] < .05 && rawB[,-c(1)][j,..i] < .05) {
        indices <- rbind(indices, data.frame(sample = colnames(rawA[,-c(1)])[i] , chr = rawA$bin_id[j], stringsAsFactors = FALSE))
      }
    }
  }
}
indices$chr <- as.character(indices$chr)

for (i in 1:nrow(indices)) {
  value <- as.numeric(chr.TE[which(chr.TE$bin_id == indices$chr[i]), indices$sample[i]])
  print(value)
  if (value < 1.10) {
    print("pass change parameter")
    print(chr.TE[[indices$sample[i]]][which(chr.TE$bin_id == indices$chr[i])])
    chr.TE[[indices$sample[i]]] = replace(chr.TE[[indices$sample[i]]], which(chr.TE$bin_id == indices$chr[i]), .1)
    #chr.TE[[indices$sample[i]]][which(chr.TE$bin_id == indices$chr[i])] <- .1 #:= .1
    print(chr.TE[[indices$sample[i]]][which(chr.TE$bin_id == indices$chr[i])])
    print("changed")
  }
}

#########################Data for Pval Plot###########################

#calculate pvals
#source(sprintf("%s/scripts/calculating_pvals_byarm_nikos.R", scriptsdir))
source(sprintf("%s/scripts/pval_grouped_samples.R", scriptsdir))
source(sprintf("%s/scripts/calculating_control_pvals_bychr_v3.R", scriptsdir))

#Add arm level annotations to ourpreds
preds <- data.table(ourpreds[["cns"]])
ganno2 <- ganno[,c(1,6:9)]
preds2 <- setkey(preds, seqnames, start, end, id)
ganno3 <- setkey(ganno2, seqnames, start, end, id)
preds3 <- merge(ganno3, preds2)
OLR_preds_warm <- cbind(preds3[,c("arm")], preds3[,-c(1:5)])

#total expression from OLR binned by arm
OLR_preds_byarm <- OLR_preds_warm %>% 
                     group_by(arm) %>%
                       summarise_all(mean, na.rm = TRUE)

#reorder OLR preds columns like pval matrix
#OLR_preds_byarm <- OLR_preds_byarm[,colnames(pval_matrix_control_byarm)]


########################Data for Cluster Plot##########################

#TPM object normalized with chr bins
adt.bychr <- readRDS(file=sprintf("%s/aggregated_results/adt.bychr.rds", dirpath))
#absolute difference between allele cnts per chr
abs_allele_diff_mat.bychr <- readRDS(file=sprintf("%s/aggregated_results/abs_allele_diff_mat.bychr.rds", dirpath))
#total SNP counts over both alleles taken after separate allele normalization. ~same as var_mat.bychr. 
var_mat_sep.bychr <- readRDS(file=sprintf("%s/aggregated_results/var_mat_sep.bychr.rds", dirpath))

#add in quadrature and scale by sqrt(2) so that the point (1,1) is 1 in aggregated dim
matrix_adt.bychr <- as.matrix(adt.bychr[,-c(1)]); matrix_var.bychr <- as.matrix(var_mat_sep.bychr[,-c(1)]);
matrix_CN_Ratio <- sqrt( ( (matrix_adt.bychr^2) + (matrix_var.bychr^2) ) ) / sqrt(2)
CN_ratio.bychr <- data.table(cbind(adt.bychr[,c(1)], matrix_CN_Ratio))

message("Loading completed!")

#################################################################################
# Sanity check plots:
#################################################################################

#read in pvals
non_diploid_chrs_w_vals <- readRDS(file = sprintf("%s/aggregated_results/non_diploid_chromosomes.rds", dirpath))

#plot distribution of pvals across chromosomes
#tmp <- as.numeric(non_diploid_chrs_w_vals[which(non_diploid_chrs_w_vals[,c("sample")] %in% controlSampleIDs), c("chr")])
#hist(tmp, breaks=seq(min(tmp)-0.5, max(tmp)+0.5, by=1))
#length(test[,1])/length(unique(test[,c("sample")]))
#hist(as.numeric(non_diploid_chrs_w_vals[which(non_diploid_chrs_w_vals[,c("sample")] %in% unlist(unique(golden_samples[["normal"]][,c("sample_id")]))), c("chr")]), breaks=seq(min((as.numeric(non_diploid_chrs_w_vals[which(non_diploid_chrs_w_vals[,c("sample")] %in% unlist(unique(golden_samples[["normal"]][,c("sample_id")]))), c("chr")]))-0.5, max((as.numeric(non_diploid_chrs_w_vals[which(non_diploid_chrs_w_vals[,c("sample")] %in% unlist(unique(golden_samples[["normal"]][,c("sample_id")]))), c("chr")]))+0.5, by=1))))

#group and count number of genes in each chromosome
#adt_byarm <- adt[,c(1:2)] %>% 
#               group_by(seqnames) %>%
#                 summarise_all(length)

#Barplot of genes per chr
#barplot(height = adt_byarm$id, names = adt_byarm$seqnames, cex.names = .75, main = "Number of Genes per chr", ylab = "Number of Genes", xlab = "Chromosomes")


#################################################################################
# plots:
#################################################################################

require(doParallel)
registerDoParallel(25)

#Pull information for a given sample from ML prediction object
get_flat_format_tbl <- function(mysampid) {
  probs <- data.table(sample_id = mysampid, cbind(ourpreds$rowdata,  ourpreds$preds[[mysampid]]))
  setkeyv(probs, c("seqnames", "start", "end"))
  setkeyv(arms, c("seqnames", "start", "end"))
  tmp <- foverlaps(probs, arms)
  tmp2 <- tmp[, c("sample_id", "seqnames", "arm", "i.start", "i.end", "id", "1", "2", "3"), with=F]
  setnames(x = tmp2, old = c("i.start", "i.end", "1", "2", "3"), new = c("start", "end", "c1", "c2", "c3"))
  return(tmp2)
}

#single cpu plot function for a sample
plot_pdf <- function(myid = "170223_A5a", pdfFile = "./", chrs = chrs) {
  pdf(file = pdfFile, width = 10, height = 14)
  message("Doing ", pdfFile)
  #Get probs for specific sample
  probs <- get_flat_format_tbl(myid)
  #print table of pvals in pdf
  if (length(non_diploid_chrs_w_vals[which(non_diploid_chrs_w_vals[,1] == myid),c("sample")]) >= 1 ) {
    grid.table(non_diploid_chrs_w_vals[which(non_diploid_chrs_w_vals[,1] == myid),])
  }
  #Plot Bar plot based of allelic and OLR information
  plot_barplots_of_AllelicAndExpBiasPerSamples(dfs = chr.TE, wt = rep(1:nrow(chr.TE)), 
                                               fracs = list(Af = chr.Af, Bf = chr.Bf), nogeno = c(),
                                               ids = c(myid), chr = "")
  colnames(chr.TE[,order()])
  for(chr in chrs) {
    #Chr information
    message(chr)
    print(dim(na.omit(probs[seqnames == chr])))
    #Plot Machine learning visual function
    plot_cn_and_intermediate_prediction(probs = na.omit(probs[seqnames == chr]), controlIDs = controlIDs, plot_pvals = F)
    #Plot raw data visuals function
    plot_raw_data_and_prediction_boxplots2(myid = myid, chr = chr, adt = adt,
                                           nonzeros.zs = nonzeros.zs, coding = coding,
                                           preds = adt_fake, controlSampleIDs = controlSampleIDs)
    #Plot p-vals, tpm ratio, and OLR prediction boxplots for both arms
    plot_pvals_and_tpm_distributions(myid=myid, chr=chr, destDir=destDir)
    
  }
  dev.off()
}

#Variables
destDir <- sprintf("%s/visual_results", dirpath)
system(sprintf("mkdir -p %s", destDir))
#chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr11", "chr12", "chr16", "chr17", "chr19")
chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
#arms_to_exclude <- c("10q", "13q", 18", "20", "21", "22", "X")
# 14, 15 have had no p arm values and the rest have intentionally excluded arms 

#Get the samples that can be visualized
samples_to_visualize <- intersect(col_anno$WTA.plate, samples_to_use)

#make plots for every sample in anno list "col_anno"
#This is experimental parallel computation
foreach( i = c( 1:length(samples_to_visualize) ) ) %dopar% {
#foreach( i = c( 1:10 ) ) %dopar% {
    
#for (i in 1:length(samples_to_visualize)) {}
#source(sprintf("%s/plots/plot_raw_data_and_prediction_boxplots_w_2nd_allelic.R", scriptsdir))
#for (i in 1:2) {
  
  myid <- col_anno$WTA.plate[col_anno$WTA.plate %in% samples_to_use][i]
  #myid <- "190628_3D"
  pairid <- col_anno$Pairs[col_anno$Pairs %in% samples_to_use][i]
  if(!is.na(pairid)) {
    pdfFile <- sprintf("%s/%s.pair_%s.All_Plots.pdf", destDir, myid, pairid)
  } else {
    pdfFile <- sprintf("%s/%s.All_Plots.pdf", destDir, myid)
  }
  myid <- as.character(myid)
  plot_pdf(myid = myid, pdfFile = pdfFile, chrs = chrs)
  
}

#Make barcharts by family
setkey(anno, WTA.plate)
families <- sort(table(anno[samples_to_visualize][Pairs != "NA" & Pairs != "mother"]$Pairs), decreasing = T)
#exclude_chrs <- c("chr13", "chr18", "chr21", "chr22", "chrX") #for cluster plot, first instated in "normalize_tpm_and_vars_byfamily.R"
exclude_chrs <- c()

make_path <- sprintf("mkdir -p %s/byfamily/", destDir)
system(make_path)

foreach(i = c(1:length(names(families)))) %dopar% {
#for(i in 3:4) {
  myfamily = names(families)[i]
  message(myfamily)
  myids <- anno[Pairs %in% myfamily]$WTA.plate
  pdf(file = sprintf("%s/byfamily/%s.family.pdf", destDir, myfamily), width = 18, height = 12)
  #plot Barcharts by family
  p <- plot_barplots_of_AllelicAndExpBiasPerSamples(dfs = chr.TE, 
                                                    wt = rep(1:nrow(chr.TE)), fracs = list(Af = chr.Af, Bf = chr.Bf), nogeno = c(),
                                                    ids = myids, chr = "")
  #Plot cluster map output by family
  ClusterMap3D_byfamily(adt.bychr, var_mat_sep.bychr, abs_allele_diff_mat.bychr, myids, myfamily, exclude_chrs)
  #print(p)
  dev.off()
}

#Just a Sanity check that the code ran to completion
print("Done with making visuals")
