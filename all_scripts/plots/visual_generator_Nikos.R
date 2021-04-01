#################################################################################
# Passed arguments:
#################################################################################
#args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/data/processed_bam", "/pellmanlab/nikos/Stam_Etai_Data")
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
source(sprintf("%s/plots/class_prob_to_cn_level.R", scriptsdir))
source(sprintf("%s/plots/plot_cn_prediction_probabilities.R", scriptsdir))
source(sprintf("%s/plots/calc_and_plot_intermediate_pvalue.R", scriptsdir))
source(sprintf("%s/plots/plot_raw_data_and_prediction_boxplots2.R", scriptsdir))
source(sprintf("%s/plots/plot_barplots_of_AllelicAndExpBiasPerSamples.R", scriptsdir))
source(sprintf('%s/plots/plot_cn_and_intermediate_prediction_latest_model_intermediate_bars.R', scriptsdir))

#############################Parameters##################################

#Sample families to exclude
exclude <- c("113015_T12", "061516_T12", "070516_LookRNAseq", "072516_LookRNAseq", 
             "071816_LookRNAseq", "071816_LookRNAseq", "161005_LookRNAseq", 
             "061316_LookRNAseq", "161005_LookRNAseq_control", "Jinyu", "170726_MCF10A_control",
             "040416_ATCC", "161017_LookRNAseq", "062716_LookRNAseq", "080816LookRNAseq", "080816_LookRNAseq",
             "071416_LookRNAseq", "161003_LookRNAseq")

#Window size used in models
NUM_OF_FEATURES <- 50 

##########################Data for Raw Plots###########################
message("Loading data...")

#log-space centered tpm object
adt <- readRDS(sprintf("%s/aggregated_results/adt.rds", dirpath))
adt.na <- readRDS(sprintf("%s/aggregated_results/adt.na.rds", dirpath))
adt <- adt.default <- adt[order(seqnames, start, end)]
adt <- cbind( adt[,c(1:4)], setcolorder(adt[,-c(1:4)], order(colnames(adt[,-c(1:4)]))) )

#TMP: Create a fake object to feed as ML preds. I COULD JUST FEED ADT AGAIN.
adt_fake <- adt[,-c(1:4)]
adt_fake[,] <- 0
adt_fake <- cbind(adt[,c(1:4)], adt_fake)

#raw tpm object
rsemtpm <- readRDS(sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", dirpath))
high_qc_ids <- names(which(colSums(rsemtpm>5)>4000))

#columns annotations used in file naming convention for visuals
col_anno <- data.table(readRDS( sprintf("%s/work_in_progress/Annotation_list_long_vnikos_210129.rds", datadir) ))
dim(col_anno)
col_anno <- col_anno[ WTA.plate %in% high_qc_ids]
dim(col_anno)
col_anno[Pairs == "NA", Pairs := NA]

anno <- data.table(readRDS(sprintf("%s/work_in_progress/Annotation_list_long_vnikos_210129.rds", datadir) ))
samples_to_use <- anno[!LookRNAseq.Exp %in% exclude]$WTA.plate
columns <- colnames(adt)[-c(1:4)]
samples_to_use <- c(intersect(columns, samples_to_use))

#Variant matrix in the coding and UTR regions
coding <- readRDS(sprintf("%s/aggregated_results/ASE.coding.rds", dirpath))

#fraction of nonzero tpm (binned) object
nonzeros.zs <- readRDS(sprintf("%s/aggregated_results/nonzeros.zs.bin50.rds", dirpath))

#Gene and cell annotaions
geneRanges <- readRDS(sprintf("%s/aggregated_results/geneRanges.rds", dirpath))
controlSampleIDs <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))
controlSampleIDs2 <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs2.rds", dirpath))
controlIDs <- readRDS(file = sprintf("%s/aggregated_results/reduced_controlIDs.rds", dirpath))
arms <- ganno <- data.table(readRDS(file = sprintf("%s/aggregated_results/ganno.rds", dirpath))) #called arms in visualization script
centromeres <- readRDS( sprintf("%s/centromeres.rds", datadir) )

#Why???
controlSampleIDs <- controlSampleIDs[-which(controlSampleIDs %in% "170425_A3")]
controlSampleIDs2 <- controlSampleIDs2[-which(controlSampleIDs2 %in% "170425_A3")]

#triIDs <- col_anno[control_group %in% c("Tri12_LCM", "Tri12_facs")]$WTA.plate[col_anno[control_group %in% c("Tri12_LCM", "Tri12_facs")]$WTA.plate %in% names(which(colSums(rsemtpm>0)>6000))]

#We can resolve this later with the annotation stuff
#anno1 <- data.table(readRDS("/homes10/ejacob/WORK/secondaryanalysis/Stamatis/reports/Stamatis_list_v12_180404.QCRNAv3.annoqc.rds"))
#anno2 <- data.table(read_excel("/singlecellcenter/etai/ExperimentsData/Stamatis/May072018/Stamatis_list_v12_180518.xlsx"))
#anno3 <- merge(anno1, anno2[, c("WTA.plate", "Pairs", "control_group")], by = "WTA.plate")
#anno3$Pairs <- anno3$Pairs.y
#anno3$control_group <- anno3$control_group.y
#anno3$control_group.x <- anno3$control_group.y <- NULL
#anno3$Pairs.y <- anno3$Pairs.x <- NULL
#tri12ids <- anno3[control_group == "Tri12_facs" & th1 > 6000]$WTA.plate
#tri8ids <- anno3[control_group == "Tri8_facs" & th1 > 6000]$WTA.plate
#tri21ids <- anno3[control_group == "Tri21_facs" & th1 > 6000]$WTA.plate
#tris <- c(tri12ids, tri21ids, tri8ids)

##########################Data for OLR Plots###########################
message("Loading more data...")

#Main ML Predictions object
ourpreds <- readRDS(file = sprintf("%s/ML_data/preds.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES))
predicted_samples <- names(ourpreds$preds)
samples_to_use <- c(intersect(predicted_samples, samples_to_use))

#for (i in 1:length(samples_to_use)) {
#  myid <- samples_to_use[10]
#  probs <- ourpreds[["preds"]][[myid]]
#  preds <- round(probs)
#  for (i in 1:ncol(preds)) {preds[,i] <- preds[,i]*i}
#  CN <- rowSums(preds)
#  avg_CN <- mean(CN)
#  entry <- cbind(myid, avg_CN)
#  OLR_preds <- rbind(OLR_preds, entry)
#}

#Auxiliary ML Data Objects
interstat.chr <- data.table(readRDS(file = sprintf("%s/ML_data/interstat.chrlevel.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES)))
interstat.arm <- data.table(readRDS(file = sprintf("%s/ML_data/interstat.armlevel.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES)))
fracstat <- data.table(readRDS(file = sprintf("%s/ML_data/fracstat.chrlevel.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES)))

##########################Data for Boxcharts###########################
message("Loading more data...")

#SSM, MA, and Allelic CN information
Ms.chr <- readRDS(file = sprintf("%s/CN_data/CN_predictions.bychr.rds", dirpath))
#dirpath <- "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2/"
#Ms.chr <- readRDS(file = sprintf("%s/CN_predictions.bychr.rds", dirpath))

#Total expression SSM:
#chr.TE <- copy(Ms.chr$SSM.TE)
#setnames(chr.TE, old = "seqnames", new = "bin_id")
#chr.TE <- cbind(chr.TE[,c("bin_id")], 2^(chr.TE[,-c("bin_id")]))
#chr.TE <- as.data.frame(chr.TE)

#total expression by OLR
OLR_preds <- cbind(ourpreds[["cns"]][,c(2)], ourpreds[["cns"]][,-c(1:4)])
OLR_preds_bychr <- OLR_preds %>% 
  group_by(seqnames) %>%
  summarise_all(mean, na.rm = TRUE)
chr.TE <- setnames(OLR_preds_bychr, old = "seqnames", new = "bin_id")
chr.TE <- as.data.frame(chr.TE)

#Total expression from Normalized TPM
#chr.TE <- as.data.frame(readRDS(sprintf("%s/aggregated_results/normalized_rsemtpm_bychr.rds", dirpath)))
#colnames(chr.TE)[1] <- "bin_id"

#Total Expression LGBM
#chr.TE <- copy(Ms.chr$ML)
#setnames(chr.TE, old = "seqnames", new = "bin_id")
#chr.TE <- cbind(chr.TE[,c("bin_id")], (chr.TE[,-c("bin_id")]))
#chr.TE <- as.data.frame(chr.TE)

#Allele expression (option B (raw aggs by gene)):
chr.Af <- copy(Ms.chr$raw.A)
setnames(chr.Af, old = "seqnames", new = "bin_id")
chr.Bf <- copy(Ms.chr$raw.B)
setnames(chr.Bf, old = "seqnames", new = "bin_id")
chr.Af <- cbind(chr.Af[,1], chr.Af[,-1]/(chr.Af[, -1] + chr.Bf[, -1]))
chr.Bf <- cbind(chr.Bf[,1], chr.Bf[,-1]/(chr.Af[, -1] + chr.Bf[, -1]))
IDs <- intersect(colnames(chr.TE)[-1], colnames(chr.Af[-1]))
samples_to_use <- c(intersect(samples_to_use, IDs))

#######################################################################

message("Loading completed!")

#################################################################################
# plots:
#################################################################################

require(doParallel)
registerDoParallel(20)

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
  #Plot Bar plot based of allelic and SSM information
  plot_barplots_of_AllelicAndExpBiasPerSamples(dfs = chr.TE, wt = rep(1:nrow(chr.TE)), 
                                               fracs = list(Af = chr.Af, Bf = chr.Bf), nogeno = c(),
                                               ids = c(myid), chr = "")
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
  }
  dev.off()
}

#Variables
destDir <- sprintf("%s/visual_results", dirpath)
system(sprintf("mkdir -p %s", destDir))
chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr11", "chr12", "chr16", "chr17", "chr19")
# 14, 15 have had no p arm values and the rest have intentionally excluded arms 

#Get the samples that can be visualized
samples_to_visualize <- intersect(col_anno$WTA.plate, samples_to_use)

#make plots for every sample in anno list "col_anno"
#This is experimental parallel computation
foreach( i = c( 1:length(samples_to_visualize) ) ) %dopar% {
#foreach( i = c( 1:10 ) ) %dopar% {
    
#for (i in 1:length(samples_to_visualize)) {}
#source(sprintf("%s/plots/plot_raw_data_and_prediction_boxplots2.R", scriptsdir))
#for (i in 1:2) {
  
  myid <- col_anno$WTA.plate[col_anno$WTA.plate %in% samples_to_use][i]
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

make_path <- sprintf("mkdir -p %s/byfamily/", destDir)
system(make_path)

foreach(i = c(1:length(names(families)))) %dopar% {
#for(i in 3:4) {
  myfamily = names(families)[i]
  message(myfamily)
  myids <- anno[Pairs %in% myfamily]$WTA.plate
  pdf(file = sprintf("%s/byfamily/%s.barplots.pdf", destDir, myfamily), width = 18, height = 12)
  p <- plot_barplots_of_AllelicAndExpBiasPerSamples(dfs = chr.TE, 
                                                    wt = rep(1:nrow(chr.TE)), fracs = list(Af = chr.Af, Bf = chr.Bf), nogeno = c(),
                                                    ids = myids, chr = "")
  #print(p)
  dev.off()
}
#Just a Sanity check that the code ran to completion
print("Done with making visuals")
