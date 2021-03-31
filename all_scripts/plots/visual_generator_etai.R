#################################################################################
# Passed arguments:
#################################################################################
#args <- c("/pellmanlab/nikos/Stam_Etai_Scripts", "/pellmanlab/stam_niko/data/processed_bam", "/pellmanlab/nikos/Stam_Etai_Data")
args <- commandArgs(trailingOnly = TRUE)
print(args)
scriptsdir <- args[1]
dirpath <- args[2]
dirpath <- sprintf("%s/test", dirpath)
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
source( sprintf("%s/plots/class_prob_to_cn_level.R", scriptsdir) )
source( sprintf("%s/scripts/CNML.R", scriptsdir) )
source( sprintf("%s/plots/plot_cn_prediction_probabilities.R", scriptsdir) )
source( sprintf("%s/plots/calc_and_plot_intermediate_pvalue.R", scriptsdir) )
source( sprintf("%s/plots/plot_raw_data_and_prediction_boxplots2.R", scriptsdir) )

#Sample families to exclude
exclude <- c("113015_T12", "061516_T12", "070516_LookRNAseq", "072516_LookRNAseq", 
             "071816_LookRNAseq", "071816_LookRNAseq", "161005_LookRNAseq", 
             "061316_LookRNAseq", "161005_LookRNAseq_control", "Jinyu", "170726_MCF10A_control",
             "040416_ATCC", "161017_LookRNAseq", "062716_LookRNAseq", "080816LookRNAseq", "080816_LookRNAseq",
             "071416_LookRNAseq", "161003_LookRNAseq")

NUM_OF_FEATURES <- 50 

###################
# START Etai data #
###################
#Must be run after etai version of Run_ML
adt_fake <- adt[,-c(1:4)]
adt_fake[,] <- 0
adt_fake <- cbind(adt[,c(1:4)], adt_fake)
arms <- ganno
#################
# END Etai data #
#################

#columns annotations used in file naming convention for visuals
#col_anno <- data.table(read.csv( sprintf("%s/work_in_progress/Annotation_list_long_vnikos_210129.csv", datadir) ))
col_anno <- data.table(read_xlsx(path = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2/data/samples_annotation.xlsx"))
dim(col_anno)
col_anno <- col_anno[ WTA.plate %in% high_qc_ids]
dim(col_anno)
col_anno[Pairs == "NA", Pairs := NA]

#anno <- data.table(read.csv(sprintf("%s/work_in_progress/Annotation_list_long_vnikos_210129.csv", datadir) ))
#samples_to_use <- anno[!LookRNAseq.Exp %in% exclude]$WTA.plate
#columns <- colnames(adt)[-c(1:4)]
#samples_to_use <- c(intersect(columns, samples_to_use))
anno <- data.table(read_xlsx(path = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2/data/samples_annotation.xlsx"))
samples_to_use <- anno[!LookRNAseq.Exp %in% exclude]$WTA.plate

#Why???
controlSampleIDs <- controlSampleIDs[-which(controlSampleIDs %in% "170425_A3")]
controlSampleIDs2 <- controlSampleIDs2[-which(controlSampleIDs2 %in% "170425_A3")]

controlIDs <- readRDS(file = sprintf("%s/aggregated_results/reduced_controlIDs.rds", dirpath))

triIDs <- col_anno[control_group %in% c("Tri12_LCM", "Tri12_facs")]$WTA.plate[col_anno[control_group %in% c("Tri12_LCM", "Tri12_facs")]$WTA.plate %in% names(which(colSums(rsemtpm>0)>6000))]

#We can resolve this later with the annotation stuff
anno1 <- data.table(readRDS("/homes10/ejacob/WORK/secondaryanalysis/Stamatis/reports/Stamatis_list_v12_180404.QCRNAv3.annoqc.rds"))
anno2 <- data.table(read_excel("/singlecellcenter/etai/ExperimentsData/Stamatis/May072018/Stamatis_list_v12_180518.xlsx"))
anno3 <- merge(anno1, anno2[, c("WTA.plate", "Pairs", "control_group")], by = "WTA.plate")
anno3$Pairs <- anno3$Pairs.y
anno3$control_group <- anno3$control_group.y
anno3$control_group.x <- anno3$control_group.y <- NULL
anno3$Pairs.y <- anno3$Pairs.x <- NULL

tri12ids <- anno3[control_group == "Tri12_facs" & th1 > 6000]$WTA.plate
tri8ids <- anno3[control_group == "Tri8_facs" & th1 > 6000]$WTA.plate
tri21ids <- anno3[control_group == "Tri21_facs" & th1 > 6000]$WTA.plate
tris <- c(tri12ids, tri21ids, tri8ids)

#Main ML Predictions object
ourpreds <- readRDS(file = sprintf("%s/ML_data/preds.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES))
predicted_samples <- names(ourpreds$preds)
samples_to_use <- c(intersect(predicted_samples, samples_to_use))

#Auxiliary ML Data Objects
#interstat.chr <- data.table(readRDS(file = sprintf("%s/ML_data/interstat.chrlevel.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES)))
#interstat.arm <- data.table(readRDS(file = sprintf("%s/ML_data/interstat.armlevel.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES)))
#fracstat <- data.table(readRDS(file = sprintf("%s/ML_data/fracstat.chrlevel.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES)))
interstat.chr <- data.table(readRDS(file = sprintf("%s/latest_model/interstat.chrlevel.probs_olr_model.WIN%d.rds", dstpath, NUM_OF_FEATURES)))
interstat.arm <- data.table(readRDS(file = sprintf("%s/latest_model/interstat.armlevel.probs_olr_model.WIN%d.rds", dstpath, NUM_OF_FEATURES)))
fracstat <- data.table(readRDS(file = sprintf("%s/latest_model/fracstat.chrlevel.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES)))
centromeres <- data.table(readRDS("/homes10/ejacob/WORK/secondaryanalysis/Stamatis/data/centromeres.rds"))

message("Loading completed!")

#################################################################################
# plots:
#################################################################################

#Load ML plot scripts
source( sprintf('%s/plots/plot_cn_and_intermediate_prediction_latest_model.R', scriptsdir) )
#source( sprintf('%s/plots/plot_cn_and_intermediate_prediction_latest_model_intermediate_bars.R', scriptsdir) )

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
  probs <- get_flat_format_tbl(myid)
  for(chr in chrs) {
    
    #Chr information
    message(chr)
    print(dim(na.omit(probs[seqnames == chr])))
    
    #Plot Machine learning visual function
    plot_cn_and_intermediate_prediction(probs = na.omit(probs[seqnames == chr]), controlIDs = controlIDs, plot_pvals = F)
    
    #Plot raw data visuals function
    #plot_raw_data_and_prediction_boxplots2(myid = myid, chr = chr, adt = adt,
    #                                              nonzeros.zs = nonzeros.zs, coding = coding,
    #                                              preds = adt_fake, controlSampleIDs = controlSampleIDs)
    #
  
  }
  dev.off()
}

#Variables
destDir <- sprintf("%s/visual_results_alldata3", dstpath)
system(sprintf("mkdir -p %s", destDir))
chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr11", "chr12", "chr16", "chr17", "chr19")
# 14, 15 have had no p arm values and the rest have intentionally excluded arms 

#Get the samples that can be visualized
samples_to_visualize <- intersect(col_anno$WTA.plate, samples_to_use)

samples_to_visualize <- c("190429_1D", "190701_3A", "190628_4B", "190702_6D", "190701_9F", "190627_3G")
#make plots for every sample in anno list "col_anno"
#for(i in 1:length(samples_to_visualize)) {

for(i in 1:length(samples_to_visualize)) {
#for(i in 1:2) {

  myid <- col_anno$WTA.plate[col_anno$WTA.plate %in% samples_to_visualize][i]
  pairid <- col_anno$Pairs[col_anno$Pairs %in% samples_to_visualize][i]
  if(!is.na(pairid)) {
    pdfFile <- sprintf("%s/%s.pair_%s.All_Plots.pdf", destDir, myid, pairid)
  } else {
    pdfFile <- sprintf("%s/%s.All_Plots.pdf", destDir, myid)
  }
  myid <- as.character(myid)
  plot_pdf(myid = myid, pdfFile = pdfFile, chrs = chrs)
  
}

#Just a Sanity check that the code ran to completion
print("Done with making visuals")

