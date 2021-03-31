#################################################################################
# loading data:
#################################################################################

require(data.table)
require(gtools)
#source('~/WORK/secondaryanalysis/methods_paper_results/R/CNML.R')
#source('~/WORK/Papers/MLpaper/R/utilities/plot_cn_prediction_probabilities.R')

mlinput2.na <- readRDS("/singlecellcenter/etai/ExperimentsData/Stamatis/reports.CNV12.0/data/mlinput.with_na.rds")
adt <- mlinput2.na$adt[,-c(5,6)]

rsemtpm <- readRDS("/pellmanlab/stam_niko/etai_code/experiments_etai/Stamatis_list_v15_190910.rsemtpm.rds")
#high_qc_ids <- names(which(colSums(rsemtpm>5)>5000))
high_qc_ids <- names(which(colSums(rsemtpm>5)>4000))

col_anno_old <- readRDS("/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/col_anno.Apr122018.rds")

col_anno <- data.table(read_excel("/singlecellcenter/etai/ExperimentsData/Stamatis/Stamatis_list_v15_190910.xlsx"))
dim(col_anno)
col_anno <- col_anno[ WTA.plate %in% high_qc_ids]
dim(col_anno)
col_anno[Pairs == "NA", Pairs := NA]

message("Loading data...")
coding <- readRDS("/pellmanlab/stam_niko/etai_code/analysis_result_data/ASE.coding.rds")
noncoding <- readRDS("/pellmanlab/stam_niko/etai_code/analysis_result_data/ASE.noncoding.rds")

nonzeros.zs <- readRDS("/pellmanlab/stam_niko/etai_code/analysis_result_data/nonzeros.zs.bin50.rds")

load("/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.geneRanges.RData")
names(geneRanges)[1] <- "id"

centromeres <- readRDS("/pellmanlab/stam_niko/etai_code/analysis_result_data/centromeres.rds")

MLREG <- readRDS("/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/LGBM.regression.rds")

require(readxl)
anno <- data.table(read_excel("/singlecellcenter/etai/ExperimentsData/Stamatis/May072018/Stamatis_list_v12_180508.QCRNAv3.withQCRNAv3.xlsx"))
controlSampleIDs <- anno[ (mainGroup == "A_level_control" | mainGroup == "B_level_control") & th1 > 6000 ]$WTA.plate
controlSampleIDs <- controlSampleIDs[controlSampleIDs %in% colnames(rsemtpm)]

controlSampleIDs2 <- anno[ (mainGroup == "A_level_control") & th1 > 6000 ]$WTA.plate
controlSampleIDs2 <- controlSampleIDs2[controlSampleIDs2 %in% colnames(rsemtpm)]

controlSampleIDs <- controlSampleIDs[-which(controlSampleIDs %in% "170425_A3")]
controlSampleIDs2 <- controlSampleIDs2[-which(controlSampleIDs2 %in% "170425_A3")]

controlIDs <- col_anno[control_group %in% c("control_LCM_notreatment", "control_gen2_notreatment", "noco_gen2_noMN_lcm", "noco_gen1_noMN", "p53_gen2_noMN")]$WTA.plate

triIDs <- col_anno[control_group %in% c("Tri12_LCM", "Tri12_facs")]$WTA.plate[col_anno[control_group %in% c("Tri12_LCM", "Tri12_facs")]$WTA.plate %in% names(which(colSums(rsemtpm>0)>6000))]

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

message("Loading more data...")
anno <- data.table(read_excel("/singlecellcenter/etai/ExperimentsData/Stamatis/Stamatis_list_v14_181103.xlsx"))
arms <- readRDS(file = "/pellmanlab/stam_niko/etai_code/analysis_result_data/arms.rds")

interstat.arm <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/interstat.armlevel.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")
interstat.chr <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/interstat.chrlevel.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")

lgbmpreds5 <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/tmp/lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")
fracstat <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/fracstat.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")

message("Loading completed!")
#################################################################################
# plots:
#################################################################################


#single cpu:
source('/pellmanlab/stam_niko/etai_code/DFCI.scRNAseq.workflows/plots/plot_raw_data_and_prediction_boxplots.R')
source('/pellmanlab/stam_niko/etai_code/DFCI.scRNAseq.workflows/plots/plot_cn_and_intermediate_prediction.R')

library(matrixStats)


temp <- lgbmpreds5[sample_id == "190429_1C" & seqnames == "chr1"]

plot_pdf <- function(myid = "170223_A5a", pdfFile = "./", chrs = chrs) {
  #pdfFile <- sprintf("%s/%s.CNinterPlots.CNV12.2.pdf", destdir, myid)
  pdf(file = pdfFile, width = 10, height = 14)
  message("Doing ", pdfFile)
  mychrs <- as.character(unique(lgbmpreds5$seqnames))[as.character(unique(lgbmpreds5$seqnames)) %in% chrs]
  for(chr in mychrs) {
    message(chr)
    tmp <- plot_cn_and_intermediate_prediction(lgbmpreds5[sample_id == myid & seqnames == chr], plot_pvals = F)
    tmp <- plot_raw_data_and_prediction_boxplots(myid = myid, chr = chr, 
                                                 adt = adt, nonzeros.zs = nonzeros.zs, coding = coding, 
                                                 MLREG = MLREG, controlSampleIDs = controlSampleIDs)
    ###Nikos Code###
    #Should I only pass pair IDs to this expression
    #Get from barplots_and_heatmaps_for_CN_predictions.Oct302019.R
    #mp <- plot_barplots_of_AllelicAndExpBiasPerSamples(ids = myid, chr = chr, 
    #                                                   dfs = chr.TE, wt = rep(1:nrow(chr.TE)), 
    #                                                   fracs = list(Af = chr.Af, Bf = chr.Bf), ymax = NULL,
    #                                                   plotOnlyDepth = F, pval.th = 0.01, nogeno = c(), 
    #                                                   plotAxisText = T)
    ###Nikos Code###
    grid.arrange(grobs = c(tmp), ncol=1, as.table=F)
  }
  dev.off()
}

destDir <- "/pellmanlab/stam_niko/etai_code/experiments_etai/Final_Visualizations_includingf"
chrs <- as.character(unique(MLREG$seqnames))

for(i in 1:nrow(col_anno)) {
  
  myid <- col_anno$WTA.plate[i]
  pairid <- col_anno$Pairs[i]
  if(!is.na(pairid)) {
    pdfFile <- sprintf("%s/%s.pair_%s.CNinterPlots.CNV12.2.pdf", destDir, myid, pairid)
  } else {
    pdfFile <- sprintf("%s/%s.CNinterPlots.CNV12.2.pdf", destDir, myid)
  }
  plot_pdf(myid = myid, pdfFile = pdfFile, chrs = chrs)
  
}



### Experimental parallel approach ###
#TODO: Fix read-in data so that it includes f-experiment 
#TODO: Unserialize SNOW so that it doesnt parallelize the same task

myid.list <- col_anno$WTA.plate
pdfFile.list <- list()

for(i in 1:nrow(col_anno)) {
  myid <- col_anno$WTA.plate[i]
  pairid <- col_anno$Pairs[i]
  if(!is.na(pairid)) {
    pdfFile.list <- append(pdfFile.list, sprintf("%s/%s.pair_%s.CNinterPlots.CNV12.2.pdf", destDir, myid, pairid))
  } else {
    pdfFile.list <- append(pdfFile.list, sprintf("%s/%s.CNinterPlots.CNV12.2.pdf", destDir, myid))
  }
}

cpus <- 24
require(snow)
cl <- snow::makeCluster(cpus, type = "SOCK")

clusterEvalQ(cl, library(matrixStats))
clusterEvalQ(cl, source('/pellmanlab/stam_niko/etai_code/DFCI.scRNAseq.workflows/plots/plot_raw_data_and_prediction_boxplots.R'))
clusterEvalQ(cl, source('/pellmanlab/stam_niko/etai_code/DFCI.scRNAseq.workflows/plots/plot_cn_and_intermediate_prediction.R'))

snow::clusterExport(cl = cl, list = c("interstat.arm", "interstat.chr",  "controlIDs", "centromeres"))
snow::clusterExport(cl = cl, list = c("rsemtpm"))
snow::clusterExport(cl = cl, list = c("lgbmpreds5", "adt", "nonzeros.zs", "coding", "MLREG", "controlSampleIDs"))
snow::clusterExport(cl = cl, list = c("myid.list",  "pdfFile.list", "chrs",  "sample_id"))
snow::clusterExport(cl = cl, list = c("plot_pdf", "plot_cn_and_intermediate_prediction", "plot_raw_data_and_prediction_boxplots"))

parLapply(cl, myid.list, plot_pdf, pdfFile.list, chrs)

stopCluster(cl)

### Experimental parallel approach ###

#Just a Sanity check that the code ran to completion
system("touch DONE.txt", destDir)
