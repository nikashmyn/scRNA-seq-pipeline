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

#CN <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/tmp/CN.LGBM_0.05.mono_allele_th_0.25_preds.rds")

#reinc filtering:
#require(readxl)
# reinc <- data.table(read_excel("/singlecellcenter/etai/ExperimentsData/Stamatis/Feb_13_2018_reincorporation_list_St.xlsx"))
# 
# reinc <- reinc[ID %in% high_qc_ids,]
# reinc <- reinc[!pairing %in% names(which(table(reinc$pairing) < 2))]
# anno_row <- as.data.frame(reinc[ ,1:2])
# rownames(anno_row) <- anno_row$ID
# anno_row$ID <- NULL
# anno_row$relashionship <- reinc$`family relationship (imaging)`
# colnames(anno_row)[1] <- c("family")
# anno_row$relashionship <- factor(anno_row$relashionship)
# 
# reinc <- data.table(read_excel("/singlecellcenter/etai/ExperimentsData/Stamatis/Feb_13_2018_reincorporation_list_St.xlsx"))

message("Loading data...")
coding <- readRDS("/pellmanlab/stam_niko/etai_code/analysis_result_data/ASE.coding.rds")
noncoding <- readRDS("/pellmanlab/stam_niko/etai_code/analysis_result_data/ASE.noncoding.rds")

#nonzeros <- readRDS("/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/nonzeros.bin50.rds")
nonzeros.zs <- readRDS("/pellmanlab/stam_niko/etai_code/analysis_result_data/nonzeros.zs.bin50.rds")

load("/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.geneRanges.RData")
#names(geneRanges)[1] <- "id"

#centromeres <- readRDS("/pellmanlab/stam_niko/etai_code/analysis_result_data/centromeres.rds")

#MLREG <- readRDS("/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/LGBM.regression.rds")

require(readxl)
anno <- data.table(read_excel("/singlecellcenter/etai/ExperimentsData/Stamatis/May072018/Stamatis_list_v12_180508.QCRNAv3.withQCRNAv3.xlsx"))
controlSampleIDs <- anno[ (mainGroup == "A_level_control" | mainGroup == "B_level_control") & th1 > 6000 ]$WTA.plate
controlSampleIDs <- controlSampleIDs[controlSampleIDs %in% colnames(rsemtpm)]

controlSampleIDs2 <- anno[ (mainGroup == "A_level_control") & th1 > 6000 ]$WTA.plate
controlSampleIDs2 <- controlSampleIDs2[controlSampleIDs2 %in% colnames(rsemtpm)]

controlSampleIDs <- controlSampleIDs[-which(controlSampleIDs %in% "170425_A3")]
controlSampleIDs2 <- controlSampleIDs2[-which(controlSampleIDs2 %in% "170425_A3")]

#col_anno <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/col_anno.Apr122018.rds")
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
#arms <- readRDS(file = "/pellmanlab/stam_niko/etai_code/analysis_result_data/arms.rds")
#alleles.all <- readRDS("/singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025e/alleles.all.rds")

#interstat.arm <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/interstat.armlevel.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")
#interstat.chr <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/interstat.chrlevel.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")

#lgbmpreds5 <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/tmp/lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")
#fracstat <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/fracstat.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")

#### Bar Plot Data inclusion -Nikos ###
#dirpath <- "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2/"
#anno <- data.table(read_xlsx(path = sprintf("%s/data/samples_annotation.xlsx", dirpath)))
#
#Ms.chr <- readRDS(file = sprintf("%s/CN_predictions.bychr.rds", dirpath))
#
##Total expression:
#chr.TE <- copy(Ms.chr$ML)
#setnames(chr.TE, old = "seqnames", new = "bin_id")
#chr.TE <- as.data.frame(chr.TE)
#
##Allele expression (option B (raw aggs by gene)):
#chr.Af <- copy(Ms.chr$raw.A)
#setnames(chr.Af, old = "seqnames", new = "bin_id")
#chr.Bf <- copy(Ms.chr$raw.B)
#setnames(chr.Bf, old = "seqnames", new = "bin_id")
#chr.Af <- cbind(chr.Af[,1], chr.Af[,-1]/(chr.Af[, -1] + chr.Bf[, -1]))
#chr.Bf <- cbind(chr.Bf[,1], chr.Bf[,-1]/(chr.Af[, -1] + chr.Bf[, -1]))

### End Bar Plot Data Inclusion ###

message("Loading completed!")
#################################################################################
# plots:
#################################################################################


#single cpu:
source('/pellmanlab/stam_niko/etai_code/DFCI.scRNAseq.workflows/plots/plot_raw_data_and_prediction_boxplots.R')
#source('/pellmanlab/stam_niko/etai_code/DFCI.scRNAseq.workflows/plots/plot_cn_and_intermediate_prediction.R')

library(matrixStats)


plot_pdf <- function(myid = "170223_A5a", pdfFile = "./", chrs = chrs) {
  #pdfFile <- sprintf("%s/%s.CNinterPlots.CNV12.2.pdf", destdir, myid)
  pdf(file = pdfFile, width = 10, height = 14)
  message("Doing ", pdfFile)
  #mychrs <- as.character(unique(lgbmpreds5$seqnames))[as.character(unique(lgbmpreds5$seqnames)) %in% chrs]
  for(chr in chrs) {
    message(chr)
    #tmp <- plot_cn_and_intermediate_prediction(lgbmpreds5[sample_id == myid & seqnames == chr], plot_pvals = F)
    tmp <- plot_raw_data_and_prediction_boxplots(myid = myid, chr = chr, 
                                                 adt = adt, nonzeros.zs = nonzeros.zs, coding = coding, 
                                                 MLREG = NA, controlSampleIDs = controlSampleIDs)
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
mychrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr5", "chr6", "chr7", "chr8", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")


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

#Just a Sanity check that the code ran to completion
system("touch DONE.txt", destDir)
