#mostly relevant for plots of Sep192019

args <- commandArgs(trailingOnly = TRUE)



if(length(args) < 3) {
  stop("Wrong input.\nRscript run_prediction_for_sample.R sample_id secondary_id dstdir")
}

mysample_id <- args[1]
secondary_id <- args[2]
dstdir <- args[3]

dirpath <- "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2"

require(data.table)


message("Loading data..")
adt.na <- readRDS(file = sprintf("%s/data/adt.na.rds", dirpath))
adt <- readRDS(file = sprintf("%s/data/adt.rds", dirpath))
nonzeros <- readRDS(file = sprintf("%s/data/nonzeros.rds", dirpath))
geneSpecific <- readRDS(file = sprintf("%s/data/geneSpecific.rds", dirpath))
ganno <- readRDS(file = sprintf("%s/data/ganno.rds", dirpath))

message("Loading data for plots..")
controlSampleIDs2 <- readRDS(sprintf("%s/data/controlSampleIDs2.rds", dirpath))
controlSampleIDs <- readRDS(sprintf("%s/data/controlSampleIDs.rds", dirpath))
geneSpecific <- readRDS(file = sprintf("%s/data/geneSpecific.rds", dirpath))
rsemtpm <- readRDS(file = sprintf("%s/data/rsemtpm.rds", dirpath))
#ASE <- readRDS(file = sprintf("%s/data/ASE.rds", dirpath))
coding <- readRDS(file = sprintf("%s/data/ASE.coding.rds", dirpath))

configs <- readRDS(sprintf("%s/data/param_config_list.rds", dirpath))
geneRanges <- readRDS(sprintf("%s/data/geneRanges.rds", dirpath))
#adt.old <- readRDS("/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/adt.na.rds")
high_qc_ids <- names(which(colSums(rsemtpm>5)>5000))
nonzeros.zs <- readRDS(file = sprintf("%s/data/nonzeros.zs.bin50.rds", dirpath))

col_anno <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/col_anno.Apr122018.rds")
dim(col_anno)
col_anno <- col_anno[ WTA.plate %in% high_qc_ids]
dim(col_anno)
col_anno[Pairs == "NA", Pairs := NA]

controlIDs <- col_anno[experimental_label %in% c("control_LCM_notreatment", "control_gen2_notreatment", "noco_gen2_noMN_lcm", "noco_gen1_noMN", "p53_gen2_noMN")]$WTA.plate
triIDs <- col_anno[experimental_label %in% c("Tri12_LCM", "Tri12_facs")]$WTA.plate[col_anno[experimental_label %in% c("Tri12_LCM", "Tri12_facs")]$WTA.plate %in% names(which(colSums(rsemtpm>0)>6000))]

interstat.arm <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/interstat.armlevel.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")
interstat.chr <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/interstat.chrlevel.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")
fracstat <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/fracstat.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")
centromeres <- readRDS("~ejacob/WORK/secondaryanalysis/Stamatis/data/centromeres.rds")

source('~/WORK/Papers/MLpaper/R/utilities/calc_and_plot_intermediate_pvalue.R')
source('~/WORK/Papers/MLpaper/R/utilities/plot_cn_and_intermediate_prediction.R')
source('~/WORK/Papers/MLpaper/R/utilities/plot_raw_data_and_prediction_boxplots2.R')
library(matrixStats)
require(gridExtra)

par(mfrow=c(1,1))

manual_selected_vars6 <- c("geneSpecific.MEDIAN.MEDIAN.GC",
                           "geneSpecific.MEDIAN.MEDIAN.cv.nz",
                           "geneSpecific.MEDIAN.MEDIAN.Length",
                           "geneSpecific.MEAN.MEAN.interspace",
                           "TE.all.SD",
                           "TE.onlyExpressed.MAXSTRETCH",
                           "Frac.nonzeros.MAXSTRETCH",
                           "TE.all.Q95",
                           "TE.onlyExpressed.Q95",
                           "Frac.nonzeros.SD",
                           "TE.onlyExpressed.Q25",
                           "TE.onlyExpressed.Q75",
                           "TE.all.Q05",
                           "Frac.nonzeros.MEAN",
                           "TE.onlyExpressed.MEDIAN",
                           "TE.onlyExpressed.MEAN",
                           "TE.all.Q25",
                           "TE.all.Q75",
                           "TE.all.MEDIAN",
                           "TE.all.MEAN")


#test:
#load db:
#dball <- readRDS("/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/tmp/db.all_samples_ftrs.zscored.winSize50.rds")

#chr <- "chr12"
#mysample_id = "190705_2A" #"190701_3D" #jinyu "190528_eHAP_B3" #"190627_1A" #ctrl facs "190528_RPE1_D2" #"190705_6E" #trilcm"190705_5D" #trilcm 190705_2A" #ctrl "190429_3B" #controlSampleIDs[11] #"170512_B2"

message("Calculating features for predictions..")
source('~/WORK/Papers/MLpaper/R/utilities/generate_ftrs_for_input.R')
tt <- generate_ftrs_for_input(adt = adt, adt.na = adt.na, nonzeros = nonzeros, ganno = ganno,
                              sample_id = mysample_id)
dim(tt)
m <- as.matrix((tt[,manual_selected_vars6, with=F]))
source('~/WORK/Papers/MLpaper/R/utilities/predict_using_lightgbm.R')
tmp.preds <- cbind(tt[,1:6], predict_using_lightgbm(m))
dim(tmp.preds)
#chr <- "chr12"

par(mfrow=c(1,1))
source('~/WORK/Papers/MLpaper/R/utilities/calc_and_plot_intermediate_pvalue.R') 
source('~/WORK/Papers/MLpaper/R/utilities/plot_cn_and_intermediate_prediction.R')

plot_pdf <- function(myid, dstdir) {
  pdfFile <- sprintf("%s/%s.%s.CNinterPlots.CNV12.2.pdf", dstdir, myid, secondary_id)
  message("Generating plots in file: ", pdfFile)
  pdf(file = pdfFile, width = 10, height = 14)
  message("Doing ", pdfFile)
  mychrs <- as.character(unique(tmp.preds$seqnames))
  for(chr in mychrs) {
    message(chr)
    tmp <- plot_cn_and_intermediate_prediction(tmp.preds[seqnames == chr], plot_pvals = F)
    tmp <- plot_raw_data_and_prediction_boxplots2(myid = mysample_id, chr = chr, 
                                                  adt = adt, nonzeros.zs = nonzeros.zs, coding = coding, 
                                                  preds = tmp.preds, controlSampleIDs = controlSampleIDs)
    grid.arrange(grobs = c(tmp), ncol=1, as.table=F)
  }
  dev.off()
}

plot_pdf(myid = mysample_id, dstdir = dstdir)




#TODO:
#add raw data boxplots

#old:
# tt2 <- (dball[sample_id == mysample_id, colnames(tt), with=F])
# dim(tt2)
# m2 <- as.matrix((tt2[,manual_selected_vars6, with=F]))
# tmp2.preds <- cbind(tt2[,1:6], predict_using_lightgbm(m2))
# dim(tmp2.preds)
# 
# sum(tt2$id != tt$id)
# plot(tmp.preds[[7]], tmp2.preds[[7]])
# abline(0,1, col="red", lwd=2)
# plot(tmp.preds[[8]], tmp2.preds[[8]])
# abline(0,1, col="red", lwd=2)
# plot(tmp.preds[[9]], tmp2.preds[[9]])
# abline(0,1, col="red", lwd=2)
# 


#tmp <- plot_cn_and_intermediate_prediction(tmp.preds[seqnames == chr], plot_pvals = F, main = "new")

#old
#tmp <- plot_cn_and_intermediate_prediction(tmp2.preds[seqnames == chr], plot_pvals = F, main = "old")    


