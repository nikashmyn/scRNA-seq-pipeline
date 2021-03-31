#running all type of predictions for a sample

args <- c("160627_B1", "", "/pellmanlab/stam_niko/data/processed_bam", "/pellmanlab/nikos/Stam_Etai_Data")
args <- commandArgs(trailingOnly = TRUE)


if(length(args) < 3) {
  stop("Wrong input.\nRscript run_prediction_for_sample.R sample_id secondary_id dstdir")
}

mysample_id <- args[1] #triIDs[2]  #"170512_B4"  #args[1]
secondary_id <- args[2] #"" #args[2]
dirpath <- args[3]
dstdir <- sprintf("%s/SSM_data", dirpath)
datapath <- args[4]

require(data.table)


message("Loading data..")
#Read in Data
adt.na <- readRDS(file = sprintf("%s/aggregated_results/adt.na.rds", dirpath))
adt <- readRDS(file = sprintf("%s/aggregated_results/adt.rds", dirpath))
nonzeros <- readRDS(file = sprintf("%s/aggregated_results/nonzeros.rds", dirpath))
alleles.all <- readRDS(file = sprintf("%s/aggregated_results/alleles.all.rds", dirpath)) # "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2/data/alleles.all.rds"
ganno <- readRDS(file = sprintf("%s/aggregated_results/ganno.rds", dirpath))

#Read in "given" objects. Need to figure out how to recreate these.
geneSpecific <- readRDS(file = sprintf("%s/geneSpecific.rds", datapath))
centromeres <- readRDS(file = sprintf("%s/centromeres.rds", datapath))

message("Loading more data..")
controlSampleIDs2 <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs2.rds", dirpath))
controlSampleIDs <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))
rsemtpm <- readRDS(file = sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", dirpath))
#ASE <- readRDS(file = sprintf("%s/data/ASE.rds", dirpath))
coding <- readRDS(file = sprintf("%s/aggregated_results/ASE.coding.rds", dirpath))
configs <- readRDS(sprintf("%s/aggregated_results/param_config_list.rds", dirpath))
geneRanges <- readRDS(sprintf("%s/aggregated_results/geneRanges.rds", dirpath))
high_qc_ids <- names(which(colSums(rsemtpm>5)>5000))
nonzeros.zs <- readRDS(file = sprintf("%s/aggregated_results/nonzeros.zs.bin50.rds", dirpath))

col_anno1 <- data.table(read.csv( sprintf("%s/work_in_progress/Annotation_list_long_vnikos_210129.csv", datapath) ))
dim(col_anno)
col_anno <- col_anno[ WTA.plate %in% high_qc_ids]
dim(col_anno)
col_anno[Pairs == "NA", Pairs := NA]

controlIDs <- col_anno[experimental_label %in% c("control_LCM_notreatment", "control_gen2_notreatment", "noco_gen2_noMN_lcm", "noco_gen1_noMN", "p53_gen2_noMN")]$WTA.plate
triIDs <- col_anno[experimental_label %in% c("Tri12_LCM", "Tri12_facs")]$WTA.plate[col_anno[experimental_label %in% c("Tri12_LCM", "Tri12_facs")]$WTA.plate %in% names(which(colSums(rsemtpm>0)>6000))]

#interstat.arm <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/interstat.armlevel.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")
#interstat.chr <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/interstat.chrlevel.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")
#fracstat <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/fracstat.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")

#source('~/WORK/Papers/MLpaper/R/utilities/calc_and_plot_intermediate_pvalue.R')
#source('~/WORK/Papers/MLpaper/R/utilities/plot_cn_and_intermediate_prediction.R')
#source('~/WORK/Papers/MLpaper/R/utilities/plot_raw_data_and_prediction_boxplots2.R')
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


########## RUNNING SSM (+1) for TE ###################################
source('/homes10/ejacob/WORK/Papers/MLpaper/R/utilities/GSSM_utils.R')

run_my_SSM1 <- function() {
  chrsToExcludeFromNormalization = c("chrX", "chrM", "chr10")
  SSM1 <- adt[, lapply(.SD, getMyGSSM),
            .SDcols = mysample_id, by = seqnames]
  
  dtm <- SSM1[, lapply(.SD, median),
            .SDcols = -1, by = "seqnames"]
  myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in% chrsToExcludeFromNormalization), -1]))
  SSM1 <- 2^(SSM1[[2]] - myCellMedians)

  return(SSM1)
}
preds.SSM.TE <- run_my_SSM1()

### BELOW UNFINISHED ###

#binding together all data as one table
SSM.TE <- preds.SSM.TE
preds <- SSM.TE[order(seqnames, start)]

outfname <- sprintf("%s/%s.SSM_CN_predictions.rds", dstdir, mysample_id)

saveRDS(preds, file = outfname)
message("Predictions saved to ", outfname)


