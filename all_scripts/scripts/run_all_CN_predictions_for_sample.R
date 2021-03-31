#running all type of predictions for a sample

args <- commandArgs(trailingOnly = TRUE)


if(length(args) < 3) {
  stop("Wrong input.\nRscript run_prediction_for_sample.R sample_id secondary_id dstdir")
}

mysample_id <-  "160627_B1" #args[1] #triIDs[2]  #"170512_B4"  #args[1]
secondary_id <- "" #args[2] #args[2]
dstdir <-  "/pellmanlab/stam_niko/data/processed_bam/CN_data_etai" #args[3]

dirpath <- "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2"

require(data.table)


message("Loading data..")
adt.na <- readRDS(file = sprintf("%s/data/adt.na.rds", dirpath))
adt <- readRDS(file = sprintf("%s/data/adt.rds", dirpath))
nonzeros <- readRDS(file = sprintf("%s/data/nonzeros.rds", dirpath))
geneSpecific <- readRDS(file = sprintf("%s/data/geneSpecific.rds", dirpath))
ganno <- readRDS(file = sprintf("%s/data/ganno.rds", dirpath))
alleles.all <- readRDS(file = sprintf("%s/data/alleles.all.rds", dirpath)) # "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2/data/alleles.all.rds"


message("Loading more data..")
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

#interstat.arm <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/interstat.armlevel.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")
#interstat.chr <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/interstat.chrlevel.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")
#fracstat <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/fracstat.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")
#centromeres <- readRDS("~ejacob/WORK/secondaryanalysis/Stamatis/data/centromeres.rds")

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



########## RUNNING ML ###################################
message("Calculating features for predictions..")
source('/homes10/ejacob/WORK/Papers/MLpaper/R/utilities/generate_ftrs_for_input.R')
tt <- generate_ftrs_for_input(adt = adt, adt.na = adt.na, nonzeros = nonzeros, ganno = ganno,
                              sample_id = mysample_id)
dim(tt)
m <- as.matrix((tt[,manual_selected_vars6, with=F]))
source('/homes10/ejacob/WORK/Papers/MLpaper/R/utilities/predict_using_lightgbm.R')
preds.ML <- cbind(tt[,1:6], predict_using_lightgbm(m))
dim(preds.ML)
colnames(preds.ML)[7:9] <- c("ML.CN1", "ML.CN2", "ML.CN3")
#plot(max.col(preds.ML[,7:9]))
########## RUNNING MA ###################################

run_my_MA <- function(MA_winSize = 100) {
  chrsToExcludeFromNormalization = c("chrX", "chrM", "chr10")
  MA <- adt[, lapply(.SD, rollmean, k = MA_winSize, fill = "extend", align = "center"),
            .SDcols = mysample_id, by = seqnames]
  dtm <- MA[, lapply(.SD, median),
            .SDcols = -1, by = "seqnames"]
  myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in%chrsToExcludeFromNormalization), -1]))
  
  MA2 <- 2^sweep(MA[,-1], 2, myCellMedians, "-")
  
  return(MA2)
}

preds.MA <- do.call(cbind, lapply(c(50, 75, 100), run_my_MA))
colnames(preds.MA) <- c("MA50", "MA75", "MA100")

stopifnot(sum(adt$id != preds.MA$id) == 0)


########## RUNNING SSM (+1) for TE ###################################
source('~/WORK/Papers/MLpaper/R/utilities/GSSM_utils.R')

run_my_SSM1 <- function() {
  chrsToExcludeFromNormalization = c("chrX", "chrM", "chr10")
  SSM1 <- adt[, lapply(.SD, getMyGSSM),
            .SDcols = mysample_id, by = seqnames]
  
  dtm <- SSM1[, lapply(.SD, median),
            .SDcols = -1, by = "seqnames"]
  myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in%chrsToExcludeFromNormalization), -1]))
  SSM1 <- 2^(SSM1[[2]] - myCellMedians)

  return(SSM1)
}
preds.SSM.TE <- run_my_SSM1()


########## RUNNING IMF ###################################
source('~/WORK/Papers/MLpaper/R/utilities/runIterativeMedianFilter.R')
run_my_IMF <- function() {
  chrsToExcludeFromNormalization = c("chrX", "chrM", "chr10")
  IMF <- adt[, lapply(.SD, function(x) runIterativeMedianFilter(x, min_win = 15, levels = 5)[[1]]),
              .SDcols = mysample_id, by = seqnames]

  dtm <- IMF[, lapply(.SD, median),
              .SDcols = -1, by = "seqnames"]
  myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in%chrsToExcludeFromNormalization), -1]))
  IMF <- 2^(IMF[[2]] - myCellMedians)
  
  return(IMF)
}
preds.IMF <- run_my_IMF()

########## RUNNING MV-SSM for alleles and TE  ###################################

#Calculation of MV-SSM for 3 variable alleles +0 and +1:
#-------------------------------------------------------
chrsToExcludeFromNormalization = c("chrX", "chrM", "chr10")

source('~/WORK/Papers/MLpaper/R/utilities/mergeTotalExpAndAlleleData.R')

message("Doing MV-SSM predictions..")
run_my_MVSSM3 <- function() {
  getMyAlleleLevelSSMs <- function(i, alleles.all, adt, chr = NULL) {
    myalleles <- mergeTotalExpAndAlleleData(adt, A = alleles.all$An, B = alleles.all$Bn, id = i)
    if(is.null(chr)) {
      tt <- data.table(getMy3VariateGSSM(mv3 = as.matrix(myalleles[, c("TE", "A", "B"), with=F])))
      return(cbind(myalleles, tt))
    } else {
      tt <- data.table(getMy3VariateGSSM(mv3 = as.matrix(myalleles[seqnames==chr, c("TE", "A", "B"), with=F])))
      return(cbind(myalleles[seqnames == chr], tt))
    }
  }
  
  tt <- rbindlist(lapply(as.character(unique(adt$seqnames)), 
                         function(chr) getMyAlleleLevelSSMs(i = mysample_id, alleles.all = alleles.all, adt = adt, chr = chr)))
  dtm <- tt[, lapply(.SD, median),
             .SDcols = c("level.TE", "level.A", "level.B"), by = "seqnames"]
  myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in%chrsToExcludeFromNormalization), -1]))
  tt2 <- cbind(tt[,1:7], 2^sweep(tt[, c("level.TE", "level.A", "level.B")], 2, myCellMedians, "-"))
  setnames(tt2, old = c("level.TE", "level.A", "level.B"), new = c("SSM.TEAB.TE", "SSM.TEAB.A", "SSM.TEAB.B"))
  return(tt2)
}

preds.SSM.TEAB <- run_my_MVSSM3()



########## RUNNING MV-SSM for allele A or B and TE  ###################################



run_my_MVSSM2 <- function() {
  getMyOneAlleleLevelSSMs <- function(i, alleles.all, adt, chr = NULL, allele = "A") {
    myalleles <- mergeTotalExpAndAlleleData(adt, A = alleles.all$An, B = alleles.all$Bn, id = i)
    if(is.null(chr)) {
      tt <- data.table(getMy2VariateGSSM(mv2 = as.matrix(myalleles[, c(allele, "TE"), with=F])))
      return(cbind(myalleles, tt))
    } else {
      tt <- data.table(getMy2VariateGSSM(mv2 = as.matrix(myalleles[seqnames==chr, c(allele, "TE"), with=F])))
      return(cbind(myalleles[seqnames == chr], tt))
    }
  }
  
  tt.A <- rbindlist(lapply(as.character(unique(adt$seqnames)), 
                           function(chr) getMyOneAlleleLevelSSMs(i = mysample_id, alleles.all = alleles.all, adt = adt, chr = chr, allele = "A")))
  dtm <- tt.A[, lapply(.SD, median),
            .SDcols = c("V1"), by = "seqnames"]
  myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in%chrsToExcludeFromNormalization), -1]))
  tt2.A <- 2^(as.numeric(tt.A$V1) - myCellMedians)
  
  tt.B <- rbindlist(lapply(as.character(unique(adt$seqnames)), 
                           function(chr) getMyOneAlleleLevelSSMs(i = mysample_id, alleles.all = alleles.all, adt = adt, chr = chr, allele = "B")))
  dtm <- tt.B[, lapply(.SD, median),
              .SDcols = c("V1"), by = "seqnames"]
  myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in%chrsToExcludeFromNormalization), -1]))
  tt2.B <- 2^(as.numeric(tt.B$V1) - myCellMedians)
  stopifnot(sum(tt.A$gene_id != tt.B$gene_id) == 0)
  
  tt <- cbind(tt.A[, 1:7], SSM.TEAorB.A = tt2.A, SSM.TEAorB.B = tt2.B)
  return(tt)
}

preds.SSM.TEAorB <- run_my_MVSSM2()

#binding together all data as one table

stopifnot(sum(preds.ML$id != adt$id) == 0)
preds <- cbind(preds.ML[, 1:6], TPM.LFR = adt[[mysample_id]], preds.ML[,7:9], preds.MA, SSM.TE = preds.SSM.TE, IMF = preds.IMF)
preds2 <- merge(preds, preds.SSM.TEAB[, -c(5:7)], by.y = "gene_id", by.x = "id", all = T)
preds3 <- merge(preds2, preds.SSM.TEAorB[, -c(2:7)], by.y = "gene_id", by.x = "id", all = T)
preds3 <- preds3[order(seqnames, start)]

outfname <- sprintf("%s/%s.all_CN_predictions.rds", dstdir, mysample_id)

saveRDS(preds3, file = outfname)
message("Predictions saved to ", outfname)

# 
# chr <- "chr12"
# plot(preds3[seqnames == chr]$start, max.col(preds3[seqnames == chr][, c("ML.CN1", "ML.CN2", "ML.CN3")]), ylim=c(0,4))
# points(preds3[seqnames == chr]$start, preds3[seqnames == chr]$SSM.TEAB.A, col="blue")
# points(preds3[seqnames == chr]$start, preds3[seqnames == chr]$SSM.TEAB.B, col="red")
# points(preds.SSM.TEAB[seqnames == chr]$start, preds.SSM.TEAB[seqnames == chr]$level.B, col="red")
# points(preds3[seqnames == chr]$start, preds3[seqnames == chr]$SSM.TEAB.TE*2, col="brown")
# 
# points(preds3[seqnames == chr]$start, preds3[seqnames == chr]$SSM.TEAorB.A, col="darkblue")
# points(preds3[seqnames == chr]$start, preds3[seqnames == chr]$SSM.TEAorB.B, col="darkred")
# 
# 
# 
# 
# 

