####### Global Variables ########
args <- c("/pellmanlab/nikos/Stam_Etai_Scripts", "/pellmanlab/stam_niko/data/processed_bam", "/pellmanlab/nikos/Stam_Etai_Data", "SIS1025a",  "SIS1025b",  "SIS1025d", "SIS1025e", "SIS1025f_Lane1" , "SIS1025f_Lane2")
#args <- commandArgs(trailingOnly = TRUE)
print(args)
scriptsdir <- args[1]
wkpath <- args[2]
datadir <- args[3]
experiments <- args[4:length(args)]
mk_agg_dir <- sprintf("mkdir -p %s/aggregated_results", wkpath)
system(mk_agg_dir)

require(data.table)
require(readxl)
require(data.table)
require(matrixStats)
source( sprintf('%s/scripts/runRNAseqPipelineFor_BF9VP.utils.R', scriptsdir) )
source( sprintf('%s/scripts/generate_baseline_datasets_for_ML_and_stats_analysis.R', scriptsdir) )
source( sprintf('%s/plots/class_prob_to_cn_level.R', scriptsdir) )

anno <- data.table(read.csv( sprintf("%s/work_in_progress/Annotation_list_long_vnikos_210129.csv", datadir) ))

#read in rsemtpm object
rsemtpm <- readRDS(file = sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", wkpath))

######################################################################################
######################################################################################
##############################Graph1##################################################
######################################################################################
######################################################################################


fromRawTPMtoExprsRatio <- function(rsemtpm, geneRanges, controlSampleIDs, max_pos_col = 6, plusOne = 1, min_sd = 0.0001, 
                                   zerosAsNA = F, normBySd = F,
                                   maxExp = 1200, quantileOfExp = 0.99, minDetectionLevel = 5, minNumOfSamplesToDetect = 3, 
                                   geneInterspaceLogBase = 8, interspaceCoef = 5, doNormalizeByControls = F) {

  wi <- which(rowSums(rsemtpm[, controlSampleIDs] >= minDetectionLevel) >= minNumOfSamplesToDetect)
  dt <- data.table(rsemtpm[wi,], keep.rownames = "id")
  if(length(names(which(table(colnames(dt))>1))) > 0)
    dt <- dt[, -which(colnames(dt) %in% names(which(table(colnames(dt))>1))), with=F]
  setkey(geneRanges, id)
  setkey(dt, id)
  dt <- merge(geneRanges, dt)
  require(gtools)
  dt$seqnames <- as.character(dt$seqnames)
  dt$seqnames <- factor(dt$seqnames, mixedsort(unique(dt$seqnames)))
  dt <- dt[order(seqnames, start, end)]
  dt <- dt[seqnames != "chrM" & seqnames != "chrY"]
  
  M <- copy(dt[, controlSampleIDs, with = F])
  if(zerosAsNA) 
    M[M==0] <- NA  
  
  message("NAs = ", sum(is.na(M)))
  
  q <- min(maxExp, unname(quantile(as.matrix(M), quantileOfExp, na.rm=T)))
  M[M > q] <- q
  message("Quantile q = ", q)
  dt2go <- copy(dt[, -c(1:max_pos_col), with=F])
  dt2go[dt2go>q] <- q
  
  #log before mean:
  dt2go <- log2(dt2go + plusOne)
  M <- log2(M + plusOne)
  
  message("NAs = ", sum(is.na(M)))
  #stats if doRatioBeforeScaleAndZscore == F
  means <- rowMeans(M, na.rm = T)
  sds <- rowSds(as.matrix(M), na.rm = T)
  
  #stats if doRatioBeforeScaleAndZscore == T
  M2 <- sweep(M, 1, means, "-")
  if(normBySd)
    M2 <- sweep(M2, 1, sds, "/") #* sqrt(mean(sds^2, na.rm=T))
  
  means2 <- rowMeans(M2, na.rm = T)
  sds2 <- rowSds(as.matrix(M2), na.rm = T)
  sds2 <- ifelse(sds2 < min_sd, min_sd, sds2)
  
  #gene interspace:
  interspace <- unlist(lapply(unique(dt$seqnames), 
                              function(chr) c(dt[seqnames == chr]$start[-1], dplyr::last(dt[seqnames == chr]$end)+1) - dt[seqnames == chr]$end + geneInterspaceLogBase))
  interspace[interspace <= 0] <- 1
  interspace <- round(log(x = interspace, base = geneInterspaceLogBase), digits = 2)
  interspace <- exp(-interspaceCoef*interspace/max(interspace))
  
  #calculating final dt:
  #subtracting the mean centers the data. Or you can divide by the std to standardize.
  if(doNormalizeByControls) {
    dt <- cbind(dt[, 1:max_pos_col, with = F], sweep(dt2go, 1, means, "-") )
  }
  if(normBySd) {
    dt <- cbind(dt[, 1:max_pos_col, with = F], sweep(dt2go, 1, sds, "/")) #* sqrt(mean(sds^2, na.rm=T)) )
  }
  if(doNormalizeByControls) {
    if(normBySd) {
      dt <- cbind(dt[, 1:max_pos_col, with = F], sweep(dt2go, 1, means, "-") )
      dt <- cbind(dt[, 1:max_pos_col, with = F], sweep(dt[, -c(1:max_pos_col), with = F], 1, sds, "/")) 
    }
  }
  if(!doNormalizeByControls) {
    if(!normBySd) {
      dt <- cbind(dt[, 1:max_pos_col, with = F], dt2go)
    }
  }
  
  res <- list(adt = dt, means = means, means2 = means2, sds = sds, sds2 = sds2, interspace = interspace, controlSampleIDs = controlSampleIDs)
  return(res)
}

generate_adt_adt0_adt.na_and_nonzeros_data <- function(rsemtpm, th = 4, dirpath = "/pellmanlab/stam_niko/data/processed_bam", scriptsdir = "/pellmanlab/nikos/Stam_Etai_Scripts", normBySd = F) {
  
  geneRanges <- readRDS(sprintf("%s/aggregated_results/geneRanges.rds", dirpath))
  controlSampleIDs2 <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs2.rds", dirpath))
  controlSampleIDs <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))
  
  #adt.na:
  #-------
  mlinput2.na <- with(configs, fromRawTPMtoExprsRatio(rsemtpm = rsemtpm, geneRanges = geneRanges, 
                                                      controlSampleIDs = controlSampleIDs2, max_pos_col = 6, plusOne = 1, 
                                                      maxExp = maxExp, quantileOfExp = quantileOfExp, 
                                                      minDetectionLevel = minDetectionLevel,
                                                      zerosAsNA = T, minNumOfSamplesToDetect = 3, #minNumOfSamplesToDetect, 
                                                      doNormalizeByControls = T, normBySd = normBySd))
  
  adt <- mlinput2.na$adt[,-c(5,6)]
  dim(adt)
  #controls:
  columns <- c("id", "seqnames", "start", "end", controlSampleIDs2)
  tpmc <- rsemtpm[names(which(rowSums(rsemtpm[, controlSampleIDs2]>0)> th)) ,controlSampleIDs2] 
  tpmc[tpmc>0] <- 1
  
  #samples:
  tpm <- adt[id %in% names(which(rowSums(rsemtpm[, controlSampleIDs2]>0)> th))]
  dim(tpm)
  mynonzeros <- rsemtpm[tpm$id, colnames(tpm)[-c(1:4)]]
  mytpm <- as.matrix(tpm[,-c(1:4)])
  stopifnot(sum(colnames(mytpm) != colnames(mynonzeros)) == 0)
  mytpm[mynonzeros == 0] <- NA
  adt.na <- cbind(tpm[,1:4], mytpm)
  adt.na <- adt.na[order(seqnames, start)]
  setkey(adt.na, id)
  dim(adt.na)
  #saveRDS(adt.na, file = sprintf("%s/aggregated_results/adt.na.rds", dirpath))
  
  #adt:
  #----
  mlinput <- with(configs, fromRawTPMtoExprsRatio(rsemtpm = rsemtpm, geneRanges = geneRanges, 
                                                  controlSampleIDs = controlSampleIDs2, max_pos_col = 6, plusOne = 1, 
                                                  maxExp = maxExp, quantileOfExp = quantileOfExp, 
                                                  minDetectionLevel = minDetectionLevel,
                                                  zerosAsNA = F, minNumOfSamplesToDetect = 3, #minNumOfSamplesToDetect, 
                                                  doNormalizeByControls = T, normBySd = normBySd))
  
  
  #remove unneed annotations
  adt <- mlinput$adt[,-c(5,6)]
  
  dim(adt)
  stopifnot(sum(!adt.na$id %in% adt$id) == 0)
  setkey(adt, id)
  
  adt <- adt[adt.na$id]
  dim(adt)
  
  return( list(adt, adt.na) )
}

#creates geneRanges object
load("/pellmanlab/stam_niko/refgenomes/Gencode/v25/gencode.v25.primary_assembly.ERCC92.annotation.geneRanges.RData")
names(geneRanges)[1] <- "id"

#creates adt rds objects for ML
tmppath <- "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2"
configs <- readRDS(sprintf("%s/data/param_config_list.rds", tmppath))

#To me it looks like etai's df is standardized as well
adt.noSd <- generate_adt_adt0_adt.na_and_nonzeros_data(rsemtpm, dirpath = wkpath, th = 4, scriptsdir = scriptsdir, normBySd = F) 

adt.noSd1 <- adt.noSd[[1]]
adt.na.noSd <- adt.noSd[[2]]

saveRDS(adt.noSd1, "/pellmanlab/stam_niko/data/processed_bam/test/adt.noSd.rds")
saveRDS(adt.na.noSd, "/pellmanlab/stam_niko/data/processed_bam/test/adt.na.noSd.rds")


######################################################################################
######################################################################################
##############################Graph2##################################################
######################################################################################
######################################################################################

generate_adt_adt0_adt.na_and_nonzeros_data_nolog_nocent <- function(rsemtpm, th = 4, dirpath = "/pellmanlab/stam_niko/data/processed_bam", scriptsdir = "/pellmanlab/nikos/Stam_Etai_Scripts", normBySd = F) {
  
  geneRanges <- readRDS(sprintf("%s/aggregated_results/geneRanges.rds", dirpath))
  controlSampleIDs2 <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs2.rds", dirpath))
  controlSampleIDs <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))
  
  #adt.na:
  #-------
  mlinput2.na <- with(configs, fromRawTPMtoExprsRatio_nolog(rsemtpm = rsemtpm, geneRanges = geneRanges, 
                                                            controlSampleIDs = controlSampleIDs2, max_pos_col = 6, plusOne = 1, 
                                                            maxExp = maxExp, quantileOfExp = quantileOfExp, 
                                                            minDetectionLevel = minDetectionLevel,
                                                            zerosAsNA = T, minNumOfSamplesToDetect = 3, #minNumOfSamplesToDetect, 
                                                            doNormalizeByControls = F, normBySd = normBySd))
  
  adt <- mlinput2.na$adt[,-c(5,6)]
  dim(adt)
  #controls:
  columns <- c("id", "seqnames", "start", "end", controlSampleIDs2)
  tpmc <- rsemtpm[names(which(rowSums(rsemtpm[, controlSampleIDs2]>0)> th)) ,controlSampleIDs2] 
  tpmc[tpmc>0] <- 1
  
  #samples:
  tpm <- adt[id %in% names(which(rowSums(rsemtpm[, controlSampleIDs2]>0)> th))]
  dim(tpm)
  mynonzeros <- rsemtpm[tpm$id, colnames(tpm)[-c(1:4)]]
  mytpm <- as.matrix(tpm[,-c(1:4)])
  stopifnot(sum(colnames(mytpm) != colnames(mynonzeros)) == 0)
  mytpm[mynonzeros == 0] <- NA
  adt.na <- cbind(tpm[,1:4], mytpm)
  adt.na <- adt.na[order(seqnames, start)]
  setkey(adt.na, id)
  dim(adt.na)
  #saveRDS(adt.na, file = sprintf("%s/aggregated_results/adt.na.rds", dirpath))
  
  #adt:
  #----
  mlinput <- with(configs, fromRawTPMtoExprsRatio_nolog(rsemtpm = rsemtpm, geneRanges = geneRanges, 
                                                        controlSampleIDs = controlSampleIDs2, max_pos_col = 6, plusOne = 1, 
                                                        maxExp = maxExp, quantileOfExp = quantileOfExp, 
                                                        minDetectionLevel = minDetectionLevel,
                                                        zerosAsNA = F, minNumOfSamplesToDetect = 3, #minNumOfSamplesToDetect, 
                                                        doNormalizeByControls = F, normBySd = normBySd))
  
  
  #remove unneed annotations
  adt <- mlinput$adt[,-c(5,6)]
  
  dim(adt)
  stopifnot(sum(!adt.na$id %in% adt$id) == 0)
  setkey(adt, id)
  
  adt <- adt[adt.na$id]
  dim(adt)
  return(adt)
}

#creates geneRanges object
load("/pellmanlab/stam_niko/refgenomes/Gencode/v25/gencode.v25.primary_assembly.ERCC92.annotation.geneRanges.RData")
names(geneRanges)[1] <- "id"
saveRDS(geneRanges, sprintf("%s/aggregated_results/geneRanges.rds", wkpath))

#creates adt rds objects for ML
tmppath <- "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2"
configs <- readRDS(sprintf("%s/data/param_config_list.rds", tmppath))

#To me it looks like etai's df is standardized as well
adt_nolog_nocent <- generate_adt_adt0_adt.na_and_nonzeros_data_nolog_nocent(rsemtpm, dirpath = wkpath, th = 4, scriptsdir = scriptsdir, normBySd = F) 

saveRDS(adt_nolog_nocent, "/pellmanlab/stam_niko/data/comparable_datasets/nikos_dataset/aggregated_results2/adt_nolog_nocent.rds")


######################################################################################
######################################################################################
##############################Graph3##################################################
######################################################################################
######################################################################################


fromRawTPMtoExprsRatio_nolog_meandiv <- function(rsemtpm, geneRanges, controlSampleIDs, max_pos_col = 6, plusOne = 1, min_sd = 0.0001, 
                                                 zerosAsNA = F, normBySd = F,
                                                 maxExp = 1200, quantileOfExp = 0.99, minDetectionLevel = 5, minNumOfSamplesToDetect = 3, 
                                                 geneInterspaceLogBase = 8, interspaceCoef = 5, doNormalizeByControls = F) {
  
  wi <- which(rowSums(rsemtpm[, controlSampleIDs] >= minDetectionLevel) >= minNumOfSamplesToDetect)
  dt <- data.table(rsemtpm[wi,], keep.rownames = "id")
  if(length(names(which(table(colnames(dt))>1))) > 0)
    dt <- dt[, -which(colnames(dt) %in% names(which(table(colnames(dt))>1))), with=F]
  setkey(geneRanges, id)
  setkey(dt, id)
  dt <- merge(geneRanges, dt)
  require(gtools)
  dt$seqnames <- as.character(dt$seqnames)
  dt$seqnames <- factor(dt$seqnames, mixedsort(unique(dt$seqnames)))
  dt <- dt[order(seqnames, start, end)]
  dt <- dt[seqnames != "chrM" & seqnames != "chrY"]
  
  M <- copy(dt[, controlSampleIDs, with = F])
  if(zerosAsNA) 
    M[M==0] <- NA  
  
  message("NAs = ", sum(is.na(M)))
  
  q <- min(maxExp, unname(quantile(as.matrix(M), quantileOfExp, na.rm=T)))
  M[M > q] <- q
  message("Quantile q = ", q)
  dt2go <- copy(dt[, -c(1:max_pos_col), with=F])
  dt2go[dt2go>q] <- q
  
  #log before mean:
  #dt2go <- log2(dt2go + plusOne)
  #M <- log2(M + plusOne)
  
  message("NAs = ", sum(is.na(M)))
  #stats if doRatioBeforeScaleAndZscore == F
  means <- rowMeans(M, na.rm = T)
  sds <- rowSds(as.matrix(M), na.rm = T)
  
  #stats if doRatioBeforeScaleAndZscore == T
  M2 <- sweep(M, 1, means, "-")
  if(normBySd)
    M2 <- sweep(M2, 1, sds, "/") #* sqrt(mean(sds^2, na.rm=T))
  
  means2 <- rowMeans(M2, na.rm = T)
  sds2 <- rowSds(as.matrix(M2), na.rm = T)
  sds2 <- ifelse(sds2 < min_sd, min_sd, sds2)
  
  #gene interspace:
  interspace <- unlist(lapply(unique(dt$seqnames), 
                              function(chr) c(dt[seqnames == chr]$start[-1], dplyr::last(dt[seqnames == chr]$end)+1) - dt[seqnames == chr]$end + geneInterspaceLogBase))
  interspace[interspace <= 0] <- 1
  interspace <- round(log(x = interspace, base = geneInterspaceLogBase), digits = 2)
  interspace <- exp(-interspaceCoef*interspace/max(interspace))
  
  #calculating final dt:
  #subtracting the mean centers the data. Or you can divide by the std to standardize.
  if(doNormalizeByControls) {
    ddt <- cbind(dt[, 1:max_pos_col, with = F], sweep(dt2go, 1, means, "/") )
  }
  if(normBySd) {
    ddt <- cbind(dt[, 1:max_pos_col, with = F], sweep(dt2go, 1, sds, "/")) #* sqrt(mean(sds^2, na.rm=T)) )
  }
  if(doNormalizeByControls) {
    if(normBySd) {
      ddt <- cbind(dt[, 1:max_pos_col, with = F], sweep(dt2go, 1, means, "/") )
      ddt <- cbind(dt[, 1:max_pos_col, with = F], sweep(ddt[, -c(1:max_pos_col), with = F], 1, sds, "/")) 
    }
  }
  if(!doNormalizeByControls) {
    if(!normBySd) {
      ddt <- cbind(dt[, 1:max_pos_col, with = F], dt2go)
    }
  }
  
  res <- list(adt = ddt, means = means, means2 = means2, sds = sds, sds2 = sds2, interspace = interspace, controlSampleIDs = controlSampleIDs)
  return(res)
}

generate_adt_adt0_adt.na_and_nonzeros_data_nolog_meandiv_sddiv <- function(rsemtpm, th = 4, dirpath = "/pellmanlab/stam_niko/data/processed_bam", scriptsdir = "/pellmanlab/nikos/Stam_Etai_Scripts", normBySd = F) {
  
  geneRanges <- readRDS(sprintf("%s/aggregated_results/geneRanges.rds", dirpath))
  controlSampleIDs2 <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs2.rds", dirpath))
  controlSampleIDs <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))
  
  #adt.na:
  #-------
  mlinput2.na <- with(configs, fromRawTPMtoExprsRatio_nolog_meandiv(rsemtpm = rsemtpm, geneRanges = geneRanges, 
                                                                    controlSampleIDs = controlSampleIDs2, max_pos_col = 6, plusOne = 1, 
                                                                    maxExp = maxExp, quantileOfExp = quantileOfExp, 
                                                                    minDetectionLevel = minDetectionLevel,
                                                                    zerosAsNA = T, minNumOfSamplesToDetect = 3, #minNumOfSamplesToDetect, 
                                                                    doNormalizeByControls = T, normBySd = normBySd))
  
  adt <- mlinput2.na$adt[,-c(5,6)]
  dim(adt)
  #controls:
  columns <- c("id", "seqnames", "start", "end", controlSampleIDs2)
  tpmc <- rsemtpm[names(which(rowSums(rsemtpm[, controlSampleIDs2]>0)> th)) ,controlSampleIDs2] 
  tpmc[tpmc>0] <- 1
  
  #samples:
  tpm <- adt[id %in% names(which(rowSums(rsemtpm[, controlSampleIDs2]>0)> th))]
  dim(tpm)
  mynonzeros <- rsemtpm[tpm$id, colnames(tpm)[-c(1:4)]]
  mytpm <- as.matrix(tpm[,-c(1:4)])
  stopifnot(sum(colnames(mytpm) != colnames(mynonzeros)) == 0)
  mytpm[mynonzeros == 0] <- NA
  adt.na <- cbind(tpm[,1:4], mytpm)
  adt.na <- adt.na[order(seqnames, start)]
  setkey(adt.na, id)
  dim(adt.na)
  #saveRDS(adt.na, file = sprintf("%s/aggregated_results/adt.na.rds", dirpath))
  
  #adt:
  #----
  mlinput <- with(configs, fromRawTPMtoExprsRatio_nolog_meandiv(rsemtpm = rsemtpm, geneRanges = geneRanges, 
                                                                controlSampleIDs = controlSampleIDs2, max_pos_col = 6, plusOne = 1, 
                                                                maxExp = maxExp, quantileOfExp = quantileOfExp, 
                                                                minDetectionLevel = minDetectionLevel,
                                                                zerosAsNA = F, minNumOfSamplesToDetect = 3, #minNumOfSamplesToDetect, 
                                                                doNormalizeByControls = T, normBySd = normBySd))
  
  
  #remove unneed annotations
  adt <- mlinput$adt[,-c(5,6)]
  
  dim(adt)
  stopifnot(sum(!adt.na$id %in% adt$id) == 0)
  setkey(adt, id)
  
  adt <- adt[adt.na$id]
  dim(adt)
  return(adt)
}

adt_nolog_meandiv <- generate_adt_adt0_adt.na_and_nonzeros_data_nolog_meandiv_sddiv(rsemtpm, dirpath = wkpath, th = 4, scriptsdir = scriptsdir, normBySd = F) 

saveRDS(adt_nolog_meandiv, "/pellmanlab/stam_niko/data/comparable_datasets/nikos_dataset/aggregated_results2/adt_nolog_meandiv.rds")


######################################################################################
######################################################################################
##############################Graph4##################################################
######################################################################################
######################################################################################


fromRawTPMtoExprsRatio_nolog_meandiv <- function(rsemtpm, geneRanges, controlSampleIDs, max_pos_col = 6, plusOne = 1, min_sd = 0.0001, 
                                                 zerosAsNA = F, normBySd = F,
                                                 maxExp = 1200, quantileOfExp = 0.99, minDetectionLevel = 5, minNumOfSamplesToDetect = 3, 
                                                 geneInterspaceLogBase = 8, interspaceCoef = 5, doNormalizeByControls = F) {
  
  wi <- which(rowSums(rsemtpm[, controlSampleIDs] >= minDetectionLevel) >= minNumOfSamplesToDetect)
  dt <- data.table(rsemtpm[wi,], keep.rownames = "id")
  if(length(names(which(table(colnames(dt))>1))) > 0)
    dt <- dt[, -which(colnames(dt) %in% names(which(table(colnames(dt))>1))), with=F]
  setkey(geneRanges, id)
  setkey(dt, id)
  dt <- merge(geneRanges, dt)
  require(gtools)
  dt$seqnames <- as.character(dt$seqnames)
  dt$seqnames <- factor(dt$seqnames, mixedsort(unique(dt$seqnames)))
  dt <- dt[order(seqnames, start, end)]
  dt <- dt[seqnames != "chrM" & seqnames != "chrY"]
  
  M <- copy(dt[, controlSampleIDs, with = F])
  if(zerosAsNA) 
    M[M==0] <- NA  
  
  message("NAs = ", sum(is.na(M)))
  
  q <- min(maxExp, unname(quantile(as.matrix(M), quantileOfExp, na.rm=T)))
  M[M > q] <- q
  message("Quantile q = ", q)
  dt2go <- copy(dt[, -c(1:max_pos_col), with=F])
  dt2go[dt2go>q] <- q
  
  #log before mean:
  #dt2go <- log2(dt2go + plusOne)
  #M <- log2(M + plusOne)
  
  message("NAs = ", sum(is.na(M)))
  #stats if doRatioBeforeScaleAndZscore == F
  means <- rowMeans(M, na.rm = T)
  sds <- rowSds(as.matrix(M), na.rm = T)
  
  #stats if doRatioBeforeScaleAndZscore == T
  M2 <- sweep(M, 1, means, "-")
  if(normBySd)
    M2 <- sweep(M2, 1, sds, "/") #* sqrt(mean(sds^2, na.rm=T))
  
  means2 <- rowMeans(M2, na.rm = T)
  sds2 <- rowSds(as.matrix(M2), na.rm = T)
  sds2 <- ifelse(sds2 < min_sd, min_sd, sds2)
  
  #gene interspace:
  interspace <- unlist(lapply(unique(dt$seqnames), 
                              function(chr) c(dt[seqnames == chr]$start[-1], dplyr::last(dt[seqnames == chr]$end)+1) - dt[seqnames == chr]$end + geneInterspaceLogBase))
  interspace[interspace <= 0] <- 1
  interspace <- round(log(x = interspace, base = geneInterspaceLogBase), digits = 2)
  interspace <- exp(-interspaceCoef*interspace/max(interspace))
  
  #calculating final dt:
  #subtracting the mean centers the data. Or you can divide by the std to standardize.
  if(doNormalizeByControls) {
    ddt <- cbind(dt[, 1:max_pos_col, with = F], sweep(dt2go, 1, means, "/") )
  }
  if(normBySd) {
    ddt <- cbind(dt[, 1:max_pos_col, with = F], sweep(dt2go, 1, sds, "/")) #* sqrt(mean(sds^2, na.rm=T)) )
  }
  if(doNormalizeByControls) {
    if(normBySd) {
      ddt <- cbind(dt[, 1:max_pos_col, with = F], sweep(dt2go, 1, means, "/") )
      ddt <- cbind(dt[, 1:max_pos_col, with = F], sweep(ddt[, -c(1:max_pos_col), with = F], 1, sds, "/")) 
    }
  }
  if(!doNormalizeByControls) {
    if(!normBySd) {
      ddt <- cbind(dt[, 1:max_pos_col, with = F], dt2go)
    }
  }
  
  res <- list(adt = ddt, means = means, means2 = means2, sds = sds, sds2 = sds2, interspace = interspace, controlSampleIDs = controlSampleIDs)
  return(res)
}

generate_adt_adt0_adt.na_and_nonzeros_data_nolog_meandiv_sddiv <- function(rsemtpm, th = 4, dirpath = "/pellmanlab/stam_niko/data/processed_bam", scriptsdir = "/pellmanlab/nikos/Stam_Etai_Scripts", normBySd = F) {
  
  geneRanges <- readRDS(sprintf("%s/aggregated_results/geneRanges.rds", dirpath))
  controlSampleIDs2 <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs2.rds", dirpath))
  controlSampleIDs <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))
  
  #adt.na:
  #-------
  mlinput2.na <- with(configs, fromRawTPMtoExprsRatio_nolog_meandiv(rsemtpm = rsemtpm, geneRanges = geneRanges, 
                                                                    controlSampleIDs = controlSampleIDs2, max_pos_col = 6, plusOne = 1, 
                                                                    maxExp = maxExp, quantileOfExp = quantileOfExp, 
                                                                    minDetectionLevel = minDetectionLevel,
                                                                    zerosAsNA = T, minNumOfSamplesToDetect = 3, #minNumOfSamplesToDetect, 
                                                                    doNormalizeByControls = T, normBySd = normBySd))
  
  adt <- mlinput2.na$adt[,-c(5,6)]
  dim(adt)
  #controls:
  columns <- c("id", "seqnames", "start", "end", controlSampleIDs2)
  tpmc <- rsemtpm[names(which(rowSums(rsemtpm[, controlSampleIDs2]>0)> th)) ,controlSampleIDs2] 
  tpmc[tpmc>0] <- 1
  
  #samples:
  tpm <- adt[id %in% names(which(rowSums(rsemtpm[, controlSampleIDs2]>0)> th))]
  dim(tpm)
  mynonzeros <- rsemtpm[tpm$id, colnames(tpm)[-c(1:4)]]
  mytpm <- as.matrix(tpm[,-c(1:4)])
  stopifnot(sum(colnames(mytpm) != colnames(mynonzeros)) == 0)
  mytpm[mynonzeros == 0] <- NA
  adt.na <- cbind(tpm[,1:4], mytpm)
  adt.na <- adt.na[order(seqnames, start)]
  setkey(adt.na, id)
  dim(adt.na)
  #saveRDS(adt.na, file = sprintf("%s/aggregated_results/adt.na.rds", dirpath))
  
  #adt:
  #----
  mlinput <- with(configs, fromRawTPMtoExprsRatio_nolog_meandiv(rsemtpm = rsemtpm, geneRanges = geneRanges, 
                                                                controlSampleIDs = controlSampleIDs2, max_pos_col = 6, plusOne = 1, 
                                                                maxExp = maxExp, quantileOfExp = quantileOfExp, 
                                                                minDetectionLevel = minDetectionLevel,
                                                                zerosAsNA = F, minNumOfSamplesToDetect = 3, #minNumOfSamplesToDetect, 
                                                                doNormalizeByControls = T, normBySd = normBySd))
  
  
  #remove unneed annotations
  adt <- mlinput$adt[,-c(5,6)]
  
  dim(adt)
  stopifnot(sum(!adt.na$id %in% adt$id) == 0)
  setkey(adt, id)
  
  adt <- adt[adt.na$id]
  dim(adt)
  return(adt)
}

adt_nolog_meandiv_sddiv <- generate_adt_adt0_adt.na_and_nonzeros_data_nolog_meandiv_sddiv(rsemtpm, dirpath = wkpath, th = 4, scriptsdir = scriptsdir, normBySd = T) 

saveRDS(adt_nolog_meandiv_sddiv, "/pellmanlab/stam_niko/data/comparable_datasets/nikos_dataset/aggregated_results2/adt_nolog_meandiv_sddiv.rds")

