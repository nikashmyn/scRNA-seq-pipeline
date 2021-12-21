#generate_baseline_datasets_for_ML_and_stats_analysis.R
#for CNV12.2

#Loading data:
require(data.table)

generate_all_genes_annotation <- function() {
  gtf <- rtracklayer::import("/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.annotation.gtf")
  tmp <- as.data.table(gtf)
  genesanno <- tmp[type=="gene"]
  saveRDS(genesanno, file = sprintf("%s/data/all_genes_anno.rds", dirpath))
}

generate_coding_snps_data <- function(generateAlleleData = F, dirpath = "/pellmanlab/stam_niko/data/processed_bam") {
  #Input Data
  adt.na <- readRDS(sprintf("%s/aggregated_results/adt.na.rds", dirpath))
  nonzeros <- readRDS(sprintf("%s/aggregated_results/nonzeros.rds", dirpath))
  adt <- readRDS(sprintf("%s/aggregated_results/adt.rds", dirpath))
  ASE <- readRDS(file = sprintf("%s/aggregated_results/ASE.rds", dirpath))
    
  m3 <- ASE
  #strict filters
  sum(m3$features$allele==0)
  wi <- with(m3$features, which(biAllelicFilter == T & (!is.na(coding) | (!is.na(fiveUTR)) | (!is.na(threeUTR)))))
  stopifnot(length(wi) > 0)
  coding <- list(features = m3$features[wi,], A = m3$A[wi,], B = m3$B[wi, ])
  wi2 <- with(m3$features, which(!(!is.na(coding) | (!is.na(fiveUTR)) | (!is.na(threeUTR)))))
  noncoding <- list(features = m3$features[wi2,], A = m3$A[wi2,], B = m3$B[wi2, ])
  
  #coding
  coding$A[is.na(coding$A)] <- 0
  coding$B[is.na(coding$B)] <- 0
  coding$A[coding$A > 0] <- 1
  coding$B[coding$B > 0] <- 1
  
  wi <- which(rowSums(coding$A + coding$B)==0)
  
  if(length(wi) == 0) {
    coding$A <- cbind(coding$features[,1:4], coding$A)
    coding$B <- cbind(coding$features[,1:4], coding$B)
  } else {
    coding$A <- cbind(coding$features[-wi,1:4], coding$A[-wi,])
    coding$B <- cbind(coding$features[-wi,1:4], coding$B[-wi,])
    
  }
  coding$A <- coding$A[order(seqnames, start, decreasing = F)]
  coding$B <- coding$B[order(seqnames, start, decreasing = F)]
  unique(coding$B$seqnames)
  stopifnot(sum(coding$A$id != coding$B$id)==0)
  
  mygenes <- adt.na[,1:4]
  
  setkeyv(mygenes, c("seqnames", "start", "end"))
  
  #A:
  setnames(coding$A, old = "id", new = "snp_id")
  setkeyv(coding$A, c("seqnames", "start", "end"))
  tmp <- foverlaps(mygenes, coding$A)
  tmp2 <- tmp[,c(c("id", "snp_id", "seqnames", "i.start", "i.end"), colnames(tmp)[-c(1:4, ncol(tmp)-c(0:2))]), with=F]
  setnames(tmp2, old = c("i.start", "i.end"), new = c("start", "end"))
  A <- tmp2
  
  #B:
  setnames(coding$B, old = "id", new = "snp_id")
  setkeyv(coding$B, c("seqnames", "start", "end"))
  tmp <- foverlaps(mygenes, coding$B)
  tmp2 <- tmp[,c(c("id", "snp_id", "seqnames", "i.start", "i.end"), colnames(tmp)[-c(1:4, ncol(tmp)-c(0:2))]), with=F]
  setnames(tmp2, old = c("i.start", "i.end"), new = c("start", "end"))
  B <- tmp2
  
  stopifnot(sum(A$snp_id != B$snp_id, na.rm = T)==0)
  coding_snps <- list(A = A, B = B)
  saveRDS(coding_snps, file = sprintf("%s/aggregated_results/coding_snps.rds", dirpath))
  
  nonzeros <- nonzeros[order(seqnames, start)]
  
  bins <- lapply(unique(nonzeros$seqnames), function(chr) unlist(lapply(1:ceiling(nrow(nonzeros[seqnames==chr])/50), rep, 50))[1:nrow(nonzeros[seqnames==chr])])
  names(bins) <- unique(nonzeros$seqnames)
  bins2 <- rbindlist(lapply(names(bins), function(chr) data.table(seqnames = chr, bin = bins[[chr]])))
  stopifnot(sum(bins2$seqnames != nonzeros$seqnames)==0)
  nonzeros2 <- cbind(bin = bins2$bin, nonzeros[,1:4])
  mybins <- nonzeros2
  mybins2 <- merge(mybins[, lapply(.SD, min), 
                          .SDcols = c("start"), by = c("seqnames", "bin")], 
                   mybins[, lapply(.SD, max), 
                          .SDcols = c("end"), by = c("seqnames", "bin")], 
                   by = c("seqnames", "bin"))
  
  mybins <- mybins2
  
  
  #merging:
  setkeyv(mybins, c("seqnames", "start", "end"))
  setkeyv(coding$A, c("seqnames", "start", "end"))
  setkeyv(coding$B, c("seqnames", "start", "end"))
  
  coding$bins.A <- foverlaps(x = coding$A, y = mybins, nomatch = NULL)
  coding$bins.B <- foverlaps(x = coding$B, y = mybins, nomatch = NULL)
  stopifnot(sum(coding$bins.A$id != coding$bins.B$id) == 0)
  
  #coding:
  coding$cnts.A <- coding$bins.A[, lapply(.SD, function(x) sum(x, na.rm = T)), .SDcols = colnames(coding$A)[-c(1:4)], by = c("seqnames", "bin")]
  coding$cnts.B <- coding$bins.B[, lapply(.SD, function(x) sum(x, na.rm = T)), .SDcols = colnames(coding$B)[-c(1:4)], by = c("seqnames", "bin")]
  
  coding$aggs.A <- coding$bins.A[, lapply(.SD, function(x) mean(x, na.rm = T)), .SDcols = colnames(coding$A)[-c(1:4)], by = c("seqnames", "bin")]
  coding$aggs.B <- coding$bins.B[, lapply(.SD, function(x) mean(x, na.rm = T)), .SDcols = colnames(coding$B)[-c(1:4)], by = c("seqnames", "bin")]
  
  coding$aggs.A.zs <- cbind(coding$aggs.A[, c(1:2)], scale(coding$aggs.A[, -c(1:2)], center = T, scale = T))
  coding$aggs.B.zs <- cbind(coding$aggs.B[, c(1:2)], scale(coding$aggs.B[, -c(1:2)], center = T, scale = T))
  
  
}

get_snps_fraction_bins_zscore <- function(th = 10, saveMe = F, generateAlleleData = F, dirpath = "/pellmanlab/stam_niko/data/processed_bam") {
  #prep:
  m3 <- readRDS(file = sprintf("%s/aggregated_results/ASE.rds", dirpath)) 
  adt.na <- readRDS(sprintf("%s/aggregated_results/adt.na.rds", dirpath))
  nonzeros <- readRDS(sprintf("%s/aggregated_results/nonzeros.rds", dirpath))
  adt <- readRDS(sprintf("%s/aggregated_results/adt.rds", dirpath))
  
  #strict filters
  sum(m3$features$allele==0)
  wi <- with(m3$features, which(biAllelicFilter == T & (!is.na(coding) | (!is.na(fiveUTR)) | (!is.na(threeUTR)))))
  
  coding <- list(features = m3$features[wi,], A = m3$A[wi,], B = m3$B[wi, ])
  wi2 <- with(m3$features, which(!(!is.na(coding) | (!is.na(fiveUTR)) | (!is.na(threeUTR)))))
  noncoding <- list(features = m3$features[wi2,], A = m3$A[wi2,], B = m3$B[wi2, ])
  
  #coding
  coding$A[is.na(coding$A)] <- 0
  coding$B[is.na(coding$B)] <- 0
  coding$A[coding$A > 0] <- 1
  coding$B[coding$B > 0] <- 1
  
  wi <- which(rowSums(coding$A + coding$B)==0)
  
  if(length(wi) == 0) {
    coding$A <- cbind(coding$features[,1:4], coding$A)
    coding$B <- cbind(coding$features[,1:4], coding$B)
  } else {
    coding$A <- cbind(coding$features[-wi,1:4], coding$A[-wi,])
    coding$B <- cbind(coding$features[-wi,1:4], coding$B[-wi,])
    
  }
  coding$A <- coding$A[order(seqnames, start, decreasing = F)]
  coding$B <- coding$B[order(seqnames, start, decreasing = F)]
  unique(coding$B$seqnames)
  stopifnot(sum(coding$A$id != coding$B$id)==0)
  #saveRDS(coding, file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/ASE.coding.rds")
  
  #noncoding:
  noncoding$A[is.na(noncoding$A)] <- 0
  noncoding$B[is.na(noncoding$B)] <- 0
  noncoding$A[noncoding$A > 0] <- 1
  noncoding$B[noncoding$B > 0] <- 1
  
  wi <- which(rowSums(noncoding$A + noncoding$B)==0)
  
  noncoding$A <- cbind(noncoding$features[-wi,1:4], noncoding$A[-wi,])
  noncoding$B <- cbind(noncoding$features[-wi,1:4], noncoding$B[-wi,])
  
  noncoding$A <- noncoding$A[order(seqnames, start, decreasing = F)]
  noncoding$B <- noncoding$B[order(seqnames, start, decreasing = F)]
  unique(noncoding$B$seqnames)
  
  stopifnot(sum(noncoding$A$id != noncoding$B$id)==0)
  #saveRDS(noncoding, file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/ASE.noncoding.rds")
  
  #calculating bins from genes:
  # columns <- c("id", "seqnames", "start", "end")
  # ml <- MLREG[,columns, with=F]
  # nonzeros <- rsemtpm[names(which(rowSums(rsemtpm[, controlSampleIDs]>0)> th)), ]
  # nonzeros[nonzeros > 0] <- 1
  # nonzeros <- data.table(nonzeros, keep.rownames = "id")
  # nonzeros <- merge(ml, nonzeros, by = "id")
  nonzeros <- nonzeros[order(seqnames, start)]
  
  bins <- lapply(unique(nonzeros$seqnames), function(chr) unlist(lapply(1:ceiling(nrow(nonzeros[seqnames==chr])/50), rep, 50))[1:nrow(nonzeros[seqnames==chr])])
  names(bins) <- unique(nonzeros$seqnames)
  bins2 <- rbindlist(lapply(names(bins), function(chr) data.table(seqnames = chr, bin = bins[[chr]])))
  stopifnot(sum(bins2$seqnames != nonzeros$seqnames)==0)
  nonzeros2 <- cbind(bin = bins2$bin, nonzeros[,1:4])
  mybins <- nonzeros2
  mybins2 <- merge(mybins[, lapply(.SD, min), 
                          .SDcols = c("start"), by = c("seqnames", "bin")], 
                   mybins[, lapply(.SD, max), 
                          .SDcols = c("end"), by = c("seqnames", "bin")], 
                   by = c("seqnames", "bin"))
  
  mybins <- mybins2
  
  
  #merging:
  setkeyv(mybins, c("seqnames", "start", "end"))
  setkeyv(coding$A, c("seqnames", "start", "end"))
  setkeyv(coding$B, c("seqnames", "start", "end"))
  setkeyv(noncoding$A, c("seqnames", "start", "end"))
  setkeyv(noncoding$B, c("seqnames", "start", "end"))
  
  coding$bins.A <- foverlaps(x = coding$A, y = mybins, nomatch = NULL)
  coding$bins.B <- foverlaps(x = coding$B, y = mybins, nomatch = NULL)
  stopifnot(sum(coding$bins.A$id != coding$bins.B$id) == 0)
  
  noncoding$bins.A <- foverlaps(x = noncoding$A, y = mybins, nomatch = NULL)
  noncoding$bins.B <- foverlaps(x = noncoding$B, y = mybins, nomatch = NULL)
  stopifnot(sum(noncoding$bins.A$id != noncoding$bins.B$id) == 0)
  
  #coding:
  coding$cnts.A <- coding$bins.A[, lapply(.SD, function(x) sum(x, na.rm = T)), .SDcols = colnames(coding$A)[-c(1:4)], by = c("seqnames", "bin")]
  coding$cnts.B <- coding$bins.B[, lapply(.SD, function(x) sum(x, na.rm = T)), .SDcols = colnames(coding$B)[-c(1:4)], by = c("seqnames", "bin")]
  
  coding$aggs.A <- coding$bins.A[, lapply(.SD, function(x) mean(x, na.rm = T)), .SDcols = colnames(coding$A)[-c(1:4)], by = c("seqnames", "bin")]
  coding$aggs.B <- coding$bins.B[, lapply(.SD, function(x) mean(x, na.rm = T)), .SDcols = colnames(coding$B)[-c(1:4)], by = c("seqnames", "bin")]
  
  coding$aggs.A.zs <- cbind(coding$aggs.A[, c(1:2)], scale(coding$aggs.A[, -c(1:2)], center = T, scale = T))
  coding$aggs.B.zs <- cbind(coding$aggs.B[, c(1:2)], scale(coding$aggs.B[, -c(1:2)], center = T, scale = T))
  
  #Original cmd: hist((as.matrix(coding$aggs.B.zs[, controlSampleIDs, with=F])), breaks=33, col="seagreen")
  hist((as.matrix(coding$aggs.B.zs[, ..controlSampleIDs])), breaks=33, col="seagreen")
  
  if(saveMe)
    saveRDS(coding, file = sprintf("%s/aggregated_results/ASE.coding.rds", dirpath))
  
  #noncoding:
  noncoding$cnts.A <- noncoding$bins.A[, lapply(.SD, function(x) sum(x, na.rm = T)), .SDcols = colnames(noncoding$A)[-c(1:4)], by = c("seqnames", "bin")]
  noncoding$cnts.B <- noncoding$bins.B[, lapply(.SD, function(x) sum(x, na.rm = T)), .SDcols = colnames(noncoding$B)[-c(1:4)], by = c("seqnames", "bin")]
  
  
  noncoding$aggs.A <- noncoding$bins.A[, lapply(.SD, function(x) mean(x, na.rm = T)), .SDcols = colnames(noncoding$A)[-c(1:4)], by = c("seqnames", "bin")]
  noncoding$aggs.B <- noncoding$bins.B[, lapply(.SD, function(x) mean(x, na.rm = T)), .SDcols = colnames(noncoding$B)[-c(1:4)], by = c("seqnames", "bin")]
  
  noncoding$aggs.A.zs <- cbind(noncoding$aggs.A[, c(1:2)], scale(noncoding$aggs.A[, -c(1:2)], center = T, scale = T))
  noncoding$aggs.B.zs <- cbind(noncoding$aggs.B[, c(1:2)], scale(noncoding$aggs.B[, -c(1:2)], center = T, scale = T))
  
  if(saveMe)
    saveRDS(noncoding, file = sprintf("%s/aggregated_results/ASE.noncoding.rds", dirpath))
  
}

generate_adt_adt0_adt.na_and_nonzeros_data <- function(th = 4, dirpath = "/pellmanlab/stam_niko/data/processed_bam", scriptsdir = "/pellmanlab/nikos/Stam_Etai_Scripts", normBySd = F, datadir = "/pellmanlab/nikos/Stam_Etai_Data") {
  
  #source('~/WORK/secondaryanalysis/methods_paper_results/R/CNML.R')
  #load_control_data()
  source(sprintf('%s/scripts/fromRawTPMtoExprsRatio.R', scriptsdir))
  #source(sprintf('%s/scripts/fromRawTPMtoExprsRatio_by_limma.R', scriptsdir))
  rsemtpm <- readRDS(file = sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", dirpath))
  geneRanges <- readRDS(sprintf("%s/geneRanges.rds", datadir))
  controlSampleIDs2 <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs2.rds", dirpath))
  controlSampleIDs <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))
  
  #adt.na:
  #-------
  mlinput2.na <- with(configs, fromRawTPMtoExprsRatio(rsemtpm = rsemtpm, geneRanges = geneRanges, 
                                                      controlSampleIDs = controlSampleIDs, max_pos_col = 6, plusOne = 1, 
                                                      maxExp = maxExp, quantileOfExp = quantileOfExp, 
                                                      minDetectionLevel = minDetectionLevel,
                                                      zerosAsNA = T, minNumOfSamplesToDetect = minNumOfSamplesToDetect,  
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
  #stopifnot(sum(!adt.old$id %in% adt.na$id) == 0)
  #sum(!adt.old$id %in% adt.na$id)
  setkey(adt.na, id)
  #adt.na <- adt.na[adt.old$id]
  dim(adt.na)
  saveRDS(adt.na, file = sprintf("%s/aggregated_results/adt.na.rds", dirpath))
  
  #adt:
  #----
  mlinput <- with(configs, fromRawTPMtoExprsRatio(rsemtpm = rsemtpm, geneRanges = geneRanges, 
                                                  controlSampleIDs = controlSampleIDs, max_pos_col = 6, plusOne = 1, 
                                                  maxExp = maxExp, quantileOfExp = quantileOfExp, 
                                                  minDetectionLevel = minDetectionLevel,
                                                  zerosAsNA = F, minNumOfSamplesToDetect = minNumOfSamplesToDetect, #minNumOfSamplesToDetect, 
                                                  doNormalizeByControls = T, normBySd = normBySd))
  
  
  #remove unneed annotations
  adt <- mlinput$adt[,-c(5,6)]
  
  dim(adt)
  stopifnot(sum(!adt.na$id %in% adt$id) == 0)
  setkey(adt, id)
  
  adt <- adt[adt.na$id]
  dim(adt)
  saveRDS(adt, file = sprintf("%s/aggregated_results/adt.rds", dirpath))
  
  #adt0: - should exclude from analysis
  #-----
  
  # mlinput <- with(configs, fromRawTPMtoExprsRatio(rsemtpm = rsemtpm, geneRanges = geneRanges, 
  #                                                 controlSampleIDs = controlSampleIDs2, max_pos_col = 6, plusOne = 1, 
  #                                                 maxExp = maxExp, quantileOfExp = quantileOfExp, 
  #                                                 minDetectionLevel = minDetectionLevel, 
  #                                                 zerosAsNA = F, minNumOfSamplesToDetect = 3, #minNumOfSamplesToDetect, 
  #                                                 doNormalizeByControls = F))
  # 
  # 
  # adt <- mlinput$adt[,-c(5,6)]
  # dim(adt)
  # tpmdata <- fromRawTPMtoMLInput(rsemtpm = rsemtpm, geneRanges = geneRanges[,1:4], controlSampleIDs = controlSampleIDs2, max_pos_col = 4, plusOne = 1, min_sd = 0.0001,
  #                                maxExp = 1200, quantileOfExp = 0.99, minDetectionLevel = 5, minNumOfSamplesToDetect = 5, 
  #                                geneInterspaceLogBase = 8, doNormalizeByControls = F)
  # adt <- tpmdata$adt
  #IDs <- names(which(colSums(rsemtpm>0)> 2000))
  #saveRDS(adt, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/reports.CNV12.0/data/adt.rds")
  
  #Preparing data for CN calculations:
  # tmp <- log2(adt[,-(1:4)] + 1) #adt[,-(1:4)] #
  # M <- copy(tmp[, controlSampleIDs2, with = F] + 1)
  # q <- min(log2(1200), unname(quantile(as.matrix(M), 0.99)))
  # M[M > q] <- q
  # 
  # myMeans <- rowMeans(M, na.rm = T)
  # adt2 <- cbind(adt[,1:4], sweep(tmp, 1, myMeans, "-"))
  # adt <- adt2
  # adt <- adt[order(seqnames, start)]
  # require(gtools)
  # adt <- adt[mixedorder(seqnames)]
  # 
  # 
  #Total expression with +0:
  #-------------------------
  #tpmdata is the same for +1 or +0 since we do not normalize
  #now we add to log2 input + 1 works good when allele level is A and B and not An and Bn
  # adt0 <- tpmdata$adt
  # M <- copy(adt0[, controlSampleIDs2, with = F])
  # q <- min(1200, unname(quantile(as.matrix(M), 0.99)))
  # M[M > q] <- q
  # 
  # myMeans <- rowMeans(M, na.rm = T)
  # adt0 <- cbind(adt0[,1:4], log2(sweep(adt0[,-(1:4)], 1, myMeans, "-") + 1))
  # adt0 <- adt0[order(seqnames, start)]
  # 
  
  #nonzeros:
  #---------
  nonzeros <- adt.na[,-c(1:4)]
  nonzeros[!is.na(nonzeros)] <- 1
  nonzeros[is.na(nonzeros)] <- 0
  nonzeros <- cbind(adt.na[, 1:4], nonzeros)
  saveRDS(nonzeros, file = sprintf("%s/aggregated_results/nonzeros.rds", dirpath))

}

compute_nonZero_bins_zscore <- function(saveMe = F, dirpath = "/pellmanlab/stam_niko/data/processed_bam") {
  
  #input data
  nonzeros <- readRDS(sprintf("%s/aggregated_results/nonzeros.rds", dirpath))
  controlSampleIDs2 <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs2.rds", dirpath))
  controlSampleIDs <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))
  rsemtpm <- readRDS(file = sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", dirpath))
  
  #sample
  nonzeros <- nonzeros[order(seqnames, start)]
  
  #creates a column in nonzeros with the bin number. putting genes in groups of 50.
  bins <- lapply(unique(nonzeros$seqnames), function(chr) unlist(lapply(1:ceiling(nrow(nonzeros[seqnames==chr])/50), rep, 50))[1:nrow(nonzeros[seqnames==chr])])
  names(bins) <- unique(nonzeros$seqnames)
  bins2 <- rbindlist(lapply(names(bins), function(chr) data.table(seqnames = chr, bin = bins[[chr]])))
  stopifnot(sum(bins2$seqnames != nonzeros$seqnames)==0)
  nonzeros2 <- cbind(bin = bins2$bin, nonzeros)
  
  #Takes the mean fraction of expression across the 50 genes for each bin in each sample.
  nonzeros3 <- nonzeros2[, lapply(.SD, function(x) mean(x, na.rm = T)), .SDcols = c(colnames(rsemtpm)), by = c("seqnames", "bin")]
  nonzeros <- nonzeros3
  
  #scale sweeps through df and replaces fraction of expression with z-score. But this scales by column/sample.
  nonzeros.zs <- scale(nonzeros[, -c(1:2)], center = T, scale = T)
  nonzeros.zs <- cbind(nonzeros[, 1:2], nonzeros.zs)
  
  #subtract median of chr median fraction of expressed genes option.
  #nonzeros.zs <- t(t(nonzeros[,-c(1:2)]) - colMedians(as.matrix(nonzeros[, lapply(.SD, median, na.rm = T), .SDcols=-c(1:2), by = "seqnames"][,-c(1)]), na.rm = T))
  #nonzeros.zs <- cbind(nonzeros[,c(1:2)], nonzeros.zs)
  
  if(saveMe) {
    saveRDS(nonzeros, file = sprintf("%s/aggregated_results/nonzeros.bin50.rds", dirpath))
    saveRDS(nonzeros.zs, file = sprintf("%s/aggregated_results/nonzeros.zs.bin50.rds", dirpath))
  }
  
  require(matrixStats)
  hist(rowMedians(as.matrix(nonzeros[, controlSampleIDs, with=F])), breaks=23, col="seagreen")
  hist(rowMedians(as.matrix(nonzeros.zs[, controlSampleIDs, with=F])), breaks=23, col="seagreen")
  return(c(nonzeros, nonzeros.zs))
}


generate_RPE1_GeneSpecificData <- function(myMeans, controlSampleIDs, mySds, chr = NULL, base = 2, digits = 2, GC_length_rds_fname = NULL, dirpath = "/pellmanlab/stam_niko/data/processed_bam") {
  
  adt.na <- readRDS("/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/adt.na.rds")
  adt.raw.na <- readRDS("/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/adt.raw.na.rds")
  controlSampleIDs2 <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs2.rds", dirpath))
  controlSampleIDs <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))
  myids <- intersect(adt.na$id, adt.raw.na$id)
  setkeyv(adt.raw.na, cols = "id")
  adt <- adt.raw.na[myids, ]
  
  if(is.null(GC_length_rds_fname)) {
    glen <- fread("/singlecellcenter/etai/ReferenceData/Gencode/GC_lengths.tsv") #produced by: GTF2LengthGC.R in annotations dir - the sum of all exons' length
    setnames(glen, old = "V1", new = "id")
  } else {
    message("Reading GC and length from input file..")
    glen <- readRDS(GC_length_rds_fname)
    names(glen)[1] <- "id"
  }
  
  M <- as.matrix(adt[, controlSampleIDs, with=F])
  #mean without including zeros:
  mean.nz <- rowMeans(M, na.rm = T)
  sd.nz <- rowSds(M, na.rm = T)
  cv.nz <- sd.nz/mean.nz
  #frac of zeros:
  mean.z <- rowMeans(is.na(M))
  #mean including zeros:
  Mwz <- t(apply(M, 1, function(m) ifelse(is.na(m), min(m, na.rm = T), m)))
  mean.wz <- rowMeans(Mwz)
  sd.wz <- rowSds(Mwz)
  cv.wz <- sd.wz/mean.wz
  
  #interspace:
  interspace <- unlist(lapply(unique(adt$seqnames), function(chr) c(adt[seqnames == chr]$start[-1], data.table::last(adt[seqnames == chr]$end)+1) - adt[seqnames == chr]$end + 1.5))
  interspace[interspace <= 0] <- 1
  interspace <- round(log(x = interspace, base = base), digits = digits)
  
  message("Gene specific features: ", nrow(glen), " from tpm data: ", nrow(adt))
  setkey(glen, "id")
  glen <- glen[adt$id]
  message("Following reduction - Gene specific features: ", nrow(glen), " from tpm data: ", nrow(adt))
  
  geneSpecific <- data.table(adt[,1:4], 
                             glen[, 2:3],
                             mean.nz = mean.nz, sd.nz = sd.nz, cv.nz = cv.nz, 
                             zeroFrac = mean.z, 
                             mean.wz = mean.wz, sd.wz = sd.wz, cv.wz = cv.wz, 
                             interspace = interspace)
  # myTCV <- log2((mySds + 0.1)/(myMeans))
  # myTCV[which(myTCV < -5)] <- -5
  # myTCV[which(myTCV > 5)] <- 5
  # 
  # myTCV <- round(scale(myTCV, center = T, scale = T), digits = digits)
  # geneSpecific <- myTCV #myMeans => v7
  # #TODO: maybe use the new probability based 1-exp transformation function
  # 
  # 
  # geneSpecific <- list(TCV = myTCV, interspace = interspace, exonsLenSum = glen$Length, exonsGC = glen$GC)
  saveRDS(geneSpecific, file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/geneSpecific.rds")
  
  geneSpecific <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/geneSpecific.rds")
  
  myfuncs <- list(MEAN = function(x) mean(x, na.rm=TRUE), 
                  MEDIAN = function(x) median(x, na.rm=TRUE),
                  Q05 = function(x) quantile(x, probs = c(0.05), na.rm = TRUE),
                  Q25 = function(x) quantile(x, probs = c(0.25), na.rm = TRUE),
                  Q75 = function(x) quantile(x, probs = c(0.75), na.rm = TRUE),
                  Q95 = function(x) quantile(x, probs = c(0.95), na.rm = TRUE),
                  SD = function(x) sd(x, na.rm = TRUE)
  )
  
  winSize <- 50
  myfuncs$MAXSTRETCH <- NULL
  
  
  tmp <- lapply(names(myfuncs), function(fc) geneSpecific[, lapply(.SD, rollapply, width = winSize, FUN = myfuncs[[fc]], fill = NA, align = "center"),
                                                            .SDcols = names(geneSpecific)[-c(1:4)], by = seqnames])
  names(tmp) <- names(myfuncs)
  saveRDS(tmp, file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/geneSpecific.ftrs.winSize50.raw.rds")
  
  return(geneSpecific)
}

generate_ganno <- function(dirpath = "/pellmanlab/stam_niko/data/processed_bam") {
  ganno <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/adt.na.rds")[,1:4]
  arms <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/arms.rds")
  
  setkeyv(ganno, c("seqnames", "start", "end"))
  setkeyv(arms, c("seqnames", "start", "end"))
  
  ganno <- foverlaps(ganno, arms)
  setnames(ganno, old = c("start", "end", "i.start", "i.end"), new = c("arm_start", "arm_end", "start", "end"))
  saveRDS(ganno, file = sprintf("%s/data/ganno.rds", dirpath))
}

generate_ganno2 <- function(dirpath = "/pellmanlab/stam_niko/data/processed_bam", datapath = "/pellmanlab/nikos/Stam_Etai_Data") {
  ### revised by nikos ###
  adt.na_anno <- readRDS(sprintf("%s/aggregated_results/adt.na.rds", dirpath))[,1:4]
  armRanges <- readRDS( sprintf("%s/armRanges.rds", datapath) )
  
  setkeyv(adt.na_anno, c("seqnames", "start", "end"))
  setkeyv(armRanges, c("seqnames", "start", "end"))
  
  ganno <- foverlaps(adt.na_anno, armRanges)
  setnames(ganno, old = c("start", "end", "i.start", "i.end"), new = c("arm_start", "arm_end", "start", "end"))
  saveRDS(ganno, file = sprintf("%s/aggregated_results/ganno.rds", dirpath))
}
