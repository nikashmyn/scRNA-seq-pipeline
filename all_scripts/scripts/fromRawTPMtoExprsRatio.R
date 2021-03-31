#the same as "fromRawTPMtoMLInput" but mean on log space and not mean before log transform
fromRawTPMtoExprsRatio <- function(rsemtpm, geneRanges, controlSampleIDs, max_pos_col = 6, plusOne = 1, min_sd = 0.0001, 
                                 zerosAsNA = F, normBySd = F,
                                 maxExp = 1200, quantileOfExp = 0.99, minDetectionLevel = 5, minNumOfSamplesToDetect = 3, 
                                 geneInterspaceLogBase = 8, interspaceCoef = 5, doNormalizeByControls = F) {
  require(data.table)
  require(matrixStats)
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