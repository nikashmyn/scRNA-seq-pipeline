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
  
  #cap tpm expression 
  dt2go <- copy(dt[, -c(1:max_pos_col), with=F])
  q <- min(maxExp, unname(quantile(as.matrix(dt2go), quantileOfExp, na.rm=T)))
  message("Quantile q = ", q)
  dt2go[dt2go>q] <- q
  
  #logCPM by limma::voom
  voom_tpm <- limma::voom(dt2go) #, normalize.method = "cyclicloess")

  #log before mean:
  M <- voom_tpm[["E"]]
  
  M2 <- copy(M[, c(controlSampleIDs)])
  if(zerosAsNA) {
    M2[M2<0] <- NA 
  } 

  message("NAs = ", sum(is.na(M2)))
  #stats if doRatioBeforeScaleAndZscore == F
  means <- rowMeans(M2, na.rm = T)
  sds <- rowSds(as.matrix(M2), na.rm = T)
  
  #calculating final dt:
  #subtracting the mean centers the data. Or you can divide by the std to standardize.
  if(doNormalizeByControls) {
    dt <- cbind(dt[, 1:max_pos_col, with = F], sweep(M, 1, means, "-") )
  }
  if(normBySd) {
    dt <- cbind(dt[, 1:max_pos_col, with = F], sweep(M, 1, sds, "/")) #* sqrt(mean(sds^2, na.rm=T)) )
  }
  if(doNormalizeByControls) {
    if(normBySd) {
      dt <- cbind(dt[, 1:max_pos_col, with = F], sweep(M, 1, means, "-") )
      dt <- cbind(dt[, 1:max_pos_col, with = F], sweep(dt[, -c(1:max_pos_col), with = F], 1, sds, "/")) 
    }
  }
  if(!doNormalizeByControls) {
    if(!normBySd) {
      dt <- cbind(dt[, 1:max_pos_col, with = F], M)
    }
  }
  
  #chrs <- as.character(dt$seqnames)
  #chrs_red <- chrs[!duplicated(chrs)]
  #dt_norm <- c()
  #for (i in 1:length(chrs_red)) {
  #  chr <- chrs_red[i]
  #  dt_chr <- dt[which(chrs %in% chr),]
  #  dt_chrnorm <- scale(dt_chr[, -c(1:max_pos_col), with = F], center = T, scale = F)
  #  dt_chrnorm <- cbind(dt_chr[, c(1:max_pos_col), with = F], dt_chrnorm)
  #  dt_norm <- rbind(dt_norm, dt_chrnorm)
  #}
  #dt <- dt_norm
  res <- list(adt = dt, means = means, sds = sds, controlSampleIDs = controlSampleIDs)
  return(res)
}