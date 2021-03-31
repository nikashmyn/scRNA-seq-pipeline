

getMyGenomicRangesClusters <- function(id, myclusts) {
  require(GenomicRanges)
  getMyDTClusters <- function(chr, x) {
    rr <- rle(myclusts[seqnames==chr][[x]])
    dtr <- (data.table(seqnames = chr, 
                       start = myclusts[seqnames == chr]$start[c(1,cumsum(rr$lengths[-length(rr$length)]))],
                       end = myclusts[seqnames == chr]$end[cumsum(rr$lengths)],
                       CN = rr$values,
                       diff.fwd = c(0,diff(rr$values)),
                       diff.rev = c(0,diff(rev(rr$values)))))
    return(dtr)
  }
  
  gr <- GRanges(rbindlist(lapply(unique(myclusts$seqnames), getMyDTClusters, id)))
  return(gr)
}

produceSegmentsFromRndLocations <- function(mydt, IDs, chrs = NA, win_sizes = seq(3, 100, 10), N = 30) {
  
  getMyRnd <- function(chr, win_size) {
    rng <- range(which(mydt$seqnames == chr))
    
    wi <- which(mydt$seqnames == chr)
    if(win_size > length(wi)) {
      message("Setting win size to ", length(wi) - 2, " since it is too big..")
      win_size <- length(wi) - 2 
    }
    z <- ceiling(((win_size - 1)/2))
    
    wi <- wi[-c(1:(z), (length(wi)-z + 1):length(wi))]
    
    L <- length(wi)
    if(N > L) {
      message("Sampling cannot be bigger than chromosome length.\nAdjusting to chromosome size*0.8 (", L, ")")
      N <- round(L*0.8)
    }
    
    rnd <- sample(x = wi, size = N, replace = F)  
    rnd <- data.table(from = rnd - z, to = rnd + z)
    return(rnd)
  }
  
  getMyRndSegMeans <- function(i, win_size) {
    if(length(chrs) == 1) {
      rnds <- getMyRnd(chr = chrs, win_size = win_size) 
    } else {
      rnds <- rbindlist(lapply(chrs, getMyRnd, win_size))  
    }
    return(rowMeans(t(apply(rnds, 1, function(x) mydt[[i]][x[1]:x[2]]))))
  }
  
  getMyRndVariousSegMeans <- function(win_size) {
    message(win_size)
    res <- unlist(lapply(IDs, getMyRndSegMeans, win_size))
    return(res)
  }
  
  if(is.na(chrs[1])) {
    chrs <- as.character(unique(mydt$seqnames))
  } else {
    chrs <- chrs[chrs %in% as.character(unique(mydt$seqnames))]
  }
  
  tt <- lapply(win_sizes, getMyRndVariousSegMeans)
  names(tt) <- win_sizes
  
  return(tt)
  
}


getMyChrPreds <- function(MAT, method_name, maxn = 4, IDs = tri12ids, excludeChrs = c("chrX", "chr10")) {
  M <- MAT[[method_name]]
  dtm <- M[, lapply(.SD, median),
           .SDcols = colnames(M)[-c(1:maxn)], by = "seqnames"]
  
  data.table(seqnames = dtm$seqnames, 2^sweep(as.matrix(dtm[,IDs, with=F]), 2, colMedians(as.matrix(dtm[!(seqnames %in% excludeChrs),IDs, with=F])), "-"))
}


getMyGenePreds <- function(MAT, method_name, maxn = 4, IDs = tri12ids, excludeChrs = c("chrX", "chr10"), dfc = NA, eps = 0.01) {
  M <- MAT[[method_name]]
  dtm <- M[, lapply(.SD, median),
           .SDcols = colnames(M)[-c(1:maxn)], by = "seqnames"]
  
  mydt <- data.table(2^sweep(as.matrix(M[,IDs, with=F]), 2, colMeans(as.matrix(dtm[!(seqnames %in% excludeChrs),IDs, with=F])), "-"))
  
  
  if(!is.na(dfc[1]))
    mydt <- sweep(mydt, 1, rowMeans(as.matrix(dfc[[method_name]][,-c(1:maxn), with=F])) + eps, "/")
  mydt[mydt>4] <- 4
  return(data.table(M[, 1:maxn, with=F], mydt))
}



runPredOnIDsMultiCores <- function(IDs, dt2, alleles.all, CPUs = 26) {
  require(data.table)
  require(KFAS)
  require(signal)
  require(snow)
  cl <- makeCluster(CPUs, type = "SOCK")
  #snow::clusterExport(cl = cl, list = c(names(environment())), envir = environment())
  snow::clusterExport(cl,c("dt2", "alleles.all"), envir = environment())
  clusterEvalQ(cl, expr = source("~ejacob/WORK/secondaryanalysis/Stamatis/G2/runRNAseqPipelineFor_BF9VP.utils.R"))
  myres <- parLapply(cl, IDs, 
                     function(i) {
                       message(i)
                       myalleles <- mergeTotalExpAndAlleleData(dt2, A = alleles.all$An, B = alleles.all$Bn, id = i)
                       chrs <- as.character(unique(myalleles$seqnames))
                       mychrres <- lapply(chrs, 
                                          function(chr) {
                                            message(chr)
                                            result <- runMyMedianAndSSMFiltering(x = dt2[seqnames == chr][[i]], 
                                                                                 x2 = dt2[seqnames == chr]$interspace, 
                                                                                 A = alleles.all$An[seqnames == chr][[i]], B = alleles.all$Bn[seqnames == chr][[i]], abpos = alleles.all$An[seqnames == chr]$start, 
                                                                                 myalleles = myalleles[seqnames == chr],
                                                                                 t = dt2[seqnames == chr]$start, 
                                                                                 main = i, plotme = F)
                                            return(result)
                                          })
                       names(mychrres) <- chrs
                       return(mychrres)
                     })
  
  names(myres) <- IDs
  stopCluster(cl)
  return(myres)
}



getBiallelicPval.binomial <- function(an, bn) {
  if(an + bn == 0)
    return(1)
  return(binom.test(x = an, n = an + bn , p = 0.5, alternative = "two.sided")$p.value)
}

generateExpressionBasedSNPFeatures <- function(m, myIDs, cpus = 44) {
  require(data.table)
  message("Preparing data..")
  A <- copy(m$A)
  B <- copy(m$B)
  features <- copy(m$features)
  
  A2 <- A[, myIDs, with=F]
  B2 <- B[, myIDs, with=F]
  A2[is.na(A2)] <- 0
  B2[is.na(B2)] <- 0
  message("Calculating features based on expression..")
  AN <- rowSums(A2 > B2)
  BN <- rowSums(B2 > A2)
  
  features$above0.A <- rowSums(A2>0)
  features$above0.B <- rowSums(B2>0)
  features$AgB <- AN
  features$BgA <- BN
  
  am <- rowMeans(A[, myIDs, with=F], na.rm = T)
  bm <- rowMeans(B[, myIDs, with=F], na.rm = T)
  features$avg.A <- am
  features$avg.B <- bm
  
  return(features)
  
  message("Executing binomial test on ", cpus, " clusters..")
  require(snow)
  cl <- makeCluster(3, type = "SOCK")
  #snow::clusterExport(cl = cl, list = c(names(environment())), envir = environment())
  snow::clusterExport(cl,c("getBiallelicPval.binomial", "AN", "BN"), envir = environment())
  mypvals0.binom0 <- parSapply(cl, 1:1000, function(i) { 
    getBiallelicPval.binomial(a = AN[i], b = BN[i]) })
  
  # mypvals0.binom0 <- parSapply(cl, 1:length(AN), function(i) { 
  #   getBiallelicPval.binomial(a = AN[i], b = BN[i]) })
  features$expBiallelicPvals <- mypvals0.binom0
  stopCluster(cl)
  
  message("Task completed.")
  return(features)
}

phaseMyAllelesBasedOnHaplotypes <- function(m, th = 0) {
  require(data.table)
  require(pryr)
  a <- copy(m$A)
  b <- copy(m$B)
  c <- copy(m$A)
  
  wi <- which(m$features$allele < th)
  a[wi, ] <- b[wi, ]
  b[wi, ] <- c[wi, ]
  a <- a[-which(m$features$allele == 0)]
  b <- b[-which(m$features$allele == 0)]
  features <- m$features[-which(m$features$allele == 0)]
  return(list(A = a, B = b, features = features))
}

collectAllVCFdata <- function(patternFile = "/singlecellcenter/etai/ExperimentsData/Stamatis/May072018/vcfs/alleles.ADs.v5.chr%s.rds",
                              mychrs = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22, "X")) {
  M <- lapply(sprintf(patternFile, mychrs), function(f) { message("Working on file: ", f); return(readRDS(f)) })
  names(M) <- sprintf("chr%s", mychrs)
  A <- rbindlist(lapply(M, function(x) data.table(x$A)))
  B <- rbindlist(lapply(M, function(x) data.table(x$B)))
  features <- rbindlist(lapply(M, function(x) data.table(x$features)))
  
  return(list(A = A, B = B, features = features))
}

loadAndPhaseSNPdataFromRPE1ExperimentsBasedOnHaplotypeAndExpression <- function(out_fname = "/singlecellcenter/etai/ExperimentsData/Stamatis/180803_M01209_0114_000000000-BF9VP/reports/haplotypePhasedADs.rds") {
  require(readxl)
  anno <- readRDS("~ejacob/WORK/secondaryanalysis/Stamatis/reports/Stamatis_list_v12_180404.QCRNAv3.annoqc.rds")
  anno1 <- read_excel(path = "/singlecellcenter/etai/ExperimentsData/Stamatis/May072018/Stamatis_list_v12_180508.QCRNAv3.withQCRNAv3.xlsx", sheet = "main", na = c("", "NA"))
  anno2 <- data.table(read_excel(path = "/singlecellcenter/etai/ExperimentsData/Stamatis/May072018/Stamatis_list_v12_180508.QCRNAv3.withQCRNAv3.xlsx", sheet = "all", na = c("", "NA")))
  myIDs <- anno2[cellType1 == "RPE1" & th1 > 6000 & totalReads > 200000 & (control_group != "bunch_of_grapes" | is.na(control_group))]$WTA.plate
  
  #collect SNP data from all chrs except Y and M: 
  m <- collectAllVCFdata(patternFile = "/singlecellcenter/etai/ExperimentsData/Stamatis/May072018/vcfs/alleles.ADs.v5.chr%s.rds")
  A <- m$A
  B <- m$B
  features <- m$features
  
  
  #modifying column names of Hiseq run (the names of columns 1-382 are already modified):
  colnames(A)[-c(1:382)] <- anno$WTA.plate[sapply(colnames(A)[-c(1:382)], grep, x = anno$Fastq_files3, value = F)]
  colnames(B)[-c(1:382)] <- anno$WTA.plate[sapply(colnames(B)[-c(1:382)], grep, x = anno$Fastq_files3, value = F)]
  #modifying column names of Miseq run
  colnames(A)[grep("-", colnames(A), value = F)] <- gsub("-", "_", gsub(pattern = "_S\\d+_L001_R1_001", "", grep("-", colnames(A), value = T)))
  colnames(B)[grep("-", colnames(B), value = F)] <- gsub("-", "_", gsub(pattern = "_S\\d+_L001_R1_001", "", grep("-", colnames(B), value = T)))
  
  #from this point A is the ref allele read counts, B is B allele read count and features is the SNP information
  
  #Phasing - DNA based haplotypes:
  m <- phaseMyAllelesBasedOnHaplotypes(m = list(A = A, B = B, features = features))
  #saveRDS(m, file = out_fname)
  #Phasing - RNA based expression:
  expFeatures <- generateExpressionBasedSNPFeatures(m = m, myIDs = myIDs)
  
  fracPerChr <- rbindlist(lapply(unique(expFeatures$seqnames), function(chr) {
                        data.table(seqnames = chr, 
                                   frac1.5 = round(100*sum(with(expFeatures[ seqnames == chr & above0.A > 2 & above0.B > 2], abs(log2((AgB + 1) / 
                                                 (BgA + 1))) > log2(1.5)))/nrow(expFeatures[ seqnames == chr & above0.A > 2 & above0.B > 2]), digits = 1),
                                   frac2 = round(100*sum(with(expFeatures[ seqnames == chr & above0.A > 2 & above0.B > 2], abs(log2((AgB + 1) / 
                                                 (BgA + 1))) > log2(2)))/nrow(expFeatures[ seqnames == chr & above0.A > 2 & above0.B > 2]), digits = 1),
                                   frac3 = round(100*sum(with(expFeatures[ seqnames == chr & above0.A > 2 & above0.B > 2], abs(log2((AgB + 1) / 
                                                 (BgA + 1))) > log2(3)))/nrow(expFeatures[ seqnames == chr & above0.A > 2 & above0.B > 2]), digits = 1),
                                   frac4 = round(100*sum(with(expFeatures[ seqnames == chr & above0.A > 2 & above0.B > 2], abs(log2((AgB + 1) / 
                                                 (BgA + 1))) > log2(4)))/nrow(expFeatures[ seqnames == chr & above0.A > 2 & above0.B > 2]), digits = 1),
                                   frac5 = round(100*sum(with(expFeatures[ seqnames == chr & above0.A > 2 & above0.B > 2], abs(log2((AgB + 1) / 
                                                 (BgA + 1))) > log2(5)))/nrow(expFeatures[ seqnames == chr & above0.A > 2 & above0.B > 2]), digits = 1))
                      }))
  
  #The largest threshold (i.e. accounting for as much alleles as possible in the analysis) of imblanced bi-allelic expression that can still detect the 10q gain is frac3 based on the "fracPerChr" analysis
  m$imbalancedFractionOfSNPs <- fracPerChr
  m$refIDs <- myIDs
  m$features <- expFeatures
  #saveRDS(m, file = out_fname)
  
  ########################################################################
  #filtering out all non biallelic expressed genes except in 10q and chrX
  ########################################################################
  #We take any monoallelic expressed SNP from chrX:
  wi.chrX <- which(m$features$seqnames == "chrX" & (m$features$above0.A > 2 | m$features$above0.B > 2))
  #we take all 10q SNP with very high allelic expression imbalance (frac5) since we know this arm has 3 copies in RPE-1
  start10q <- max(centromeres[seqnames == "chr10"]$end)
  #The fraction of 10q and 10p above th frac6 (3.5% vs. 2.1%)
  length(with(m$features, which(seqnames == "chr10" & start > start10q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(6))))/length(with(m$features, which(seqnames == "chr10" & start > start10q & (above0.A > 2 & above0.B > 2) )))
  length(with(m$features, which(seqnames == "chr10" & start < start10q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(6))))/length(with(m$features, which(seqnames == "chr10" & start < start10q & (above0.A > 2 & above0.B > 2) )))
  #The fraction of 10q and 10p above th frac5 (4.5% vs. 2.4%)
  length(with(m$features, which(seqnames == "chr10" & start > start10q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(5))))/length(with(m$features, which(seqnames == "chr10" & start > start10q & (above0.A > 2 & above0.B > 2) )))
  length(with(m$features, which(seqnames == "chr10" & start < start10q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(5))))/length(with(m$features, which(seqnames == "chr10" & start < start10q & (above0.A > 2 & above0.B > 2) )))
  #The fraction of 10q and 10p above th frac3 (15.5% vs. 6%)
  length(with(m$features, which(seqnames == "chr10" & start > start10q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(3))))/length(with(m$features, which(seqnames == "chr10" & start > start10q & (above0.A > 2 & above0.B > 2) )))
  length(with(m$features, which(seqnames == "chr10" & start < start10q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(3))))/length(with(m$features, which(seqnames == "chr10" & start < start10q & (above0.A > 2 & above0.B > 2) )))
  
  #For instance, in chr1 (1p and 1q) yo do not see this bias as in chr10
  start1q <- max(centromeres[seqnames == "chr1"]$end)
  #The fraction of 10q and 10p above th frac6 (2.5% vs. 4.4%)
  length(with(m$features, which(seqnames == "chr1" & start > start1q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(6))))/length(with(m$features, which(seqnames == "chr1" & start > start1q & (above0.A > 2 & above0.B > 2) )))
  length(with(m$features, which(seqnames == "chr1" & start < start1q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(6))))/length(with(m$features, which(seqnames == "chr1" & start < start1q & (above0.A > 2 & above0.B > 2) )))
  #The fraction of 10q and 10p above th frac5 (3.2% vs. 5.1%)
  length(with(m$features, which(seqnames == "chr1" & start > start1q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(5))))/length(with(m$features, which(seqnames == "chr1" & start > start1q & (above0.A > 2 & above0.B > 2) )))
  length(with(m$features, which(seqnames == "chr1" & start < start1q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(5))))/length(with(m$features, which(seqnames == "chr1" & start < start1q & (above0.A > 2 & above0.B > 2) )))
  #The fraction of 10q and 10p above th frac3 (7.5% vs. 10%)
  length(with(m$features, which(seqnames == "chr1" & start > start1q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(3))))/length(with(m$features, which(seqnames == "chr1" & start > start1q & (above0.A > 2 & above0.B > 2) )))
  length(with(m$features, which(seqnames == "chr1" & start < start1q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(3))))/length(with(m$features, which(seqnames == "chr1" & start < start1q & (above0.A > 2 & above0.B > 2) )))
  
  #Therefore, we will use frac5 as a threshold for 10q
  wi.chr10q <- with(m$features, which(seqnames == "chr10" & start > start10q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) <= log2(5)))
  
  #chr12 has a high frequency of trisomies, therefore, we will use a higher threshold of frac4 for filtering out those snps:
  wi.chr12 <- with(m$features, which(seqnames == "chr12" & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) <= log2(4)))
  wi.all <- with(m$features, which((above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) <= log2(3)))
  
  #and now we union all SNPs
  wi <- unique(c(wi.all, wi.chr10q, wi.chr12, wi.chrX))
  m$features$biAllelicFilter <- FALSE
  m$features$biAllelicFilter[wi] <- TRUE
  
  #Coding SNPs include also 5` and 3` UTRs:
  with(m$features, which(biAllelicFilter == T & (!is.na(coding) | (!is.na(fiveUTR)) | (!is.na(threeUTR)))))
  saveRDS(m, file = out_fname)
  return(m)
}

#TODO: edit - take relevant info for SSM
generateGeneSpecificFtrs <- function(adt, myMeans, mySds, chr = NULL, base = 2, digits = 2) {
  glen <- fread("/singlecellcenter/etai/ReferenceData/Gencode/GC_lengths.tsv") #produced by: GTF2LengthGC.R in annotations dir - the sum of all exons' length
  setnames(glen, old = "V1", new = "id")
  setkey(glen, "id")
  glen <- glen[adt$id]
  
  myTCV <- log2((mySds + 0.1)/(myMeans))
  myTCV[which(myTCV < -5)] <- -5
  myTCV[which(myTCV > 5)] <- 5
  
  myTCV <- round(scale(myTCV, center = T, scale = T), digits = digits)
  geneSpecific <- myTCV #myMeans => v7
  
  interspace <- unlist(lapply(unique(adt$seqnames), function(chr) c(adt[seqnames == chr]$start[-1], dplyr::last(adt[seqnames == chr]$end)+1) - adt[seqnames == chr]$end + 1.5))
  interspace[interspace <= 0] <- 1
  interspace <- round(log(x = interspace, base = base), digits = digits)
  
  geneSpecific <- list(TCV = myTCV, interspace = interspace, exonsLenSum = glen$Length, exonsGC = glen$GC)
  return(geneSpecific)
}

#m - allele level for all backgeound data, #mym - allele level for lib specific data
prepareAlleleData <- function(m, mym = NA, max.th = 120, useFraction = F, controlSampleIDs = NA,
                              geneInterspaceLogBase = 8, interspaceCoef = 5) {
  
  if(!is.na(mym)) {
    stopifnot(sum(m$features$id != mym$features$id) == 0)
    mym2 <- list()
    mym2$features <- m$features
    mym2$A  <- cbind(m$A, mym$A)
    mym2$B  <- cbind(m$B, mym$B)
  } else {
    mym2 <- m
  }
  wi <- which(with(mym2$features, (fiveUTR == T | threeUTR == T | coding == T) & biAllelicFilter == T))
  mym2$A <- mym2$A[wi,]
  mym2$B <- mym2$B[wi,]
  mym2$features <- mym2$features[wi,]
  #checks:
  colSums(mym2$A>0, na.rm = T)
  colSums(mym2$B>0, na.rm = T)
  sum(mym2$A, na.rm = T)
  sum(mym2$B, na.rm = T)

  mym2$A[is.na(mym2$A)] <- 0
  mym2$B[is.na(mym2$B)] <- 0
  mym2$A[mym2$A > max.th] <- max.th
  mym2$B[mym2$B > max.th] <- max.th
  
  A <- mym2$A
  B <- mym2$B
  
  if(useFraction) {
    Ac <- round(A/(A+B), digits = 2)
    Ac[is.na(Ac)] <- 0
    Bc <- round(B/(A+B), digits = 2)
    Bc[is.na(Bc)] <- 0
  } else {
    #(1) TPM like
    Ac <- scale(A, center=FALSE, scale=(colSums(A + B)+1))*1e6 + 1
    Bc <- scale(B, center=FALSE, scale=(colSums(B + B)+1))*1e6 + 1
    #(2) and log:
    #Ac <- log2(Ac + 1)
    #Bc <- log2(Bc + 1)
  }
  
  #interspace:
  interspace <- unlist(lapply(unique(mym2$features$seqnames), 
                              function(chr) c(mym2$features[seqnames == chr]$start[-1], 
                                              dplyr::last(mym2$features[seqnames == chr]$end)+1) - mym2$features[seqnames == chr]$end + geneInterspaceLogBase))
  interspace[interspace <= 0] <- 1
  interspace <- round(log(x = interspace, base = geneInterspaceLogBase), digits = 2)
  interspace <- exp(-interspaceCoef*interspace/max(interspace))
  
  A <- cbind(mym2$features[, 1:4], interspace, Ac)
  B <- cbind(mym2$features[, 1:4], interspace, Bc)
  setnames(x = A, old = "id", "snp_id")
  setnames(x = B, old = "id", "snp_id")
  
  An <- Bn <- NA
  if(!is.na(controlSampleIDs[1]) & !useFraction) {
    
    alleleNorm.A <- rowMeans(as.matrix(A[, controlSampleIDs, with=F]))
    alleleNorm.B <- rowMeans(as.matrix(B[, controlSampleIDs, with=F]))
    An <- cbind(A[, 1:5], log2(sweep(A[, -c(1:5)], 1, alleleNorm.A, "/") ))
    Bn <- cbind(B[, 1:5], log2(sweep(B[, -c(1:5)], 1, alleleNorm.B, "/") ))
    A <- cbind(A[, 1:5], log2(A[, -c(1:5)]))
    B <- cbind(B[, 1:5], log2(B[, -c(1:5)]))
  }
  
  return(list(A = A, B = B, An = An, Bn = Bn))
}

prepareTotalExpData <- function(rsemtpm, geneRanges, controlSampleIDs, max_pos_col = 6,
                                maxExp = 1200, quantileOfExp = 0.99, minDetectionLevel = 5, minNumOfSamplesToDetect = 3, 
                                geneInterspaceLogBase = 8, interspaceCoef = 5, doNormalizeByControls = T) {
  require(data.table)
  require(matrixStats)
  
  wi <- which(rowSums(rsemtpm[, controlSampleIDs] >= minDetectionLevel) >= minNumOfSamplesToDetect)
  dt <- data.table(rsemtpm[wi,], keep.rownames = "id")
  if(length(names(which(table(colnames(dt))>1))) > 0)
    dt <- dt[, -which(colnames(dt) %in% names(which(table(colnames(dt))>1))), with=F]
  setkey(geneRanges, id)
  setkey(dt, id)
  dt <- merge(geneRanges, dt)
  dt <- dt[order(seqnames, start, end)]
  dt <- dt[seqnames != "chrM" & seqnames != "chrY"]
  
  M <- copy(dt[, controlSampleIDs, with = F] + 1)
  q <- min(maxExp, unname(quantile(as.matrix(M), quantileOfExp)))
  M[M > q] <- q
  dt2go <- copy(dt[, -c(1:max_pos_col), with=F])
  dt2go[dt2go>q] <- q
  
  normalizationFactor <- rowMeans(M)
  
  #gene interspace:
  interspace <- unlist(lapply(unique(dt$seqnames), 
                              function(chr) c(dt[seqnames == chr]$start[-1], dplyr::last(dt[seqnames == chr]$end)+1) - dt[seqnames == chr]$end + geneInterspaceLogBase))
  interspace[interspace <= 0] <- 1
  interspace <- round(log(x = interspace, base = geneInterspaceLogBase), digits = 2)
  interspace <- exp(-interspaceCoef*interspace/max(interspace))
  #completing:
  if(doNormalizeByControls) {
    dt <- cbind(dt[, 1:max_pos_col, with = F], interspace, log2(sweep(dt2go + 1, 1, normalizationFactor, "/") ))
  } else {
    dt <- cbind(dt[, 1:max_pos_col, with = F], interspace, dt2go)
  }
  
  return(dt)
}

mergeTotalExpAndAlleleData <- function(dt, A, B, id = "170320_A7") {
  require(data.table)
  if("id" %in% names(dt))
    setnames(dt, old = "id", new = "gene_id")
  TE <- data.table(dt[, 1:4], TE = dt[[id]])
  setkey(TE, seqnames, start, end)
  AB <- data.table(A[, 1:4], A = A[[id]], B = B[[id]])
  
  setkey(AB, seqnames, start, end)
  abt <- foverlaps(TE, AB, nomatch = NA)
  abt$start <- NULL
  abt$end <- NULL
  setnames(abt, old = c("i.start", "i.end"), c("start", "end"))
  abt <- abt[order(seqnames, start, end)]
  aggabt <- abt[, lapply(.SD, function(x) round(mean(x, na.rm=TRUE), digits = 2)), by=c("gene_id"), .SDcols=c("A", "B", "TE") ] 
  aggabt <- merge(na.omit(aggabt), TE[,1:4], by = "gene_id")
  aggabt <- aggabt[order(seqnames, start, end)]
  return(aggabt)
}

getAllMatricesForItrMedFilter <- function(myres.all) {
  params <- names(myres.all[[1]][[1]])
  
  myfunc3 <- function(param_name) { 
    myfunc2 <- function(sample_id) {
      myfunc <- function(y) lapply(myres.all[[sample_id]], function(x) data.table(x[[y]]))
      xx <- myfunc(param_name)
      xxx <- rbindlist(lapply(names(xx), function(x) data.table(seqnames = x, xx[[x]])))
      setnames(xxx, old = "V1", new = sample_id)
      return(xxx[,2])
    }
    
    
    message(param_name)
    result <- do.call(cbind, lapply(names(myres.all), myfunc2))
    return(result)
  }
  mats <- lapply(params, myfunc3) 
  names(mats) <- params
  return(c(mats))
}

runIterativeMedianFilterOnAllSamples <- function(dt2, maxn = 7, CPUs = 26, mylevels = 3:12, win_size = 15) {
  require(data.table)
  require(KFAS)
  require(signal)
  require(snow)
  cl <- makeCluster(CPUs, type = "SOCK")
  IDs <- colnames(dt2)[-c(1:maxn)]
  chrs <- as.character(unique(dt2$seqnames))
  #snow::clusterExport(cl = cl, list = c(names(environment())), envir = environment())
  clusterEvalQ(cl, library(data.table))
  clusterEvalQ(cl, library(signal)) 
  snow::clusterExport(cl,c("dt2", "IDs", "chrs", "mylevels", "win_size"), envir = environment())
  #clusterEvalQ(cl, expr = source("library(data.table)"))
  snow::clusterEvalQ(cl, expr = source("~ejacob/WORK/secondaryanalysis/Stamatis/G2/runRNAseqPipelineFor_BF9VP.utils.R"))
  myres <- parLapply(cl, IDs, 
                     function(i) {
                       message(i)
                       mychrres <- lapply(chrs, 
                                          function(chr) {
                                            message(chr)
                                            result <- runIterativeMedianFilter(x = dt2[seqnames == chr][[i]], min_win = win_size, levels = mylevels)
                                            return(result)
                                          })
                       names(mychrres) <- chrs
                       return(mychrres)
                     })
  
  names(myres) <- IDs
  stopCluster(cl)
  return(myres)
}

runIterativeMedianFilter <- function(x, t = NA, min_win = 15, levels = 3, plotMe = F) {
  require(signal)
  m <- 0
  if(plotMe)
    plot(t, x)
  xx <- list()
  for(i in 1:max(levels)) {
    if(plotMe)
      message(min_win*i + m)
    x <- signal::filter(filt = MedianFilter(min_win*i + m), x = x)
    if(i %in% levels)
      xx[[sprintf("level%d", i)]] <- x
    m <- 1
  }
  if(plotMe) {
    lines(t, x, lwd = 3, col="red")
    abline(v = 65600208, lty = "dashed", col = "green", lwd=2)
  }
  
  return(xx)
}
runMyMedianAndSSMFiltering <- function(x = NA, x2 = NA, A = NA, B = NA, abpos = NA, myalleles = NA, t = NA, 
                                       chr = "chr5", main = "", xlab = "", ylim = c(0,3), plotme = T) {
  require(matrixStats)
  require(signal)
  require(ggplot2)
  require(KFAS)
  
  res <- list() 
  res$te.t <- t
  res$te.raw <- x
  if(plotme) {
    par(mfrow=c(1,1))
    plot(t, x, type = "l", col = alpha(colour = "black", alpha = 0.5), main = main, xlab = xlab, ylab = "Signal output", ylim = ylim)
  }
  # 7-point filter
  #lines(t, filter(MedianFilter(11), x), col = "blue", lwd=2) # another way to call it
  # 7-point and 15 points recursive filter
  res$te.rec_median <- signal::filter(MedianFilter(71), signal::filter(MedianFilter(41), signal::filter(MedianFilter(15), x)))
  if(plotme) {
    lines(t, res$te.rec_median, col = "darkgreen", lwd=4)
  }
  #SSM:
  
  if(is.na(x2[1])) {
    res$te.ssm <- getMyGSSM(x)
    if(plotme) {
      lines(t, res$te.ssm, col = "green", lwd=3)  
    }
    # legend("topright", inset=.05, title=NULL, cex = 0.8,
    #        legend = c("TPM ratio signal", "Median Smooth", "Recursive median Smooth", "Simple SSM"),
    #        fill=c(alpha(colour = "black", alpha = 0.5), "blue", "green", "red"), horiz=F)
  } else {
    res$te.ssm <- getMyGSSM(x)
    res$te.inter_mvssm <- getMy2VariateGSSM(mv2 = cbind(x, x2))
    if(plotme) {
      lines(t, res$te.ssm, col = "green", lwd=3)  
      lines(t, res$te.inter_mvssm, col = "magenta", lwd=3, lty = "dashed")
    }
  }
  if(!is.na(A[1])) {
    #lines(abpos, getMyGSSM(A), col = alpha(colour = "blue", alpha = 0.3), lwd=3)
    res$A.pos <- abpos
    res$A.rec_median <- signal::filter(MedianFilter(71), signal::filter(MedianFilter(41), signal::filter(MedianFilter(25), A)))
    if(plotme) {
      lines(abpos, res$A.rec_median, col = alpha(colour = "blue", alpha = 0.5), lwd=3)
    }
  }
  if(!is.na(B[1])) {
    res$B.pos <- abpos
    res$B.rec_median <- signal::filter(MedianFilter(71), signal::filter(MedianFilter(41), signal::filter(MedianFilter(25), B)))
    #lines(abpos, getMyGSSM(B), col = alpha(colour = "red", alpha = 0.3), lwd=2)
    if(plotme) {
      lines(abpos, res$B.rec_median, col = alpha(colour = "red", alpha = 0.55), lwd=2)
    }
  }
  if(plotme) {
    if(chr == "chr5")
      abline(v = 65600208, lty = "dashed", col = "green")
  
    legend("topright", inset=.05, title=NULL, cex = 0.8,
           legend = c("TPM ratio signal", "Recursive median Smooth", "Simple SSM", "MV interspace SSM", "Raw A", "Raw B", "TE MV-SSM combined"),
           fill=c(alpha(colour = "black", alpha = 0.5), "darkgreen", "green", "magenta", "blue", "red", "purple"), horiz=F)
  }
  
  if(!is.na(myalleles[1])) {
    tt <- getMy3VariateGSSM(mv3 = as.matrix(myalleles[, c("TE", "A", "B"), with=F]))
    res$te.teAB_mvssm <- tt
    res$teAB.pos <- myalleles$start
    if(plotme) {
      lines(myalleles$start, tt[,1], col = alpha(colour = "purple", alpha = 0.85), lwd=4)
      rng <- range(tt[,2:3])
      #plot(myalleles$start, tt[,1], type = "p", col = alpha(colour = "black", alpha = 0.5), main = "Combined analysis", xlab = xlab, ylab = "Signal output", ylim = rng)
      plot(myalleles$start, tt[,2], col = alpha(colour = "blue", alpha = 0.5), main = "Combined analysis", xlab = xlab, ylab = "Signal output", ylim = rng)
      #points(myalleles$start, tt[,2], col = alpha(colour = "blue", alpha = 0.5), lwd=3)
      points(myalleles$start, tt[,3], col = alpha(colour = "red", alpha = 0.7), lwd=2)
      legend("topright", inset=.05, title=NULL, cex = 0.8,
             legend = c("MV-SSM", "A", "B"),
             fill=c(alpha(colour = "black", alpha = 0.5), "blue", "red"), horiz=F)
    }
    
  }
  if(plotme) {
    if(chr == "chr5")
      abline(v = 65600208, lty = "dashed", col = "green")
  }
  return(res)
}

wapply <- function(Y, width, by = NULL)
{
  if (is.null(by)) by <- width
  
  lenX <- length(Y)
  SEQ1 <- seq(1, lenX - width + 1, by = by)
  SEQ2 <- lapply(SEQ1, function(x) Y[x:(x + width - 1)])
  return(SEQ2)
}

getMyLocalGSSM <- function(mv1) {
  res <- sapply(wapply(dt2[seqnames == "chr5"]$`170512_B1`, width = 30, by=5), function(x) median(getMyGSSM(x)))
  return(res)
}
getMyGSSM <- function(mv1) {
  require(KFAS)
  Zt <- matrix(c(1, 0), 1, 2)
  Ht <- matrix(NA)
  Tt <- matrix(c(1, 0, 1, 1), 2, 2)
  Rt <- matrix(c(1, 0), 2, 1)
  Qt <- matrix(NA)
  a1 <- matrix(c(1, 0), 2, 1)
  P1 <- matrix(0, 2, 2)
  P1inf <- diag(2)
  
  if(sum(mv1, na.rm = T) == 0)
    mv1[sample(length(mv1), 2)] <- 0.00001
  
  model_gaussian <- SSModel(mv1 ~ -1 +
                              SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf), 
                            H = Ht)
  fit_gaussian <- fitSSM(model_gaussian, inits = c(0, 0), method = "BFGS")
  out_gaussian <- KFS(fit_gaussian$model)
  return((coef(out_gaussian)[,1]))
}


getMy3VariateGSSM <- function(mv3, polyDegree = 1) {
  N <- 3
  #In cases where all signal is zero we seed some random values to avoid crash
  if(sum(mv3[,1], na.rm = T) == 0)
    mv3[sample(nrow(mv3), 2),1] <- 0.00001
  if(sum(mv3[,2], na.rm = T) == 0)
    mv3[sample(nrow(mv3), 2),2] <- 0.00001
  if(sum(mv3[,3], na.rm = T) == 0)
    mv3[sample(nrow(mv3), 2),3] <- 0.00001
  
  if(polyDegree == 2) {
    model <- SSModel(mv3 ~ -1 + SSMtrend(degree = 2, Q = list(matrix(NA, N, N), matrix(0, N, N))) +
                       SSMcustom(Z = diag(1, N), T = diag(0, N), Q = matrix(NA, N, N), P1 = matrix(NA, N, N)))
  } else if(polyDegree == 1) {
    model <- SSModel(mv3 ~ -1 + SSMtrend(degree = 1, Q = list(matrix(NA, N, N))) +
                       SSMcustom(Z = diag(1, N), T = diag(0, N), Q = matrix(NA, N, N), P1 = matrix(NA, N, N)))
  }
  
  init <- chol(cov(mv3), pivot = T)
  inits <- rep(c(log(diag(init)), init[upper.tri(init)]), 2)
  inits[is.nan(inits)] <- mean(inits, na.rm=T)
  inits[is.infinite(inits)] <- mean(inits[!is.infinite(inits)], na.rm=T)
  
  fitinit1 <- fitSSM(model, updatefn = updatefn3,
                     inits = inits,
                     method = "BFGS")
  fit1 <- fitSSM(model, updatefn = updatefn3, inits = fitinit1$optim.out$par,
                 method = "BFGS")
  
  out <- KFS(fit1$model) 
  
  return(coef(out))
  
}

updatefn3 <- function(pars, model, ...) {
  Q <- diag(exp(pars[1:3]))
  Q[upper.tri(Q)] <- pars[4:6]  
  model["Q", etas = "level"] <- crossprod(Q)
  Q <- diag(exp(pars[7:9]))
  Q[upper.tri(Q)] <- pars[10:12]
  model["Q", etas = "custom"] <- model["P1", states = "custom"] <- crossprod(Q)
  model
}

getMy2VariateGSSM <- function(mv2, polyDegree = 1) {
  require(KFAS)
  
  N <- 2
  #In cases where all signal is zero we seed some random values to avoid crash
  if(sum(mv2[,1], na.rm = T) == 0)
    mv2[sample(nrow(mv2), 2),1] <- 0.00001
  if(sum(mv2[,2], na.rm = T) == 0)
    mv2[sample(nrow(mv2), 2),2] <- 0.00001
  
  if(polyDegree == 2) { #poly = 2
    model <- SSModel(mv2 ~ -1 + SSMtrend(degree = 2, Q = list(matrix(NA, N, N), matrix(0, N, N))) +
                       SSMcustom(Z = diag(1, N), T = diag(0, N), Q = matrix(NA, N, N), P1 = matrix(NA, N, N)))
  } else {
    model <- SSModel(mv2 ~ -1 + SSMtrend(degree = 1, Q = list(matrix(NA, N, N))) +
                       SSMcustom(Z = diag(1, N), T = diag(0, N), Q = matrix(NA, N, N), P1 = matrix(NA, N, N)))
  }
  
  init <- chol(cov(mv2), pivot = T)
  inits <- rep(c(log(diag(init)), init[upper.tri(init)]), 2)
  inits[is.nan(inits)] <- mean(inits, na.rm=T)
  inits[is.na(inits)] <- mean(inits, na.rm=T)
  inits[is.infinite(inits)] <- mean(inits[!is.infinite(inits)], na.rm=T)
  
  
  fitinit <- fitSSM(model, updatefn = updatefn2,
                    inits = inits,
                    method = "BFGS")
  fit <- fitSSM(model, updatefn = updatefn2, inits = fitinit$optim.out$par,
                method = "BFGS")
  out <- KFS(fit$model) 
  return(coef(out)[,1])
}

getMy2VariateGSSM.raw <- function(mv2, polyDegree = 1) {
  require(KFAS)
  
  N <- 2
  #In cases where all signal is zero we seed some random values to avoid crash
  if(sum(mv2[,1], na.rm = T) == 0)
    mv2[sample(nrow(mv2), 2),1] <- 0.00001
  if(sum(mv2[,2], na.rm = T) == 0)
    mv2[sample(nrow(mv2), 2),2] <- 0.00001
  
  if(polyDegree == 2) { #poly = 2
    model <- SSModel(mv2 ~ -1 + SSMtrend(degree = 2, Q = list(matrix(NA, N, N), matrix(0, N, N))) +
                       SSMcustom(Z = diag(1, N), T = diag(0, N), Q = matrix(NA, N, N), P1 = matrix(NA, N, N)))
  } else {
    model <- SSModel(mv2 ~ -1 + SSMtrend(degree = 1, Q = list(matrix(NA, N, N))) +
                       SSMcustom(Z = diag(1, N), T = diag(0, N), Q = matrix(NA, N, N), P1 = matrix(NA, N, N)))
  }
  
  init <- chol(cov(mv2), pivot = T)
  inits <- rep(c(log(diag(init)), init[upper.tri(init)]), 2)
  inits[is.nan(inits)] <- mean(inits, na.rm=T)
  inits[is.na(inits)] <- mean(inits, na.rm=T)
  inits[is.infinite(inits)] <- mean(inits[!is.infinite(inits)], na.rm=T)
  
  
  fitinit <- fitSSM(model, updatefn = updatefn2,
                    inits = inits,
                    method = "BFGS")
  fit <- fitSSM(model, updatefn = updatefn2, inits = fitinit$optim.out$par,
                method = "BFGS")
  out <- KFS(fit$model) 
  return(coef(out))
}

updatefn2 <- function(pars, model, ...) {
  Q <- diag(exp(pars[1:2]))
  Q[upper.tri(Q)] <- pars[3]  
  model["Q", etas = "level"] <- crossprod(Q)
  Q <- diag(exp(pars[4:5]))
  Q[upper.tri(Q)] <- pars[6]
  model["Q", etas = "custom"] <- model["P1", states = "custom"] <- crossprod(Q)
  model
}

AddFeaturesandPhasing <- function(m, anno, centromeres ) {
# BY NIKOS MYNHIER
  require(readxl)
  #anno <- data.table(read_excel("/singlecellcenter/etai/ExperimentsData/Stamatis/Stamatis_list_v14_181103.xlsx"))
  #anno <- data.table(read_excel("/pellmanlab/nikos/Stam_Etai_Data/Stamatis_list_v15_190910.xlsx"))
  #old_anno <- data.table(read_excel(path = "/pellmanlab/nikos/Stam_Etai_Data/Stamatis_list_v12_180508.QCRNAv3.withQCRNAv3.xlsx", sheet = "all", na = c("", "NA")))
  #old_anno2 <- readRDS("/homes10/ejacob/WORK/secondaryanalysis/Stamatis/reports/Stamatis_list_v12_180404.QCRNAv3.annoqc.rds")
  #centromeres <- readRDS("/pellmanlab/stam_niko/etai_code/analysis_result_data/centromeres.rds")
  #myIDs <- anno[QC_th1 > 6000 & (control_group != "bunch_of_grapes" | is.na(control_group))]$WTA.plate

  #modifying column names of Hiseq run (the names of columns 1-382 are already modified):
  #colnames(m$A)[-c(1:382)] <- old_anno2$WTA.plate[sapply(colnames(m$A)[-c(1:382)], grep, x = old_anno2$Fastq_files3, value = F)]
  #colnames(m$B)[-c(1:382)] <- old_anno2$WTA.plate[sapply(colnames(m$B)[-c(1:382)], grep, x = old_anno2$Fastq_files3, value = F)]
  #modifying column names of Miseq run
  #colnames(m$A)[grep("-", colnames(m$A), value = F)] <- gsub("-", "_", gsub(pattern = "_S\\d+_L001_R1_001", "", grep("-", colnames(m$A), value = T)))
  #colnames(m$B)[grep("-", colnames(m$B), value = F)] <- gsub("-", "_", gsub(pattern = "_S\\d+_L001_R1_001", "", grep("-", colnames(m$B), value = T)))
  
  #Phasing - RNA based expression:
  expFeatures <- generateExpressionBasedSNPFeatures_simple(m = m)
  print("expression done")
  
  fracPerChr <- rbindlist(lapply(unique(expFeatures$seqnames), function(chr) {
    data.table(seqnames = chr, 
               frac1.5 = round(100*sum(with(expFeatures[ seqnames == chr & above0.A > 2 & above0.B > 2], abs(log2((AgB + 1) / 
                                                                                                                    (BgA + 1))) > log2(1.5)))/nrow(expFeatures[ seqnames == chr & above0.A > 2 & above0.B > 2]), digits = 1),
               frac2 = round(100*sum(with(expFeatures[ seqnames == chr & above0.A > 2 & above0.B > 2], abs(log2((AgB + 1) / 
                                                                                                                  (BgA + 1))) > log2(2)))/nrow(expFeatures[ seqnames == chr & above0.A > 2 & above0.B > 2]), digits = 1),
               frac3 = round(100*sum(with(expFeatures[ seqnames == chr & above0.A > 2 & above0.B > 2], abs(log2((AgB + 1) / 
                                                                                                                  (BgA + 1))) > log2(3)))/nrow(expFeatures[ seqnames == chr & above0.A > 2 & above0.B > 2]), digits = 1),
               frac4 = round(100*sum(with(expFeatures[ seqnames == chr & above0.A > 2 & above0.B > 2], abs(log2((AgB + 1) / 
                                                                                                                  (BgA + 1))) > log2(4)))/nrow(expFeatures[ seqnames == chr & above0.A > 2 & above0.B > 2]), digits = 1),
               frac5 = round(100*sum(with(expFeatures[ seqnames == chr & above0.A > 2 & above0.B > 2], abs(log2((AgB + 1) / 
                                                                                                                  (BgA + 1))) > log2(5)))/nrow(expFeatures[ seqnames == chr & above0.A > 2 & above0.B > 2]), digits = 1))
  }))
  
  #The largest threshold (i.e. accounting for as much alleles as possible in the analysis) of imblanced bi-allelic expression that can still detect the 10q gain is frac3 based on the "fracPerChr" analysis
  m$imbalancedFractionOfSNPs <- fracPerChr
  #m$refIDs <- myIDs
  m$features <- expFeatures
  #saveRDS(m, file = out_fname)
  
  ########################################################################
  #filtering out all non biallelic expressed genes except in 10q and chrX
  ########################################################################
  #We take any monoallelic expressed SNP from chrX:
  wi.chrX <- which(m$features$seqnames == "chrX" & (m$features$above0.A > 2 | m$features$above0.B > 2))
  #we take all 10q SNP with very high allelic expression imbalance (frac5) since we know this arm has 3 copies in RPE-1
  start10q <- max(centromeres[seqnames == "chr10"]$end)
  #The fraction of 10q and 10p above th frac6 (3.5% vs. 2.1%)
  length(with(m$features, which(seqnames == "chr10" & start > start10q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(6))))/length(with(m$features, which(seqnames == "chr10" & start > start10q & (above0.A > 2 & above0.B > 2) )))
  length(with(m$features, which(seqnames == "chr10" & start < start10q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(6))))/length(with(m$features, which(seqnames == "chr10" & start < start10q & (above0.A > 2 & above0.B > 2) )))
  #The fraction of 10q and 10p above th frac5 (4.5% vs. 2.4%)
  length(with(m$features, which(seqnames == "chr10" & start > start10q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(5))))/length(with(m$features, which(seqnames == "chr10" & start > start10q & (above0.A > 2 & above0.B > 2) )))
  length(with(m$features, which(seqnames == "chr10" & start < start10q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(5))))/length(with(m$features, which(seqnames == "chr10" & start < start10q & (above0.A > 2 & above0.B > 2) )))
  #The fraction of 10q and 10p above th frac3 (15.5% vs. 6%)
  length(with(m$features, which(seqnames == "chr10" & start > start10q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(3))))/length(with(m$features, which(seqnames == "chr10" & start > start10q & (above0.A > 2 & above0.B > 2) )))
  length(with(m$features, which(seqnames == "chr10" & start < start10q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(3))))/length(with(m$features, which(seqnames == "chr10" & start < start10q & (above0.A > 2 & above0.B > 2) )))
  
  #For instance, in chr1 (1p and 1q) yo do not see this bias as in chr10
  start1q <- max(centromeres[seqnames == "chr1"]$end)
  #The fraction of 10q and 10p above th frac6 (2.5% vs. 4.4%)
  length(with(m$features, which(seqnames == "chr1" & start > start1q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(6))))/length(with(m$features, which(seqnames == "chr1" & start > start1q & (above0.A > 2 & above0.B > 2) )))
  length(with(m$features, which(seqnames == "chr1" & start < start1q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(6))))/length(with(m$features, which(seqnames == "chr1" & start < start1q & (above0.A > 2 & above0.B > 2) )))
  #The fraction of 10q and 10p above th frac5 (3.2% vs. 5.1%)
  length(with(m$features, which(seqnames == "chr1" & start > start1q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(5))))/length(with(m$features, which(seqnames == "chr1" & start > start1q & (above0.A > 2 & above0.B > 2) )))
  length(with(m$features, which(seqnames == "chr1" & start < start1q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(5))))/length(with(m$features, which(seqnames == "chr1" & start < start1q & (above0.A > 2 & above0.B > 2) )))
  #The fraction of 10q and 10p above th frac3 (7.5% vs. 10%)
  length(with(m$features, which(seqnames == "chr1" & start > start1q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(3))))/length(with(m$features, which(seqnames == "chr1" & start > start1q & (above0.A > 2 & above0.B > 2) )))
  length(with(m$features, which(seqnames == "chr1" & start < start1q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) > log2(3))))/length(with(m$features, which(seqnames == "chr1" & start < start1q & (above0.A > 2 & above0.B > 2) )))
  
  #Therefore, we will use frac5 as a threshold for 10q
  wi.chr10q <- with(m$features, which(seqnames == "chr10" & start > start10q & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) <= log2(5)))
  
  #chr12 has a high frequency of trisomies, therefore, we will use a higher threshold of frac4 for filtering out those snps:
  wi.chr12 <- with(m$features, which(seqnames == "chr12" & (above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) <= log2(4)))
  wi.all <- with(m$features, which((above0.A > 2 & above0.B > 2) & abs(log2((AgB + 1) / (BgA + 1))) <= log2(3)))
  
  #and now we union all SNPs
  wi <- unique(c(wi.all, wi.chr10q, wi.chr12, wi.chrX))
  m$features$biAllelicFilter <- FALSE
  m$features$biAllelicFilter[wi] <- TRUE
  
  #Coding SNPs include also 5` and 3` UTRs:
  with(m$features, which(biAllelicFilter == T & (!is.na(coding) | (!is.na(fiveUTR)) | (!is.na(threeUTR)))))
  return(m)
}

generateExpressionBasedSNPFeatures_simple <- function(m) {
  # BY NIKOS MYNHIER
  require(data.table)
  message("Preparing data..")
  A <- copy(m$A)
  B <- copy(m$B)
  features <- copy(m$features)
  
  A[is.na(A)] <- 0
  B[is.na(B)] <- 0
  message("Calculating features based on expression..")
  AN <- rowSums(A > B)
  BN <- rowSums(B > A)
  
  features$above0.A <- rowSums(A>0)
  features$above0.B <- rowSums(B>0)
  features$AgB <- AN
  features$BgA <- BN
  
  am <- rowMeans(A, na.rm = T)
  bm <- rowMeans(B, na.rm = T)
  features$avg.A <- am
  features$avg.B <- bm
  
  return(features)
}
