
#SSM +1 estimates:
#-------------------
IDs <- names(which(colSums(rsemtpm>0)> 2000))
cpus <- 24
require(snow)
cl <- snow::makeCluster(cpus, type = "SOCK")
clusterEvalQ(cl, library(data.table))
clusterEvalQ(cl, source('/homes10/ejacob/WORK/secondaryanalysis/methods_paper_results/R/CNML.R'))

snow::clusterExport(cl = cl, list = c("getMyGSSM", "adt", "IDs"))

tt <- do.call(cbind, parLapply(cl, IDs, function(i) unlist(lapply(unique(adt$seqnames), function(chr, i) getMyGSSM(adt[seqnames == chr][[i]]), i))))
colnames(tt) <- IDs
SSM <- data.table(adt[, 1:4], tt)

stopCluster(cl)

dtm <- SSM[, lapply(.SD, median),
           .SDcols = -(1:4), by = "seqnames"]
myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in% chrsToExcludeFromNormalization), -1]))
SSM2 <- 2^sweep(SSM[,-(1:4)], 2, myCellMedians, "-")
SSM <- cbind(adt[, 1:4], SSM2)
hist(as.matrix(SSM[, -(1:4)]), breaks=111, col="darkblue", xlim=c(0,4))
saveRDS(SSM, file = sprintf("%s/../Stamatis_list_v14_181025.rsemtpm.SSM.GM.rds", wkdir_e))
saveRDS(SSM, file = sprintf("%s/../SSM1.rds", wkdir_e))

#SSM +0 estimates:  FAILED EXECUTION: issues with snow an cluster connection (dont close cluster above)
#-------------------
snow::clusterExport(cl = cl, list = c("getMyGSSM", "adt0", "IDs"))

tt <- do.call(cbind, parLapply(cl, IDs, function(i) unlist(lapply(unique(adt0$seqnames), function(chr, i) getMyGSSM(adt0[seqnames == chr][[i]]), i))))
colnames(tt) <- IDs
SSM0<- data.table(adt0[, 1:4], tt)


dtm <- SSM0[, lapply(.SD, median),
           .SDcols = -(1:4), by = "seqnames"]
myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in% chrsToExcludeFromNormalization), -1]))
SSM02 <- sweep(SSM0[,-(1:4)], 2, myCellMedians, "/")
hist(as.matrix(SSM02[SSM02<4.1]), breaks=111, col="darkblue")
SSM0 <- cbind(adt[, 1:4], SSM02)

saveRDS(SSM0, file = sprintf("%s/../SSM0.rds", wkdir_e))

#MV-SSM calculations +1
#----------------------
cpus <- 24
require(snow)
cl <- snow::makeCluster(cpus, type = "SOCK")
clusterEvalQ(cl, library(data.table))
clusterEvalQ(cl, source('/homes10/ejacob/WORK/secondaryanalysis/methods_paper_results/R/CNML.R'))
clusterEvalQ(cl, source('/homes10/ejacob/WORK/secondaryanalysis/Stamatis/G2/runRNAseqPipelineFor_BF9VP.utils.R'))

snow::clusterExport(cl = cl, list = c("getMy2VariateGSSM", "adt", "IDs", "tpmdata"))

tt <- do.call(cbind, parLapply(cl, IDs, function(i) unlist(lapply(unique(adt$seqnames), 
                                                                    function(chr, i) getMy2VariateGSSM(mv2 = cbind(adt[seqnames == chr][[i]], 
                                                                                                                   tpmdata$interspace[adt$seqnames == chr])), i))))
colnames(tt) <- IDs
MVSSM <- data.table(adt[, 1:4], tt)

stopCluster(cl)

dtm <- MVSSM[, lapply(.SD, median),
           .SDcols = -(1:4), by = "seqnames"]
myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in% chrsToExcludeFromNormalization), -1]))
MVSSM2 <- 2^sweep(MVSSM[,-(1:4)], 2, myCellMedians, "-")
MVSSM <- cbind(adt[, 1:4], MVSSM2)

saveRDS(MVSSM, file = sprintf("%s/../Stamatis_list_v14_181025.rsemtpm.MVSSM.GM.rds", wkdir_e))
saveRDS(MVSSM, file = sprintf("%s/../MVSSM1.TEIS.rds", wkdir_e))

#######################################################################
#save.image("/pellmanlab/stam_niko/etai_code/experiments/Save_states/post_mv-ssm1.RData")
#load("/pellmanlab/stam_niko/etai_code/experiments/Save_states/post_mv-ssm1.RData")
#######################################################################

#calculation of MA:
#-------------------
MA_winSize <- 100
MA <- adt[, lapply(.SD, rollmean, k = MA_winSize, fill = "extend", align = "center"),
          .SDcols = names(adt)[-c(1:4)], by = seqnames]

dtm <- MA[, lapply(.SD, median),
          .SDcols = -1, by = "seqnames"]

myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in%chrsToExcludeFromNormalization), -1]))
MA2 <- 2^sweep(MA[,-1], 2, myCellMedians, "-")
MA <- cbind(adt[, 1:4], MA2)
hist(as.matrix(MA[, -c(1:4)]))

saveRDS(MA, file = sprintf("%s/../MA.rds", wkdir_e))
#######################################################################
#save.image("/pellmanlab/stam_niko/etai_code/experiments/Save_states/post_ma.RData")
#load("/pellmanlab/stam_niko/etai_code/experiments/Save_states/post_ma.RData")
#######################################################################

#calculation of IMF:
#-------------------
medfilterRes <- runIterativeMedianFilterOnAllSamples(dt2 = adt, CPUs = 42)
mats.imf <- getAllMatricesForItrMedFilter(myres.all = medfilterRes)
mats.imf2 <- lapply(mats.imf, function(x) cbind(adt[,1:4], x))
IMF <- mats.imf2
saveRDS(IMF, file = sprintf("%s/../Stamatis_list_v14_181025.rsemtpm.IMF.GM.rds", wkdir_e))
dtm <- IMF$level6[, lapply(.SD, median),
          .SDcols = -c(1:4), by = "seqnames"]

myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in%chrsToExcludeFromNormalization), -1]))
IMF6 <- cbind(IMF$level6[, 1:4], 2^sweep(IMF$level6[,-(1:4)], 2, myCellMedians, "-"))
saveRDS(IMF6, file = sprintf("%s/../IMF6.rds", wkdir_e))

#######################################################################
save.image("/pellmanlab/stam_niko/etai_code/experiments/Save_states/post_imf.RData")
#load("/pellmanlab/stam_niko/etai_code/experiments/Save_states/post_imf.RData")
#######################################################################

#calculation of ML:
#------------------
mlIDs <- gsub(pattern = "_predicted_by_id17.rds", replacement = "", x = list.files(path = "/singlecellcenter/etai/tmp/CNML/pointOfChange/run_id17/output/", pattern = "*predicted_by_id17.rds"))
tmp <- lapply(mlIDs, function(id) 
  readRDS(sprintf("/singlecellcenter/etai/tmp/CNML/pointOfChange/run_id17/output/%s_predicted_by_id17reg.rds", id)))
names(tmp) <- mlIDs

tmp2 <- do.call(cbind, lapply(tmp, function(x) x$prediction))
tmp3 <- cbind(tmp[[1]][,1:4], tmp2)
MLREG <- tmp3
saveRDS(MLREG, file = sprintf("%s/../MLREG.rds", wkdir_e))


#calculation of ML regression:
#------------------
mlIDs <- gsub(pattern = "_predicted_by_id17reg.rds", replacement = "", x = list.files(path = "/singlecellcenter/etai/tmp/CNML/pointOfChange/run_id17/output/", pattern = "*predicted_by_id17reg.rds"))
tmp <- lapply(mlIDs, function(id) 
  readRDS(sprintf("/singlecellcenter/etai/tmp/CNML/pointOfChange/run_id17/output/%s_predicted_by_id17.rds", id)))
names(tmp) <- mlIDs

tmp2 <- do.call(cbind, lapply(tmp, function(x) x$prediction))
tmp3 <- cbind(tmp[[1]][,1:4], tmp2)
ML <- tmp3
saveRDS(ML, file = sprintf("%s/../ML.rds", wkdir_e))

#Calculation of MV-SSM for 3 variable alleles +0 and +1:
#-------------------------------------------------------
chrsToExcludeFromNormalization = c("chrX", "chrM", "chr10")

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

getMyAlleleLevelSSMs0 <- function(i, alleles.all, adt0, chr = NULL) {
  myalleles <- mergeTotalExpAndAlleleData(adt0, A = alleles.all$A, B = alleles.all$B, id = i)
  if(is.null(chr)) {
    tt <- data.table(getMy3VariateGSSM(mv3 = as.matrix(myalleles[, c("TE", "A", "B"), with=F])))
    return(cbind(myalleles, tt))
  } else {
    tt <- data.table(getMy3VariateGSSM(mv3 = as.matrix(myalleles[seqnames==chr, c("TE", "A", "B"), with=F])))
    return(cbind(myalleles[seqnames == chr], tt))
  }
}

cpus <- 24
require(snow)
cl <- snow::makeCluster(cpus, type = "SOCK")
clusterEvalQ(cl, library(data.table))
clusterEvalQ(cl, library(KFAS))
clusterEvalQ(cl, source('/homes10/ejacob/WORK/secondaryanalysis/methods_paper_results/R/CNML.R'))
clusterEvalQ(cl, source('/homes10/ejacob/WORK/secondaryanalysis/Stamatis/G2/runRNAseqPipelineFor_BF9VP.utils.R'))

#choosing IDs:
IDs <- names(which(colSums(rsemtpm>0)> 2000))
IDs <- IDs[IDs %in% colnames(alleles.all$A)]


snow::clusterExport(cl = cl, list = c("getMyAlleleLevelSSMs", "adt", "IDs", "alleles.all"))
snow::clusterExport(cl = cl, list = c("getMyAlleleLevelSSMs0", "adt0"))

#running for +0:
tt <- parLapply(cl, IDs, function(i) rbindlist(lapply(as.character(unique(adt0$seqnames)), 
                                                      function(chr, i) getMyAlleleLevelSSMs0(i = i, alleles.all = alleles.all, adt0 = adt0, chr = chr), i)))
names(tt) <- IDs
MVSSM0.AB <- tt
saveRDS(MVSSM0.AB, file = sprintf("%s/../Stamatis_list_v14_181025.rsemtpm.MVSSM0.AB.rds", wkdir_e))

#collecting per signal:
#TE
tmp <- data.table(MVSSM0.AB[[1]][, c("gene_id", "seqnames", "start", "end")], do.call(cbind, lapply(MVSSM0.AB, function(x) x$level.TE)))
dtm <- tmp[, lapply(.SD, median),
             .SDcols = -(1:4), by = "seqnames"]
myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in% chrsToExcludeFromNormalization), -1]))
tmp2 <- sweep(tmp[,-(1:4)], 2, myCellMedians, "/")
par(mfrow=c(1,1))
hist(as.matrix(tmp2[tmp2<4.1]), breaks=111, col="darkblue")
MVSSM0.TE <- cbind(tmp[, 1:4], tmp2)
#A
tmp <- data.table(MVSSM0.AB[[1]][, c("gene_id", "seqnames", "start", "end")], do.call(cbind, lapply(MVSSM0.AB, function(x) x$level.A)))
dtm <- tmp[, lapply(.SD, median),
           .SDcols = -(1:4), by = "seqnames"]
myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in% chrsToExcludeFromNormalization), -1]))
tmp2 <- sweep(tmp[,-(1:4)], 2, myCellMedians, "/")
par(mfrow=c(1,1))
hist(as.matrix(tmp2[tmp2<4.1]), breaks=111, col="darkblue")
MVSSM0.A <- cbind(tmp[, 1:4], tmp2)
#B
tmp <- data.table(MVSSM0.AB[[1]][, c("gene_id", "seqnames", "start", "end")], do.call(cbind, lapply(MVSSM0.AB, function(x) x$level.B)))
dtm <- tmp[, lapply(.SD, median),
           .SDcols = -(1:4), by = "seqnames"]
myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in% chrsToExcludeFromNormalization), -1]))
tmp2 <- sweep(tmp[,-(1:4)], 2, myCellMedians, "/")
par(mfrow=c(1,1))
hist(as.matrix(tmp2[tmp2<4.1]), breaks=111, col="darkblue")
MVSSM0.B <- cbind(tmp[, 1:4], tmp2)
#saving:
saveRDS(MVSSM0.TE, file = sprintf("%s/../MVSSM0.TEAB.TE.rds", wkdir_e))
saveRDS(MVSSM0.A, file  = sprintf("%s/../MVSSM0.TEAB.A.rds" , wkdir_e))
saveRDS(MVSSM0.B, file  = sprintf("%s/../MVSSM0.TEAB.B.rds" , wkdir_e))


#running for +1
tt1 <- parLapply(cl, IDs, function(i) rbindlist(lapply(as.character(unique(adt$seqnames)), 
                                                      function(chr, i) getMyAlleleLevelSSMs(i = i, alleles.all = alleles.all, adt = adt, chr = chr), i)))
names(tt1) <- IDs
MVSSM.AB <- tt1
saveRDS(MVSSM.AB, file = sprinf("%s/../Stamatis_list_v14_181025.rsemtpm.MVSSM.AB.rds", wkdir_e))

#collecting per signal:
#TE
tmp <- data.table(MVSSM.AB[[1]][, c("gene_id", "seqnames", "start", "end")], do.call(cbind, lapply(MVSSM.AB, function(x) x$level.TE)))
dtm <- tmp[, lapply(.SD, median),
           .SDcols = -(1:4), by = "seqnames"]
myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in% chrsToExcludeFromNormalization), -1]))
tmp2 <- 2^sweep(tmp[,-(1:4)], 2, myCellMedians, "-")
par(mfrow=c(1,1))
hist(as.matrix(tmp2[tmp2<4.1]), breaks=111, col="darkblue")
MVSSM1.TE <- cbind(tmp[, 1:4], tmp2)
#A
tmp <- data.table(MVSSM.AB[[1]][, c("gene_id", "seqnames", "start", "end")], do.call(cbind, lapply(MVSSM.AB, function(x) x$level.A)))
dtm <- tmp[, lapply(.SD, median),
           .SDcols = -(1:4), by = "seqnames"]
myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in% chrsToExcludeFromNormalization), -1]))
tmp2 <- 2^sweep(tmp[,-(1:4)], 2, myCellMedians, "-")
par(mfrow=c(1,1))
hist(as.matrix(tmp2[tmp2<4.1]), breaks=111, col="darkblue")
MVSSM1.A <- cbind(tmp[, 1:4], tmp2)
#B
tmp <- data.table(MVSSM.AB[[1]][, c("gene_id", "seqnames", "start", "end")], do.call(cbind, lapply(MVSSM.AB, function(x) x$level.B)))
dtm <- tmp[, lapply(.SD, median),
           .SDcols = -(1:4), by = "seqnames"]
myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in% chrsToExcludeFromNormalization), -1]))
tmp2 <- 2^sweep(tmp[,-(1:4)], 2, myCellMedians, "-")
par(mfrow=c(1,1))
hist(as.matrix(tmp2[tmp2<4.1]), breaks=111, col="darkblue")
MVSSM1.B <- cbind(tmp[, 1:4], tmp2)
#saving:
saveRDS(MVSSM1.TE, file = sprintf("%s/../MVSSM1.TEAB.TE.rds", wkdir_e))
saveRDS(MVSSM1.A, file  = sprintf("%s/../MVSSM1.TEAB.A.rds" , wkdir_e))
saveRDS(MVSSM1.B, file  = sprintf("%s/../MVSSM1.TEAB.B.rds" , wkdir_e))


stopCluster(cl)


#########################################
#Calculating independently allele levels:
#########################################

#Calculation of MV-SSM for 2 variable alleles +0 and +1:
#-------------------------------------------------------
chrsToExcludeFromNormalization = c("chrX", "chrM", "chr10")

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


getMyOneAlleleLevelSSMs0 <- function(i, alleles.all, adt0, chr = NULL, allele = "A") {
  myalleles <- mergeTotalExpAndAlleleData(adt0, A = alleles.all$A, B = alleles.all$B, id = i)
  if(is.null(chr)) {
    tt <- data.table(getMy2VariateGSSM(mv2 = as.matrix(myalleles[, c(allele, "TE"), with=F])))
    return(cbind(myalleles, tt))
  } else {
    tt <- data.table(getMy2VariateGSSM(mv2 = as.matrix(myalleles[seqnames==chr, c(allele, "TE"), with=F])))
    return(cbind(myalleles[seqnames == chr], tt))
  }
}

cpus <- 24
require(snow)
cl <- snow::makeCluster(cpus, type = "SOCK")
clusterEvalQ(cl, library(data.table))
clusterEvalQ(cl, library(KFAS))
clusterEvalQ(cl, source('/homes10/ejacob/WORK/secondaryanalysis/methods_paper_results/R/CNML.R'))
clusterEvalQ(cl, source('/homes10/ejacob/WORK/secondaryanalysis/Stamatis/G2/runRNAseqPipelineFor_BF9VP.utils.R'))

#choosing IDs:
IDs <- names(which(colSums(rsemtpm>0)> 2000))
IDs <- IDs[IDs %in% colnames(alleles.all$A)]


snow::clusterExport(cl = cl, list = c("getMyOneAlleleLevelSSMs", "adt", "IDs", "alleles.all"))
snow::clusterExport(cl = cl, list = c("getMyOneAlleleLevelSSMs0", "adt0"))

#running for +0:
#---------------
tt.A <- parLapply(cl, IDs, function(i) rbindlist(lapply(as.character(unique(adt0$seqnames)), 
                                                      function(chr, i) getMyOneAlleleLevelSSMs0(i = i, alleles.all = alleles.all, adt0 = adt0, chr = chr, allele = "A"), i)))
names(tt.A) <- IDs
MVSSM0.TEA <- tt.A
saveRDS(MVSSM0.TEA, file = sprintf("%s/../Stamatis_list_v14_181025.rsemtpm.MVSSM0.TEA.rds", wkdir_e))

tt.B <- parLapply(cl, IDs, function(i) rbindlist(lapply(as.character(unique(adt0$seqnames)), 
                                                        function(chr, i) getMyOneAlleleLevelSSMs0(i = i, alleles.all = alleles.all, adt0 = adt0, chr = chr, allele = "B"), i)))
names(tt.B) <- IDs
MVSSM0.TEB <- tt.B
saveRDS(MVSSM0.TEB, file = sprintf("%s/../Stamatis_list_v14_181025.rsemtpm.MVSSM0.TEB.rds", wkdir_e))


#collecting per signal:

#A
tmp <- data.table(MVSSM0.TEA[[1]][, c("gene_id", "seqnames", "start", "end")], do.call(cbind, lapply(MVSSM0.TEA, function(x) as.numeric(x$V1))))
dtm <- tmp[, lapply(.SD, median),
           .SDcols = -(1:4), by = "seqnames"]
myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in% chrsToExcludeFromNormalization), -1]))
tmp2 <- sweep(tmp[,-(1:4)], 2, myCellMedians, "/")
par(mfrow=c(1,1))
hist(as.matrix(tmp2[tmp2<4.1]), breaks=111, col="darkblue")
MVSSM0.TEA <- cbind(tmp[, 1:4], tmp2)
#B
tmp <- data.table(MVSSM0.TEB[[1]][, c("gene_id", "seqnames", "start", "end")], do.call(cbind, lapply(MVSSM0.TEB, function(x) as.numeric(x$V1))))
dtm <- tmp[, lapply(.SD, median),
           .SDcols = -(1:4), by = "seqnames"]
myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in% chrsToExcludeFromNormalization), -1]))
tmp2 <- sweep(tmp[,-(1:4)], 2, myCellMedians, "/")
par(mfrow=c(1,1))
hist(as.matrix(tmp2[tmp2<4.1]), breaks=111, col="darkblue")
MVSSM0.TEB <- cbind(tmp[, 1:4], tmp2)

#saving:
saveRDS(MVSSM0.TEA, file = sprintf("%s/../MVSSM0.TEA.rds", wkdir_e))
saveRDS(MVSSM0.TEB, file = sprintf("%s/../MVSSM0.TEB.rds", wkdir_e))



#running for +1:
#---------------
tt.A <- parLapply(cl, IDs, function(i) rbindlist(lapply(as.character(unique(adt$seqnames)), 
                                                        function(chr, i) getMyOneAlleleLevelSSMs(i = i, alleles.all = alleles.all, adt = adt, chr = chr, allele = "A"), i)))
names(tt.A) <- IDs
MVSSM1.TEA <- tt.A
saveRDS(MVSSM1.TEA, file = sprintf("%s/../Stamatis_list_v14_181025.rsemtpm.MVSSM1.TEA.rds", wkdir_e))

tt.B <- parLapply(cl, IDs, function(i) rbindlist(lapply(as.character(unique(adt$seqnames)), 
                                                        function(chr, i) getMyOneAlleleLevelSSMs(i = i, alleles.all = alleles.all, adt = adt, chr = chr, allele = "B"), i)))
names(tt.B) <- IDs
MVSSM1.TEB <- tt.B
saveRDS(MVSSM1.TEB, file = sprintf("%s/../Stamatis_list_v14_181025.rsemtpm.MVSSM1.TEB.rds", wkdir_e))


#collecting per signal:

#A
tmp <- data.table(MVSSM1.TEA[[1]][, c("gene_id", "seqnames", "start", "end")], do.call(cbind, lapply(MVSSM1.TEA, function(x) as.numeric(x$V1))))
dtm <- tmp[, lapply(.SD, median),
           .SDcols = -(1:4), by = "seqnames"]
myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in% chrsToExcludeFromNormalization), -1]))
tmp2 <- sweep(tmp[,-(1:4)], 2, myCellMedians, "-")
par(mfrow=c(1,1))
hist(as.matrix(tmp2[tmp2<4.1]), breaks=111, col="darkblue")
MVSSM1.TEA <- cbind(tmp[, 1:4], tmp2)
#B
tmp <- data.table(MVSSM1.TEB[[1]][, c("gene_id", "seqnames", "start", "end")], do.call(cbind, lapply(MVSSM1.TEB, function(x) as.numeric(x$V1))))
dtm <- tmp[, lapply(.SD, median),
           .SDcols = -(1:4), by = "seqnames"]
myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in% chrsToExcludeFromNormalization), -1]))
tmp2 <- sweep(tmp[,-(1:4)], 2, myCellMedians, "-")
par(mfrow=c(1,1))
hist(as.matrix(tmp2[tmp2<4.1]), breaks=111, col="darkblue")
MVSSM1.TEB <- cbind(tmp[, 1:4], tmp2)

#saving:
saveRDS(MVSSM1.TEA, file = sprintf("%s/../MVSSM1.TEA.rds", wkdir_e))
saveRDS(MVSSM1.TEB, file = sprintf("%s/../MVSSM1.TEB.rds", wkdir_e))

stopCluster(cl)


