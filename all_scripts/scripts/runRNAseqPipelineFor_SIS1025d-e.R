
scriptsdir <- "/pellmanlab/stam_niko/etai_code/DFCI.scRNAseq.workflows/scripts"
wkdir_d <- "/pellmanlab/stam_niko/etai_code/experiments_etai/SIS1025d"
wkdir_e <- "/pellmanlab/stam_niko/etai_code/experiments_etai/SIS1025e"
mk_wkdir_d <- sprintf("mkdir -p /pellmanlab/stam_niko/etai_code/experiments_etai/SIS1025d", wkdir_d)
mk_wkdir_e <- sprintf("mkdir -p /pellmanlab/stam_niko/etai_code/experiments_etai/SIS1025e", wkdir_e)
system(mk_wkdir_d)
system(mk_wkdir_e)


#TODO: MOVE HARD DATA TO ACCESSABLE AND MODULAR LOCATION AND RUN THE FIRST BAM CMD

#SIS1025d

##################################################################
#Bash basic genotye pipeline cmds (from fastq to TPM and SNP RDs):
##################################################################

#trasforming bam files to fastqs:
#--------------------------------
#first option
#setwd(wkdir_d)
#system("mkdir -p bams")
#setwd("bams")
#system("cp /singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025d/bams/*.bam  ${PWD}/")
#system("parallel -j 16 samtools sort -n {} -o ${PWD}/{/}.qsort ::: /papathanasiou/SIS1025d/*.bam &")
#setwd(wkdir_d)
#system("mkdir -p fastqs")
#setwd("fastqs")
#bamtofastq <- sprintf("parallel -j 22 bedtools bamtofastq -i {} -fq ${PWD}/{/}.R1.fastq -fq2 ${PWD}/{/}.R2.fastq ::: %s/bams/*.bam &", wkdir_d)
#system(bamtofastq)
#system("parallel -j 22 gzip {} ::: ${PWD}/*.fastq")
#second option: if you are given fastqs directly
#mydir <- "/singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025d"
#system("cp -s /papathanasiou/SIS1025d/*.bam /singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025d/bams")
#bamFiles <- list.files(sprintf("%s/bams", mydir), pattern = "*.bam$", full.names = T)
#cmds <- sprintf("bash ~ejacob/WORK/RNA/Workflows/bam2fastq.pairedend.sh %s", bamFiles)
#writeLines(text = cmds, con = sprintf("%s/bam2fastq.cmds", mydir))
#system("parallel --jobs 45 < bam2fastq.cmds &")

#running Pipeline:
#-----------------
#execute_workflow <- sprintf("dash %s/runBasic_RNA_pipeline_on_fastqs.sh %s 10 %s &> SIS1025d.basicSCRNAseq.stdouterr &", scriptsdir, wkdir_d, scriptsdir)
#system(execute_workflow)
#setwd(wkdir_d)
#system("mkdir -p vcf")
#setwd("vcf")
#make_bamlist <- sprintf("ls %s/gatk/*.bam > bamfiles.list", wkdir_d)
#system(make_bamlist)
#gen_genotype_cmds <- sprintf("bash %s/RPE-1_GRCh38_Genotype_etai.sh ${PWD}/bamfiles.list %s/vcf/SIS1025d.gatkBamFiles 1", scriptsdir, wkdir_d)
#system(gen_genotype_cmds)
#run_genotype_cmds <- sprintf("parallel --jobs 8 < %s/vcf/SIS1025d.gatkBamFiles_RPE_hets_GT.UG_jobs.list &> %s/vcf/SIS1025d.gatkBamFiles_RPE_hets_GT.UG_jobs.stdouterr", wkdir_d, wkdir_d)
#system(run_genotype_cmds)
#RNA_genotype <- sprintf("Rscript %s/prepareGenotypeDataFromRNAvcfs2.R 14 %s/vcf %s/vcf/SIS1025d.gatkBamFiles_RPE_hets.GT_UG.%s.vcf", scriptsdir, wkdir_d, wkdir_d, '%s')
#system(RNA_genotype)

#SIS1025e

##################################################################
#Bash basic genotye pipeline cmds (from fastq to TPM and SNP RDs):
##################################################################

#trasforming bam files to fastqs:
#--------------------------------
#first option
#setwd(wkdir_e)
#system("mkdir -p bams")
#setwd("bams")
#system("cp /singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025e/bams/*.bam  ${PWD}/")
#system("parallel -j 16 samtools sort -n {} -o ${PWD}/{/}.qsort.bam ::: /papathanasiou/SIS1025e/*.bam &")
#setwd(wkdir_e)
#system("mkdir -p fastqs")
#setwd("fastqs")
#bamtofastq <- sprintf("parallel -j 22 bedtools bamtofastq -i {} -fq ${PWD}/{/}.R1.fastq -fq2 ${PWD}/{/}.R2.fastq  ::: %s/bams/*.bam &", wkdir_e)
#system(bamtofastq)
#system("parallel -j 22 gzip {}  :::  ${PWD}/*.fastq")
#second option:
#mydir <- "/singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025e"
#system("cp -s /papathanasiou/SIS1025e/*.bam /singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025e/bams")
#bamFiles <- list.files(sprintf("%s/bams", mydir), pattern = "*.bam$", full.names = T)
#cmds <- sprintf("bash ~ejacob/WORK/RNA/Workflows/bam2fastq.pairedend.sh %s", bamFiles)
#writeLines(text = cmds, con = sprintf("%s/bam2fastq.cmds", mydir))
#system("parallel --jobs 45 < bam2fastq.cmds &")

#running Pipeline:
#-----------------
#execute_workflow <- sprintf("dash %s/runBasic_RNA_pipeline_on_fastqs.sh %s 10 %s &> SIS1025e.basicSCRNAseq.stdouterr &", scriptsdir, wkdir_e, scriptsdir)
#system(execute_workflow)
#setwd(wkdir_e)
#system("mkdir -p vcf")
#setwd("vcf")
#make_bamlist <- sprintf("ls %s/gatk/*.bam > bamfiles.list", wkdir_e)
#system(make_bamlist)
#gen_genotype_cmds <- sprintf("bash %s/RPE-1_GRCh38_Genotype_etai.sh ${PWD}/bamfiles.list %s/vcf/SIS1025e.gatkBamFiles 1", scriptsdir, wkdir_e)
#system(gen_genotype_cmds)
#run_genotype_cmds <- sprintf("parallel --jobs 8 < %s/vcf/SIS1025e.gatkBamFiles_RPE_hets_GT.UG_jobs.list &> %s/vcf/SIS1025e.gatkBamFiles_RPE_hets_GT.UG_jobs.stdouterr", wkdir_e, wkdir_e)
#system(run_genotype_cmds)
#RNA_genotype <- sprintf("Rscript %s/prepareGenotypeDataFromRNAvcfs2.R 14 %s/vcf %s/vcf/SIS1025e.gatkBamFiles_RPE_hets.GT_UG.%s.vcf", scriptsdir, wkdir_e, wkdir_e, '%s')
#system(RNA_genotype)


########################################################################
########################################################################
########################################################################
########################################################################
########################################################################





##########
#basic QC:
##########
source("/homes10/ejacob/WORK/Misc/mine_old/SEQ2MAT/R/RNAseqUtils.R")
require(readxl)
library(data.table)
anno <- data.table(read_excel("/singlecellcenter/etai/ExperimentsData/Stamatis/Stamatis_list_v14_181025.xlsx"))

rnaseqOutDir <- sprintf("%s/star/", wkdir_d)
expr1 <- collectAllSTARCountsFromFileDir(workdir = rnaseqOutDir)
rnaseqOutDir <- sprintf("%s/star/", wkdir_e)
expr2 <- collectAllSTARCountsFromFileDir(workdir = rnaseqOutDir)

qc <- cbind(expr1$qc, expr2$qc)
counts <- cbind(expr1$counts, expr2$counts)
tmp <- rbindlist(lapply(1:length(anno$Fastq_files), function(x)  {
  res <- grep(pattern = sprintf("_%s_", anno$Fastq_file[x]), ignore.case = T, x = colnames(qc))
  data.table(id = anno$WTA.plate[x], 
             idx = ifelse(length(res) > 0, res, -1)) }))

colnames(qc)[tmp[idx>0]$idx] <- tmp[idx>0]$id
colnames(counts)[tmp[idx>0]$idx] <- tmp[idx>0]$id
qc2 <- round(sweep(qc, 2, colSums(qc), FUN = "/")*100)
qcbygene <- rbind(colSums(counts>0), colSums(counts>4), colSums(counts>9))
rownames(qcbygene) <- c("th1", "th5", "th10")

qc3 <- t(rbind(qc2, qcbygene, totalReads=colSums(qc)))
write.csv(x = qc3, file = sprintf("%s/QC.csv", wkdir_d))
write.csv(x = qc3, file = sprintf("%s/QC.csv", wkdir_e))
qcs <- data.table(qc3, keep.rownames = "id")

#Total Gene Counts with Jinyu cells
#hist(qcs$geneCounts, main="", xlab="Gene Counts", cex.lab=1.5, cex.axis=1.4, col="darkgreen", xlim=c(0,100), breaks=24)

#Gene Counts of just Jinyu cells
#hist(qcs[grep("Jinyu", qcs$id),]$geneCounts, main="", xlab="Gene Counts", cex.lab=1.5, cex.axis=1.4, col="darkgreen", xlim=c(0,100), breaks=24)

#Histogram with out Jinyu cells
hist(qcs[-grep("Jinyu", qcs$id),]$geneCounts, main="", xlab="Gene Counts", cex.lab=1.2, cex.axis=1.4, col="darkgreen", xlim=c(0,100), breaks=24)

#######################################################################
#save.image("/pellmanlab/stam_niko/etai_code/experiments/Save_states/post_qc_state.RData")
#load("/pellmanlab/stam_niko/etai_code/experiments/Save_states/post_qc_state.RData")
#######################################################################

#######################################################################
#Collection of total expression data (controls and samples to predict):
#######################################################################

#New lib data collection:
#------------------------
#SIS1025d
myrsem <- collectAllRSEMresultsFromDir(dir = sprintf("%s/rsem", wkdir_d))
saveRDS(myrsem, sprintf("%s/rsemoutput.rds", wkdir_d))
myrsemtpm1 <- myrsem$genes$TPM

#SIS1025e
myrsem <- collectAllRSEMresultsFromDir(dir = sprintf("%s/rsem", wkdir_e))
saveRDS(myrsem, sprintf("%s/rsemoutput.rds", wkdir_e))
myrsemtpm2 <- myrsem$genes$TPM

stopifnot(sum(rownames(myrsemtpm1) != rownames(myrsemtpm2)) == 0)

myrsemtpm <- cbind(myrsemtpm1, myrsemtpm2)
mytmpids <- do.call(rbind, lapply(rownames(myrsemtpm), function(x) str_split(string = x, pattern = "_", n = 2)[[1]]))[,1]
mytmpids[grep("_PAR_Y", x = rownames(myrsemtpm))] <- paste(mytmpids[grep("_PAR_Y", x = rownames(myrsemtpm))], "_PAR_Y", sep = "")
rownames(myrsemtpm) <- mytmpids

#Old lib data Binding:
#------------------------

rawTPMfile <- "/singlecellcenter/etai/ExperimentsData/Stamatis/Stamatis_list_v12_180404.rsemtpm.rds"
prev_rsemtpm <- readRDS(rawTPMfile)

stopifnot(sum(rownames(myrsemtpm) != rownames(prev_rsemtpm)) == 0)

rsemtpm <- cbind(prev_rsemtpm, myrsemtpm)

tmp <- rbindlist(lapply(1:length(anno$Fastq_files), function(x)  {
  res <- grep(pattern = sprintf("_%s_", anno$Fastq_file[x]), ignore.case = T, x = colnames(rsemtpm))
  data.table(id = anno$WTA.plate[x], 
             idx = ifelse(length(res) > 0, res, -1)) }))

colnames(rsemtpm)[tmp[idx>0]$idx] <- tmp[idx>0]$id

#all tpm we have so far:
saveRDS(rsemtpm, file = sprintf("%s/../Stamatis_list_v14_181025.rsemtpm.rds", wkdir_d))

############################################
# Control Samples Collection for Annotations !!!!!! INCLUDE IN EXPERIMENT COMBINING STEP
############################################

require(readxl)
anno <- data.table(read_excel("/singlecellcenter/etai/ExperimentsData/Stamatis/May072018/Stamatis_list_v12_180508.QCRNAv3.withQCRNAv3.xlsx"))
controlSampleIDs <- anno[ (mainGroup == "A_level_control" | mainGroup == "B_level_control") & th1 > 6000 ]$WTA.plate
controlSampleIDs <- controlSampleIDs[controlSampleIDs %in% colnames(rsemtpm)]

controlSampleIDs2 <- anno[ (mainGroup == "A_level_control") & th1 > 6000 ]$WTA.plate
controlSampleIDs2 <- controlSampleIDs2[controlSampleIDs2 %in% colnames(rsemtpm)]

anno <- data.table(read_excel("/singlecellcenter/etai/ExperimentsData/Stamatis/Stamatis_list_v14_181103.xlsx"))

#############################################################################################################################
#Collection of SNP data (controls and samples to predict) and phasing them based on genotype data and expression (biallelic):
#############################################################################################################################

source('/homes10/ejacob/WORK/secondaryanalysis/Stamatis/G2/runRNAseqPipelineFor_BF9VP.utils.R')

m <- readRDS("/singlecellcenter/etai/ExperimentsData/Stamatis/180803_M01209_0114_000000000-BF9VP/reports/haplotypePhasedADs.rds")

#Load and phase data from specific lib:
#--------------------------------------
mym1 <- collectAllVCFdata(patternFile = sprintf("%s/vcf/alleles.ADs.v5.chr%s.rds", wkdir_e, '%s'))
mym2 <- collectAllVCFdata(patternFile = sprintf("%s/vcf/alleles.ADs.v5.chr%s.rds", wkdir_d, '%s'))

stopifnot(sum(mym2$features$id != mym1$features$id)==0)
mym <- list(A = cbind(mym1$A, mym2$A), B = cbind(mym1$B, mym2$B), features = mym1$features)

tmp <- rbindlist(lapply(1:length(anno$Fastq_files), function(x)  {
  res <- grep(pattern = sprintf("_%s_", anno$Fastq_file[x]), ignore.case = T, x = colnames(mym$A))
  data.table(id = anno$WTA.plate[x], 
             idx = ifelse(length(res) > 0, res, -1)) }))

colnames(mym$A)[tmp[idx>0]$idx] <- tmp[idx>0]$id

tmp <- rbindlist(lapply(1:length(anno$Fastq_files), function(x)  {
  res <- grep(pattern = sprintf("_%s_", anno$Fastq_file[x]), ignore.case = T, x = colnames(mym$B))
  data.table(id = anno$WTA.plate[x], 
             idx = ifelse(length(res) > 0, res, -1)) }))

colnames(mym$B)[tmp[idx>0]$idx] <- tmp[idx>0]$id

mym <- phaseMyAllelesBasedOnHaplotypes(m = mym)

saveRDS(object = mym, file = sprintf("%s/haplotypePhasedADs.SIS1025d-e.rds", wkdir_e))


############################################################################
# Data Restructuring of Variants and Gene Expr (TPM) into format for Models:
############################################################################

#Prepare data for GSSM: features totalExp, A, B, minor allele freq, gene interdist
#------------------------------------------------------------------------------------

#Allele level data:
alleles.all <- prepareAlleleData(m = m, mym = mym, useFraction = F, controlSampleIDs = controlSampleIDs)
saveRDS(alleles.all, file = sprintf("%s/alleles.all.rds", wkdir_e))

alleles.all.frac <- prepareAlleleData(m = m, mym = mym, useFraction = T, controlSampleIDs = controlSampleIDs)

#data preparation for CNV estimates:
#----------------------------------
source('/homes10/ejacob/WORK/secondaryanalysis/methods_paper_results/R/CNML.R')

#Total expression with +1:
#-------------------------
load("/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.geneRanges.RData")
names(geneRanges)[1] <- "id"

tpmdata <- fromRawTPMtoMLInput(rsemtpm = rsemtpm, geneRanges = geneRanges[,1:4], controlSampleIDs = controlSampleIDs2, max_pos_col = 4, plusOne = 1, min_sd = 0.0001,
                               maxExp = 1200, quantileOfExp = 0.99, minDetectionLevel = 5, minNumOfSamplesToDetect = 5, 
                               geneInterspaceLogBase = 8, doNormalizeByControls = F)
adt <- tpmdata$adt
IDs <- names(which(colSums(rsemtpm>0)> 2000))
saveRDS(adt, file = sprintf("%s/../adt.rds", wkdir_e))

#Preparing data for CN calculations:
tmp <- log2(adt[,-(1:4)] + 1) #adt[,-(1:4)] #
M <- copy(tmp[, controlSampleIDs2, with = F] + 1)
q <- min(log2(1200), unname(quantile(as.matrix(M), 0.99)))
M[M > q] <- q

myMeans <- rowMeans(M, na.rm = T)
adt2 <- cbind(adt[,1:4], sweep(tmp, 1, myMeans, "-"))
adt <- adt2
adt <- adt[order(seqnames, start)]
require(gtools)
adt <- adt[mixedorder(seqnames)]

#Total expression with +0:
#-------------------------
#tpmdata is the same for +1 or +0 since we do not normalize
#now we add to log2 input + 1 works good when allele level is A and B and not An and Bn
adt0 <- tpmdata$adt
M <- copy(adt0[, controlSampleIDs2, with = F])
q <- min(1200, unname(quantile(as.matrix(M), 0.99)))
M[M > q] <- q

myMeans <- rowMeans(M, na.rm = T)
adt0 <- cbind(adt0[,1:4], interspace = tpmdata$interspace, log2(sweep(adt0[,-(1:4)], 1, myMeans, "/") + 1))
adt0 <- adt0[order(seqnames, start)]

#####################
#Temporary? QC plots:
#####################
#mid tmp plots:

require(signal)
runNewLibQuick <- function(IDs, fname = "qc_SIS1025e_5q.pdf", chr = "chr5") {
  pdf(file = fname)
  for(i in IDs) {
    message(i)
    myalleles <- mergeTotalExpAndAlleleData(adt, A = alleles.all$An, B = alleles.all$Bn, id = i)
    runMyMedianAndSSMFiltering(x = adt[seqnames == chr][[i]], 
                               x2 = tpmdata$interspace[adt$seqnames == chr], 
                               #A = alleles.all.frac$A[seqnames == chr][[i]], B = alleles.all.frac$B[seqnames == chr][[i]], abpos = alleles.all.frac$A[seqnames == chr]$start, 
                               A = alleles.all$An[seqnames == chr][[i]], B = alleles.all$Bn[seqnames == chr][[i]], abpos = alleles.all$An[seqnames == chr]$start, 
                               myalleles = myalleles[seqnames == chr],
                               t = adt[seqnames == chr]$start, 
                               main = i, xlab = "Genomic Position")  
    
  }
  
  dev.off()
}

runNewLibQuick0 <- function(IDs, fname = "qc_SIS1025e_5q.pdf", chr = "chr5") {
  pdf(file = fname)
  for(i in IDs) {
    message(i)
    myalleles <- mergeTotalExpAndAlleleData(adt0, A = alleles.all$A, B = alleles.all$B, id = i)
    runMyMedianAndSSMFiltering(x = adt0[seqnames == chr][[i]], 
                               x2 = tpmdata$interspace[adt0$seqnames == chr], 
                               #A = alleles.all.frac$A[seqnames == chr][[i]], B = alleles.all.frac$B[seqnames == chr][[i]], abpos = alleles.all.frac$A[seqnames == chr]$start, 
                               A = alleles.all$A[seqnames == chr][[i]], B = alleles.all$B[seqnames == chr][[i]], abpos = alleles.all$A[seqnames == chr]$start, 
                               myalleles = myalleles[seqnames == chr],
                               t = adt0[seqnames == chr]$start, 
                               main = i, xlab = "Genomic Position")  
    
  }
  
  dev.off()
}

chrsToExcludeFromNormalization = c("chrX") #, "chrM", "chr10")

#____________ BROKEN DOWN __________________

##########################################################
# CN Predictions
##########################################################

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



#---------------------------------------------- plots -----------------------------------------------#
#targ:
table(anno$treatment)

plotdir_d = sprintf("mkdir -p %s/../plots/SIS1025d", wkdir_d)
plotdir_e = sprintf("mkdir -p %s/../plots/SIS1025e", wkdir_e)
system(plotdir_d)
system(plotdir_e)

targMN5q <- anno[treatment == "targMN"]$WTA.plate
targMN5q <- targMN5q[targMN5q %in% colnames(SSM)]
pdf(file = sprintf("%s/../plots/SIS1025d/targMN_5q.pdf", wkdir_d))
for(i in targMN5q) {
  message(i)
  plotGenomeCN.SSM(id = i, adt = SSM[seqnames == "chr5"])
  abline(v = 65600208, lty = "dashed", col = "green")

}
dev.off()

targMN5q <- anno[treatment == "targMN"]$WTA.plate
targMN5q <- targMN5q[targMN5q %in% colnames(SSM)]
pdf(file = sprintf("%s/../plots/SIS1025d/targMN_5q.WGT.pdf", wkdir_d), width = 28, height = 38)  
for(i in targMN5q) {
  
par(mfrow=c(5,1))
myid <- i
message(mydt)
plotGenomeCN.SSM(id = myid, adt = adt, controlSampleIDs = controlSampleIDs2, main = sprintf("%s - TPM Ratio", myid), ylim = c(-10,10))
plotGenomeCN.SSM(id = myid, adt = MA, controlSampleIDs = controlSampleIDs2, main = sprintf("%s - MA", myid))
plotGenomeCN.SSM(id = myid, adt = SSM, controlSampleIDs = controlSampleIDs2, main = sprintf("%s - SSM", myid))
plotGenomeCN.SSM(id = myid, adt = MVSSM, controlSampleIDs = controlSampleIDs2, main = sprintf("%s - MV-SSM", myid))

mypreds <- readRDS(sprintf("/singlecellcenter/etai/tmp/CNML/pointOfChange/run_id17/output/%s_predicted_by_id17.rds", myid))
plotGenomeCN(mypreds, main = sprintf("%s - ML", myid))
}
dev.off()



targMN4q <- anno[treatment == "targMN, 4q_710"]$WTA.plate
targMN4q <- targMN4q[targMN4q %in% colnames(SSM)]
pdf(file = sprintf("%s/../plots/SIS1025d/targMN_4q.pdf", wkdir_d))
for(i in targMN4q) {
  message(i)
  plotGenomeCN.SSM(id = i, adt = SSM[seqnames == "chr4"])
  #abline(v = 65600208, lty = "dashed", col = "green")
  
}
dev.off()


pdf(file = sprintf("%s/../plots/SIS1025d/chr4_controls.pdf", wkdir_d))
for(i in controlSampleIDs2) {
  message(i)
  plotGenomeCN.SSM(id = i, adt = SSM[seqnames == "chr4"])
  #abline(v = 65600208, lty = "dashed", col = "green")
  
}
dev.off()


targMNXa <- anno[treatment == "targMN, Xa_721"]$WTA.plate
targMNXa <- targMNXa[targMNXa %in% colnames(SSM)]
pdf(file = "/singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025d/targMN_Xa.pdf")  
for(i in targMNXa) {
  message(i)
  plotGenomeCN.SSM(id = i, adt = SSM[seqnames == "chrX"])
  #abline(v = 65600208, lty = "dashed", col = "green")
  
}
dev.off()


targMNX <- anno[treatment == "targMN, X_955"]$WTA.plate
targMNX
targMNX <- targMNX[targMNX %in% colnames(SSM)]

pdf(file = sprintf("%s/../plots/SIS1025d/targMN_X.pdf", wkdir_d)) 
for(i in targMNX) {
  message(i)
  plotGenomeCN.SSM(id = i, adt = SSM[seqnames == "chrX"])
  #abline(v = 65600208, lty = "dashed", col = "green")
  
}
dev.off()




