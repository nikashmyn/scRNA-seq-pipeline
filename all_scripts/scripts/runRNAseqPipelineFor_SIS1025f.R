scriptsdir <- "/pellmanlab/stam_niko/etai_code/DFCI.scRNAseq.workflows/scripts"
wkdir_f <- "/pellmanlab/stam_niko/etai_code/experiments_etai/SIS1025f"
mk_wkdir_f <- sprintf("mkdir -p /pellmanlab/stam_niko/etai_code/experiments_etai/SIS1025f", wkdir_f)
system(mk_wkdir_f)

#SIS1025f

##################################################################
#Bash basic genotye pipeline cmds (from fastq to TPM and SNP RDs):
##################################################################

#trasforming bam files to fastqs:
#--------------------------------
#first option
#system("mkdir -p /singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025f/bams")
#system("cd /singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025f/bams")
#system("cp /papathanasiou/SIS1025f/* ./")
#system("parallel -j 26 samtools sort -n {} ${PWD}/{/}.qsort ::: *.bam &")
#system("rm -f *.demult.bam")
#system("rm -f *.md5")
#system("cd ..")
#system("mkdir fastqs")
#system("cd fastqs")
#system("parallel -j 18 bedtools bamtofastq -i {} -fq ${PWD}/{/}.R1.fastq -fq2 ${PWD}/{/}.R2.fastq  ::: /singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025f/bams/*.bam &")
#system("parallel -j 22 gzip {}  ::: ${PWD}/*.fastq")


#running Pipeline:
#-----------------
#system("bash ~/WORK/RNA/Workflows/runBasic_RNA_pipeline_on_fastqs.sh /singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025f/ 8 *bam.R1.fastq.gz >& SIS1025f.basicSCRNAseq.stdouterr &")
### till here ###
#system("parallel -j 10 bash ~ejacob/WORK/RNA/Workflows/CollectMultipleMetricsForMultipleBams.sh {}  ::: ${PWD}/star/*.Aligned.sortedByCoord.out.bam &") #needed here because of error

#system("cd /singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025f")
#system("mkdir /singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025f/vcf")
#system("cd /singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025f/vcf")
system("ls /pellmanlab/stam_niko/data/processed_bam/SIS1025f/VarCall_BAMs/*.bam > /pellmanlab/stam_niko/data/processed_bam/SIS1025f/bamfiles.list")
#Using CNV12.0: genotype_file=/czlab/Data/RPE-1_Genotyping_GRCh38/v.1/RPE-1.hets.vcf.gz:
#system("bash ~/WORK/RNA/Workflows/RPE-1_GRCh38_Genotype_etai.sh ${PWD}/bamfiles.list /singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025f/vcf/SIS1025f.gatkBamFiles 1")
#if command line is < 65535:
#system("parallel --jobs 8 < SIS1025f.gatkBamFiles_RPE_hets_GT.UG_jobs.list >& SIS1025f.gatkBamFiles_RPE_hets_GT.UG_jobs.stdouterr &")
#else if > 65535
#system("cat SIS1025f.gatkBamFiles_RPE_hets_GT.UG_jobs.list | parallel --jobs 8 --pipe -N1 bash >& SIS1025f.gatkBamFiles_RPE_hets_GT.UG_jobs.stdouterr &")

#system("Rscript ~/WORK/RNA/Workflows/prepareGenotypeDataFromRNAvcfs.R 6 /singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025f/vcf /singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025f/vcf/SIS1025f.gatkBamFiles_RPE_hets.GT_UG.%s.vcf")



##########
#basic QC:
##########

#on my Mac:
#rsync -Pvr ~/Data/ExperimentalData/Stamatis/RNA/Stamatis_list_v15_190910.xlsx ejacob@lynx.dfci.harvard.edu:/singlecellcenter/etai/ExperimentsData/Stamatis
  
source("/homes10/ejacob/WORK/Misc/mine_old/SEQ2MAT/R/RNAseqUtils.R")
require(readxl)
anno <- data.table(read_excel("/singlecellcenter/etai/ExperimentsData/Stamatis/Stamatis_list_v15_190910.xlsx"))

rnaseqOutDir <- sprintf("%s/star/", wkdir_f)
expr1 <- collectAllSTARCountsFromFileDir(workdir = rnaseqOutDir)

#rnaseqOutDir <- "/singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025f/star/"
#expr2 <- collectAllSTARCountsFromFileDir(workdir = rnaseqOutDir)

qc <- expr1$qc #cbind(expr1$qc, expr2$qc)
counts <- expr1$counts #cbind(expr1$counts, expr2$counts)
#two lanes:
qc.l1 <- qc[, grep("Lane1", colnames(qc))]
qc.l2 <- qc[, grep("Lane2", colnames(qc))]
counts.l1 <- counts[, grep("Lane1", colnames(counts))]
counts.l2 <- counts[, grep("Lane2", colnames(counts))]


#lane1
tmp.l1 <- rbindlist(lapply(1:length(anno$Fastq_files), function(x)  {
  res <- grep(pattern = sprintf("_%s_", anno$Fastq_file[x]), ignore.case = T, x = colnames(qc.l1))
  data.table(id = anno$WTA.plate[x], 
             idx = ifelse(length(res) > 0, res, -1)) }))
colnames(qc.l1)[tmp.l1[idx>0]$idx] <- tmp.l1[idx>0]$id
colnames(counts.l1)[tmp.l1[idx>0]$idx] <- tmp.l1[idx>0]$id

#lane2
tmp.l2 <- rbindlist(lapply(1:length(anno$Fastq_files), function(x)  {
  res <- grep(pattern = sprintf("_%s_", anno$Fastq_file[x]), ignore.case = T, x = colnames(qc.l2))
  data.table(id = anno$WTA.plate[x], 
             idx = ifelse(length(res) > 0, res, -1)) }))
colnames(qc.l2)[tmp.l2[idx>0]$idx] <- tmp.l2[idx>0]$id
colnames(counts.l2)[tmp.l2[idx>0]$idx] <- tmp.l2[idx>0]$id

qc <- as.matrix(qc.l1) + as.matrix(qc.l2)
qc2 <- round(sweep(qc, 2, colSums(qc), FUN = "/")*100)

#sanity check
sum(colnames(qc2.l1) != colnames(qc2.l2))
sum(colnames(qc.l1) != colnames(qc.l2))
sum(colnames(counts.l1) != colnames(counts.l2))

counts <- as.matrix(counts.l1) + as.matrix(counts.l2)
saveRDS(counts, sprintf("%s/counts.rds", wkdir_f))

qcbygene <- rbind(colSums(counts>0), colSums(counts>4), colSums(counts>9))
rownames(qcbygene) <- c("th1", "th5", "th10")

qc3 <- t(rbind(qc2, qcbygene, totalReads=colSums(qc)))
qcs <- data.table(qc3, keep.rownames = "id")
write.csv(x = qcs, file = sprintf("%s/QC.csv", wkdir_f))
hist(qcs$geneCounts, main="", xlab="Gene Counts", cex.lab=1.5, cex.axis=1.4, col="darkgreen", xlim=c(0,100), breaks=34)


#on my mac:
#system("rsync -Pvr ejacob@lynx.dfci.harvard.edu:/singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025f/QC.csv ~/Data/ExperimentalData/Stamatis/RNA/SIS1025f.QC.csv")



# hist(qcs[grep("Jinyu", qcs$id),]$geneCounts, main="", xlab="Gene Counts", cex.lab=1.5, cex.axis=1.4, col="darkgreen", xlim=c(0,100), breaks=24)
# hist(qcs[-grep("Jinyu", qcs$id),]$geneCounts, main="", xlab="Gene Counts", cex.lab=1.5, cex.axis=1.4, col="darkgreen", xlim=c(0,100), breaks=24)
# 
# require(readxl)
# anno <- data.table(read_excel("/singlecellcenter/etai/ExperimentsData/Stamatis/May072018/Stamatis_list_v12_180508.QCRNAv3.withQCRNAv3.xlsx"))
# controlSampleIDs <- anno[ (mainGroup == "A_level_control" | mainGroup == "B_level_control") & th1 > 6000 ]$WTA.plate
# controlSampleIDs <- controlSampleIDs[controlSampleIDs %in% colnames(rsemtpm)]
# 
# controlSampleIDs2 <- anno[ (mainGroup == "A_level_control") & th1 > 6000 ]$WTA.plate
# controlSampleIDs2 <- controlSampleIDs2[controlSampleIDs2 %in% colnames(rsemtpm)]
# 
# anno <- data.table(read_excel("/singlecellcenter/etai/ExperimentsData/Stamatis/Stamatis_list_v14_181103.xlsx"))




#######################################################################
#Collection of total expression data (controls and samples to predict):
#######################################################################

#New lib data collection:
#------------------------
#SIS1025f
myrsem <- collectAllRSEMresultsFromDir(dir = sprintf("%s/rsem", wkdir_f))
saveRDS(myrsem, sprintf("%s/rsemoutput.rds", wkdir_f))
myrsemtpm.l1 <- myrsem$genes$TPM[, grep("Lane1", colnames(myrsem$genes$TPM))]
myrsemtpm.l2 <- myrsem$genes$TPM[, grep("Lane2", colnames(myrsem$genes$TPM))]
stopifnot(sum(rownames(myrsemtpm.l1) != rownames(myrsemtpm.l2)) == 0)
#lane1
tmp <- rbindlist(lapply(1:length(anno$Fastq_files), function(x)  {
  res <- grep(pattern = sprintf("_%s_", anno$Fastq_file[x]), ignore.case = T, x = colnames(myrsemtpm.l1))
  data.table(id = anno$WTA.plate[x], 
             idx = ifelse(length(res) > 0, res, -1)) }))
colnames(myrsemtpm.l1)[tmp[idx>0]$idx] <- tmp[idx>0]$id
saveRDS(myrsemtpm.l1, sprintf("%s/rsemtpm.lane1.rds", wkdir_f))
#lane2
tmp <- rbindlist(lapply(1:length(anno$Fastq_files), function(x)  {
  res <- grep(pattern = sprintf("_%s_", anno$Fastq_file[x]), ignore.case = T, x = colnames(myrsemtpm.l2))
  data.table(id = anno$WTA.plate[x], 
             idx = ifelse(length(res) > 0, res, -1)) }))
colnames(myrsemtpm.l2)[tmp[idx>0]$idx] <- tmp[idx>0]$id
saveRDS(myrsemtpm.l2, sprintf("%s/rsemtpm.lane2.rds", wkdir_f))

stopifnot(sum(colnames(myrsemtpm.l1) != colnames(myrsemtpm.l2)) == 0)
#merge two lanes:
myrsemtpm <- (myrsemtpm.l1 + myrsemtpm.l2)/2
saveRDS(myrsemtpm, sprintf("%s/rsemtpm.rds", wkdir_f))

#TEMP
myrsemtpm <- readRDS(sprintf("%s/rsemtpm.rds", wkdir_f))
#TEMP

mytmpids <- do.call(rbind, lapply(rownames(myrsemtpm), function(x) str_split(string = x, pattern = "_", n = 2)[[1]]))[,1]
mytmpids[grep("_PAR_Y", x = rownames(myrsemtpm))] <- paste(mytmpids[grep("_PAR_Y", x = rownames(myrsemtpm))], "_PAR_Y", sep = "")
rownames(myrsemtpm) <- mytmpids

#merging it with previous rsem data: "%s/../Stamatis_list_v14_181025.rsemtpm.rds
prev_rsemtpm <- readRDS(sprintf("%s/../Stamatis_list_v14_181025.rsemtpm.rds", wkdir_f))
#stopifnot(sum(rownames(myrsemtpm) != rownames(prev_rsemtpm)) == 0)
rsemtpm <- cbind(prev_rsemtpm, myrsemtpm)

#Save all tpm we have so far:
saveRDS(rsemtpm, file = sprintf("%s/../Stamatis_list_v15_190910.rsemtpm.rds", wkdir_f))




#############################################################################################################################
#Collection of SNP data (controls and samples to predict) and phasing them based on genotype data and expression (biallelic):
#############################################################################################################################

source('/homes10/ejacob/WORK/secondaryanalysis/Stamatis/G2/runRNAseqPipelineFor_BF9VP.utils.R')

#Load and phase data from specific lib:
#--------------------------------------
source('/homes10/ejacob/WORK/RNA/Workflows/utils/collectAllVCFdata.R')
mym <- collectAllVCFdata(patternFile = "/singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025f/vcf/alleles.ADs.v5.chr%s.rds")
#phase SNPs:
source('/homes10/ejacob/WORK/RNA/Workflows/utils/phaseMyAllelesBasedOnHaplotypes.R')
mym <- phaseMyAllelesBasedOnHaplotypes(m = mym)
stopifnot(sum(prev_m$features$id != mym$features$id)==0)

#aggregating two lanes together:
#-------------------------------
mym.l1 <- list(A = mym$A[, grep("Lane1", colnames(mym$A)), with=F],
               B = mym$B[, grep("Lane1", colnames(mym$B)), with=F])
stopifnot(sum(colnames(mym.l1$A) != colnames(mym.l1$B)) == 0)
mym.l2 <- list(A = mym$A[, grep("Lane2", colnames(mym$A)), with=F],
               B = mym$B[, grep("Lane2", colnames(mym$B)), with=F])
stopifnot(sum(colnames(mym.l2$A) != colnames(mym.l2$B)) == 0)


#lane1
tmp <- rbindlist(lapply(1:length(anno$Fastq_files), function(x)  {
  res <- grep(pattern = sprintf("_%s_", anno$Fastq_file[x]), ignore.case = T, x = colnames(mym.l1$A))
  data.table(id = anno$WTA.plate[x], 
             idx = ifelse(length(res) > 0, res, -1)) }))
colnames(mym.l1$A)[tmp[idx>0]$idx] <- tmp[idx>0]$id
colnames(mym.l1$B)[tmp[idx>0]$idx] <- tmp[idx>0]$id

#lane2
tmp <- rbindlist(lapply(1:length(anno$Fastq_files), function(x)  {
  res <- grep(pattern = sprintf("_%s_", anno$Fastq_file[x]), ignore.case = T, x = colnames(mym.l2$A))
  data.table(id = anno$WTA.plate[x], 
             idx = ifelse(length(res) > 0, res, -1)) }))
colnames(mym.l2$A)[tmp[idx>0]$idx] <- tmp[idx>0]$id
colnames(mym.l2$B)[tmp[idx>0]$idx] <- tmp[idx>0]$id

stopifnot(sum(colnames(mym.l1$A) != colnames(mym.l2$A)) == 0)
stopifnot(sum(colnames(mym.l1$B) != colnames(mym.l2$B)) == 0)

#replacing NAs by zeros:
mym.l1$A[is.na(mym.l1$A)] <- 0
mym.l1$B[is.na(mym.l1$B)] <- 0
mym.l2$A[is.na(mym.l2$A)] <- 0
mym.l2$B[is.na(mym.l2$B)] <- 0

mym2 <- list(A = mym.l1$A + mym.l2$A, B = mym.l1$B + mym.l2$B, features = mym$features)
#going back to NA representation:
mym2$A[mym2$A == 0] <- NA
mym2$B[mym2$B == 0] <- NA
saveRDS(object = mym2, file = sprintf("%s/haplotypePhasedADs.SIS1025f.rds", wkdir_f))

#Merge current run with previous data:
#-------------------------------------
#loading previous allele data (already phased):
#These data also includes biallelicFilter based on the function: loadAndPhaseSNPdataFromRPE1ExperimentsBasedOnHaplotypeAndExpression.R in workflows
prev_m <- readRDS("/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/ASE.rds")
stopifnot(sum(mym2$features$id != prev_m$features$id) == 0)

m3 <- list()
m3$features <- prev_m$features #using the features from prev_m since they include biallelicFilter info
m3$A <- cbind(prev_m$A, mym2$A)
m3$B <- cbind(prev_m$B, mym2$B)

saveRDS(object = m3, file = sprintf("%s/../Stamatis_list_v15_190910.ASE.rds", wkdir_f))

#2. prepare data for GSSM: features totalExp, A, B, minor allele freq, gene interdist
#------------------------------------------------------------------------------------

dirpath <- "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2"

controlSampleIDs2 <- readRDS(sprintf("%s/data/controlSampleIDs2.rds", dirpath))
controlSampleIDs <- readRDS(sprintf("%s/data/controlSampleIDs.rds", dirpath))
geneSpecific <- readRDS(file = sprintf("%s/data/geneSpecific.rds", dirpath))
rsemtpm <- readRDS(file = sprintf("%s/data/rsemtpm.rds", dirpath))
ASE <- readRDS(file = sprintf("%s/data/ASE.rds", dirpath))
coding <- readRDS(sprintf("%s/data/coding_snps.rds", dirpath))
configs <- readRDS(sprintf("%s/data/param_config_list.rds", dirpath))
geneRanges <- readRDS(sprintf("%s/data/geneRanges.rds", dirpath))
adt.old <- readRDS("/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/adt.na.rds")


m <- readRDS("/singlecellcenter/etai/ExperimentsData/Stamatis/Stamatis_list_v15_190910.ASE.rds")

#Allele level data:
alleles.all <- prepareAlleleData(m = m, useFraction = F, controlSampleIDs = controlSampleIDs)
saveRDS(alleles.all, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025f/alleles.all.rds")

alleles.all.frac <- prepareAlleleData(m = m, useFraction = T, controlSampleIDs = controlSampleIDs)
saveRDS(alleles.all.frac, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025f/alleles.all.frac.rds")

####### till here ##########


#data preparation for CNV estimates:
#----------------------------------
source('~/WORK/secondaryanalysis/methods_paper_results/R/CNML.R')


#Total expression with +1:
#-------------------------
load("/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.geneRanges.RData")
names(geneRanges)[1] <- "id"

tpmdata <- fromRawTPMtoMLInput(rsemtpm = rsemtpm, geneRanges = geneRanges[,1:4], controlSampleIDs = controlSampleIDs2, max_pos_col = 4, plusOne = 1, min_sd = 0.0001,
                               maxExp = 1200, quantileOfExp = 0.99, minDetectionLevel = 5, minNumOfSamplesToDetect = 5, 
                               geneInterspaceLogBase = 8, doNormalizeByControls = F)
adt <- tpmdata$adt
IDs <- names(which(colSums(rsemtpm>0)> 2000))
saveRDS(adt, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/reports.CNV12.0/data/adt.rds")

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


#QC plots:
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

runNewLibQuick(IDs = targ5qids, fname = "QC_SIS1025e_5q_3.pdf", chr = "chr5")
runNewLibQuick0(IDs = targ5qids, fname = "QC_SIS1025e_5q_4.pdf", chr = "chr5")


chrsToExcludeFromNormalization = c("chrX") #, "chrM", "chr10")

##########################################################
# CN Predictions
##########################################################


#SSM +1 estimates:
#-------------------
IDs <- names(which(colSums(rsemtpm>0)> 2000))
cpus <- 42
require(snow)
cl <- snow::makeCluster(cpus, type = "SOCK")
clusterEvalQ(cl, library(data.table))
clusterEvalQ(cl, source('~/WORK/secondaryanalysis/methods_paper_results/R/CNML.R'))

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
saveRDS(SSM, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/Stamatis_list_v14_181025.rsemtpm.SSM.GM.rds")
saveRDS(SSM, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/reports.CNV12.0/data/SSM1.rds")

#SSM +0 estimates:
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

saveRDS(SSM0, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/reports.CNV12.0/data/SSM0.rds")

#MV-SSM calculations +1
#----------------------
cpus <- 42
require(snow)
cl <- snow::makeCluster(cpus, type = "SOCK")
clusterEvalQ(cl, library(data.table))
clusterEvalQ(cl, source('~/WORK/secondaryanalysis/methods_paper_results/R/CNML.R'))
clusterEvalQ(cl, source('~/WORK/secondaryanalysis/Stamatis/G2/runRNAseqPipelineFor_BF9VP.utils.R'))

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

saveRDS(MVSSM, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/Stamatis_list_v14_181025.rsemtpm.MVSSM.GM.rds")
saveRDS(MVSSM, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/reports.CNV12.0/data/MVSSM1.TEIS.rds")

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

saveRDS(MA, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/reports.CNV12.0/data/MA.rds")

#calculation of IMF:
#-------------------
medfilterRes <- runIterativeMedianFilterOnAllSamples(dt2 = adt, CPUs = 42)
mats.imf <- getAllMatricesForItrMedFilter(myres.all = medfilterRes)
mats.imf2 <- lapply(mats.imf, function(x) cbind(adt[,1:4], x))
IMF <- mats.imf2
saveRDS(IMF, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/Stamatis_list_v14_181025.rsemtpm.IMF.GM.rds")
dtm <- IMF$level6[, lapply(.SD, median),
          .SDcols = -c(1:4), by = "seqnames"]

myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in%chrsToExcludeFromNormalization), -1]))
IMF6 <- cbind(IMF$level6[, 1:4], 2^sweep(IMF$level6[,-(1:4)], 2, myCellMedians, "-"))
saveRDS(IMF6, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/reports.CNV12.0/data/IMF6.rds")


#calculation of ML:
#------------------
mlIDs <- gsub(pattern = "_predicted_by_id17.rds", replacement = "", x = list.files(path = "/singlecellcenter/etai/tmp/CNML/pointOfChange/run_id17/output/", pattern = "*predicted_by_id17.rds"))
tmp <- lapply(mlIDs, function(id) 
  readRDS(sprintf("/singlecellcenter/etai/tmp/CNML/pointOfChange/run_id17/output/%s_predicted_by_id17reg.rds", id)))
names(tmp) <- mlIDs

tmp2 <- do.call(cbind, lapply(tmp, function(x) x$prediction))
tmp3 <- cbind(tmp[[1]][,1:4], tmp2)
MLREG <- tmp3
saveRDS(MLREG, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/reports.CNV12.0/data/MLREG.rds")


#calculation of ML regression:
#------------------
mlIDs <- gsub(pattern = "_predicted_by_id17reg.rds", replacement = "", x = list.files(path = "/singlecellcenter/etai/tmp/CNML/pointOfChange/run_id17/output/", pattern = "*predicted_by_id17reg.rds"))
tmp <- lapply(mlIDs, function(id) 
  readRDS(sprintf("/singlecellcenter/etai/tmp/CNML/pointOfChange/run_id17/output/%s_predicted_by_id17.rds", id)))
names(tmp) <- mlIDs

tmp2 <- do.call(cbind, lapply(tmp, function(x) x$prediction))
tmp3 <- cbind(tmp[[1]][,1:4], tmp2)
ML <- tmp3
saveRDS(ML, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/reports.CNV12.0/data/ML.rds")

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

cpus <- 45
require(snow)
cl <- snow::makeCluster(cpus, type = "SOCK")
clusterEvalQ(cl, library(data.table))
clusterEvalQ(cl, library(KFAS))
clusterEvalQ(cl, source('~/WORK/secondaryanalysis/methods_paper_results/R/CNML.R'))
clusterEvalQ(cl, source('~/WORK/secondaryanalysis/Stamatis/G2/runRNAseqPipelineFor_BF9VP.utils.R'))

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
saveRDS(MVSSM0.AB, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/Stamatis_list_v14_181025.rsemtpm.MVSSM0.AB.rds")

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
saveRDS(MVSSM0.TE, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/reports.CNV12.0/data/MVSSM0.TEAB.TE.rds")
saveRDS(MVSSM0.A, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/reports.CNV12.0/data/MVSSM0.TEAB.A.rds")
saveRDS(MVSSM0.B, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/reports.CNV12.0/data/MVSSM0.TEAB.B.rds")


#running for +1
tt1 <- parLapply(cl, IDs, function(i) rbindlist(lapply(as.character(unique(adt$seqnames)), 
                                                      function(chr, i) getMyAlleleLevelSSMs(i = i, alleles.all = alleles.all, adt = adt, chr = chr), i)))
names(tt1) <- IDs
MVSSM.AB <- tt1
saveRDS(MVSSM.AB, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/Stamatis_list_v14_181025.rsemtpm.MVSSM.AB.rds")

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
saveRDS(MVSSM1.TE, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/reports.CNV12.0/data/MVSSM1.TEAB.TE.rds")
saveRDS(MVSSM1.A, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/reports.CNV12.0/data/MVSSM1.TEAB.A.rds")
saveRDS(MVSSM1.B, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/reports.CNV12.0/data/MVSSM1.TEAB.B.rds")




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

cpus <- 25
require(snow)
cl <- snow::makeCluster(cpus, type = "SOCK")
clusterEvalQ(cl, library(data.table))
clusterEvalQ(cl, library(KFAS))
clusterEvalQ(cl, source('~/WORK/secondaryanalysis/methods_paper_results/R/CNML.R'))
clusterEvalQ(cl, source('~/WORK/secondaryanalysis/Stamatis/G2/runRNAseqPipelineFor_BF9VP.utils.R'))

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
saveRDS(MVSSM0.TEA, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/Stamatis_list_v14_181025.rsemtpm.MVSSM0.TEA.rds")

tt.B <- parLapply(cl, IDs, function(i) rbindlist(lapply(as.character(unique(adt0$seqnames)), 
                                                        function(chr, i) getMyOneAlleleLevelSSMs0(i = i, alleles.all = alleles.all, adt0 = adt0, chr = chr, allele = "B"), i)))
names(tt.B) <- IDs
MVSSM0.TEB <- tt.B
saveRDS(MVSSM0.TEB, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/Stamatis_list_v14_181025.rsemtpm.MVSSM0.TEB.rds")


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
saveRDS(MVSSM0.TEA, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/reports.CNV12.0/data/MVSSM0.TEA.rds")
saveRDS(MVSSM0.TEB, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/reports.CNV12.0/data/MVSSM0.TEB.rds")



#running for +1:
#---------------
tt.A <- parLapply(cl, IDs, function(i) rbindlist(lapply(as.character(unique(adt$seqnames)), 
                                                        function(chr, i) getMyOneAlleleLevelSSMs(i = i, alleles.all = alleles.all, adt = adt, chr = chr, allele = "A"), i)))
names(tt.A) <- IDs
MVSSM1.TEA <- tt.A
saveRDS(MVSSM1.TEA, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/Stamatis_list_v14_181025.rsemtpm.MVSSM1.TEA.rds")

tt.B <- parLapply(cl, IDs, function(i) rbindlist(lapply(as.character(unique(adt$seqnames)), 
                                                        function(chr, i) getMyOneAlleleLevelSSMs(i = i, alleles.all = alleles.all, adt = adt, chr = chr, allele = "B"), i)))
names(tt.B) <- IDs
MVSSM1.TEB <- tt.B
saveRDS(MVSSM1.TEB, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/Stamatis_list_v14_181025.rsemtpm.MVSSM1.TEB.rds")


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
saveRDS(MVSSM1.TEA, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/reports.CNV12.0/data/MVSSM1.TEA.rds")
saveRDS(MVSSM1.TEB, file = "/singlecellcenter/etai/ExperimentsData/Stamatis/reports.CNV12.0/data/MVSSM1.TEB.rds")



stopCluster(cl)





#---------------------------------------------- plots -----------------------------------------------#
#targ:
table(anno$treatment)



targMN5q <- anno[treatment == "targMN"]$WTA.plate
targMN5q <- targMN5q[targMN5q %in% colnames(SSM)]
pdf(file = "/singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025d/targMN_5q.pdf")  
for(i in targMN5q) {
  message(i)
  plotGenomeCN.SSM(id = i, adt = SSM[seqnames == "chr5"])
  abline(v = 65600208, lty = "dashed", col = "green")

}
dev.off()

targMN5q <- anno[treatment == "targMN"]$WTA.plate
targMN5q <- targMN5q[targMN5q %in% colnames(SSM)]
pdf(file = "/singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025d/targMN_5q.WGT.pdf", width = 28, height = 38)  
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
pdf(file = "/singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025d/targMN_4q.pdf")  
for(i in targMN4q) {
  message(i)
  plotGenomeCN.SSM(id = i, adt = SSM[seqnames == "chr4"])
  #abline(v = 65600208, lty = "dashed", col = "green")
  
}
dev.off()


pdf(file = "/singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025d/chr4_controls.pdf")  
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

pdf(file = "/singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025d/targMN_X.pdf")  
for(i in targMNX) {
  message(i)
  plotGenomeCN.SSM(id = i, adt = SSM[seqnames == "chrX"])
  #abline(v = 65600208, lty = "dashed", col = "green")
  
}
dev.off()

