##########################
### Script Explanation ###
##########################

#Run this script as a part of Data Aggregation script to run the SSM and MA models.
#running all type of predictions for a sample

#-------------------------------------------------------
#This script was written by Nikos Mynhier and Etai Jacob.

##############################
### Load Packages and Data ###
##############################

require(zoo)
library(matrixStats)
require(gridExtra)

message("Loading data..")
#Read in Data
adt.na <- readRDS(file = sprintf("%s/aggregated_results/adt.na.rds", dirpath))
adt <- readRDS(file = sprintf("%s/aggregated_results/adt.rds", dirpath))
nonzeros <- readRDS(file = sprintf("%s/aggregated_results/nonzeros.rds", dirpath))
alleles.all <- readRDS(file = sprintf("%s/aggregated_results/alleles.all.rds", dirpath)) # "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2/data/alleles.all.rds"
ganno <- readRDS(file = sprintf("%s/aggregated_results/ganno.rds", dirpath))

#Variant matrix in the coding and UTR regions
coding <- readRDS(sprintf("%s/aggregated_results/ASE.coding.rds", dirpath))

#Standard Annotations and parameters for analysis
configs <- readRDS(sprintf("%s/param_config_list.rds", datapath))
geneRanges <- readRDS(sprintf("%s/geneRanges.rds", datapath))
centromeres <- readRDS(file = sprintf("%s/centromeres.rds", datapath))

message("Loading more data..")
controlSampleIDs2 <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs2.rds", dirpath))
controlSampleIDs <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))
rsemtpm <- readRDS(file = sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", dirpath))
#ASE <- readRDS(file = sprintf("%s/data/ASE.rds", dirpath))
nonzeros.zs <- readRDS(file = sprintf("%s/aggregated_results/nonzeros.zs.bin50.rds", dirpath))

#get samples to run based on samples in df
anno <- data.table(read.csv( sprintf("%s/work_in_progress/annotation_list.csv", datadir)))
columns <- colnames(adt)[-c(1:4)]
samples_to_use <- c(intersect(columns, anno$WTA.plate))

#get high quality samples from raw tpm values
all_QC <- readRDS(file = sprintf("%s/aggregated_results/all_QC.rds", dirpath))

#Get annotation list of all samples that are in adt (tpm ratio) object
col_anno <- data.table(read.csv( sprintf("%s/work_in_progress/annotation_list.csv", datadir)))
dim(col_anno)
col_anno <- col_anno[ WTA.plate %in% samples_to_use]
dim(col_anno)
col_anno[Pairs == "NA", Pairs := NA]

controlIDs <- col_anno[experimental_label %in% c("control_LCM_notreatment", "control_gen2_notreatment", "noco_gen2_noMN_lcm", "noco_gen1_noMN", "p53_gen2_noMN")]$WTA.plate
triIDs <- col_anno[experimental_label %in% c("Tri12_LCM", "Tri12_facs")]$WTA.plate[col_anno[experimental_label %in% c("Tri12_LCM", "Tri12_facs")]$WTA.plate %in% names(which(colSums(rsemtpm>0)>6000))]
saveRDS(controlIDs, file = sprintf("%s/aggregated_results/controlIDs.rds", dirpath))

#This scripts requires these objects ordered
adt <- adt[order(seqnames, start, end)]
adt.na <- adt.na[order(seqnames, start, end)]
nonzeros <- nonzeros[order(seqnames, start, end)]
ganno <- ganno[order(seqnames, start, end)]
stopifnot(sum(ganno$id != adt$id) == 0)

#Manual input data
par(mfrow=c(1,1))
manual_selected_vars6 <- c("geneSpecific.MEDIAN.MEDIAN.GC", "geneSpecific.MEDIAN.MEDIAN.cv.nz", "geneSpecific.MEDIAN.MEDIAN.Length", "geneSpecific.MEAN.MEAN.interspace", "TE.all.SD", "TE.onlyExpressed.MAXSTRETCH", "Frac.nonzeros.MAXSTRETCH", "TE.all.Q95", "TE.onlyExpressed.Q95", "Frac.nonzeros.SD", "TE.onlyExpressed.Q25", "TE.onlyExpressed.Q75", "TE.all.Q05", "Frac.nonzeros.MEAN", "TE.onlyExpressed.MEDIAN", "TE.onlyExpressed.MEAN", "TE.all.Q25", "TE.all.Q75", "TE.all.MEDIAN", "TE.all.MEAN")

#################################
### Compute MA and SSM Models ###
#################################

#Function for running the models
predict_CN <- function(mysample_id, outfname) {

  ########## RUNNING MA ###################################
  run_my_MA <- function(MA_winSize = 100) {
    chrsToExcludeFromNormalization = c("chrX", "chrM", "chr10")
    MA <- adt[, lapply(.SD, rollapply, width = MA_winSize, FUN=mean, partial = T, align = "center"), # fill = NA_real_ 
              .SDcols = mysample_id, by = seqnames]
    dtm <- MA[, lapply(.SD, median),
              .SDcols = -1, by = "seqnames"]
    myCellMedians <- colMedians(as.matrix(dtm[!(seqnames %in% chrsToExcludeFromNormalization), -1]))
    
    MA2 <- 2^sweep(MA[,-1], 2, myCellMedians, "-")
    
    return(MA2)
  }
  
  preds.MA <- do.call(cbind, lapply(c(50, 75, 100), run_my_MA))
  colnames(preds.MA) <- c("MA50", "MA75", "MA100")
  
  stopifnot(sum(adt$id != preds.MA$id) == 0)
  print("Done with MA")
  
  ########## RUNNING SSM (+1) for TE ###################################
  source(sprintf('%s/scripts/GSSM_utils.R', scriptsdir))

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
  print("Done with SSM")
  
  ############ binding together all data as one table ############
  
  preds <- cbind(ganno, preds.MA, SSM.TE = preds.SSM.TE)
  preds3 <- preds[order(seqnames, start, end)]
  
  saveRDS(preds3, file = outfname)
  message("Predictions saved to ", outfname)
  
}

#Run function in loop for each sample
print("Number of samples: ")
print(nrow(col_anno))
for(i in 1:nrow(col_anno)) {
  
  myid <- col_anno$WTA.plate[i]
  pairid <- col_anno$Pairs[i]
  if(!is.na(pairid)) {
    File <- sprintf("%s/%s.pair_%s.all_CN_predictions.rds", CNdir, myid, pairid)
  } else {
    File <- sprintf("%s/%s.all_CN_predictions.rds", CNdir, myid)
  }
  myid <- as.character(myid)
  predict_CN(mysample_id = myid, outfname = File)
  
}

#Print sanity check that the code ran to completion
print("Done with SSM & MA CN predictions")

#####################################
### Aggregate CN from all samples ###
#####################################

inputdir <- CNdir
#list of matrices of all predictions at the gene level
Ms <- list()

files <- list.files(path = inputdir, pattern = "*all_CN_predictions.rds", full.names = T)
files_names <- strsplit(basename(files), ".", fixed = T)
file_names <- unlist(lapply(files_names, `[[`, 1))

wta.names <- c(col_anno$WTA.plate) #this was c(col_anno[,1][[1]])
row_anno_col <- c("id", "seqnames", "arm", "start", "end")
col <- append(c("id", "seqnames", "start", "end"),  as.character(col_anno$WTA.plate)) #this col_anno was called same as above 
col <- as.character(col)
col_n <- which(col %in% colnames(coding$A))

#Select SNPs to take if they have at least 10 control samples with hits
snpstotake <- which(rowSums(coding$A[, controlSampleIDs, with=F])>=10 | rowSums(coding$B[, controlSampleIDs, with=F])>=10)

#Reduce to SNPs that are represented in the control samples
A <- as.data.frame(coding$A)[snpstotake, col_n] 
B <- as.data.frame(coding$B)[snpstotake, col_n] 
A <- data.table(A)
B <- data.table(B)

setkeyv(ganno, cols = c("seqnames", "start", "end"))
setkeyv(A, cols = c("seqnames", "start", "end"))
setkeyv(B, cols = c("seqnames", "start", "end"))
setnames(x = A, old = "id", new = "snp_id")
setnames(x = B, old = "id", new = "snp_id")

#aggregate by gene:
A2 <- foverlaps(x = ganno[, c("id", "seqnames", "start", "end", "arm")], y = A)
setcolorder(A2, neworder = c("id", "seqnames", "arm", "i.start", "i.end", "snp_id", "start", "end"))
setnames(A2, old = c("i.start", "i.end", "start", "end"), new = c("start", "end", "snp_start", "snp_end"))
A3 <- A2[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), by="id", .SDcols = -colnames(A2)[1:8]] 
A4 <- merge(ganno[, c("id", "seqnames", "arm", "start", "end")], A3, by = "id")
A4 <- A4[order(seqnames,start, end)]

B2 <- foverlaps(x = ganno[, c("id", "seqnames", "start", "end", "arm")], y = B)
setcolorder(B2, neworder = c("id", "seqnames", "arm", "i.start", "i.end", "snp_id", "start", "end"))
setnames(B2, old = c("i.start", "i.end", "start", "end"), new = c("start", "end", "snp_start", "snp_end"))
B3 <- B2[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), by="id", .SDcols = -colnames(B2)[1:8]] 
B4 <- merge(ganno[, c("id", "seqnames", "arm", "start", "end")], B3, by = "id")
B4 <- B4[order(seqnames,start, end)]

#Add allelic information to Ms (model aggregate) object
Ms$raw.A <- copy(A4)
Ms$raw.B <- copy(B4)

#note! These are not the same dimensions go back and do all cells for models if needed.
do.call(rbind, lapply(Ms, function(x) dim(x)))

#aggregating all other preds:
#----------------------------
preds_to_agg <- c( "MA50", "MA75", "MA100", "SSM.TE") 

agg_samples <- function(file, col_name = "MA50") {
  preds <- readRDS(file)
  mypreds <- preds[, col_name, with = F]
  name <- file_names[which(files == file)]
  setnames(mypreds, old = col_name, new = name)
  return(mypreds)
}

tt <- lapply(preds_to_agg, function(x)  cbind(B4[, row_anno_col, with = F], do.call(cbind, lapply(files, agg_samples, col_name = x))))
names(tt) <- preds_to_agg
Ms <- c(Ms, tt)

#verifing everything is with the same dimension:
do.call(rbind, lapply(Ms, function(x) dim(x)))
saveRDS(Ms, sprintf("%s/CN_predictions.bygene.rds", CNdir))

#now let's aggregate data by arm and chromosome and save data for further analysis:
agg_by_chr <- function(M) {
  agg <- M[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), by="seqnames", .SDcols = -row_anno_col] 
  return(agg)
}

agg_by_arm <- function(M) {
  agg <- M[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), by="arm", .SDcols = -row_anno_col] 
  return(agg)
}

Ms.chr <- lapply(Ms, agg_by_chr)
Ms.arm <- lapply(Ms, agg_by_arm)

do.call(rbind, lapply(Ms.chr, function(x) dim(x)))
saveRDS(Ms.chr, file = sprintf("%s/CN_predictions.bychr.rds", CNdir))

do.call(rbind, lapply(Ms.arm, function(x) dim(x)))
saveRDS(Ms.arm, file = sprintf("%s/CN_predictions.byarm.rds", CNdir))

