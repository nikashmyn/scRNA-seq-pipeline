#running all type of predictions for a sample

require(data.table)
require(zoo)

message("Loading data..")
#Read in Data
adt.na <- readRDS(file = sprintf("%s/aggregated_results/adt.na.rds", dirpath))
adt <- readRDS(file = sprintf("%s/aggregated_results/adt.rds", dirpath))
nonzeros <- readRDS(file = sprintf("%s/aggregated_results/nonzeros.rds", dirpath))
alleles.all <- readRDS(file = sprintf("%s/aggregated_results/alleles.all.rds", dirpath)) # "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2/data/alleles.all.rds"
ganno <- readRDS(file = sprintf("%s/aggregated_results/ganno.rds", dirpath))

#Standard Annotations and parameters for analysis
configs <- readRDS(sprintf("%s/param_config_list.rds", datapath))
geneRanges <- readRDS(sprintf("%s/geneRanges.rds", datapath))

#Read in "given" objects. Need to figure out how to recreate these.
#geneSpecific <- readRDS(file = sprintf("%s/geneSpecific.rds", datapath))
centromeres <- readRDS(file = sprintf("%s/centromeres.rds", datapath))

message("Loading more data..")
controlSampleIDs2 <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs2.rds", dirpath))
controlSampleIDs <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))
rsemtpm <- readRDS(file = sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", dirpath))
#ASE <- readRDS(file = sprintf("%s/data/ASE.rds", dirpath))
coding <- readRDS(file = sprintf("%s/aggregated_results/ASE.coding.rds", dirpath))
high_qc_ids <- names(which(colSums(rsemtpm>5)>4000))
nonzeros.zs <- readRDS(file = sprintf("%s/aggregated_results/nonzeros.zs.bin50.rds", dirpath))

col_anno <- data.table(readRDS( sprintf("%s/work_in_progress/Annotation_list_long_vnikos_210129.rds", datapath) ))
dim(col_anno)
col_anno <- col_anno[ WTA.plate %in% high_qc_ids]
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


predict_CN <- function(mysample_id, outfname) {

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

print("Number of samples: ")
print(nrow(col_anno))
for(i in 1:nrow(col_anno)) {
#for(i in 1:5) {
  
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

#Just a Sanity check that the code ran to completion
print("Done with SSM & MA CN predictions")


#####################################
### Aggregate CN from all samples ###
#####################################

inputdir <- CNdir
#list of matrices of all predictions at the gene level
Ms <- list()

files <- list.files(path = inputdir, pattern = "*all_CN_predictions.rds", full.names = T)

row_anno_col <- c("id", "seqnames", "arm", "start", "end")
col <- append(c("id", "seqnames", "start", "end"),  col_anno$WTA.plate)

preds <- readRDS(files[1])

#Raw allele aggregates:
#TODO: add allele fraction from coding$A and B at the chr and arm level

snpstotake <- which(rowSums(coding$A[, controlSampleIDs, with=F])>=10 | rowSums(coding$B[, controlSampleIDs, with=F])>=10)
length(snpstotake)

A <- copy(coding$A[snpstotake, ..col])
B <- copy(coding$B[snpstotake, ..col])

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
  setnames(mypreds, old = col_name, new = file)
  return(mypreds)
}

tt <- lapply(preds_to_agg, function(x)  cbind(B4[, row_anno_col, with = F], do.call(cbind, lapply(files, agg_samples, col_name = x))))
names(tt) <- preds_to_agg
for ( i in 1:length(tt) ) {names(tt[[i]])[-c(1:5)] <- col_anno$WTA.plate}
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

