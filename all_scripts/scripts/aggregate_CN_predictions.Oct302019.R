#aggregate_CN_predictions.Oct302019.R

#this script takes all predicitions per sample and aggregate each prediction type by arm and whole chr level. 
#Then it saves everything to matrices.

#args <- c("/pellmanlab/stam_niko/data/processed_bam", "/pellmanlab/stam_niko/data/processed_bam/CN_data", "/pellmanlab/nikos/Stam_Etai_Data")
args <- commandArgs(trailingOnly = TRUE)
dirpath <- args[1]
dstdir <- args[2]
datapath <- args[3]

require(data.table)
coding <- readRDS(file = sprintf("%s/aggregated_results/ASE.coding.rds", dirpath))
controlSampleIDs2 <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs2.rds", dirpath))
controlSampleIDs <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))
ganno <- readRDS(file = sprintf("%s/aggregated_results/ganno.rds", dirpath))

inputdir <- dstdir
#list of matrices of all predictions at the gene level
Ms <- list()

files <- list.files(path = inputdir, pattern = "*.rds", full.names = T)

row_anno_col <- c("id", "seqnames", "arm", "start", "end")

preds <- readRDS(files[1])

#Raw allele aggregates:
#TODO: add allele fraction from coding$A and B at the chr and arm level

snpstotake <- which(rowSums(coding$A[, controlSampleIDs, with=F])>=10 | rowSums(coding$B[, controlSampleIDs, with=F])>=10)
length(snpstotake)

A <- copy(coding$A[snpstotake, ])
B <- copy(coding$B[snpstotake, ])

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
Ms <- c(Ms, tt)

#verifing everything is with the same dimension:
do.call(rbind, lapply(Ms, function(x) dim(x)))
saveRDS(Ms, sprintf("%s/CN_predictions.bygene.rds", dstdir))


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
saveRDS(Ms.chr, file = sprintf("%s/CN_predictions.bychr.rds", dstdir))


do.call(rbind, lapply(Ms.arm, function(x) dim(x)))
saveRDS(Ms.arm, file = sprintf("%s/CN_predictions.byarm.rds", dstdir))

