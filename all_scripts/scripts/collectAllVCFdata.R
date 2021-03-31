collectAllVCFdata <- function(patternFile = "/singlecellcenter/etai/ExperimentsData/Stamatis/May072018/vcfs/alleles.ADs.v5.chr%s.rds",
                              mychrs = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22, "X")) {
  require(data.table)
  M <- lapply(sprintf(patternFile, mychrs), function(f) { message("Working on file: ", f); return(readRDS(f)) })
  names(M) <- sprintf("chr%s", mychrs)
  A <- rbindlist(lapply(M, function(x) data.table(x$A)))
  B <- rbindlist(lapply(M, function(x) data.table(x$B)))
  features <- rbindlist(lapply(M, function(x) data.table(x$features)))
  
  return(list(A = A, B = B, features = features))
}
