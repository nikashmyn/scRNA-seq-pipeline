#preload of: lgbmpreds5
#preload of: adt for plot only
class_prob_to_cn_level <- function(probs = NA, myid = "170223_A5a", chr = "chr8", anno_col_numbers = c(1:6), plotme = F, add_anno = F) {
  require(data.table)
  require(zoo)
  
  if(is.na(probs))
    probs <- lgbmpreds5[sample_id == myid & seqnames == chr]
  
  maxprob <- apply(probs[, -anno_col_numbers, with = F], 1, max)
  midpred <- apply(probs[, -anno_col_numbers, with = F], 1, function(x) order(x)[2])
  cns <- max.col(probs[, -anno_col_numbers, with = F])
  
  p1 <- probs[[length(anno_col_numbers) + 1]]
  p2 <- probs[[length(anno_col_numbers) + 2]]
  p3 <- probs[[length(anno_col_numbers) + 3]]
  
  
  cn_level <- ifelse(cns == 1, 2 - p1, 0) + 
              ifelse(cns == 2, ifelse(midpred == 3, 3 - p2, 1 + p2), 0) +
              ifelse(cns == 3, 2 + p3, 0)
  
  
  if(!plotme) {
    if(add_anno) {
      cn_level <- cbind(lgbmpreds5[sample_id == myid & seqnames == chr, c("id", "seqnames", "start", "end")], pred = cn_level)
      setnames(cn_level, old = "pred", new = myid)
    }
    return(cn_level)
  }
  plot(MLREG[seqnames == "chr8"][["170223_A5a"]], class_prob_to_cn_level())
  
  myma <- rollmean(x = 2*2^adt[probs$id][[myid]], k = 50, fill = "extend")
  plot(cn_level, ylim=c(1,4), main = cor(myma, cn_level))
  points(myma, col="red")
              
  
}
