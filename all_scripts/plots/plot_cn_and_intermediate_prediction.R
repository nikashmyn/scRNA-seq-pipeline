#plot_cn_and_intermediate_prediction.R

require(data.table)
require(ggplot2)

#TODO: 
#squares = CN. circles = probs
#probs <- lgbmpreds5[sample_id == myid & seqnames == mychr]
#probs is a data.table
#centromeres <- readRDS("~ejacob/WORK/secondaryanalysis/Stamatis/data/centromeres.rds")
#p2 <- plot_cn_and_intermediate_prediction(lgbmpreds5[sample_id == "170126_A5" & seqnames == "chr1"])
#function requires the interstat and interstat.arm to be loaded to the workspace prior to calling
#function requires fracstat to be loaded to the workspace prior to calling
plot_cn_and_intermediate_prediction <- function(probs, anno_col_numbers = c(1:6), th = 2, main = NA, plot_pvals = T) {
  par.before <- par()
  # c(bottom, left, top, right) 
  par(mar = c(4, 4, 4, 3.5) + 0.1)
  myid <- probs$sample_id[1]
  chr <- probs$seqnames[1]
  myqarm <- grep(pattern = "q", unique(probs[seqnames == chr]$arm), value = T)
  myparm <- grep(pattern = "p", unique(probs[seqnames == chr]$arm), value = T)
  
  pvals.arm <- calc_and_plot_intermediate_pvalue.arm(plotme = F, myid = myid, controlSampleIDs = controlSampleIDs, myseqname = "1p")
  pvals.chr <- calc_and_plot_intermediate_pvalue.chr(plotme = F, myid = myid, controlSampleIDs = controlSampleIDs, myseqname = "chr1")
  #print((pvals.chr))
  rng <- range(centromeres[seqnames==chr, c("start", "end")])
  
  

  p1 <- probs[[length(anno_col_numbers) + 1]]
  p2 <- probs[[length(anno_col_numbers) + 2]]
  p3 <- probs[[length(anno_col_numbers) + 3]]
  
  pos <- rowMeans(probs[, c("start", "end")])
  
  maxprob <- apply(probs[, -anno_col_numbers, with = F], 1, max)
  midpred <- apply(probs[, -anno_col_numbers, with = F], 1, function(x) order(x)[2])
  cns <- max.col(probs[, -anno_col_numbers, with = F])
  #print(midpred)
  #print(cns)
  intermediate <- ifelse(maxprob < 0.67, 0.5, NA)
  
  mypreds <- cns
  mypreds[maxprob < 0.67] <- 0 
  reps <- rle(mypreds)
  #par(mfrow=c(2,1))
  #plot(mypreds, ylim = c(1,3.2))
  myStretches <- unlist(lapply(reps$lengths, function(x) if(x > th) { rep(1, x) } else { rep(0, x)}))
  
  
  
  cncols <- c(alpha("blue", 0.8), alpha("darkgreen", 0.8), alpha("red", 0.8))
  im <- mypreds*myStretches
  plot(im, ylim=c(0.66,3.6), xlim=c(-length(im)*.1, length(im)*1.1), col=alpha(cncols[cns], 0.5), 
       main = ifelse(is.na(main), sprintf("%s", myid), main),
       cex.axis=1.4, cex.lab = 1.5, pch = 15, 
       xlab = sprintf("Position in %s (Mbs)", chr), yaxt = "n", xaxt = "n", ylab = "", lwd=2)
  
  points(ifelse(cns == 1, 2 - p1, NA), col=ifelse(!is.na(intermediate), 
                            alpha("black", 0.8), alpha("black", 0.5)), 
         lwd = 1, pch = 20, cex = ifelse(!is.na(intermediate), 1.3, 0.8))
  #return(p2)
  points(ifelse(cns == 2, ifelse(midpred == 3, 3 - p2, 1 + p2), NA), col=ifelse(!is.na(intermediate), 
                            alpha("black", 0.8), alpha("black", 0.5)), 
         lwd = 1, pch = 20, cex = ifelse(!is.na(intermediate), 1.3, 0.8))
  points(ifelse(cns == 3, 2 + p3, NA), col=ifelse(!is.na(intermediate), 
                            alpha("black", 0.8), alpha("black", 0.5)), 
         lwd = 1, pch = 20, cex = ifelse(!is.na(intermediate), 1.3, 0.8))
  
  
  if(sum(im==2) > 0) {
    rect(xleft = which(im == 2), xright = which(im == 2)+1, ybottom = 1.68, ytop = 2.32, 
         col = alpha("darkgreen", 0.6), border = NA)
  }
  if(sum(im==3) > 0)  {
    rect(xleft = which(im == 3), xright = which(im == 3)+1, ybottom = 2.68, ytop = 3.32, 
         col = alpha("red", 0.6), border = NA)
  }
  if(sum(im==1) > 0) {
    rect(xleft = which(im == 1), xright = which(im == 1)+1, ybottom = 0.68, ytop = 1.32, 
         col = alpha("blue", 0.6), border = NA)
  }
  # if(sum(mypreds==0) > 0)
  #   rect(xleft = which(mypreds == 1), xright = which(mypreds == 1)+1, ybottom = 0.68, ytop = 1.32, 
  #        col = alpha("blue", 0.6), border = NA)
  
  
  # plot(cns, ylim=c(0,4), col=cncols[cns], main = sprintf("%s", myid), cex.axis=1.4, cex.lab = 1.3,
  #      xlab = sprintf("Gene index (%s)", chr), yaxt = "n", ylab = "", lwd=2)
  
  x_idx <- unique(c(seq(1, length(cns), 50), length(cns)))
  axis(1, at = x_idx,
       labels = c(round(pos[x_idx]/1e6)), 
       col.axis= "black", las=2, cex.axis = 1.2, font = 2)
  
  axis(2, at=c(1),
       labels=c("Loss"), 
       col.axis= "blue", las=0, cex.axis = 1.8, font = 2)
  axis(2, at=c(2),
       labels=c("Normal"), 
       col.axis= "darkgreen", las=0, cex.axis = 1.8, font = 2)
  axis(2, at=c(3),
       labels=c("Gain"), 
       col.axis= "red", las=0, cex.axis = 1.8, font = 2)
  # axis(4, at=c(0.5, 1.5, 2.5),
  #      labels=c(rep("Mid", 3)), 
  #      col.axis= "black", las=0, cex.axis = 1.2, font = 2)
  
  # rect(xleft = -10, xright = 1e6, ybottom = 0.34, ytop = 0.66, 
  #      col = alpha("gray60", 0.4), border = NA)
  rect(xleft = 0, xright = length(im)+1, ybottom = 1 + 0.34, ytop = 1 + 0.66, 
       col = alpha("gray60", 0.4), border = NA)
  rect(xleft = 0, xright = length(im)+1, ybottom = 2 + 0.34, ytop = 2 + 0.66, 
       col = alpha("gray60", 0.4), border = NA)
  
  segments(x0 = 1, x1 = length(cns), y0 = 1.5, y1 = 1.5, col = "gray20", lty = "dotted", lwd = 2)
  segments(x0 = 1, x1 = length(cns), y0 = 2.5, y1 = 2.5, col = "gray20", lty = "dotted", lwd = 2)
  
  segments(x0 = 1, x1 = length(cns), y0 = 1, y1 = 1, col = "gray40", lty = "dotted", lwd = 3)
  segments(x0 = 1, x1 = length(cns), y0 = 2, y1 = 2, col = "gray40", lty = "dotted", lwd = 3)
  segments(x0 = 1, x1 = length(cns), y0 = 3, y1 = 3, col = "gray40", lty = "dotted", lwd = 3)
  
  
  
  #calculating intermediates and CNs:
  im2 <- cns
  im2[midpred == 2 & cns == 3 & maxprob < 0.67 ] <- 2.5
  im2[midpred == 3 & cns == 2 & maxprob < 0.67 ] <- 2.5
  im2[midpred == 1 & cns == 2 & maxprob < 0.67 ] <- 1.5
  im2[midpred == 2 & cns == 1 & maxprob < 0.67 ] <- 1.5
  fracs <- round(100*table(factor(im2, levels = c(1,1.5,2,2.5,3)))/length(im))
  
  #plot centomere:
  rect(xleft = max(which(probs$start < rng[1])), xright = min(which(probs$start > rng[2])+4), ybottom = -10, ytop = 10,
       col = alpha("black", 0.6), border = NA)
  #plot pvals on figure:
  if(plot_pvals) {
    if(length(myparm)> 0 & nrow(pvals.arm)>0)
      text(x = 0, y = 3.5, labels = sprintf("P-val(%%%d) < %.2f", 
                                            round(pvals.arm[seqnames == myparm]$frac*100), pvals.arm[seqnames == myparm]$pval), adj = 0, cex = 1.35)
    if(length(myqarm)> 0 & nrow(pvals.arm)>0)
      text(x = length(cns), y = 3.5, labels = sprintf("P-val(%%%d) < %.2f", 
                                                      round(pvals.arm[seqnames == myqarm]$frac*100), pvals.arm[seqnames == myqarm]$pval), adj = 1, cex = 1.35)
    if(nrow(pvals.chr)>0)
      text(x = length(cns)/2, y = 3.5, labels = sprintf("P-val(%%%d) < %.2f", 
                                                        round(pvals.chr[seqnames == chr]$frac*100), pvals.chr[seqnames == chr]$pval), adj = 1, cex = 1.35)
    
    pval.2CNseg <- sum((1-fracstat[seqnames == chr][["2"]]) >= 1-fracs[["2"]]/100)/length(fracstat[seqnames == chr][["2"]])  
    text(x = c(length(cns)) + 1, y = c(1.8), labels = c(sprintf("P-val < %.3f", pval.2CNseg)), col = "darkgreen", pos = 4, font = 1, cex = 1)
    pval.3CNseg <- sum((1-fracstat[seqnames == chr][["2"]]) >= 1-fracs[["3"]]/100)/length(fracstat[seqnames == chr][["2"]])  
    text(x = c(length(cns) + 1), y = c(2.8), labels = c(sprintf("P-val < %.3f", pval.3CNseg)), col = "red", pos = 4, font = 1, cex = 1)
    pval.1CNseg <- sum((1-fracstat[seqnames == chr][["2"]]) >= 1-fracs[["1"]]/100)/length(fracstat[seqnames == chr][["2"]])  
    text(x = c(length(cns) + 1), y = c(0.8), labels = c(sprintf("P-val < %.3f", pval.1CNseg)), col = "blue", pos = 4, font = 1, cex = 1)
    
  }
  text(x = c(0, 0), y = c(1.5, 2.5), labels = c("Mid"), pos = 2, font = 2, cex = 1.8, col = "gray20")
  
  
  
  text(x = c(length(cns) + 1), y = c(1.5), labels = c(sprintf("%d", fracs["1.5"])), col = "gray20", pos = 4, font = 2, cex = 1.6)
  text(x = c(length(cns) + 1), y = c(2.5), labels = c(sprintf("%d", fracs["2.5"])), col = "gray20", pos = 4, font = 2, cex = 1.6)
  text(x = c(length(cns) + 1), y = c(1), labels = c(sprintf("%d", fracs["1"])), col = "blue", pos = 4, font = 2, cex = 1.6)
  text(x = c(length(cns) + 1), y = c(2), labels = c(sprintf("%d", fracs["2"])), col = "darkgreen", pos = 4, font = 2, cex = 1.6)
  text(x = c(length(cns) + 1), y = c(3), labels = c(sprintf("%d", fracs["3"])), col = "red", pos = 4, font = 2, cex = 1.6)
  
  par(mar = par.before$mar)
  return(list(cn_percents = fracs, pvals = list(arms = pvals.arm, chrs = pvals.chr)))
  
}

calculate_fraction_of_cns_for_chr <- function(probs, anno_col_numbers = c(1:6)) {
  maxprob <- apply(probs[, -anno_col_numbers, with = F], 1, max)
  midpred <- apply(probs[, -anno_col_numbers, with = F], 1, function(x) order(x)[2])
  cns <- max.col(probs[, -anno_col_numbers, with = F])
  
  
  im2 <- cns
  im2[midpred == 2 & cns == 3 & maxprob < 0.67 ] <- 2.5
  im2[midpred == 3 & cns == 2 & maxprob < 0.67 ] <- 2.5
  im2[midpred == 1 & cns == 2 & maxprob < 0.67 ] <- 1.5
  im2[midpred == 2 & cns == 1 & maxprob < 0.67 ] <- 1.5
  fracs <- table(factor(im2, levels = c(1,1.5,2,2.5,3)))/length(im2)
  return(fracs)
  
}


