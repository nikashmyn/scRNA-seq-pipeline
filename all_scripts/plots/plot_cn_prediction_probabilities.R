#plot_cn_prediction_probabilities.R

require(data.table)
require(ggplot2)

#source('~/WORK/Papers/MNpaper/R/calc_and_plot_intermediate_pvalue.R')

#probs <- lgbmpreds5[sample_id == myid & seqnames == mychr]
#probs is a data.table
#centromeres <- readRDS("~ejacob/WORK/secondaryanalysis/Stamatis/data/centromeres.rds")
plot_cn_prediction_probabilities <- function(probs, anno_col_numbers = c(1:6), th = 0.2) {
  myid <- probs$sample_id[1]
  chr <- probs$seqnames[1]
  myqarm <- grep(pattern = "q", unique(probs[seqnames == chr]$arm), value = T)
  myparm <- grep(pattern = "p", unique(probs[seqnames == chr]$arm), value = T)
  
  pvals.arm <- calc_and_plot_intermediate_pvalue.arm(plotme = F, myid = myid, controlSampleIDs = controlSampleIDs, myseqname = "1p")
  pvals.chr <- calc_and_plot_intermediate_pvalue.arm(plotme = F, myid = myid, controlSampleIDs = controlSampleIDs, myseqname = "chr1")
  print((pvals.arm))
  rng <- range(centromeres[seqnames==chr, c("start", "end")])
  
  

  p1 <- probs[,length(anno_col_numbers) + 1, with = F]
  p2 <- probs[,length(anno_col_numbers) + 2, with = F]
  p3 <- probs[,length(anno_col_numbers) + 3, with = F]
  
  intermediate <- ifelse((abs(p2 - p1) < th & p1 > 0.2 & p2 > 0.2) | 
                         (abs(p3 - p2) < th & p3 > 0.2 & p2 > 0.2), 0.5, NA)
  intermediate <- as.vector(intermediate)
  
  cns <- max.col(probs[, -anno_col_numbers, with = F])
  #cns[!is.na(intermediate)] <- NA
  
  cncols <- c(alpha("blue", 0.8), alpha("darkgreen", 0.8), alpha("red", 0.8))

  plot(0.7 + cns*0.5, ylim=c(0,2.6), col=cncols[cns], main = sprintf("%s", myid), cex.axis=1.4, cex.lab = 1.3,
       xlab = sprintf("Gene index (%s)", chr), yaxt = "n", ylab = "", lwd=2)
  axis(2, at=c(0,0.5,1,1.2,1.7,2.2),labels=c("P = 0","P=0.5","P=1.0","CN<1", "CN=2", "CN>3"), col.axis="black", las=2, cex.axis = 1.4)
  
  if(length(which(!is.na(intermediate))) > 0)
    rect(xleft = which(!is.na(intermediate)), xright = which(!is.na(intermediate))+1, ybottom = -10, ytop = 10, 
         col = alpha("pink", 0.6), border = NA)
  
  rect(xleft = max(which(probs$start < rng[1])), xright = min(which(probs$start > rng[2])+3), ybottom = -10, ytop = 10, 
       col = alpha("black", 0.6), border = NA)
  if(length(myparm)> 0 & nrow(pvals.arm)>0)
    text(x = 0, y = 2.5, labels = sprintf("P-val(%%%d) < %.2f", round(pvals.arm[seqnames == myparm]$frac*100), pvals.arm[seqnames == myparm]$pval), adj = 0, cex = 1.25)
  if(length(myqarm)> 0 & nrow(pvals.arm)>0)
    text(x = length(cns), y = 2.5, labels = sprintf("P-val(%%%d) < %.2f", round(pvals.arm[seqnames == myqarm]$frac*100), pvals.arm[seqnames == myqarm]$pval), adj = 1, cex = 1.25)
  
  # rect(xleft = which(!is.na(intermediate))[diff(which(!is.na(intermediate))) == 1], 
  #      xright = which(!is.na(intermediate))[diff(which(!is.na(intermediate))) == 1] + 1, 
  #      ybottom = -10, ytop = 10, 
  #      col = alpha("pink", 0.6), border = NA)
  abline(h=c(1.2, 1.7, 2.2), col="gray40", lty = "dashed", lwd = 2)
  abline(h=c(0,1.0), col="gray60", lwd = 4)
  points(0.7 + cns*0.5, col=cncols[cns], main = sprintf("%s - %s", myid, chr))
  
  abline(h=0.5, col="red", lty="dotted", lwd=2)
  #points(p1, col=alpha("blue", 0.4), lwd = 1, pch = 20, cex = 0.7)
  points(p1, col=ifelse(!is.na(intermediate) & p1 > 0.2, alpha("black", 0.9), alpha("blue", 0.5)), lwd = 1, pch = 20, cex = ifelse(!is.na(intermediate) & p1 > 0.2, 0.9, 0.7))
  points(p2, col=ifelse(!is.na(intermediate) & p2 > 0.2, alpha("black", 0.9), alpha("darkgreen", 0.5)), lwd = 1, pch = 20, cex = ifelse(!is.na(intermediate) & p2 > 0.2, 0.9, 0.7))
  points(p3, col=ifelse(!is.na(intermediate) & p3 > 0.2, alpha("black", 0.9), alpha("red", 0.5)), lwd = 1, pch = 20, cex = ifelse(!is.na(intermediate) & p3 > 0.2, 0.9, 0.7))
  #points(p2, col=alpha("darkgreen", 0.4), lwd = 1, pch = 20, cex = 0.7)
  #points(p3, col=alpha("red", 0.4), lwd = 1, pch = 20, cex = 0.7)
  
  #points(intermediate, col = alpha("black", 0.4), pch=1, cex = 2)
  
  
  return(intermediate)
  
  # points(ifelse(abs(p3 - p2) < th & p3 > 0.2 & p2 > 0.2, 0.5, NA), col = alpha("black", 0.6), pch=15, lwd=5)
  # print(abs(p2 - p1))
  # print(abs(p1))
  # # points(ifelse(abs(probs[,2]-probs[,1]) < 0.2, 0.5, NA), col = alpha("magenta", 0.3), pch=15)
  # points(ifelse(abs(probs[,3]-probs[,2]) < 0.2, 0.5, NA), col = alpha("magenta", 0.3), pch=15)
  # points(predict(object = fit.lm, newdata = nd), ylim=c(0,4), main = "LM", col="seagreen")
  # points(max.col(keraspreds[sample_id == myid & seqnames == mychr][,-c(1:6)]), ylim=c(0,4), main = "DNN", col="orange")
  # #points(max.col(keraspreds2[sample_id == myid & seqnames == mychr][,-c(1:6)]), ylim=c(0,4), main = "DNN", col="brown")
  # point(max.col(lgbmpreds[sample_id == myid & seqnames == mychr][,-c(1:6)]), ylim=c(0,4), main = "DNN", col="seagreen")
  # points(max.col(lgbmpreds2[sample_id == myid & seqnames == mychr][,-c(1:6)]), ylim=c(0,4), main = "DNN", col="green")
  # points(max.col(lgbmpreds5[sample_id == myid & seqnames == mychr][,-c(1:6)]), ylim=c(0,4), main = "DNN", col="blue")
  # points(max.col(lgbmpreds6[sample_id == myid & seqnames == mychr][,-c(1:6)]), ylim=c(0,4), main = "DNN", col="magenta")
  
}
