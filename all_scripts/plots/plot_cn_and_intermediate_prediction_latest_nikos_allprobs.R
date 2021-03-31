#plot_cn_and_intermediate_prediction.R

require(data.table)
require(ggplot2)

plot_cn_and_intermediate_prediction <- function(probs, anno_col_numbers = c(1:6), th = 2, main = NA, controlIDs, plot_pvals = T) {
  par.before <- par()
  par(mar = c(4, 4, 4, 3.5) + 0.1)
  
  #The sample we are visualizing from func input
  myid <- probs$sample_id[1]
  #The chr we are visualizing from func input
  chr <- probs$seqnames[1]
  #Subset the predictions into a group for each arm
  myqarm <- grep(pattern = "q", unique(probs[seqnames == chr]$arm), value = T)
  myparm <- grep(pattern = "p", unique(probs[seqnames == chr]$arm), value = T)

  #get the range of the chr in genomic coordinates
  rng <- range(centromeres[seqnames==chr, c("start", "end")])
  
  #Get vectors with the predictions of each CN state
  p1 <- probs[[length(anno_col_numbers) + 1]]
  p2 <- probs[[length(anno_col_numbers) + 2]]
  p3 <- probs[[length(anno_col_numbers) + 3]]
  
  #Get the mean position of the gene 
  pos <- rowMeans(probs[, c("start", "end")])
  
  #Get a vector of the highest probability of all the CN states
  maxprob <- apply(probs[, -anno_col_numbers, with = F], 1, max)
  #Get a vector of the middle probability of all the CN states
  midpred <- apply(probs[, -anno_col_numbers, with = F], 1, function(x) order(x)[2])
  #Get vector with the CN state with the highest probability for each gene
  cns <- max.col(probs[, -anno_col_numbers, with = F])

  #Get vector with the transparency of each rectangle based on the probability
  #of the CN state using max-min normalization.
  #Transparency_p1 <- c(); Transparency_p2 <- c(); Transparency_p3 <- c();
  #for ( i in 1:length(p1) ) {
  #  x1 = p1[i]; x2 = p2[i]; x3 = p3[i]
  #  Transparency_p1 <- append( Transparency_p1, ((x1 - min(p1)) / (max(p1) - min(p1))) )
  #  Transparency_p2 <- append( Transparency_p2, ((x2 - min(p2)) / (max(p2) - min(p2))) )
  #  Transparency_p3 <- append( Transparency_p3, ((x3 - min(p3)) / (max(p3) - min(p3))) )
  #}
  
  cncols <- c(alpha("blue", 0.8), alpha("darkgreen", 0.8), alpha("red", 0.8))
  cn1col <- c(); cn2col <- c(); cn3col <- c();
  for (i in 1:length(cns)) {
    cn1col <- append(cn1col, alpha("blue", p1[i]))
    cn2col <- append(cn2col, alpha("darkgreen", p2[i]))
    cn3col <- append(cn3col, alpha("red", p3[i]))
  }
  
  plot(cns, ylim=c(0.66,3.6), xlim=c(-length(cns)*.1, length(cns)*1.1), col=alpha(cncols[cns], 0.5), 
       main = ifelse(is.na(main), sprintf("%s", myid), main),
       cex.axis=1.4, cex.lab = 1.5, pch = 15, 
       xlab = sprintf("Position in %s (Mbs)", chr), yaxt = "n", xaxt = "n", ylab = "", lwd=2)
  
  rect(xleft = c(1:length(cns)), xright = c(1:length(cns))+1, ybottom = 1.5, ytop = 2.5, 
       col = cn2col, border = NA)
  rect(xleft = c(1:length(cns)), xright = c(1:length(cns))+1, ybottom = 2.5, ytop = 3.5, 
       col = cn3col, border = NA)
  rect(xleft = c(1:length(cns)), xright = c(1:length(cns))+1, ybottom = .5, ytop = 1.5, 
       col = cn1col, border = NA)
  
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
  
  segments(x0 = 1, x1 = length(cns), y0 = 1.5, y1 = 1.5, col = "gray20", lty = "dotted", lwd = 2)
  segments(x0 = 1, x1 = length(cns), y0 = 2.5, y1 = 2.5, col = "gray20", lty = "dotted", lwd = 2)
  
  segments(x0 = 1, x1 = length(cns), y0 = 1, y1 = 1, col = "gray40", lty = "dotted", lwd = 3)
  segments(x0 = 1, x1 = length(cns), y0 = 2, y1 = 2, col = "gray40", lty = "dotted", lwd = 3)
  segments(x0 = 1, x1 = length(cns), y0 = 3, y1 = 3, col = "gray40", lty = "dotted", lwd = 3)
  
  #plot centomere:
  rect(xleft = max(which(probs$start < rng[1])), xright = min(which(probs$start > rng[2])+4), ybottom = -10, ytop = 10,
       col = alpha("black", 0.6), border = NA)

  text(x = c(0, 0), y = c(1.5, 2.5), labels = c("Mid"), pos = 2, font = 2, cex = 1.8, col = "gray20")
  
  text(x = c(length(cns) + 1), y = c(1), labels = c(sprintf("%d", round(100*length(cns[cns==1])/length(cns)) )), col = "blue", pos = 4, font = 2, cex = 1.6)
  text(x = c(length(cns) + 1), y = c(2), labels = c(sprintf("%d", round(100*length(cns[cns==2])/length(cns)) )), col = "darkgreen", pos = 4, font = 2, cex = 1.6)
  text(x = c(length(cns) + 1), y = c(3), labels = c(sprintf("%d", round(100*length(cns[cns==3])/length(cns)) )), col = "red", pos = 4, font = 2, cex = 1.6)
  
  par(mar = par.before$mar)
  #return(NA)
  #return(list(cn_percents = fracs, pvals = list(arms = pvals.arm, chrs = pvals.chr)))
  
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


