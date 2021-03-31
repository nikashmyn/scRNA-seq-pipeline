#plot_raw_data_and_prediction_boxplots.R

#need to load source('~/WORK/Papers/MLpaper/R/utilities/class_prob_to_cn_level.R')

plot_raw_data_and_prediction_boxplots2 <- function(myid = "170512_B4", chr = "chr5", th = 10, 
                                                  adt, nonzeros.zs, coding, preds, eps = 0, controlSampleIDs, doPlot = T) {
  require(ggplot2)
  require(matrixStats)
  
  
  #controls:
  columns <- c("id", "seqnames", "start", "end", controlSampleIDs)
  tpmc <- rsemtpm[names(which(rowSums(rsemtpm[, controlSampleIDs]>0)> th)) ,controlSampleIDs] 
  tpmc[tpmc>0] <- 1
  cnonzeros <- round(rowMeans(sweep(tpmc, 2, colMeans(tpmc), "/")), digits = 3) #round(rowMeans(tpmc), digits = 3)
  
  #sample
  columns <- c("id", "seqnames", "start", "end", myid)
  #using adjusted probabilities:
  ml <- cbind(preds[, c("id", "seqnames", "start", "end")], MLreg = class_prob_to_cn_level(myid = myid, probs = preds, chr = chr, plotme = F, add_anno = F))
  myYlab.p1 <- "Adjusted CN prob."
  
  #return(ml)  
  
  
  message("All tpm = ", nrow(adt[seqnames==chr]))
  tpm <- adt[id %in% names(which(rowSums(rsemtpm[, controlSampleIDs]>0)> th)) ,columns, with=F] #mlinput$adt[ ,columns, with=F]
  message("tpm > ", th, " = ", nrow(tpm[seqnames==chr]))
  
  #setnames(ml, myid, "MLreg")
  setnames(tpm, myid, "tpm")
  tpm[id %in% names(which(rsemtpm[, myid] == 0))]$tpm <- NA
  
  message("median of tpm ratio: ", median(tpm[, lapply(.SD, median, na.rm=T), .SDcols="tpm", by = "seqnames"]$tpm))
  tpm$tpm <- tpm$tpm - median(tpm[, lapply(.SD, median, na.rm = T), .SDcols="tpm", by = "seqnames"]$tpm)
  #tmp <- tpm$tpm - median(tpm[, lapply(.SD, median, na.rm = T), .SDcols="tpm", by = "seqnames"]$tpm)
  #plot(tmp, tpm$tpm)
  mlt <- merge(ml, tpm, by = c("id", "seqnames", "start", "end"))
  mlt <- mlt[order(seqnames, start)]
  
  bns = unlist(lapply(1:floor(nrow(mlt[seqnames==chr])/50+0.5), rep, 50))[1:nrow(mlt[seqnames==chr])]
  bns[is.na(bns)] <- last(bns[!is.na(bns)])
  mlt <- cbind(mlt[seqnames==chr], bin = bns)
  
  #return(mlt)
  binRanges <- Reduce(merge, 
                      list(mlt[, lapply(.SD, min), 
                               .SDcols = c("start"), by = bin], 
                           mlt[, lapply(.SD, max), 
                               .SDcols = c("end"), by = bin],
                           mlt[, lapply(.SD, function(x) mean(x, na.rm = T)), 
                               .SDcols = c("tpm"), by = bin],
                           mlt[, lapply(.SD, function(x) sd(x, na.rm = T)), 
                               .SDcols = c("tpm"), by = bin]))
  setnames(binRanges, c("start", "end", "tpm.x", "tpm.y"), c("binStart", "binEnd", "mean", "sd"))
  binRanges$pos <- binRanges$binStart + (round((binRanges$binEnd - binRanges$binStart)/2))
  binRanges$proportion <- (binRanges$binEnd - binRanges$binStart)/sum(binRanges$binEnd - binRanges$binStart)
  mlt1 <- merge(mlt[seqnames==chr], binRanges, by="bin")
  #return(mlt1)
  #plot(mlt1$pos, mlt1$MLreg)
  
  tt1 <- mlt1[, c("bin", "pos", "MLreg", "proportion")]
  tt1$method="MLreg"
  posvalues <- tt1$pos
  tt1$pos <- factor(tt1$pos)
  setnames(tt1, "MLreg", "value")
  # if(chr == "chrX") { no need for that
  #   message("changing chrX values.")
  #   #tt1$value <- tt1$value/2
  # }
  tt2 <- mlt1[, c("bin", "pos", "tpm", "proportion")]
  tt2$method="TPM Ratio"
  tt2$pos <- factor(tt2$pos)
  #return(tt2)
  #tt2$pos <- tt2$pos + 1e5
  setnames(tt2, "tpm", "value")
  print((mean(colSums(rsemtpm[, controlSampleIDs]>0))/sum(rsemtpm[, myid]>0)))
  tt2$value <- 2^(tt2$value)
  
  
  #alleles analysis:
  #-----------------
  
  
  #coding SNPs - allele A:
  aggs.A.zs <- coding$aggs.A
  
  #normalizing by sample
  aggs.A.zs <- cbind(aggs.A.zs[,1:2], sweep(aggs.A.zs[, -c(1:2)], 2, colMeans(as.matrix(aggs.A.zs[, -c(1:2)]))+eps, FUN = "/"))
  #return(aggs.A.zs)
  #normalizing by bin
  aggs.A.zs <- cbind(aggs.A.zs[,1:2], sweep(aggs.A.zs[, -c(1:2)], 1, rowMedians(as.matrix(aggs.A.zs[, controlSampleIDs, with=F]))+eps, FUN = "/"))
  
  totalCntsA <- cbind(coding$cnts.A[, c(1:2)], meanSites = rowMeans(coding$cnts.A[, controlSampleIDs, with=F]))
  wi <- totalCntsA[seqnames == chr & meanSites < 3]$bin
  
  aggsA <- merge(unique(tt1[method == "MLreg", c("bin", "pos", "proportion")]), aggs.A.zs[seqnames == chr][, c("bin", myid), with=F])
  aggsAc <- (merge(unique(tt1[method == "MLreg", c("bin", "pos", "proportion")]), coding$cnts.A[seqnames == chr][, c("bin", myid), with=F]))
  names(aggsA)[4] <- "fracNonZero"
  names(aggsAc)[4] <- "cntsNonZero"
  if(length(wi)>0) {
    aggsA$fracNonZero[aggsA$bin == wi] = NA
    aggsAc$cntsNonZero[aggsA$bin == wi] = NA
  }
  aggsA <- merge(aggsA, aggsAc, by = c("bin", "pos"))
  
  aggsA$method = "Frac Allele A"
  
  
  tmp <- aggs.A.zs[seqnames == chr, controlSampleIDs, with=F]
  tmp <- melt(t(tmp))[,2:3]
  colnames(tmp) <- c("bin", "value")
  if(length(wi)>0)
    tmp$value[tmp$bin == wi] = NA
  
  tmp2A <- merge(unique(tt1[method == "MLreg", c("bin", "pos", "proportion")]), tmp)
  
  tmp2A$method = "Ctrl. Frac A"
  
  #coding SNPs - allele B:
  aggs.B.zs <- coding$aggs.B
  #normalizing by sample
  aggs.B.zs <- cbind(aggs.B.zs[,1:2], sweep(aggs.B.zs[, -c(1:2)], 2, colMeans(as.matrix(aggs.B.zs[, -c(1:2)]))+eps, FUN = "/"))
  #normalizing by bin
  aggs.B.zs <- cbind(aggs.B.zs[,1:2], sweep(aggs.B.zs[, -c(1:2)], 1, rowMedians(as.matrix(aggs.B.zs[, controlSampleIDs, with=F]))+eps, FUN = "/"))
  
  totalCntsB <- cbind(coding$cnts.B[, c(1:2)], meanSites = rowMeans(coding$cnts.B[, controlSampleIDs, with=F]))
  wi <- totalCntsB[seqnames == chr & meanSites < 3]$bin
  
  aggsB <- merge(unique(tt1[method == "MLreg", c("bin", "pos", "proportion")]), aggs.B.zs[seqnames == chr][, c("bin", myid), with=F])
  aggsBc <- merge(unique(tt1[method == "MLreg", c("bin", "pos", "proportion")]), coding$cnts.B[seqnames == chr][, c("bin", myid), with=F])
  names(aggsB)[4] <- "fracNonZero"
  names(aggsBc)[4] <- "cntsNonZero"
  if(length(wi)>0) {
    aggsB$fracNonZero[aggsB$bin == wi] = NA
    aggsBc$cntsNonZero[aggsB$bin == wi] = NA
  }
  aggsB <- merge(aggsB, aggsBc, by = c("bin", "pos"))
  aggsB$method = "Frac Allele B"
  
  tmp <- aggs.B.zs[seqnames == chr, controlSampleIDs, with=F]
  tmp <- melt(t(tmp))[,2:3]
  colnames(tmp) <- c("bin", "value")
  if(length(wi)>0)
    tmp$value[tmp$bin == wi] = NA
  tmp2B <- merge(unique(tt1[method == "MLreg", c("bin", "pos", "proportion")]), tmp)
  tmp2B$method = "Ctrl. Frac B"
  
  
  aggsAB <- rbind(aggsA, aggsB)
  #return(aggsAB)
  tmp2AB <- rbind(tmp2A, tmp2B)
  
  aggsAB$pos <- factor(aggsAB$pos)
  tmp2AB$pos <- factor(tmp2AB$pos)
  
  aggsAB$method <- factor(aggsAB$method)
  tmp2AB$method <- factor(tmp2AB$method)
  
  #frac non zeros
  #--------------
  
  
  #normalizing by bin
  mynonzeros <- nonzeros.zs
  mynonzeros0 <- cbind(mynonzeros[,1:2], sweep(mynonzeros[, -c(1:2)], 1, rowMedians(as.matrix(mynonzeros[, controlSampleIDs, with=F]))+eps, FUN = "-"))
  
  
  mynonzeros <- merge(unique(tt1[method == "MLreg", c("bin", "pos", "proportion")]), mynonzeros0[seqnames == chr][, c("bin", myid), with=F])
  names(mynonzeros)[4] <- "fracNonZero"
  mynonzeros$method = "Sample Non Zeros"
  mynonzeros$pos <- factor(mynonzeros$pos)
  
  tmp <- mynonzeros0[seqnames == chr, controlSampleIDs, with=F]
  tmp <- melt(t(tmp))[,2:3]
  colnames(tmp) <- c("bin", "value")
  tmp2 <- merge(unique(tt1[method == "MLreg", c("bin", "pos", "proportion")]), tmp)
  
  if(chr == "chrX") {
    tmp2AB[method == "Ctrl. Frac A", value := NA]
    aggsAB[method == "Frac Allele A", fracNonZero := NA]
    aggsAB[method == "Frac Allele A", cntsNonZero := NA]
  }
  #return(list(tmp2AB = tmp2AB, aggsAB = aggsAB, posvalues = posvalues, tt1=tt1))
  
  require(gridExtra)
  
  rng <- range(centromeres[seqnames==chr, c("start", "end")])
  
  myf <- function(x) {
    r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  
  require(dplyr)
  
  tt2b = tt2 %>% 
    group_by(bin) %>% 
    summarise(pos = unique(pos), #.groups = "drop_last", #Use this to kill warning message
              method = unique(method),
              proportion = mean(proportion),
              median = median(value, na.rm=T), 
              se.min = quantile(value, probs = 0.25, na.rm=T), 
              se.max = quantile(value, probs = 0.75, na.rm=T))
  
  if(!doPlot) {
    allelebins <- data.table(aggsAB[method == "Frac Allele A", c("bin", "pos")],
                             A = aggsAB[method == "Frac Allele A"]$fracNonZero,
                             B = aggsAB[method == "Frac Allele B"]$fracNonZero)
    predbins <- tt1[, lapply(.SD, median, na.rm=T), .SDcols="value", by = "bin"]
    setnames(predbins, "value", "pred")
    tpmbins <- data.table(tt2b)[, c("bin", "proportion", "median")]
    setnames(tpmbins, "median", "tpm")
    nonzerosbins <- mynonzeros[, c("bin", "fracNonZero")]
    mydata <- list(predbins = predbins, tpmbins = tpmbins, nonzerosbins = nonzerosbins, allelebins = allelebins)
    mydata <- Reduce(merge, mydata)
    setcolorder(mydata, neworder = c("bin", "pos", "proportion"))
    return(mydata)
  }
  
  p1 <- ggplot(tt1, aes(x=pos, y=value, fill = method)) + ylim(c(0,3.2)) + ylab(myYlab.p1) + theme_bw() +
    geom_vline(xintercept=as.numeric(unique(tmp2AB$pos)), linetype="solid", color = alpha("gray", alpha = 0.4), size = 4) +
    geom_vline(xintercept=min(which(as.numeric(as.character(unique(tmp2AB$pos))) >= rng[1])), linetype="solid", color = alpha("black", alpha = 0.4), size = 5) +
    
    # geom_boxplot(position=position_dodge(0), fill = NA, colour = "darkred", size = 1., outlier.shape = NA,
    #              show.legend = F, varwidth = TRUE, alpha=0.6, aes(weight=proportion)) + 
    stat_summary(fun.data = myf, geom="boxplot", position = position_dodge(2), fill = NA, colour = "darkred", size = 1., outlier.shape = NA,
                 show.legend = F, varwidth = TRUE, alpha=0.6, aes(width=tt1$proportion/max(tt1$proportion)/1.7 + 0.15)) +
    
    ggtitle(sprintf("%s - %s\nMachine Learning Predictions", myid, chr)) +
    geom_hline(yintercept=c(0, 0.5, 1.0, 1.5, 2, 3), linetype="dashed", color = "black", size = c(0.4, 0.4, 0.4, 0.4, 1.05, 0.4)) +
    scale_x_discrete(breaks=posvalues, labels = as.character(round(posvalues/1e6))) + 
    #xlab(label = sprintf("Position (Mbp) %s", chr)) +
    theme(text = element_text(size=15), axis.text.x = element_text(angle=90, hjust=1), axis.title.x=element_blank())
  
  #return(tt2)
  #return(tt2b) 
  
  p2 <- ggplot(tt2b, aes(x = pos, y = median, ymin = se.min, ymax = se.max, fatten = 3+proportion/max(proportion)*15)) + #, fatten = 6+proportion/max(proportion)*15)) + #, fatten = .75 for little to no circles
    geom_vline(xintercept=tt2b$bin, linetype="solid", color = alpha("gray", alpha = 0.4), size = 4) +   
    geom_vline(xintercept=min(which(as.numeric(as.character(unique(tmp2AB$pos))) >= rng[1])), 
               linetype="solid", color = alpha("black", alpha = 0.4), size = 5) +
    geom_hline(yintercept=c(0, 0.5, 1.0, 1.5, 2), linetype="dashed", color = "black", size = c(0.7, 0.7, 1.5, 0.7, 0.7)) +
    #geom_errorbar(color = "darkred", aes(x = pos, ymin = median,  ymax = median, width = proportion/max(proportion)), size = 1) + #for horizontal error lines
    geom_pointrange(color = "darkred", aes(y = median, ymin = se.min,  ymax = se.max), size = 1) + 
    geom_linerange(color = alpha("darkred", alpha = 0.3)) +
    #geom_errorbarh(color = "darkred", aes(x = pos, y = median, xmin = pos - (proportion/max(proportion)/1.7 + 0.15)/2,  xmax = pos + (proportion/max(proportion)/1.7 + 0.15)/2), size = 1) + 
    ggtitle(sprintf("Normalized Ratio of Expressed Genes")) + ylim(c(0, min(5, max(tt2b$se.max)))) + 
    ylab("Normalized ratio") +  theme_gray(base_size=16) +
    scale_x_discrete(breaks=posvalues, labels = as.character(round(posvalues/1e6))) + 
    #xlab(label = sprintf("Position (Mbp) %s", chr)) +
    theme_bw() +
    theme(text = element_text(size=15), axis.text.x = element_text(angle=90, hjust=1), axis.title.x=element_blank())
    
  p3 <- ggplot(data = mynonzeros, aes(x=pos, y=fracNonZero)) + ylab("Standarized fraction (sdv)") +  theme_gray(base_size=16) +
    geom_vline(xintercept=as.numeric(unique(tmp2AB$pos)), linetype="solid", color = alpha("gray", alpha = 0.4), size = 4) +
    geom_vline(xintercept=min(which(as.numeric(as.character(unique(tmp2AB$pos))) >= rng[1])), linetype="solid", color = alpha("black", alpha = 0.4), size = 5) +
    ggtitle(sprintf("Fraction of Expressed Genes")) +
    stat_summary(fun.data = myf, geom="boxplot", position = position_dodge(2), fill = alpha("limegreen", alpha = 0.6), 
                 mapping = aes(x = pos, y = value, width=proportion/max(proportion)/1.7 + 0.15), data = tmp2,
                 size = 1., outlier.shape = NA,
                 show.legend = F, varwidth = TRUE, alpha=0.6)+
    
    # geom_boxplot(position=position_dodge(0), fill = alpha("limegreen", alpha = 0.6), outlier.shape = NA, 
    #              show.legend = F, varwidth = TRUE, alpha=0.6, data = tmp2, mapping = aes(x = pos, y = value, weight=proportion)) + 
    geom_point(show.legend = F, fill = alpha("darkred", alpha = 0.95), stat='identity', size=4.7, colour = "black", pch=21) +
    theme_bw() +
    scale_x_discrete(breaks=posvalues, labels = as.character(round(posvalues/1e6))) + 
    #xlab(label = sprintf("Position (Mbp) %s", chr)) +
    geom_hline(yintercept=c(-2, -1, 0, 1, 2), linetype="dashed", color = "black", size = c(0.5,0.5,1,0.5,0.5)) +
    theme(text = element_text(size=15), axis.text.x = element_text(angle=90, hjust=1), axis.title.x=element_blank())
  
  p4 <- ggplot(data = tmp2AB, mapping = aes(x = pos, y= value, fill=method)) + ylim(c(-0.3,4))+
    geom_vline(xintercept=as.numeric(unique(tmp2AB$pos)), linetype="solid", color = alpha("gray", alpha = 0.4), size = 4) +
    geom_vline(xintercept=min(which(as.numeric(as.character(unique(tmp2AB$pos))) >= rng[1])), linetype="solid", color = alpha("black", alpha = 0.4), size = 5) +
    ggtitle(sprintf("Fraction of allele specific expressed SNPs")) +
    
    # geom_boxplot(position = position_dodge(0.7), show.legend=F, varwidth=T, alpha=0.6,
    #              mapping=aes(weight=proportion+0.05, color = method, fill = NA), 
    #              outlier.shape = NA, notch=F, coef = 0, lwd = 1.1) +
    # 
    stat_summary(fun.data = myf, geom="boxplot", position = position_dodge(0.7), fill = NA, 
                 mapping = aes(x = pos, y = value, width=proportion/max(proportion)/1.7 + 0.15, color = method), 
                 size = 1., outlier.shape = NA,
                 show.legend = F, varwidth = TRUE, alpha=0.6)+
    
    geom_hline(yintercept=c(0,1,2), linetype="dashed", color = "black", size = c(0.8, 1.5, 0.8)) +
    geom_point(position = position_dodge(0.7), show.legend = F, stat='identity', size=5.5, data = aggsAB, 
               mapping = aes(x=pos, y=fracNonZero, color = method), colour = "black", pch=21) +
    geom_text(position = position_dodge(0.7), show.legend = F, stat='identity', size=2.4, data = aggsAB, colour = "white",
              mapping = aes(x=pos, y=fracNonZero, label = cntsNonZero)) +
    theme_gray(base_size=16) + ylab("Normalized Fraction") + 
    scale_fill_manual(values=c(alpha("dodgerblue", alpha = 0.8), alpha("red3", alpha = 0.8),
                               alpha("dodgerblue", alpha = 0.8), alpha("red3", alpha = 0.8))) + 
    scale_colour_manual(values=c(alpha("royalblue3", alpha = 0.95), 
                                 alpha("firebrick1", alpha = 0.95))) +
    theme_bw() +
    theme(text = element_text(size=15), axis.text.x = element_text(angle=90, hjust=1)) +
    scale_x_discrete(breaks=tt1$pos, labels = as.character(round(as.numeric(as.character(tt1$pos))/1e6))) + 
    xlab(label = sprintf("Position (Mbp) %s", chr)) 
  
  #return(list(p2, p3, p4))
  #return(p2)
  #print(p4)
  #return(grid.arrange(p2, p3, p4, nrow = 3))
  grid.arrange(p2, p3, p4, nrow = 3)
  
  #return(1)
  #return(mynonzeros0)
  #return(p4)
  
  #print(grid.arrange(p1, p2, p3, p4, nrow = 4))
  
  #dev.off()
  
  #return(p)
  
}


