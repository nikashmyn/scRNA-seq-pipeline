#plot_raw_data_and_prediction_boxplots.R

#need to load source('~/WORK/Papers/MLpaper/R/utilities/class_prob_to_cn_level.R')

plot_raw_data_and_prediction_boxplots4 <- function(myfamily_file = "F258", myid_file = "F258.3", myid = "210503_9C", chr = "chr5", destDir = "/pellmanlab/stam_niko/rerun_6_9_2021/data/visual_results", th = 10, 
                                                   adt, coding, preds, eps = 0, controlSampleIDs, doPlot = T, save.rds = F) { #nonzeros.zs,
  require(ggplot2)
  require(matrixStats)
  require(ggthemes)
  
  #get tpm objects for control cells:
  columns <- c("id", "seqnames", "start", "end", controlSampleIDs)
  tpmc <- rsemtpm[names(which(rowSums(rsemtpm[, controlSampleIDs]>0)> th)) ,controlSampleIDs] 
  tpmc[tpmc>0] <- 1
  cnonzeros <- round(rowMeans(sweep(tpmc, 2, colMeans(tpmc), "/")), digits = 3) #round(rowMeans(tpmc), digits = 3)
  
  #List of columns to pull later. Just annotations and the specified cell. 
  columns <- c("id", "seqnames", "start", "end", myid)
  
  #If using ML data get predictions for specified sample
  #using adjusted probabilities:
  ml <- cbind(preds[, c("id", "seqnames", "start", "end")], MLreg = class_prob_to_cn_level(myid = myid, probs = preds, chr = chr, plotme = F, add_anno = F))
  myYlab.p1 <- "Adjusted CN prob."
  
  #get normalized tpm (adt) for specified cell (singular) w/ only genes that have control sample expression greater than specified threshold (th)
  message("All tpm = ", nrow(adt[seqnames==chr]))
  tpm <- adt[id %in% names(which(rowSums(rsemtpm[, controlSampleIDs]>0)> th)) ,columns, with=F] #mlinput$adt[ ,columns, with=F]
  message("tpm > ", th, " = ", nrow(tpm[seqnames==chr]))
  
  #Change 0 tpm values to NA for specified cell in object from just above 
  #setnames(ml, myid, "MLreg")
  setnames(tpm, myid, "tpm")
  tpm[id %in% names(which(rsemtpm[, myid] == 0))]$tpm <- NA
  
  #bin by chr taking the median tpm value. subtract median of bin values (which are medians themselves) from all values. 
  message("median of tpm ratio: ", median(tpm[, lapply(.SD, median, na.rm=T), .SDcols="tpm", by = "seqnames"]$tpm))
  tpm$tpm <- tpm$tpm - median(tpm[, lapply(.SD, median, na.rm = T), .SDcols="tpm", by = "seqnames"]$tpm)
  
  #attach annotations and order rows
  mlt <- merge(ml, tpm, by = c("id", "seqnames", "start", "end"))
  mlt <- mlt[order(seqnames, start)]
  
  #create bins for 50 consecutive genes
  bns = unlist(lapply(1:floor(nrow(mlt[seqnames==chr])/50+0.5), rep, 50))[1:nrow(mlt[seqnames==chr])]
  bns[is.na(bns)] <- last(bns[!is.na(bns)])
  mlt <- cbind(mlt[seqnames==chr], bin = bns)
  
  #Adjust annotations columns for new bins and add proportion column representing bin size in genomic coordinates.
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
  
  #make binned object for values from ML (this is used as a place holder if not plotting OLR visual)
  tt1 <- mlt1[, c("bin", "pos", "MLreg", "proportion")]
  tt1$method="MLreg"
  posvalues <- tt1$pos
  tt1$pos <- factor(tt1$pos)
  setnames(tt1, "MLreg", "value")
  
  #Make the primary data object with normalized tpm values for specific cell binned by chr. Also unlog the tpm values that were logged before this script. 
  tt2 <- mlt1[, c("bin", "pos", "tpm", "proportion")]
  tt2$method="TPM Ratio"
  tt2$pos <- factor(tt2$pos)
  setnames(tt2, "tpm", "value")
  tt2$value <- 2^(tt2$value)
  
  #-----------------------------------------------------------------
  #allelic analysis: (fraction of nonzero SNPs and total SNP counts)
  #-----------------------------------------------------------------
  
  #coding SNPs - allele A:
  aggs.A.zs <- coding$aggs.A
  
  #normalizing by sample
  aggs.A.zs <- cbind(aggs.A.zs[,1:2], sweep(aggs.A.zs[, -c(1:2)], 2, colMeans(as.matrix(aggs.A.zs[, -c(1:2)]))+eps, FUN = "/"))
  
  #normalizing by bin
  aggs.A.zs <- cbind(aggs.A.zs[,1:2], sweep(aggs.A.zs[, -c(1:2)], 1, rowMedians(as.matrix(aggs.A.zs[, controlSampleIDs, with=F]))+eps, FUN = "/"))
  
  #Get avg of 50 gene bin total counts for all control samples and record rows with < 3 average 50 gene bin total cnts
  totalCntsA <- cbind(coding$cnts.A[, c(1:2)], meanSites = rowMeans(coding$cnts.A[, controlSampleIDs, with=F]))
  wi <- totalCntsA[seqnames == chr & meanSites < 3]$bin
  
  #create object fracNonZero (fraction of nonzero coding SNPs out of mean number of expressed coding SNPs)
  #create object cntsNonZero (total number of nonzero coding SNPs)
  aggsA <- merge(unique(tt1[method == "MLreg", c("bin", "pos", "proportion")]), aggs.A.zs[seqnames == chr][, c("bin", myid), with=F])
  aggsAc <- (merge(unique(tt1[method == "MLreg", c("bin", "pos", "proportion")]), coding$cnts.A[seqnames == chr][, c("bin", myid), with=F]))
  names(aggsA)[4] <- "fracNonZero"
  names(aggsAc)[4] <- "cntsNonZero"
  #make sites with less than 3 total nonzero SNPs into NA
  if(length(wi)>0) {
    aggsA$fracNonZero[aggsA$bin == wi] = NA
    aggsAc$cntsNonZero[aggsA$bin == wi] = NA
  }
  #combine into one dataframe
  aggsA <- merge(aggsA, aggsAc, by = c("bin", "pos"))
  #label this object as allele A
  aggsA$method = "Frac Allele A"
  
  #get same value but only for control samples  for allele A
  tmp <- aggs.A.zs[seqnames == chr, controlSampleIDs, with=F]
  tmp <- melt(t(tmp))[,2:3]
  colnames(tmp) <- c("bin", "value")
  if(length(wi)>0)
    tmp$value[tmp$bin == wi] = NA
  tmp2A <- merge(unique(tt1[method == "MLreg", c("bin", "pos", "proportion")]), tmp)
  tmp2A$method = "Ctrl. Frac A"
  
  #Do same process for Allele B
  
  #coding SNPs - allele B:
  aggs.B.zs <- coding$aggs.B
  #normalizing by sample
  aggs.B.zs <- cbind(aggs.B.zs[,1:2], sweep(aggs.B.zs[, -c(1:2)], 2, colMeans(as.matrix(aggs.B.zs[, -c(1:2)]))+eps, FUN = "/"))
  #normalizing by bin
  aggs.B.zs <- cbind(aggs.B.zs[,1:2], sweep(aggs.B.zs[, -c(1:2)], 1, rowMedians(as.matrix(aggs.B.zs[, controlSampleIDs, with=F]))+eps, FUN = "/"))
  
  #Get avg of 50 gene bin total counts for all control samples and record rows with < 3 average 50 gene bin total tpm
  totalCntsB <- cbind(coding$cnts.B[, c(1:2)], meanSites = rowMeans(coding$cnts.B[, controlSampleIDs, with=F]))
  wi <- totalCntsB[seqnames == chr & meanSites < 3]$bin
  
  #create object fracNonZero (fraction of nonzero coding SNPs out of mean number of expressed coding SNPs)
  #create object cntsNonZero (total number of nonzero coding SNPs)
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
  
  #get same values but only for control samples for allele B
  tmp <- aggs.B.zs[seqnames == chr, controlSampleIDs, with=F]
  tmp <- melt(t(tmp))[,2:3]
  colnames(tmp) <- c("bin", "value")
  if(length(wi)>0)
    tmp$value[tmp$bin == wi] = NA
  tmp2B <- merge(unique(tt1[method == "MLreg", c("bin", "pos", "proportion")]), tmp)
  tmp2B$method = "Ctrl. Frac B"
  
  
  #Combine control and regular values into two objects
  aggsAB <- rbind(aggsA, aggsB)
  tmp2AB <- rbind(tmp2A, tmp2B)
  
  aggsAB$pos <- factor(aggsAB$pos)
  tmp2AB$pos <- factor(tmp2AB$pos)
  
  aggsAB$method <- factor(aggsAB$method)
  tmp2AB$method <- factor(tmp2AB$method)
  
  #--------------------------------
  #fraction of non zeros TPM values
  #--------------------------------
  
  #normalizing nonzeros.zs (imported from previous script) by bin 
  #mynonzeros <- nonzeros.zs #nonzero.zs is the z-score of the fraction of nonzero tpm values.
  #mynonzeros0 <- cbind(mynonzeros[,1:2], sweep(mynonzeros[, -c(1:2)], 1, rowMedians(as.matrix(mynonzeros[, controlSampleIDs, with=F]))+eps, FUN = "-"))
  #
  ##merge annotations and relabel
  #mynonzeros <- merge(unique(tt1[method == "MLreg", c("bin", "pos", "proportion")]), mynonzeros0[seqnames == chr][, c("bin", myid), with=F])
  #names(mynonzeros)[4] <- "fracNonZero"
  #mynonzeros$method = "Sample Non Zeros"
  #mynonzeros$pos <- factor(mynonzeros$pos)
  #
  ##Get control version of mynonzeros object
  #tmp <- mynonzeros0[seqnames == chr, controlSampleIDs, with=F]
  #tmp <- melt(t(tmp))[,2:3]
  #colnames(tmp) <- c("bin", "value")
  #tmp2 <- merge(unique(tt1[method == "MLreg", c("bin", "pos", "proportion")]), tmp)
  
  #make chrX values to NA
  if(chr == "chrX") {
    tmp2AB[method == "Ctrl. Frac A", value := NA]
    aggsAB[method == "Frac Allele A", fracNonZero := NA]
    aggsAB[method == "Frac Allele A", cntsNonZero := NA]
  }
  
  
  #---------
  # Visuals:
  #---------
  require(gridExtra)
  require(dplyr)
  
  #get centromere ranges for each chr
  rng <- range(centromeres[seqnames==chr, c("start", "end")])
  
  #Function defining probability percentiles for boxplots later in visuals
  myf <- function(x) {
    r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  
  #Use dplyr to bin values and get summary statistics for each bin
  tt2b = tt2 %>% 
    group_by(bin) %>% 
    summarise(pos = unique(pos), #.groups = "drop_last", #Use this to kill warning message
              method = unique(method),
              proportion = mean(proportion),
              median = median(value, na.rm=T), 
              se.min = quantile(value, probs = 0.25, na.rm=T), 
              se.max = quantile(value, probs = 0.75, na.rm=T))
  
  #If not plotting return certain useful data objects
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
  
  #CZ version of p2 with discrete x-axis
  #p2 <- ggplot(tt2b, aes(x = pos, y = median)) + 
  #  geom_vline(xintercept=min(which(as.numeric(as.character(unique(tmp2AB$pos))) >= rng[1])), linetype="solid", color = alpha("black", alpha = 0.4), size = 5) +
  #  geom_hline(yintercept=c(0, 0.5, 1.0, 1.5, 2), linetype="dashed", color = "black", size = c(0.7, 0.7, 1.5, 0.7, 0.7)) +
  #  geom_point(color = "black", aes(x = pos, y = median), size = 3) + 
  #  ylim(c(0, 2.5)) + xlim(c(0, max(as.numeric(tt2b$pos)))) +
  #  scale_y_continuous(breaks = c(0,.5,1,1.5,2,2.5), labels = c("0","","1","","2",""), limits = c(0,2.5)) +
  #  scale_x_discrete(breaks=posvalues, labels = as.character(round(posvalues/1e6))) + 
  #  coord_cartesian(clip="off") +
  #  geom_rangeframe( y=seq(0, 2.5, along.with = tt2b$median)) + 
  #  theme_tufte() +
  #  theme(aspect.ratio = 1/5, text=element_text(size=18), axis.ticks.length=unit(.25, "cm"),
  #        axis.text = element_text(color="black"), axis.title = element_blank()) 
  
  #change position back to numeric 
  tt2b$pos <- as.numeric(as.character(tt2b$pos))
  #create continuous ticks
  all_ticks <- seq(0,250000000,5000000)
  major_ticks <- seq(0,250000000,25000000)
  all_ticks_breaks <- all_ticks[all_ticks < max(tt2b$pos)]; #all_ticks_breaks <- all_ticks_breaks[all_ticks_breaks > min(tt2b$pos)]
  all_ticks_labs <- as.character(round(all_ticks_breaks/1e6))
  all_ticks_labs[which(!all_ticks_breaks %in% major_ticks)] <- ""
  #CZ version of p2 with continuous x-axis
  p2 <- ggplot(tt2b, aes(x = pos, y = median)) + 
    geom_vline(xintercept=min(tt2b$pos[which(tt2b$pos >= rng[1])]), linetype="solid", color = alpha("black", alpha = 0.4), size = 5) +
    geom_hline(yintercept=c(0, 0.5, 1.0, 1.5, 2), linetype="dashed", color = "black", size = c(0.7, 0.7, 1.5, 0.7, 0.7)) +
    geom_point(color = "black", aes(x = pos, y = median), size = 3) + 
    scale_y_continuous(breaks = c(0,.5,1,1.5,2,2.5), labels = c("0","","1","","2",""), limits = c(0,2.5)) +
    scale_x_continuous(breaks=all_ticks_breaks, labels = all_ticks_labs, limits = c(0, 250000000)) + 
    coord_cartesian(clip="off") +
    geom_rangeframe(x=seq(0, max(all_ticks_breaks), along.with = tt2b$pos), y=seq(0, 2.5, along.with = tt2b$median)) + 
    theme_tufte() +
    theme(aspect.ratio = 1/5, text=element_text(size=18), axis.ticks.length=unit(.25, "cm"),
          axis.text = element_text(color="black"), axis.title = element_blank()) 

  #CZ version of p4 with discrete x-axis
  #p4 <- ggplot(data = tmp2AB, mapping = aes(x = pos, y= value, fill=method)) + 
  #  geom_vline(xintercept=min(which(as.numeric(as.character(unique(tmp2AB$pos))) >= rng[1])), linetype="solid", color = alpha("black", alpha = 0.4), size = 5) +
  #  geom_hline(yintercept=c(0, 1.0, 2), linetype="dashed", color = "black", size = c(0.7, 1.5, 0.7)) +
  #  geom_point(position = position_dodge(0.7), show.legend = F, stat='identity', size=3, data = aggsAB, 
  #             mapping = aes(x=pos, y=fracNonZero, color = method)) +
  #  scale_colour_manual(values = c("dodgerblue3", "firebrick2")) +
  #  scale_y_continuous(breaks = c(0,.5,1,1.5,2,2.5), labels = c("0","","1","","2",""), limits = c(0,2.5)) +
  #  scale_x_discrete(breaks=tt1$pos, labels = as.character(round(as.numeric(as.character(tt1$pos))/1e6))) +
  #  coord_cartesian(clip="off") +
  #  geom_rangeframe( y=seq(0, 2.5, along.with = tmp2AB$value)) + 
  #  theme_tufte() +
  #  theme(aspect.ratio = 1/5, text=element_text(size=18), axis.ticks.length=unit(.25, "cm"), 
  #        axis.text = element_text(color="black"), axis.title = element_blank()) 
  
  #change position back to numeric 
  aggsAB$pos <- as.numeric(as.character(aggsAB$pos))
  #create continuous ticks
  all_ticks <- seq(0,250000000,5000000)
  major_ticks <- seq(0,250000000,25000000)
  all_ticks_breaks <- all_ticks[all_ticks < max(aggsAB$pos)]; #all_ticks_breaks <- all_ticks_breaks[all_ticks_breaks > min(aggsAB$pos)]
  all_ticks_labs <- as.character(round(all_ticks_breaks/1e6))
  all_ticks_labs[which(!all_ticks_breaks %in% major_ticks)] <- ""
  #CZ version of p4 with continuous x-axis
  p4 <- ggplot(data = aggsAB, mapping = aes(x = pos, y = value, fill=method)) + 
    geom_hline(yintercept=c(0, 1.0, 2), linetype="dashed", color = "black", size = c(0.7, 1.5, 0.7)) +
    geom_vline(xintercept=min(aggsAB$pos[which(aggsAB$pos >= rng[1])]), linetype="solid", color = alpha("black", alpha = 0.4), size = 5) +
    geom_point(position = position_dodge(0.7), show.legend = F, stat='identity', size=3, data = aggsAB, 
               mapping = aes(x=pos, y=fracNonZero, color = method)) + 
    scale_colour_manual(values = c("dodgerblue3", "firebrick2")) +
    scale_y_continuous(breaks = c(0,.5,1,1.5,2,2.5), labels = c("0","","1","","2",""), limits = c(0,2.5)) +
    scale_x_continuous(breaks=all_ticks_breaks, labels = all_ticks_labs, limits = c(0, 250000000)) +
    coord_cartesian(clip="off") +
    geom_rangeframe(x=seq(0, max(all_ticks_breaks), along.with = aggsAB$pos), y=seq(0, 2.5, along.with = aggsAB$fracNonZero)) + 
    theme_tufte() +
    theme(aspect.ratio = 1/5, text=element_text(size=18), axis.ticks.length=unit(.25, "cm"), 
          axis.text = element_text(color="black"), axis.title = element_blank()) 
  

  ggsave(plot = p2, filename = sprintf("%s.%s.%s.TotalTPM.pdf", myfamily_file, myid_file, chr), path = destDir, device = "pdf", width = 12, height = 3, dpi = 300, units = "in")
  ggsave(plot = p4, filename = sprintf("%s.%s.%s.AllelicCov.pdf", myfamily_file, myid_file, chr), path = destDir, device = "pdf", width = 12, height = 3, dpi = 300, units = "in")
  #grid.arrange(p2, p4, nrow = 2)
  if(save.rds == T) {
    saveRDS(p2, file = sprintf("%s/%s.%s.%s.TotalTPM.rds", destDir, myfamily_file, myid_file, chr))
    saveRDS(p4, file = sprintf("%s/%s.%s.%s.AllelicCov.rds", destDir, myfamily_file, myid_file, chr))
    return( data.table(TPM_file = sprintf("%s/%s.%s.%s.TotalTPM.rds", destDir, myfamily_file, myid_file, chr), var_file = sprintf("%s/%s.%s.%s.AllelicCov.rds", destDir, myfamily_file, myid_file, chr)) )
  }
  #not needed with current macro plotting function.
  #dev.off()
  #return(p)
  
}


