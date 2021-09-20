plot_barplots_of_AllelicAndExpBiasPerSamples <- function(ids = c("170209_A5", "170209_A6"), chr = "chr12", 
                                               dfs, wt, fracs, ymax = NULL,
                                               plotOnlyDepth = F, pval.th = 0.01,
                                               nogeno = c("chr14", "chr21", "chr22", "chrY", "chrM"), plotAxisText = T) {
  
  require(gtools)
  require(ggplot2)
  require(data.table)
  require(RColorBrewer)
  
  #preparing data for plot:
  #Total hight is expBias and frac is allele A proportion in that hight
  #Four options: A, B, noGeno, mid
  ngroups <- 4
  cols <- c("dodgerblue3", "firebrick2", "green", "blueviolet") #colorRampPalette(brewer.pal(4, "Set1"))
  nmid <- 3
  
  if(sum(ids %in% colnames(fracs$Af)) < length(ids) | sum(ids %in% colnames(fracs$Bf)) < length(ids) | sum(ids %in% colnames(dfs)) < length(ids)) {
    message("Cannot plot ", ids)
    return(NA)
  }
  getMySampleDF <- function(id, chr) {
    if(chr != "" && chr != "chr") {
      chr <- sprintf("%s_", chr)
    }
    #print(chr)
    #print(id)
    A <- fracs$Af[grep(chr, fracs$Af$bin_id), c("bin_id", id), with=F]
    A$allele <- "A"
    A$sample <- id
    A <- data.table(A)
    setnames(A, old = id, new = "fraction")
    B <- fracs$Bf[grep(chr, fracs$Bf$bin_id), c("bin_id", id), with=F]
    B$allele <- "B"
    B$sample <- id
    B <- data.table(B)
    setnames(B, old = id, new = "fraction")
    
    tmp <- data.frame(bin_id = grep(chr, dfs$bin_id, value = T), amp = dfs[grep(chr, dfs$bin_id), id])
    # tmp <- merge(data.frame(amp = dfs[grep(chr, rownames(dfs)), id, drop=F]), 
    #              data.frame(wt = wt[grep(chr, names(wt))]), by="row.names")
    
    #names(tmp)[2] <- "amp"
    tmp$amp <- tmp$amp
    tmp$amp[grep("chrX", tmp$bin_id)] <- tmp$amp[grep("chrX", tmp$bin_id)]*0.5
    tmp$bias <- tmp$amp
    #tmp <- tmp[, c(-2, -3)]
    #names(tmp)[1] <- "bin_id"
    
    #allele A:
    tmp2 <- merge(tmp, A, by = "bin_id", all.x = T)
    
    wi <- which(is.na(tmp2$fraction))
    if(length(wi) > 0) {
      tmp2$fraction[wi] <- 1.0
      tmp2$allele[wi] <- "NoGeno"
      tmp2$sample[wi] <- id
    }
    #tmp2$amp <- tmp2$bias
    A <- tmp2
    
    #allele B:
    tmp2 <- merge(tmp, B, by = "bin_id", all.x = T)
    wi <- which(is.na(tmp2$fraction))
    if(length(wi) > 0) {
      tmp2$fraction[wi] <- 1.0
      tmp2$allele[wi] <- "NoGeno"
      tmp2$sample[wi] <- id
    }
    tmp2$amp <- tmp2$bias*tmp2$fraction
    B <- tmp2
    
    #allele mid:
    mid <- tmp2
    mid$allele <- "mid"
    mid$amp <- mid$bias*0.5
    mid$fraction <- 0.5
    
    #merging all 4 alternatives:
    ccp <- unique(rbindlist(list(A, B, mid)))
    ccp$bin_id <- factor(ccp$bin_id, levels = unique(mixedsort(ccp$bin_id)))
    
    return(ccp)
  }
  

  ids <- intersect(ids, colnames(dfs))
  if(length(ids) == 0)
    return(NA)
  ccps <- rbindlist(lapply(ids, getMySampleDF, chr = chr))
  if(is.null(ymax))
    ymax <- max(3.0, max(ccps$amp), na.rm = T)
  #return(ccps)
  #Plotting:
  ccps$bin_id <- factor(as.character(ccps$bin_id), levels = mixedsort(unique(as.character(ccps$bin_id))))
  #return(ccps)
  #ccps <- ccps[mixedorder(as.character(bin_id))]
  p <- ggplot(ccps[allele != "mid"], aes(x = bin_id, y = round(amp, digits = 3), 
                                         fill = allele, group = sample)) + 
    geom_bar(stat = "identity", position = "dodge", width = 0.8, size = 0.7, colour = "black")
  
  if(plotAxisText) {
    p <- p + ggtitle(sprintf("IDs: %s", paste(ids, collapse = " | "))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = +0.5, size = 15), axis.text.y = element_text(size=12)) +
      scale_fill_manual(values = cols) + #cols(ngroups)) +
      geom_errorbar(data=ccps[allele == "mid"], 
                    mapping=aes(x=bin_id, ymin=amp, ymax=amp, group = sample), 
                    width=1, size=1, position = "dodge", color= cols[3]) + #cols(ngroups)[nmid]) + 
      scale_y_continuous(name ="Expected CN", breaks = seq(0,ymax,0.25), limits = c(-0.25, ymax+0.25)) + 
      geom_text(data = ccps[ bin_id == levels(ccps$bin_id)[1] & allele == "mid"], 
                mapping = aes(label = sample, x = bin_id, y = 0, group = sample), 
                position = position_dodge(width = 0.8), hjust=1.2, color = "black", size=4, angle = 90) 
  } else {
    p <- p + ggtitle(sprintf("IDs: %s", paste(ids, collapse = " | "))) +
      theme(axis.text.x = element_blank(), axis.text.y = element_text(size=12)) +
      scale_fill_manual(values = cols(ngroups)) +
      geom_errorbar(data=ccps[allele == "mid"], 
                    mapping=aes(x=bin_id, ymin=amp, ymax=amp, group = sample), 
                    width=1, size=1, position = "dodge", color= cols[3]) + #cols(ngroups)[nmid]) + 
      scale_y_continuous(name ="Expected CN", breaks = seq(0,3.0,0.25), limits = c(-0.25, ymax+0.25)) + 
      geom_text(data = ccps[ bin_id == levels(ccps$bin_id)[1] & allele == "mid"], 
                mapping = aes(label = sample, x = bin_id, y = 0, group = sample), 
                position = position_dodge(width = 0.8), hjust=1.2, color = "black", size=4, angle = 90) 
    
  }  
  plot(p)
  return(p)
}

