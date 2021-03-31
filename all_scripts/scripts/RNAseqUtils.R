#RNAseqUtils.R

#TODO: put in rna workflow

#general RNAseq utility functions:
#1. Read a BAM file to vector of counts per gene.
#2. Read a NAM dir to a matrix of counts per gene.


#TODOS:
#0. Check why sum
#1. use findOverlaps to aggragate data per chrome and bins.
#2. plot mean bias as a function of genomic bin

#unstranded

# head 161021_M01209_0362_000000000-AWA9L_170406-A5_S5_L001_R1_001.fastq.gz.rsem.isoforms.results 
# transcript_id   gene_id length  effective_length        expected_count  TPM     FPKM    IsoPct
# ENST00000373020.8_TSPAN6-001    ENSG00000000003.14_TSPAN6       2206    1977.07 208.00  51.29   52.26   100.00
# ENST00000494424.1_TSPAN6-002    ENSG00000000003.14_TSPAN6       820     591.17  0.00    0.00    0.00    0.00
# ENST00000612152.4_TSPAN6-201    ENSG00000000003.14_TSPAN6       3796    3567.07 0.00    0.00    0.00    0.00
# ENST00000614008.4_TSPAN6-202    ENSG00000000003.14_TSPAN6       900     671.09  0.00    0.00    0.00    0.00
# ENST00000373031.4_TNMD-001      ENSG00000000005.5_TNMD  1339    1110.07 0.00    0.00    0.00    0.00
# ENST00000371582.8_DPM1-005      ENSG00000000419.12_DPM1 1161    932.07  0.00    0.00    0.00    0.00
# ENST00000371584.8_DPM1-003      ENSG00000000419.12_DPM1 1073    844.07  0.00    0.00    0.00    0.00


generate_rsems_object <- function(annofile = "/czlab/etai/ExperimentsData/Stamatis/Stamatis_list_v10_170830c_modified.csv") {
  anno <- read.csv(annofile,stringsAsFactors = F)
  anno2 <- data.frame(id = anno$WTA.plate, fname = sub("fastq", "rsem", anno$filePath_2nd))
  rsems <- collectAllRSEMresults(anno2)
  rsems$anno <- anno
  #save(rsems, file = "/czlab/etai/ExperimentsData/Stamatis/run.v2/rsems.Stamatis_list_v10_170830c.RData")
  return(rsems)
}
collectAllRSEMresultsFromDir <- function(dir = "/czlab/etai/ExperimentsData/GSE100771/rsem") {
  fnames <- list.files(path = dir, pattern = "*.genes.results", full.names = T)
  prefixs <- gsub(pattern = ".genes.results", "", fnames)
  getMyReslts <- function(filePrefix) {
    require(data.table)
    message("Getting data from file: ", filePrefix, "...")
    gfile <- sprintf("%s.genes.results", filePrefix)
    tfile <- sprintf("%s.isoforms.results", filePrefix)
    if(!file.exists(gfile) || !file.exists(tfile)) {
      message("ERROR - file ", gfile, " or file ", tfile, " do not exists.")
      return(NA)
    }
    gtbl <- fread(file = gfile, sep = "\t")
    ttbl <- fread(file = tfile, sep = "\t")
    return(list(isoforms = ttbl, genes = gtbl))
  }  
  
  require(stringr)
  lrsem <- lapply(prefixs, getMyReslts)
  names(lrsem) <- basename(prefixs)
  
  rsemListToMatrix <- function(columnName, lrsem, type) {
    tmp <- do.call(cbind, lapply(names(lrsem), function(x) lrsem[[x]][[type]][[columnName]]))
    colnames(tmp) <- names(lrsem)
    rownames(tmp) <- lrsem[[1]][[type]][[1]]
    return(tmp)
  }
  
  isoforms <- lapply(names(lrsem[[1]]$isoforms)[3:8], 
                     function(x) rsemListToMatrix(columnName = x, lrsem = lrsem, type = "isoforms"))
  names(isoforms) <- names(lrsem[[1]]$isoforms)[3:8]
  mappedIsoforms <- do.call(rbind, lapply(rownames(isoforms$TPM), function(x) str_split(string = x, pattern = "_", n = 2)[[1]]))
  genes <- lapply(names(lrsem[[1]]$genes)[3:7], 
                  function(x) rsemListToMatrix(columnName = x, lrsem = lrsem, type = "genes"))
  names(genes) <- names(lrsem[[1]]$genes)[3:7]
  mappedGenes <- do.call(rbind, lapply(rownames(genes$TPM), function(x) str_split(string = x, pattern = "_", n = 2)[[1]]))
  mappedGenes[grep("PAR_Y", rownames(genes$TPM), value = F),1] <- paste(mappedGenes[grep("PAR_Y", rownames(genes$TPM), value = F),1], "_PAR_Y", sep = "")
  mappedGenes <- gsub("PAR_Y_", replacement = "", mappedGenes[grep("PAR_Y", rownames(genes$TPM), value = F),2])
  
  return(list(isoforms = isoforms, genes = genes, anno = prefixs, 
              mappedIsoforms = mappedIsoforms, mappedGenes = mappedGenes))
  
}
#anno - data.frame with at least two fields: 
#1st column: sample id to use
#2nd column - prefix of rsem file for that sample id
#             File prefix should indicate for at least two files:
#             with the suffixes: ".rsem.genes.results" and ".rsem.isoforms.results"
#
#Example of a run:
#anno <- unname(expr$anno$yy[, c("Sample.name", "filePath_2nd")])
#anno2 <- unname(expr$anno$stam[, c("WTA.plate", "filePath_2nd")])
#names(anno) <- c("id", "filePrefix")
#names(anno2) <- c("id", "filePrefix")
#annoall <- rbind(anno, anno2)
#rsem <- collectAllRSEMresults(annoall)
collectAllRSEMresults <- function(anno) {
  require(stringr)
  lrsem <- lapply(anno[, 2], readRSEMresultsFromFile)
  names(lrsem) <- as.character(anno[, 1])
  
  rsemListToMatrix <- function(columnName, lrsem, type) {
    tmp <- do.call(cbind, lapply(names(lrsem), function(x) lrsem[[x]][[type]][[columnName]]))
    colnames(tmp) <- names(lrsem)
    rownames(tmp) <- lrsem[[1]][[type]][[1]]
    return(tmp)
  }
  isoforms <- lapply(names(lrsem[[1]]$isoforms)[3:8], 
                     function(x) rsemListToMatrix(columnName = x, lrsem = lrsem, type = "isoforms"))
  names(isoforms) <- names(lrsem[[1]]$isoforms)[3:8]
  mappedIsoforms <- do.call(rbind, lapply(rownames(isoforms$TPM), function(x) str_split(string = x, pattern = "_", n = 2)[[1]]))
  genes <- lapply(names(lrsem[[1]]$genes)[3:7], 
                     function(x) rsemListToMatrix(columnName = x, lrsem = lrsem, type = "genes"))
  names(genes) <- names(lrsem[[1]]$genes)[3:7]
  mappedGenes <- do.call(rbind, lapply(rownames(genes$TPM), function(x) str_split(string = x, pattern = "_", n = 2)[[1]]))
  mappedGenes[grep("PAR_Y", rownames(genes$TPM), value = F),1] <- paste(mappedGenes[grep("PAR_Y", rownames(genes$TPM), value = F),1], "_PAR_Y", sep = "")
  mappedGenes[grep("PAR_Y", rownames(genes$TPM), value = F),2] <- gsub("PAR_Y_", replacement = "", mappedGenes[grep("PAR_Y", rownames(genes$TPM), value = F),2])

  return(list(isoforms = isoforms, genes = genes, anno = anno, 
              mappedIsoforms = mappedIsoforms, mappedGenes = mappedGenes))
}



readRSEMresultsFromFile <- function(filePrefix = "/czlab/etai/ExperimentsData/allfastqs/161021_M01209_0362_000000000-AWA9L_170406-A5_S5_L001_R1_001.fastq.gz") {
  require(data.table)
  message("Getting data from file: ", filePrefix, "...")
  gfile <- sprintf("%s.rsem.genes.results", filePrefix)
  tfile <- sprintf("%s.rsem.isoforms.results", filePrefix)
  if(!file.exists(gfile) || !file.exists(tfile)) {
    message("ERROR - file ", gfile, " or file ", tfile, " do not exists.")
    return(NA)
  }
  gtbl <- fread(file = gfile, sep = "\t")
  ttbl <- fread(file = tfile, sep = "\t")
  return(list(isoforms = ttbl, genes = gtbl))
}



#anno - data.frame with at least two fields: 
#1st column: sample id to use
#2nd column - prefix of rsem file for that sample id
#             File prefix should indicate for at least two files:
#             with the suffixes: ".rsem.genes.results" and ".rsem.isoforms.results"
# Example:
# tmpAnno <- anno[, c("WTA.plate", "filePath_2nd")]
#
# STAR .ReadsPerGene.out.tab output:
# column 1: gene ID
# column 2: counts for unstranded RNA-seq
# column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
# column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)
collectAllSTARCountsFromFileList <- function(anno, 
                                             starOutputDir = "/czlab/etai/ExperimentsData/Stamatis/run.v2/star/",
                                             strandness = c("unstranded", "forward", "reverse")) {
  suffix=".ReadsPerGene.out.tab"
  fileidxs <- which(!is.na(anno[,2]))
  files <- sprintf("%s/%s%s", starOutputDir, basename(anno[fileidxs, 2]), suffix) 
  
  stranIdx <- match(strandness[1], c("unstranded", "forward", "reverse"))
  getMyCntTables <- function(f) {
    cat("Reading tbl:", f, "\n")
    tbl <- read.table(f, header = F, sep = "\t", 
                      stringsAsFactors = F, row.names = 1)[, stranIdx, drop=F]
    return(tbl)
  }
  htmat <- do.call(cbind, lapply(files, getMyCntTables))
  names(htmat) <- as.character(anno[fileidxs, 1])
  anno <- anno[fileidxs, ]
  
  htqc <- htmat[grep("^N_", rownames(htmat)), ]
  htmat <- htmat[-grep("^N_", rownames(htmat)), ]
  
  htqc <- rbind(htqc, geneCounts=colSums(htmat))
  return(list(counts=htmat, qc=htqc, anno=anno))
}



#anno data frame must contain at least two fields:
#Sample.name
#filePath
collectAllSTARCountsFromFileList.old <- function(anno) {
  suffix=".ReadsPerGene.out.tab"
  fileidxs <- which(!is.na(anno$filePath))
  files <- sprintf("%s%s", anno$filePath[fileidxs], suffix) 
  
  getMyCntTables <- function(f) {
    cat("Reading tbl:", f, "\n")
    tbl <- read.table(f, header = F, sep = "\t", 
                      stringsAsFactors = F, row.names = 1)[,1,drop=F]
    return(tbl)
  }
  htmat <- do.call(cbind, lapply(files, getMyCntTables))
  names(htmat) <- basename(anno$filePath[fileidxs])
  anno <- anno[fileidxs, ]
  
  htqc <- htmat[grep("^N_", rownames(htmat)), ]
  htmat <- htmat[-grep("^N_", rownames(htmat)), ]
  
  htqc <- rbind(htqc, geneCounts=colSums(htmat))
  return(list(counts=htmat, qc=htqc, anno=anno))
}


collectAllSTARCountsFromFileDir <- function(workdir="/czlab/etai/ExperimentsData/Stamatis/results_latest") {
  suffix=".ReadsPerGene.out.tab"
  
  files <- list.files(path = workdir, pattern = sprintf("*%s", suffix), full.names = T)
  
  count <- 0
  htmat <- do.call(cbind, lapply(files,
                                 function(f) {
                                   count <- count + 1
                                   message(count, ": ", f)
                                   read.table(f, header = F, sep = "\t", stringsAsFactors = F, row.names = 1)[,1,drop=F] }
                                 ))
  names(htmat) <- do.call(rbind, lapply(files, function(f) strsplit(basename(f), split = ".gz")[[1]][1]))
  
  htqc <- htmat[grep("^N_", rownames(htmat)), ]
  htmat <- htmat[-grep("^N_", rownames(htmat)), ]
  
  htqc <- rbind(htqc, geneCounts=colSums(htmat))
  return(list(counts=htmat, qc=htqc))
}

#Taken from:
#https://gist.github.com/slowkow/c6ab0348747f86e2748b
#' Convert counts to transcripts per million (TPM).
#' 
#' Convert a numeric matrix of features (rows) and conditions (columns) with
#' raw feature counts to transcripts per million.
#' 
#'    Lior Pachter. Models for transcript quantification from RNA-Seq.
#'    arXiv:1104.3889v2 
#'    
#'    Wagner, et al. Measurement of mRNA abundance using RNA-seq data:
#'    RPKM measure is inconsistent among samples. Theory Biosci. 24 July 2012.
#'    doi:10.1007/s12064-012-0162-3
#'    
#' @param counts A numeric matrix of raw feature counts i.e.
#'  fragments assigned to each gene.
#' @param featureLength A numeric vector with feature lengths.
#' @param meanFragmentLength A numeric vector with mean fragment lengths.
#' @return tpm A numeric matrix normalized by library size and feature length.
counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  featureLength <- ifelse(is.na(featureLength), 0, featureLength)
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}

#http://www.sthda.com/english/wiki/static-and-interactive-heatmap-in-r-unsupervised-machine-learning
#TODO: Add side colors sample annotation when annotation data will be available
getJaccardHeatmap <- function(M, minExprTh = 1, minSamples = 1) {
  require(reshape2)
  M <- filter.by.occurence(M, minExprTh, minSamples)
  cat("Dim in jaccard instance: ", dim(M), "\n")
  if(dim(M)[1] == 0) {
    cat("No genes passed this threshold (", minExprTh, ", ", minSamples, "\n")
    return(NA)
  }
  ll <- do.call(rbind, combn(x = 1:(dim(M)[2]), m = 2,
                             FUN = function(x) data.frame(i=names(M)[x[1]], j = names(M)[x[2]],
                                                          J=jaccard(M[,x[1]], M[,x[2]]),
                                                          stringsAsFactors = F), simplify = F))
  ll <- rbind(data.frame(i=names(M), j=names(M), J = 1, stringsAsFactors = F), ll)
  ji <- acast(ll, i~j, value.var="J")
  ji[lower.tri(ji, diag = F)] <- t(ji)[lower.tri(ji, diag = F)]


  if(length(table(ji) == 1))
    ji[1,1] <- 1.01*ji[1,1]
  require("gplots")
  par(mar = c(7, 4, 2, 2) + 0.2) 
  heatmap.2(ji, scale = "none", col = bluered(32), 
            margins = c(12,12), cexRow = 0.85, cexCol = 0.85,
            main = sprintf("Counts>%d Samples>%d\nN = %d", minExprTh, minSamples, dim(M)[1]), 
            trace = "none", density.info = "histogram",
            RowSideColors = rep(c("pink"), each = dim(ji)[1]),
            ColSideColors  = rep(c("pink"), each = dim(ji)[1]))

  return(ji)
}


downSample <- function(counts, dsread = 50000, dsfrac = NULL) {
  require(DESeq)

  if(is.null(dsfrac)) {
    counts_ds <- apply(counts, 2, function(k) rbinom(length(k), k, min(dsread/sum(k), 1.0)))
  } else {
    counts_ds <- apply(counts, 2, function(k) rbinom(length(k), k, dsfrac))
  }
  sfAll_ds <- estimateSizeFactorsForMatrix( counts_ds )
  nCounts_ds <- t( t(counts_ds) / sfAll_ds )
  
  return(nCounts_ds)

}

generateSaturationCurve <- function(counts, readBy = 50000, startat = 1000) {
  require(RColorBrewer)

  ll <- lapply(seq(startat, max(colSums(counts)), readBy), function(x) downSample(counts, dsread = x))
  ll2 <- lapply(ll, function(x) colSums(x>0))
  ll3 <- do.call(rbind, ll2)
  rownames(ll3) <- seq(startat, max(colSums(counts)), readBy)
  
  cols <- c("black", brewer.pal(n = 12, name = 'Paired'), rainbow(dim(ll3)[2]))[1:(dim(ll3)[2])]
  plot(range(as.numeric(rownames(ll3))), range(ll3), type="n", 
       xlab= "Number of reads samples", ylab = " Number of genes detected")
  
  for(i in 1:(dim(ll3)[2])) 
    lines(rownames(ll3), ll3[,i], type="b", col=cols[i], lwd=2, pch=i)
  
  legend("bottomright", inset=.05, 
         names(counts), col=cols, pch = 1:(dim(ll3)[1]), horiz=F)
  
  return(ll3)

}

filter.by.occurence <- function(
  ### Filter genes (rows) by require those kept to have atleast a min values in a min number of samples.
  ### This avoids using averages and it stable when adding samples to the study.
  df_data,
  ### Data frame to be filtered
  d_min_value,
  ### Minimum value a measurement of gene must be to be counted in a sample.
  d_min_occurence
  ### Minimum times a gene must be found in a sample (at the d_min_value) in order to not be filtered.
){
  vctr_f_keep = apply( df_data, 1, function( x ){ return( length( which( x >= d_min_value ) ) >= d_min_occurence ) } )
  return( df_data[ vctr_f_keep,] )
}

jaccard <- function(m1, m2) {
  
  J <- 1 - dist(t(cbind(unlist(m1), unlist(m2))), "binary")

  return( as.numeric( J ))
}

getOverlapingGenesBetweenSamples <- function(cnts) {
  res <- list()
  
  cnts <- as.matrix(expr$counts)
  #Take only genes that were detected in at least one experiment:
  cnts <- cnts[rowSums(cnts>0) > 0, ]
  bcnts <- (cnts>0)+0
  res[["Number of genes that were detected in at least one experiment"]] <- dim(cnts)[1]
  
  #genes in quantiles:
  cnts.means <- rowMeans(cnts)
  mq <- quantile(cnts.means, probs = c(0, seq(0.25, to = 0.75, by = 0.25), 0.95, 0.99, 1))
  
  
  getMyOLtable <- function(pf, pt) {
    bcnts.tmp <- bcnts[cnts.means > mq[pf] & cnts.means <= mq[pt], ]
    ol.detected <- combn(x = 1:(dim(bcnts.tmp)[2]), m = 2, FUN = function(x) sum((bcnts.tmp[, x[1]] + bcnts.tmp[,x[2]]) == 2))
    ol.zeros <- combn(x = 1:(dim(bcnts.tmp)[2]), m = 2, FUN = function(x) sum((bcnts.tmp[, x[1]] + bcnts.tmp[,x[2]]) == 0))
    nol <- combn(x = 1:(dim(bcnts.tmp)[2]), m = 2, FUN = function(x) sum(bcnts.tmp[, x[1]] != bcnts.tmp[,x[2]]))
    totalDetedcted <- combn(x = 1:(dim(bcnts.tmp)[2]), m = 2, FUN = function(x) sum(bcnts.tmp[, x[1]] | bcnts.tmp[,x[2]]))
    totalInGroup <- dim(bcnts.tmp)[1] #combn(x = 1:(dim(bcnts.tmp)[2]), m = 2, FUN = function(x) sum(bcnts.tmp[, x[1]] | bcnts.tmp[,x[2]]))
    df <- data.frame(t(combn(x = 1:(dim(bcnts.tmp)[2]), m = 2)), 
                     ol.detected = ol.detected, ol.zeros = ol.zeros, nol = nol, 
                     totalDetedcted = totalDetedcted, totalInGroup = totalInGroup)
    df[, 1] <- colnames(bcnts.tmp)[df[,1]]
    df[, 2] <- colnames(bcnts.tmp)[df[,2]]
    return(df)
  }
  
  #all expressed genes:
  
  ol.detected <- combn(x = 1:(dim(bcnts)[2]), m = 2, FUN = function(x) sum((bcnts[, x[1]] + bcnts[,x[2]]) == 2))
  ol.zeros <- combn(x = 1:(dim(bcnts)[2]), m = 2, FUN = function(x) sum((bcnts[, x[1]] + bcnts[,x[2]]) == 0))
  nol <- combn(x = 1:(dim(bcnts)[2]), m = 2, FUN = function(x) sum(bcnts[, x[1]] != bcnts[,x[2]]))
  totalDetedcted <- combn(x = 1:(dim(bcnts)[2]), m = 2, FUN = function(x) sum(bcnts[, x[1]] | bcnts[,x[2]]))
  totalInGroup <- dim(bcnts)[1]
  df <- data.frame(t(combn(x = 1:(dim(bcnts)[2]), m = 2)), 
                   ol.detected = ol.detected, ol.zeros = ol.zeros, nol = nol, 
                   totalDetedcted = totalDetedcted, totalInGroup = totalInGroup)
  df[, 1] <- colnames(bcnts)[df[,1]]
  df[, 2] <- colnames(bcnts)[df[,2]]
  
  res[["Overlapping genes with at least one read:"]] <- df
  
  for(pf in 1:(length(mq)-1)) {
    df <- getMyOLtable(pf, pf + 1)
    strfield <- sprintf("OVerlapping genes in percentile range %s-%s (%.2f-%.2f):", names(mq)[pf], names(mq)[pf+1], mq[pf], mq[pf+1])
    res[[strfield]] <- df
  }
  
  return(res)
  
}

plotGeneCoverage <- function(counts, mainstr="") {
  counts <- t(as.data.frame(sort(colSums(counts>0))))
  par(mar = c(7, 4, 2, 2) + 0.2) #add room for the rotated labels
  end_point = 0.5 + ncol(counts) + ncol(counts)-1
  #print(end_point)
  colnames <- colnames(counts)
  colnames(counts) <- NULL
  barplot(as.matrix(counts), beside=F, space=1, main=mainstr,
          col="lightgreen", 
          xlab="", ylab="Number of genes detected", cex.axis = 1.1)
  text(seq(1.5,end_point,by=2), par("usr")[3]-0.25, 
       srt = 60, adj= 1, xpd = TRUE,
       labels = paste(colnames), cex=0.65)  
}

plotQCBarSingleSample <- function(fname="/czlab/etai/ExperimentsData/10X/161209_M01209_0004_000000000-AY8CY/Sample_1_merged_PBMC.bam.fastq.ReadsPerGene.out.tab") {
  htmat <- read.table(fname, header = F, sep = "\t", stringsAsFactors = F, row.names = 1)
  htqc <- htmat[grep("^N_", rownames(htmat)), 1, drop=F]
  htmat <- htmat[-grep("^N_", rownames(htmat)), 1, drop=F]
 htqc <- qc <- rbind(htqc, geneCounts=colSums(htmat))

  qc <- sweep(as.matrix(qc),2,colSums(qc),`/`)*100
  barplot(as.matrix(qc), beside=T, space=1,
  col=rainbow(dim(qc)[1]),
  xlab="", ylab="Percent of all reads", cex.axis = 1.1, names.arg = rownames(qc), las=2)
  return(list(rawqc=htqc, qc=qc))  
}

plotQCBar <- function(qc) {
  par(mar = c(7, 4, 2, 2) + 0.2) #add room for the rotated labels
  end_point = 0.5 + ncol(qc) + ncol(qc)-1
  #print(end_point)
  qc <- sweep(as.matrix(qc),2,colSums(qc),`/`)*100
  colnames <- colnames(qc)
  colnames(qc) <- NULL
  barplot(as.matrix(qc), beside=F, space=1, 
          col=rainbow(dim(qc)[1]), 
          xlab="", ylab="Percent of all reads", cex.axis = 1.1)
  text(seq(1.5,end_point,by=2), par("usr")[3]-0.25, 
       srt = 60, adj= 1, xpd = TRUE,
       labels = paste(colnames), cex=0.65)
  
  legend("right", inset=.05, 
         row.names(qc), fill=rainbow(dim(qc)[1]), horiz=F)
}



#Assumes all HTseq-count files have the same order or genes and the same number of rows.
collectAllHTseqCountsFromFileList <- function(workdir="/czlab/etai/ExperimentsData/Stamatis/results_latest",
                                              suffix=".Aligned.sortedByCoord.out.bam.HTSeq-count2.txt") {
  files <- list.files(path = workdir, pattern = sprintf("*%s", suffix), full.names = T)


  htmat <- do.call(cbind, lapply(files,
                                 function(f) read.table(f, header = F, sep = "\t", stringsAsFactors = F, row.names = 1)))
  names(htmat) <- do.call(rbind, lapply(files, function(f) strsplit(basename(f), split = ".gz")[[1]][1]))

  htqc <- htmat[grep("^__", rownames(htmat)), ]
  htmat <- htmat[-grep("^__", rownames(htmat)), ]

  hist(colSums(htmat==0)/(dim(htmat)[1])*100, col="gray", xlab="% of genes with no read counts", main="", breaks=20)
  htqc <- rbind(htqc, geneCounts=colSums(htmat))
  return(list(counts=htmat, qc=htqc))
}



mergeMetaAndExpData <- function(anno) {
  #anno$basename <- sub(".gz$", "", basename(anno$filePath))
  getAllPairs <- function() {
    if(sum(anno$filePath == "NO INPUT") > 0)
      anno <- anno[ -which(anno$filePath == "NO INPUT"), ]
    tp <- table(anno$Pairs)[table(anno$Pairs)>=2]
    getpp <- function(x) {
      indexInhtrpkm = which(colnames(htrpkm) %in% anno$basename[which(anno$Pairs == names(tp[x]))])
      data.frame(event = anno$event[which(anno$Pairs == names(tp[x]))],
                 indexInhtrpkm = indexInhtrpkm,
                 fileNames = colnames(htrpkm)[indexInhtrpkm])

    }
    lapply(1:length(tp), getpp)

  }
  pr <- getAllPairs()
  pairsMat <- do.call(cbind,
                      lapply(1:length(pr),
                             function(x) combn(pr[[x]]$indexInhtrpkm, 2)))


  return(list(allPairs = pr, pairsMat = pairsMat))
}

#assumes bt file was created already using mergeMetaAndExpData as follows:
#ap <- mergeMetaAndExpData(anno)
#bt <- getGeneExpBiasPerBin(pairsMat = ap$pairsMat, txdb.genRange = txdb.genRange, chrLens = chrLens, binSize = 1e7)
plotAllPairsToPDF <- function(bt, pairsIdxs, annoNames, outputFileName="results/allPairs.binWidth1e7.pdf") {
  #ap <- mergeMetaAndExpData(anno)

  pdf(file = outputFileName, width = 15, height = 10)
  for(i in 1:(dim(pairsIdxs)[2])) {
    plotPairCrossGenomeBias(bt, 
                            pairsIdxs[1, i],
                            pairsIdxs[2, i],
                            main = annoNames[i])
                            #submain = sprintf("%s\n%s", ap$allPairs[[i]]$fileNames[1],
                             #                 ap$allPairs[[i]]$fileNames[2]))
  }
  dev.off()

}

plotPairCrossGenomeBias <- function(bt, i = 1, j = 2, main = "", submain = '') {
  bt2 <- bt[bt$i==i & bt$j==j,]
  cols <- c("blue", "red")[(bt2$pval.paired <= 0.01) + 1]

  plotMe <- function(y, ylab, main, submain) {
    barplot(y, main = main, sub = submain,
            width = rep(1, length(bt2$chr)), space = 0,
            names.arg = ifelse(duplicated(bt2$chr), "", as.character(bt2$chr)),
            col = cols,
            xlim = c(0, length(bt2$chr)+1),
            ylim = range(y[!is.nan(y)]),
            ylab = ylab, las=2)
    abline(v=which(!duplicated(bt2$chr)), col="grey", lty=6)
  }

  #y <- bt2$meanLog2Bias.sum
  #plotMe(y, ylab = "Mean log2 exp ratio", main = main, submain = submain)
  #y <- log2((bt2$Ng+1)/(bt2$nl+1))
  #plotMe(y, "Mean log2 count ratio", main = main, submain = submain)
  y <- sign(bt2$meanLog2Bias.paired)*(-log2(bt2$pval.paired))
  plotMe(y, ylab = "Signed(-log2(pval)) - paired genes", main = main, submain = submain)
  y <- (bt2$meanLog2Bias.sum)
  plotMe(y, ylab = "Bias of mean expression - grouped genes", main = main, submain = submain)


}


#TODO:
#Take care for cases where bin is empty (e.g. small bin) of genes: e.g. binSize=1e6
#see error in: bt.1e6 <- getGeneExpBiasPerBin(pairsMat = ap$pairsMat, txdb.genRange = txdb.genRange, chrLens = chrLens, binSize = 1e6)
#htrpkm is needed to be loaded before calling function
getGeneExpBiasPerBin <- function(pairsMat = combn(1:3,2),
                                 binSize = 1e24, minPairExpTh = 1,
                                 txdb.genRange, chrLens, geneLevelExp) {
  cat("preparing genomic bins by gene..\n")
  gib <- getGenomicBinsInGeneIDs(binSize = binSize,
                                 txdb.genRange = txdb.genRange, chrLens = chrLens)

  getBiasInChrBins <-function(chr, pairsMat, eps = 0.01) {
    cat("Doing chr ", chr, " with ", length(gib[[chr]]), " bins\n")
    getBiasForBin <- function(bini, chr, pairsMat) {
      mygib <- gib[[chr]][[bini]]$geneIDs
      mygib <- mygib[mygib %in% rownames(geneLevelExp)]
      
      tmp <- geneLevelExp[ mygib,]

      getBiasForPair <- function(i, j, plotMe = F) {
        #cat("i, j: ", i, " ", j, "\n")
        idxs <- (tmp[,i]+tmp[,j] > minPairExpTh)
        #cat("num of genes in bin = ", length(idxs), "\n")
        nGenesWithTh <- length(idxs)
        nGenesNoTh = dim(tmp)[1]
        meanLog2Bias.paired <- mean(log2((tmp[idxs,i]+eps)/(tmp[idxs,j]+eps)))
        meanLog2Bias.sum <- log2(mean(tmp[idxs,i]+eps)/mean(tmp[idxs,j]+eps))
        Ng <- sum(tmp[idxs,i] > tmp[idxs,j])
        Nl <- sum(tmp[idxs,i] < tmp[idxs,j])
        pval.paired <- 1
        pval.grp <- 1
        if(Ng + Nl > 5) {
          pval.paired <- wilcox.test(tmp[idxs,i], tmp[idxs,j], paired = T)$p.val
          pval.grp <- wilcox.test(tmp[idxs,i], tmp[idxs,j], paired = F)$p.val
        }
        Ng <- sum(tmp[idxs,i] > tmp[idxs,j])
        Nl <- sum(tmp[idxs,i] < tmp[idxs,j])
        if(plotMe)
          hist(log2((tmp[idxs,i]+eps)/(tmp[idxs,j]+eps)), main=sprintf("%d vs %d", i, j),
               col="gray", breaks = 22, xlab = "Log2 bias")
        return(data.frame(chr = chr,
                          binNum = bini,
                          i = i, j = j,
                          nGenesWithTh = nGenesWithTh, nGenesNoTh = nGenesNoTh,
                          meanLog2Bias.paired = meanLog2Bias.paired,
                          pval.grp = pval.grp, pval.paired = pval.paired,
                          meanLog2Bias.sum = meanLog2Bias.sum,
                          Ng = Ng, nl = Nl))

      }
      df <- do.call(rbind,
                    lapply(1:(dim(pairsMat)[2]),
                           function(x) getBiasForPair(pairsMat[1, x], pairsMat[2, x])))
      return(df)
    }
    df.chr <- do.call(rbind, lapply(1:length(gib[[chr]]), getBiasForBin, chr, pairsMat))
    return(df.chr)
  }
  df.genome <- do.call(rbind, lapply(names(gib), getBiasInChrBins, pairsMat))
  return(df.genome)
}

getGenomicBinsInGeneIDs <- function(binSize = 1e7, txdb.genRange, chrLens) {
  bins <- getGenomicbins(chrLens = chrLens, width = binSize)


  getGeneIDsForChr <- function(chr) {

    tt <- findOverlaps(query = bins[[chr]], subject = txdb.genRange)
    tt2 <- lapply(unique(queryHits(tt)), function(qh) subjectHits(tt)[queryHits(tt) == qh])
    qr <- bins[[chr]][unique(queryHits(tt))]

    tt3 <- lapply(1:length(tt2), function(i) list(range = qr[i],
                                                  geneIDs = names(txdb.genRange)[tt2[[i]]]))
    return(tt3)
  }

  genInBins <- lapply(1:length(bins), getGeneIDsForChr)
  names(genInBins) <- sapply(1:length(bins), function(x) as.character(bins[[x]][[1]]@seqnames))
  return(genInBins)
}



#This function gets the genomic coordinates in bins
getGenomicbins <- function(chrLens, width=1e7) {
  chrs <- chrLens$width
  names(chrs) <- chrLens$chr

  bins <- lapply(1:length(chrs),
                 function(i) tileGenome(chrs[i], tilewidth = min(width, chrs[i]),
                                        cut.last.tile.in.chrom = F))


  bins2 <- lapply(1:length(bins), function(i) shift(bins[[i]],
                                                    shift = chrLens[as.character(unlist(bins[[i]])@seqnames),"start"]-1))

  bins3 <- lapply(1:length(bins2), function(i) GRangesList(lapply((bins2[[i]]), trim)))
  return(bins3)



}


addRPKMtoHTobj <- function(ht, geneLen, minLibSize = 40000) {
  require(DESeq2)
  require(edgeR)
  idxs <- which(colSums(ht$counts) < minLibSize)
  if(length(idxs) > 0)
    ht$counts <- ht$counts[, -idxs]
  ht$dge <-  DGEList(ht$counts, group = names(ht$counts))
  ht$rpkm <- rpkm(x = ht$dge, normalized.lib.sizes = T, gene.length = geneLen)
  
  return(ht)
}

addDESeqToObj <- function(ht) {
  #size factor normalization: https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106
  ht$sfAll <- estimateSizeFactorsForMatrix(ht$counts)
  ht$DESeqExp <- t( t(ht$counts) / ht$sfAll )
  return(ht)
}


transformCountsMat2RPKM <- function(htmat,
                                    gtffile = "C:/Users/etai/Google Drive/DFCI/Data/gencode.v25.basic.annotation.gtf",
                                    minLibSize = 40000
) {
  require(GenomicFeatures)

  txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character())

  txdb.txLen <- transcriptLengths(txdb)
  txdb.geneLen <- max(splitAsList(txdb.txLen$tx_len, txdb.txLen$gene_id))
  ebg <- exonsBy(txdb, by="gene")
  txdb.genRange <- range(ebg)

  require(edgeR)
  htmat <- htmat[names(txdb.geneLen)[(names(txdb.geneLen) %in% rownames(htmat))],]
  htmat <- htmat[, -which(colSums(htmat) < minLibSize)]
  htdge <-  DGEList(htmat, group = names(htmat))
  htrpkm <- rpkm(x = htdge, normalized.lib.sizes = T, gene.length = txdb.geneLen)



}

bamFile2CountMat <- function(fname = "/czlab/etai/ExperimentsData/Stamatis/results_latest/10-062716-27-C5_S10_L001_R1_001.fastq.gz.Aligned.sortedByCoord.out.bam",
                             gtffile = "/czlab/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.annotation.gtf") {

  require("GenomicFeatures")
  if("txdb" %in% ls() == F) {
    txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character())

  } else {
    if(txdb$source != "gencode.v25.primary_assembly.annotation.gtf") {
      txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character())
    }
  }
}

getStrandSpecificPrimaryGeneCounts <- function(bamfilesdir = "/czlab/etai/ExperimentsData/cz/SunHur/091617/star/") {
  require(Rsamtools)
  require(AnnotationDbi)
  require(GenomicFeatures)
  
  txdb <- loadDb(file = "/czlab/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.sqlite") 
  ebg <- exonsBy(txdb, by="gene")
  
  filenames <- list.files(bamfilesdir, pattern = "*.R1.fastq.Aligned.sortedByCoord.out.bam$", full.names = T)
  bamfiles <- BamFileList(filenames, yieldSize=2000000)
  
  require("BiocParallel")
  register(BatchJobsParam(workers = 12), default = TRUE)
  
  sbf.f <- scanBamFlag(isPaired = T, isProperPair = T, isUnmappedQuery = F, 
                     hasUnmappedMate = NA, isMinusStrand = NA, isMateMinusStrand = T,
                     isFirstMateRead = NA, isSecondMateRead = NA, 
                     isSecondaryAlignment = F, isNotPassingQualityControls = F,
                     isDuplicate = F)
  sbf.r <- scanBamFlag(isPaired = T, isProperPair = T, isUnmappedQuery = F, 
                       hasUnmappedMate = NA, isMinusStrand = T, isMateMinusStrand = NA,
                       isFirstMateRead = NA, isSecondMateRead = NA, 
                       isSecondaryAlignment = F, isNotPassingQualityControls = F,
                       isDuplicate = F)
  sbf <- scanBamFlag(isPaired = T, isProperPair = T, isUnmappedQuery = F, 
                     hasUnmappedMate = NA, isMinusStrand = NA, isMateMinusStrand = NA,
                     isFirstMateRead = NA, isSecondMateRead = NA, 
                     isSecondaryAlignment = F, isNotPassingQualityControls = F,
                     isDuplicate = F)
  
  exp.f2 <- summarizeOverlaps(features=ebg, reads=bamfiles, 
                             param=ScanBamParam(flag = sbf.f, tag=c("NH"), tagFilter=list(NH=c(1))),
                             mode="Union",
                             singleEnd=FALSE,
                             ignore.strand=F,
                             fragments=TRUE )
  
  exp.r2 <- summarizeOverlaps(features=ebg, reads=bamfiles, 
                             param=ScanBamParam(flag = sbf.r, tag=c("NH"), tagFilter=list(NH=c(1))),
                             mode="Union",
                             singleEnd=FALSE,
                             ignore.strand=F,
                             fragments=TRUE )
  
  exp <- summarizeOverlaps(features=ebg, reads=bamfiles, inter.feature = F,
                             param=ScanBamParam(flag = sbf, tag=c("NH"), tagFilter=list(NH=c(1))),
                             mode="Union",
                             singleEnd=FALSE,
                             ignore.strand=FALSE,
                             fragments=TRUE )
  
  exp.all <- summarizeOverlaps(features=ebg, reads=bamfiles, inter.feature = F,
                           param=ScanBamParam(flag = sbf), #tag=c("NH"), tagFilter=list(NH=c(1))),
                           mode="Union",
                           singleEnd=FALSE,
                           ignore.strand=FALSE,
                           fragments=TRUE )
  save(exp, exp.r, exp.f, exp.all, file = "/czlab/etai/ExperimentsData/cz/SunHur/091617/exps.RData")
  load("/czlab/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.geneRanges.RData")
  gr <- geneRanges[geneRanges$seqnames != "chrM",]
  
  df <- data.frame(assay(exp))
  unstranded <- merge(gr, df, by.x="gene_id", by.y="row.names")
  df <- data.frame(assay(exp.f))
  antisense <- merge(gr, df, by.x="gene_id", by.y="row.names")
  df <- data.frame(assay(exp.f2))
  antisense2 <- merge(gr, df, by.x="gene_id", by.y="row.names")
  
  df <- data.frame(assay(exp.r))
  sense <- merge(gr, df, by.x="gene_id", by.y="row.names")
  df <- data.frame(assay(exp.r2))
  sense2 <- merge(gr, df, by.x="gene_id", by.y="row.names")
  
  df <- data.frame(assay(exp.all))
  all <- merge(gr, df, by.x="gene_id", by.y="row.names")

}

bamFileDir2RPKMmat <- function(bamfilesdir = "/cga/meyerson/Data/SingleCellSeq/RPE-1/RPE_RNA/Pool2",
                               gtffile = "/cga/meyerson/References/gtex_resources/gencode.v19.transcripts.patched_contigs.gtf",
                               doTranscriptLevel = F) {

  list.files(bamfilesdir, pattern = "*.bam$")
  filenames <- file.path(bamfilesdir, list.files(bamfilesdir, pattern = "*.bam$"))

  require("Rsamtools")
  bamfiles <- BamFileList(filenames, yieldSize=2000000)

  #Annotations:
  require("GenomicFeatures")
  if("txdb" %in% ls() == F)
    (txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character()))
  txdb.length <- transcriptLengths(txdb)
  geneLen <- max(splitAsList(txdb.length$tx_len, txdb.length$gene_id))

  (ebg <- exonsBy(txdb, by="gene"))
  (ebt <- exonsBy(txdb, by="tx", use.names = T))

  require("GenomicAlignments")
  require("BiocParallel")
  register(BatchJobsParam(workers = 24), default = TRUE)

  se.g <- summarizeOverlaps(features=ebg, reads=bamfiles,
                            mode="Union",
                            singleEnd=FALSE,
                            ignore.strand=TRUE,
                            fragments=TRUE )

  se.t <- summarizeOverlaps(features=ebt, reads=bamfiles,
                            mode="Union",
                            singleEnd=FALSE,
                            ignore.strand=TRUE,
                            fragments=TRUE )
  head(assay(se.t), 3)
  head(assay(se.g), 3)

  require(edgeR)
  se.t.edger <- DGEList(assay(se.t), group=rownames(colData(se.t)))
  se.g.edger <- DGEList(assay(se.g), group=rownames(colData(se.g)))

  #checking names - should be zero
  sum(rownames(se.t.edger) == txdb.length$tx_name)-length(rownames(se.t.edger))
  sum(rownames(se.g.edger) == names(geneLen))-length(rownames(se.g.edger))
  length(names(geneLen))

  se.t.rpkm <- rpkm(x = se.t.edger, normalized.lib.sizes = T, gene.length = txdb.length$tx_len)
  se.g.rpkm <- rpkm(x = se.g.edger, normalized.lib.sizes = T, gene.length = geneLen)
  colnames(se.g.rpkm) <- do.call(rbind, strsplit(colnames(se.g.rpkm), "\\."))[,1]

  return(list(rpkm = list(tran = se.t.rpkm, gen = se.g.rpkm),
              counts = list(tran = se.t, gen = se.g, metodMode = "union"),
              geneLen = geneLen))

}

