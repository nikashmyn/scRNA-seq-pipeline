

getKEGGgeneListByKEGGID <- function(keggid = "mmu04623") {
  require(KEGGREST)

  query <- keggGet(keggid)
  df <- data.frame(do.call(rbind, strsplit(grep(";", query[[1]]$GENE, value = T), ";")), stringsAsFactors = F)
  names(df) <- c("name", "description")
  return(df)
}


getTranscriptsDataForGenome <- function(genome = "hg38") {
  require(AnnotationDbi)
  require(GenomicFeatures)
  require(data.table)
  
  if(genome == "hg38") {
    txdb <- loadDb(file = "/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.sqlite")  
    tx_seq <-readDNAStringSet("/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.transcripts.fa")
  } else if(genome == "mm12") {
    txdb <- loadDb(file = "/singlecellcenter/etai/ReferenceData/Gencode/gencode.vM12.ERCC92.annotation.sqlite") 
    tx_seq <-readDNAStringSet("/singlecellcenter/etai/ReferenceData/Gencode/gencode.vM12.transcripts.fa")
  } else {
    stop("No other genome than hg38 or mm12 is supported at this point.")
  }
  
  lens <- transcriptLengths(txdb)
  ids <- do.call(rbind, lapply(names(tx_seq), function(x) strsplit(x, "\\|")[[1]]))
  comp <- lapply(1:length(tx_seq), function(x) letterFrequency(tx_seq[x][[1]], letters = c("ACGT"), OR = 0))
  comp2 <- do.call(rbind, comp)
  comp3 <- data.frame(ids, comp2)
  names(comp3)[1:2] <- c("tx_name", "gene_name")
  comp4 <- merge(lens, comp3, by = "tx_name")
  comp4$gene_name <- NULL

  return(comp4)
}
#genes is a matrix/data.frame/data.table with two fields: gene_name and gene_id
#genes <- read.table("/singlecellcenter/etai/ExperimentsData/10X/10Xweb/Datasets/pbmc4k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/GRCh38/genes.tsv")
#names(genes) <- c("gene_id", "gene_name")
seurat2CountData <- function(seuratObj, genes, use.raw.data = F) {
#exp(range(pbmc@data)-1)
  if(use.raw.data) {
    ys <- match(colnames(seuratObj@data), colnames(seuratObj@raw.data))
    xs <- match(rownames(seuratObj@data), rownames(seuratObj@raw.data))
    mat <- seuratObj@raw.data[,ys]
    mat <- mat[xs, ]
  } else {
    mat <- seuratObj@data
  }
  
  mat <- data.table(as.matrix(mat), keep.rownames = T)
  
  names(mat)[1] <- "gene_name"
  setkey(mat, gene_name)
  
  genes <- data.table(genes)
  setkey(genes, "gene_name")

  myident <- seuratObj@ident
  mat2 <- merge(genes, mat, all = F)
  if(sum(names(mat2)[-(1:2)] != names(myident)) > 0)
    stop("seurat ident names do not equal mat2 names")
  names(mat2)[-(1:2)] <- paste("C", as.character(myident), ".", 1:length(myident), sep = "")
  
  return(mat2)

}

getSeqLengthsFromFastaFile <- function(fa.name="/singlecellcenter/etai/ReferenceData/Yeast/GCA_000149845.2/release35/dna/Schizosaccharomyces_japonicus.GCA_000149845.2.dna_sm.toplevel.fa") {
  fa <- read.fasta(fa.name, seqtype = "DNA")
  seqLengths <- sapply(fa, length)
  return(seqLengths)

}
#get all genes in a given cytoband and chromosome:
getGenesFromCytoBand <- function(cytoband = "q", chr = "chr5", genome = "hg38") {
  require(biovizBase)
  require(GenomicRanges)
  IdeogramCyto <- biovizBase::getIdeogram(genome = genome, cytobands = TRUE)
  q <- IdeogramCyto[IdeogramCyto@seqnames == chr,]
  cytoq <- q[grep(cytoband,as.character(q@elementMetadata$name), value = F),]
  
  if((!("geneRanges" %in% ls()))) {
    cat("Loading geneRanges\n")
    load("/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.geneRanges.RData")
  }
  GR <- GenomicRanges::makeGRangesFromDataFrame(geneRanges, keep.extra.columns = T)
  tt <- GenomicRanges::findOverlaps(query = cytoq, subject = GR)
  tt2 <- lapply(unique(S4Vectors::queryHits(tt)), function(qh) GR@elementMetadata$gene_id[S4Vectors::subjectHits(tt)[S4Vectors::queryHits(tt) == qh]])
  names(tt2) <- as.character(cytoq@elementMetadata$name[unique(S4Vectors::queryHits(tt))])
  return(tt2)
}



#input can an integer vector or a character string (comma separate) of same size 
#in case of strings - splits by ',' and c
#function returns the total genomic length of all exons together excluding overlapping regions (union)
#usage example: 
#res <- do.call(rbind, lapply(1:length(mmRanges$exonStarts), 
#                             function(x) { print(x); calculateTotalLengthFromExonlist(mmRanges$exonStarts[x], mmRanges$exonEnds[x]) }))
calculateTotalLengthFromExonlist <- function(start, end) {
  if(is.na(start) | is.na(end))
    return(NA)
  require(IRanges)
  if(is.character(start)) {
    start <- as.numeric(strsplit(start, ",")[[1]])
  }
  if(is.character(end)) {
    end <- as.numeric(strsplit(end, ",")[[1]])
  }
  stopifnot(length(start) == length(end))

  ir <- IRanges(start = start, end = end)
  #return(ir)
  iru <- IRanges::union(x = ir, y=ir)
  return(data.frame(w.union = sum(width(iru)), w.original = sum(width(ir))))

}
getChrSizesFromUCSC <- function(genome = "hg38") {
  if(genome == "hg38") {
    require(data.table)
    chrSizes <- fread("https://genome.ucsc.edu/goldenpath/help/hg38.chrom.sizes")
    cs <- chrSizes$V2
    names(cs) <- chrSizes$V1
    chrSizes <- cs
    return(chrSizes)
  } else {
    return(NA)
  }
}

getGenomicArmsFromBioviz <- function(genome = "hg38") {
  require(IRanges)
  ig <- getIdeogramFromBioviz(genome = genome)
  ig$arm <- gsub("[0-9|.]+", "", as.character(ig$name))
  #q arm:
  ig.q <- IRanges::union(x = ig[ig$arm == "q"], y=ig[ig$arm == "q"])
  ig.q$arm <- paste(gsub("chr", "", seqnames(ig.q)), "q", sep = "")
  #p arm:
  ig.p <- IRanges::union(x = ig[ig$arm == "p"], y=ig[ig$arm == "p"])
  ig.p$arm <- paste(gsub("chr", "", seqnames(ig.p)), "p", sep = "")
  garms <- c(ig.p, ig.q)
  #garms$width <- width(garms)
  
  return(garms)
}

generateChrSizesFromFastaFile <- function(refFasta = "/singlecellcenter/etai/ReferenceData/Gencode/GRCm38.primary_assembly.genome.fa",
                                          ifSave = F) {
  require(seqinr)
  refFasta <- read.fasta(refFasta)
  refSizes <- sapply((refFasta), length)
  chrSizes <- refSizes
  if(ifSave)
    saveRDS(chrSizes, file = sprintf("%s.chrSizes.rds", refFasta))
  return(refSizes)

}

getCentromereLocationsForHg38 <- function() {
  require(gtools)
  url <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/centromeres.txt.gz"
  tmp <- tempfile()
  download.file(url,tmp)
  tbl <- read.csv(tmp, header = F, sep = "\t", stringsAsFactors = F)

  #names of columns are taken from http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/centromeres.sql
  names(tbl) <- c("bin", "seqnames", "start", "end", "name")
  tbl <- data.table(tbl)
  tbl <- tbl[mixedorder(seqnames, start, decreasing = F),]  
  return(tbl)
}
getIdeogramFromBioviz <- function(genome = "hg38") {
  require(biovizBase)
  require(GenomicRanges)
  require(gtools)
  
  ig <- biovizBase::getIdeogram(genome = genome, cytobands = TRUE)  
  if(genome == "hg38" | genome == "hg19" | genome == "hg18")
    ig <- ig[ seqnames(ig) %in% sprintf("chr%s", c(1:22, "X", "Y", "MT", "M")),]
  seqlevels(ig) <- (mixedsort(as.character(unique(seqlevels(ig)))))

  return(ig)
  #chr sies of hg38 can also be taken from:
  #https://genome.ucsc.edu/goldenpath/help/hg38.chrom.sizes
}

getIdeogramFromUCSC <- function(genome = "hg38") {
  require(GenomicRanges)
  require(rtracklayer)
  
  mySession <- try(browserSession("UCSC"), silent=TRUE)
  # In case it fails, use another mirror
  if(inherits(mySession, "try-error"))
    mySession <- browserSession("UCSC",
                                url="http://genome-euro.ucsc.edu/cgi-bin/")
  genome(mySession) <- genome
  #obj <- ucscTableQuery(mySession, table="gap")
  obj <- ucscTableQuery(mySession, table="cytoBand")
  tbl <- getTable(obj)
  ig <- GenomicRanges::makeGRangesFromDataFrame(tbl, keep.extra.columns = T)
  ig <- ig[ seqnames(ig) %in% sprintf("chr%s", c(1:22, "X", "Y", "MT")),]
  return(ig)
}

gene_id2HGNC <- function(gene_ids) {
  if(!("gencodeMapper" %in% ls())) {
    cat("Loading gencodeMapper.\n")
    load("~ejacob/WORK/Annotations/data/gencodeMapper.HGNC.tx_name_gene_id.RData")
  } 
  
  if(is.null(key(gencodeMapper)) || key(gencodeMapper) != "gene_id")
    setkey(gencodeMapper, gene_id)
  return(unique(gencodeMapper[gene_ids][,HGNC, gene_id]))
  #return(as.character(data.frame(na.omit(unique(gencodeMapper[gene_ids][,HGNC])))[,1]))
  
}
geneListToGencodeID <- function(gl, id_type = c("HGNC", "tx_name", "gene_id"), gencodeMapper = NULL) {
  if(is.null(gencodeMapper)) {
    cat("Loading gencodeMapper.\n")
    load("~ejacob/WORK/Annotations/data/gencodeMapper.HGNC.tx_name_gene_id.RData")
  } 
  if(!(id_type[1] %in% names(gencodeMapper)))
    stop("No such id_type")
  if(is.null(key(gencodeMapper)) || key(gencodeMapper) != id_type[1])
    setkeyv(gencodeMapper, id_type[1])
  
  return(as.character(data.frame(na.omit(unique(gencodeMapper[gl][,gene_id])))[,1]))
}

generate_HGNC_TX_GENE_map_table <- function() {
  require(AnnotationDbi)
  require(GenomicFeatures)
  require(data.table)
  txdb <- loadDb(file = "/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.sqlite")

  hgnc <- read.table(file = "/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.metadata.HGNC", 
                     sep = "\t", header = F, stringsAsFactors = F)
  names(hgnc) <- c("tx_name", "HGNC")
  df <- as.data.frame(transcripts(txdb, columns = c("tx_name", "gene_id")))
  if(sum(sapply(df$gene_id, length)>1))
    stop("multiple maps to genes!")

  df$gene_id <- unlist(df$gene_id)
  gencodeMapper <- merge(df, hgnc, all=T)
  gencodeMapper <- as.data.table(gencodeMapper)
  save(gencodeMapper, file = "~ejacob/WORK/Annotations/data/gencodeMapper.HGNC.tx_name_gene_id.RData")
  return(gencodeMapper)

}
generateGeneSetsFlatTableFromMsigDB <- function(gs.list) {
  gsl <- getGeneSetsFromMsigDB(gs.list)
  rns <- unique(unlist(gsl))
  df <- data.frame(row.names= rns)
  for(i in 1:length(gsl)) {
    df[gsl[[i]], names(gsl)[i]] <- gsl[[i]]
  }
  return(df)
}

getGeneSetsFromMsigDB <- function(gs.list) {
  if(!("msigdb" %in% ls()))
    load("~ejacob/WORK/Annotations/data/MsigDB/msigdb.v5.2.symbols.gmt.RData")
  
  gsidxs <- na.omit(match(gs.list, msigdb$geneset.names))
  if(length(gsidxs) == 0) {
    stop("No gene set was found in the list")
  }
  gsl <- msigdb$genesets[gsidxs]
  names(gsl) <- msigdb$geneset.names[gsidxs]
  return(gsl)
}
genrateMsigDBRData <- function() {
  #msigdb is downloaded from: http://software.broadinstitute.org/gsea/downloads.jsp before calling this function
  msigdb <- GSA.read.gmt(filename = "~ejacob/WORK/Annotations/data/MsigDB/msigdb.v5.2.symbols.gmt")
  save(msigdb, file = "~ejacob/WORK/Annotations/data/MsigDB/msigdb.v5.2.symbols.gmt.RData")

}
#downloaded from: https://github.com/cran/GSA/blob/master/R/GSA.read.gmt.R
GSA.read.gmt=function(filename){
  #
  ## Read in and parse a gmt file (gene set file) from the  Broad institute
  # this is tricky, because each lines (geneset) has a variable length
  #  I read the file twice, first to pick up the geneset name and description
  # in the   first two  columns, then I read it all in as a long string
  
  # The beginning and end of each gene set in the string
  # is determined by matching
  # BOTH on  geneset name and description (since geneset names sometimes
  # occur as genenames elsewhere in the file)
  
  a=scan(filename,what=list("",""),sep="\t", quote=NULL, fill=T, flush=T,multi.line=F)
  geneset.names=a[1][[1]]
  
  geneset.descriptions=a[2][[1]]
  
  dd=scan(filename,what="",sep="\t", quote=NULL)
  
  
  nn=length(geneset.names)
  n=length(dd)
  ox=rep(NA,nn)
  
  ii=1
  for(i in 1:nn){
    cat(i)
    while((dd[ii]!=geneset.names[i]) | (dd[ii+1]!=geneset.descriptions[i]) ){
      ii=ii+1
    }
    ox[i]=ii
    ii=ii+1
  }
  
  genesets=vector("list",nn)
  
  for(i in 1:(nn-1)){
    cat(i,fill=T)
    i1=ox[i]+2
    i2=ox[i+1]-1
    geneset.descriptions[i]=dd[ox[i]+1]
    genesets[[i]]=dd[i1:i2]
  }
  
  geneset.descriptions[nn]=dd[ox[nn]+1]
  genesets[[nn]]=dd[(ox[nn]+2):n]
  out=list(genesets=genesets,geneset.names=geneset.names, geneset.descriptions=geneset.descriptions)
  class(out)="GSA.genesets"
  return(out)
}

#if v is not a data.table with a field of gene_id then a data.frame or a matrix with a row.names as gene_id should be provided
mergeGeneGenomicRangesWithCountData <- function(v, exclude_suffix = F) {
  require(GenomicRanges)
  if(!("chrSizes" %in% ls())) {
    chrSizes <- getChrSizes()
  }
  if((!("geneRanges" %in% ls()))) {
    cat("Loading geneRanges\n")
    load("/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.geneRanges.RData")
  } 
  if(exclude_suffix)
    geneRanges$gene_id <- sub(pattern = '\\.\\d+$',    '', geneRanges$gene_id)
  
  if(is.null(key(geneRanges)) || key(geneRanges) != "gene_id")
    setkey(geneRanges, gene_id)
  
  #print(class(gns))
  #gns <- geneRanges #[rownames(v)]
  if(class(v)[1] != "data.table") {
    v <- data.table(v, keep.rownames = T)
    setnames(v, "rn", "gene_id")
  }
  if(is.null(key(v)) || key(v) != "gene_id")
    setkey(v, gene_id)
  AD <- merge(geneRanges, v, by.y="gene_id", by.x="gene_id", all.x=F)
  return(AD)
}
#v can be a matrix or a data.frame/table with any number of columns as long as the row.names 
#are gene_id based on gencode data
#geneID is assumed to be row.names
#Returns a list of data.tables where each entry in the list is a specific bin
#Works on any number of chrs
binDataTableByGenomicIntervalGivenGeneIDs <- function(binSize, v, exclude_suffix = F) {
  require(GenomicRanges)
  AD <- mergeGeneGenomicRangesWithCountData(v = v, exclude_suffix = exclude_suffix)
  if(!("chrSizes" %in% ls())) {
    chrSizes <- getChrSizes()
  }
  # if((!("geneRanges" %in% ls()))) {
  #   cat("Loading geneRanges\n")
  #   load("/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.geneRanges.RData")
  # } 
  # if(is.null(key(geneRanges)) || key(geneRanges) != "gene_id")
  #   setkey(geneRanges, gene_id)
  # 
  # #print(class(gns))
  # #gns <- geneRanges #[rownames(v)]
  # v <- data.table(v, keep.rownames = T)
  # setnames(v, "rn", "gene_id")
  # setkey(v, gene_id)
  # AD <- merge(geneRanges, v, by.y="gene_id", by.x="gene_id", all.x=F)
  
  # if(sum(rownames(v) != rownames(gns))) {
  #   cat("ERROR!!! - wrong number of dims ", sum(rownames(v) != rownames(gns)))
  # }
  # AD <- cbind(gns, v)
  #AD <- AD[!is.na(AD$start),]
  require(gtools)
  if(sum(is.na(AD)) > 0)
    cat("\nbinDataTableByGenomicIntervalGivenGeneIDs::ERROR!!! - NA introduced in AD matrix: ", sum(is.na(AD)), "\n")
  chrs <- mixedsort(as.character(na.omit(unique(AD$seqnames))))
  
  bins <- getGenomicbins.annoUtils(chrSizes = chrSizes[chrs, ], width = binSize)
  AD.gr <- GenomicRanges::makeGRangesFromDataFrame((AD), keep.extra.columns = T)
  require(dplyr)
  getGeneIDsForChr <- function(chr) {
    tt <- GenomicRanges::findOverlaps(query = bins[[chr]], subject = AD.gr)
    tt2 <- lapply(unique(S4Vectors::queryHits(tt)), function(qh) S4Vectors::subjectHits(tt)[S4Vectors::queryHits(tt) == qh])
    names(tt2) <- start(bins[[chr]][unique(S4Vectors::queryHits(tt))])
    mypfunc <- function(x) { 
      df <- data.frame(AD.gr[x]) 
      df <- dplyr::select(df, gene_id, everything())
      names(df)[7:(dim(df)[2])] <- names(AD)[7:(dim(AD)[2])]
      return(df)
    }
    
    tt3 <- lapply(tt2, mypfunc)

    return(tt3)
  }
  
  allbins <- lapply(1:length(bins), getGeneIDsForChr)
  names(allbins) <- sapply(1:length(bins), function(x) as.character(bins[[x]]@seqnames[1]))
  return(allbins)  
}
#returns a list of row indices for each bin
#input must have chr start end in AD
#calculation is done only for a single chr which is represented by the first row in AD
binDataTableByGenomicInterval <- function(binSize = 1e7, AD) {
  require(GenomicRanges)
  if(!("chrSizes" %in% ls())) {
    chrSizes <- getChrSizes()
  }
  chr <- strsplit((AD$ID)[1], split = ":")[[1]][1]
  
  bins <- getGenomicbins.annoUtils(chrSizes = chrSizes[chr, ], width = binSize)
  AD.gr <- makeGRangesFromDataFrame(AD)
  
  getGeneIDsForChr <- function(chr) {
    tt <- findOverlaps(query = bins[[chr]], subject = AD.gr)
    tt2 <- lapply(unique(queryHits(tt)), function(qh) subjectHits(tt)[queryHits(tt) == qh])
    names(tt2) <- start(bins[[chr]][unique(queryHits(tt))])
    return(tt2)
  }
  
  allbins <- lapply(1:length(bins), getGeneIDsForChr)
  names(allbins) <- sapply(1:length(bins), function(x) as.character(bins[[x]]@seqnames[1]))
  return(allbins)
}



#This function gets the genomic coordinates in bins
#For chr or chrs
#Assumes chr starts at genomic position of 1
getGenomicbins.annoUtils <- function(chrSizes, width=1e7, iterateChrs = T) {
  require(GenomicRanges)
  chrs <-chrSizes
  if("length" %in% names(chrSizes)) {
    chrs <- chrSizes$length
    names(chrs) <- chrSizes$chr
  }
  if(iterateChrs) {
    bins <- lapply(1:length(chrs),
                   function(i) GenomicRanges::tileGenome(seqlengths = chrs[i], tilewidth = min(width, chrs[i]),
                                          cut.last.tile.in.chrom = T))
    names(bins) <- names(chrs)
  } else {
    bins <- GenomicRanges::tileGenome(seqlengths = chrs, tilewidth = width,
                                      cut.last.tile.in.chrom = T)
  }
  return(bins)
}

getGenomicLocationOfVariantsFromGRanges <- function(gr) {
  require(VariantAnnotation)
  require(AnnotationDbi)
  txdb <- loadDb("/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.sqlite")
  loc <- locateVariants(gr, txdb, AllVariants())
  return(loc)
}

getGenomicLocationOfVariantsFromGTFfile <- function(vcf) {
  require(VariantAnnotation)
  require(AnnotationDbi)
  txdb <- loadDb("/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.sqlite")
  rd <- rowRanges(vcf)
  loc <- locateVariants(rd, txdb, AllVariants())
  return(loc)
}


getChrSizes <- function(file = "/singlecellcenter/etai/ReferenceData/UCSC/chrSizes.hg38.RData") {
  load(file)
  return(chrSizes)
}

getChrSizesFromURL <- function(url ="https://genome.ucsc.edu/goldenpath/help/hg38.chrom.sizes")  { #"https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes") {
  chrSizes <- read.table(url)
  names(chrSizes) <- c("chr", "length")
  rownames(chrSizes) <- chrSizes$chr
  return(chrSizes)
}
#fl.ref should be located in a direcotry with writing permissions
generateTabixFileFromAnnotatedVCFFile <- function(fl.ref = "/singlecellcenter/etai/ReferenceData/RPE-1_Hets_GRCh38/RPE-1_Hets_annotation/RPE_hets.recal.ann.vcf"){
  require(Rsamtools)
  fl.indexed <- sprintf("%s.tab", fl.ref)
  zipped <- bgzip(fl.ref, fl.indexed)
  idx <- indexTabix(zipped, "vcf")
  tab <- TabixFile(zipped, idx)
}

transformhg19toHg38 <- function(hg19) {
  require("AnnotationHub")
  hub <- AnnotationHub()
  chain <- query(hub, 'hg19ToHg38')[[1]]
  require("rtracklayer")
  Hg38 <- liftOver(hg19, chain)
  #all.equal(Hg38, hg19)
  return(Hg38)
}

getGrangesContainingCDSs <- function(file = "/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.sqlite") {
  require(AnnotationDbi)
  txdb <- loadDb(file)
  tt <- cds(txdb, columns = "gene_id")
  return(tt)
}

getOnlyGenesContainingCDSs <- function(file = "/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.sqlite") {
  require(AnnotationDbi)
  txdb <- loadDb(file)
  tt <- as.data.frame(cds(txdb, columns = "gene_id"))
  
  return(unique(unlist(tt$gene_id)))
}


getAllGenesInChr <- function(chr, onlyCoding=F) {
  require(AnnotationDbi)
  require(GenomicFeatures)
  txdb <- loadDb(file = "/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.sqlite")
  load("/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.metadata.HGNC.gc2hgnc.RData")
  
  if(onlyCoding) {
    df <- as.data.frame(cds(txdb, filter=list(tx_chrom=chr), columns=c("gene_id", "tx_name")))
    if(sum(sapply(df$gene_id, length) != 1) >0 || sum(sapply(df$tx_name, length) != 1))
      stop("tx_name or gene_id has more than one id")
    df$gene_id <- unlist(df$gene_id)
    df$tx_name <- unlist(df$tx_name)
  } else {
    df <- as.data.frame(transcripts(txdb, filter=list(tx_chrom=chr), columns=c("gene_id", "tx_name")))
    if(sum(sapply(df$gene_id, length) != 1) >0)
      stop("gene_id has more than one id")
    df$gene_id <- unlist(df$gene_id)
  }
  return(merge(df, gc2hgnc, by.x="tx_name", by.y="ID"))
}

generateGencode2HGNCTable <- function(fname="/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.metadata.HGNC") {
  gc2hgnc <- read.table(fname, header = F)  
  names(gc2hgnc) <- c("ID", "HGNC")
  save(gc2hgnc, file = sprintf("%s.gc2hgnc.RData", fname))
}

#gets chromosome ranges based on genes within
getChrRangesFromTxdb.genRange <- function(txdb.genRange) {
  rng <- range(unlist(txdb.genRange))
  
  getChrRange <- function(chr) {
    #return(range(rng[rng@seqnames==chr,]@ranges))
    return(data.frame(chr=chr, range(rng[rng@seqnames==chr,]@ranges)))
  }
  chrLens <- do.call(rbind, lapply(seqlevels(rng), getChrRange))
  #chrs <- chrLens$width
  rownames(chrLens) <- chrLens$chr
  
  return(chrLens)
  
}

addGene_id_to_grangeObj <- function(gr, txdb) {
  genes <- genes(txdb)
  hits <- findOverlaps(gr, genes)
  gr$gene_id <- NA
  gr$gene_id[queryHits(hits)] <- genes[subjectHits(hits)]$gene_id

  return(gr)  
}

generateTranscriptAnnotationDB <- function() {
  require(GenomicFeatures)
  require(data.table)
  #human
  txdb <- makeTxDbFromGFF(file = "/singlecellcenter/etai/ReferenceData/Gencode/v27/gencode.v27.annotation.gtf", 
                          format="gtf", circ_seqs=DEFAULT_CIRC_SEQS)
  saveDb(txdb, file = "/singlecellcenter/etai/ReferenceData/Gencode/v27/gencode.v27.annotation.sqlite")
  txs <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))
  txsRanges <- as.data.table(sort(txs))
  setcolorder(txsRanges, c("tx_name", "seqnames", "start", "end", "width", "strand", "tx_id", "gene_id"))
  setkey(txsRanges, tx_name)
  save(txsRanges, file = "/singlecellcenter/etai/ReferenceData/Gencode/v27/gencode.v27.annotation.txsRanges.RData")
  return(list(txdb = txdb, txsRanges = txsRanges))
}


generate_geneRanges_data.table <- function(db="/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.sqlite") {
  require(GenomicFeatures)
  require(data.table)
  txdb <- loadDb(db)
  txdb.txLen <- transcriptLengths(txdb)
  txdb.genLen <- max(splitAsList(txdb.txLen$tx_len, txdb.txLen$gene_id))
  ebg <- exonsBy(txdb, by="gene")
  txdb.genRange <- range(ebg)
  
  geneRanges <- as.data.table(txdb.genRange)
  names(geneRanges)[2] <- "gene_id"
  geneRanges$group <- NULL
  setkey(geneRanges, gene_id)

  #save(geneRanges, file = "/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.geneRanges.RData")
  return(geneRanges)
}

generateTxdbData <- function(gtffile = "/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.gtf") {
  require(GenomicFeatures)
  
  txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character())
  
  txdb.txLen <- transcriptLengths(txdb)
  txdb.genLen <- max(splitAsList(txdb.txLen$tx_len, txdb.txLen$gene_id))
  ebg <- exonsBy(txdb, by="gene")
  txdb.genRange <- range(ebg)
  #save(txdb.txLen, txdb.genRange, txdb.genLen, ebg, file = "/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.gtf.txdb.RData")
  #saveDb(txdb, file = "/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.sqlite")
  #txdb <- loadDb("/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.sqlite")
  tmp <- as.data.frame(txdb.genRange)
  tmp <- sort(GRanges(tmp))
  geneRanges <- as.data.frame(tmp)
  
  return(list(txdb=txdb, txLen=txdb.txLen, genRange=txdb.genRange, genLen=txdb.genLen, ebg=ebg))
  
}



generateAnnotationGTFforGenomes <- function() {
  require(GenomicFeatures)
  #Human:
  v25 <- generateTxdbData(gtffile = "/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.gtf")
  #mouse:
  mm12 <- generateTxdbData(gtffile = "/singlecellcenter/etai/ReferenceData/Gencode/gencode.vM12.ERCC92.annotation.gtf")
  #Both:
  v25_mm12 <- generateTxdbData(gtffile = "/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25_vM12_ERCC92.annotation.gtf")

  txdb <- v25$txdb
  saveDb(txdb, file = "/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.sqlite")
  
  txdb <- mm12$txdb
  saveDb(txdb, file = "/singlecellcenter/etai/ReferenceData/Gencode/gencode.vM12.ERCC92.annotation.sqlite")
  
  txdb <- v25_mm12$txdb
  saveDb(txdb, file = "/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25_vM12_ERCC92.annotation.sqlite")
  
  ebg.hs <- v25$ebg
  ebg.mm <- mm12$ebg
  ebg.hs_mm <- v25_mm12$ebg
  save(ebg.hs, ebg.mm, ebg.hs_mm, file = "/singlecellcenter/etai/ReferenceData/Gencode/ebg.hs_v25.mm_mm12.hs_mm.RData")

  #save(txdb.txLen, txdb.genRange, txdb.genLen, ebg, file = "/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.gtf.txdb.RData")
  #saveDb(txdb, file = "/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.sqlite")
  #txdb <- loadDb("/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.sqlite")
  
  
}