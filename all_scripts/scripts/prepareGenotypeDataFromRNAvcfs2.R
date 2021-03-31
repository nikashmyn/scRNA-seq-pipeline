#for >= CNV12.1
args <- commandArgs(trailingOnly = TRUE)

#Feb. 27 2019
# /singlecellcenter/RPE-1/HaplotypePhasing/RPE-1_Haplotype.txt
# 
# You only need to look at the allele_linkedReads and allele_monosomies field.
# 
# Here’re a couple of examples:
#   
# 1. All possible phased sites: allele_linkedReads ~=0 | allele_monosomies~=0
# This is equivalent to abs(linkedReads_frac)>0 | abs(ref-alt)>0
# The error rate is ~5% (estimated from the percentage of disagreement between the two sets).
# 2. Confident phase from linked reads
# allele_linkedReads~=0 & abs(linkedReads_frac)>=0.5 & linkedReads.linkedReadsCount>=100
# 3. Confident phase from monosomies (the second term uses allelic fraction among all single-cell data as a filter and needs to account for trisomy in 10q)
# allele_monosomies~=0 & covered>=50 & (abs(allelicRatio-0.5)<=0.2 | (chr==’chr10’ & pos>=60784157 & abs(min(allelicRatio,1-allelicRatio)-0.33)<=0.15)))
# 
# When imposing 2 and comparing the two sets, the percentage of disagreement is ~1% (2.08 million total)
# When imposing 3 and comparing the two sets, the percentage of disagreement is ~2% (2.10 million total)
# When imposing both 2 and 3 the percentage of disagreement is < 0.5%


if(length(args) < 3) {
  stop("Wrong input.\n Input should be:\nRscript prepareGenotypeDataFromRNAvcfs.R numOfCluster outpath patterns..\nRscript prepareGenotypeDataFromRNAvcfs.R 12 /czlab/etai/ExperimentsData/Stamatis/May072018/vcfs /czlab/etai/ExperimentsData/Stamatis/rsem/vcf/bamfilesJan232017_RPE_hets.GT_UG.%s.vcf /czlab/etai/ExperimentsData/Stamatis/SN218/vcf/SN218.gatkBamFiles_RPE_hets.GT_UG.%s.vcf")
}

require(data.table)

#source("~ejacob/WORK/DNA/R/VCFUtils.R")

genotypes <- fread("/singlecellcenter/RPE-1/HaplotypePhasing/RPE-1_Haplotype.txt", sep = "\t", header = T)

getGenomicLocationOfVariantsFromGRanges <- function(gr) {
  require(VariantAnnotation)
  require(AnnotationDbi)
  txdb <- loadDb("/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.sqlite")
  loc <- locateVariants(gr, txdb, AllVariants())
  return(loc)
}

splitNameToGenomicCoordinates <- function(AD) {
  require(GenomicRanges)
  
  df <- data.frame(do.call(rbind, strsplit(x = AD$id, split = "[_:/]" )), stringsAsFactors = F)
  names(df) <- c("chr", "start", "refGT", "altGT")
  df$start <- df$end <- as.numeric(df$start)
  df$strand <- "*"
  df <- df[,c(1,2,5,6,3,4)]
  AD2 <- cbind(id = AD$id, df, AD[,-1])
  
  return(AD2)
}

#retrieving AD matrix from vcf file: 
getAlleicDepthAndGenomicLocationMatrixFromVCFfile <- function(fls = c("/singlecellcenter/etai/ExperimentsData/Stamatis/rsem/vcf/bamfilesJan232017_RPE_hets.GT_UG.22.vcf",
                                                                      "/singlecellcenter/etai/ExperimentsData/Stamatis/SN218/vcf/SN218.gatkBamFiles_RPE_hets.GT_UG.22.vcf")) {
                                                              #na.rm = F) {
  require(VariantAnnotation)  
  require(data.table)
  source("~ejacob/WORK/Annotations/R/annoUtils.R")
  
  ###################################################################
  #1. Getting genotype read depth data for all vcfs for a single chr:
  ###################################################################
  message("Getting genotype read depth data for all vcfs for a single chr..")
  
  svp <- ScanVcfParam(geno="AD")
  ADs <- NULL
  for(fl in fls) {
    message("Working on: ", fl, "\n")
    vcfAD <- readVcf(fl, param = svp)
    AD <- data.table(geno(vcfAD)$AD, keep.rownames = "id")
    print(dim(AD))
    setkey(AD, id)
    if(is.null(ADs)) {
      ADs <- copy(AD)
      setkey(ADs, id)
    } else {
      ADs <- merge(ADs, AD, all=T)
      setkey(ADs, id)
    }
  }
  
  #####################################################################
  #2. Split id to get genomic coordinates and get location annotations:
  #####################################################################
  message("Split id to get genomic coordinates and get location annotations..")
  
  ADs2 <- splitNameToGenomicCoordinates(ADs)
  rd <- GenomicRanges::makeGRangesFromDataFrame(ADs2[,1:5], keep.extra.columns = T)
  loc <- getGenomicLocationOfVariantsFromGRanges(gr = rd)
  loc.lvls <- c("coding", "intron", "spliceSite", "fiveUTR", "threeUTR", "intergenic", "promoter")
  mylocs <- data.table(id = ADs2$id)
  mylocs[, loc.lvls] <- NA
  for(locName in names(table(loc@elementMetadata$LOCATION))) {
    mylocs[unique(mcols(loc)[mcols(loc)$LOCATION==locName, "QUERYID"]), locName] <- TRUE 
  }
  ADs2[, id := NULL]
  ADs3 <- cbind(mylocs, ADs2)
  
  #######################################################################
  #3. Merge allelic RNA data with DNA genotyping data and organize table:
  #######################################################################
  message("Merge allelic RNA data with DNA genotyping data and organize table..")
  
  ADs4 <- merge(genotypes, ADs3, by.x = c("chr", "pos"), by.y = c("chr", "start"))
  setcolorder(ADs4, neworder = c("id", "chr", "pos", "end", "strand", names(ADs4)[!names(ADs4) %in% c("id", "chr", "pos", "end", "strand")]))
  setnames(ADs4, old = c("pos", "chr"), new = c("start", "seqnames"))
  ADs4[, altGT := NULL]
  ADs4[, refGT := NULL]
  
  # if(na.rm) {
  #   cat("Reducing matrix size..\n")
  #   nona <- sapply(1:length(ADs[[1]]), 
  #                  function(y) sum(unlist(lapply(1:length(ADs), 
  #                                                function(x) !is.na(ADs[[y,x]])))))
  #   ADs <- ADs[which(nona>0),]
  # }
  
  
  message("Completed allelic depth calc..\n")
  return(ADs4)
}
#tt20 <- getAlleicDepthAndGenomicLocationMatrixFromVCFfile(fls = c("/czlab/etai/ExperimentsData/Stamatis/rsem/vcf/bamfilesJan232017_RPE_hets.GT_UG.20.vcf",
#"/czlab/etai/ExperimentsData/Stamatis/SN218/vcf/SN218.gatkBamFiles_RPE_hets.GT_UG.20.vcf"))

AD2NumericMatrix <- function(AD, prefixIdxs = 1:32) {
  require(data.table)
  message("Executing AD for ", AD$seqnames[1], "..")
  
  message("Doing allele A..")
  A <- sapply(names(AD)[-prefixIdxs], function(x) do.call(rbind, AD[[x]])[,1, drop=F] )
  message("Doing allele B..")
  B <- sapply(names(AD)[-prefixIdxs], function(x) do.call(rbind, AD[[x]])[,2, drop=F] )
  message("Completing..")
  return(list(A = A, B = B, features = AD[, prefixIdxs, with=F]))
}



################################
# Input Variables:
################################
numOfCluster <- as.numeric(args[1]) #12
outpath <- args[2] #/czlab/etai/ExperimentsData/Stamatis/May072018/vcfs
patterns <- args[-c(1, 2)] #/czlab/etai/ExperimentsData/Stamatis/rsem/vcf/bamfilesJan232017_RPE_hets.GT_UG.%s.vcf /czlab/etai/ExperimentsData/Stamatis/SN218/vcf/SN218.gatkBamFiles_RPE_hets.GT_UG.%s.vcf
maxPrefixIdx <- 32

mychrs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22, "X")
################################
# 
################################



mypatterns <- sapply(patterns, function(x) sprintf(x, mychrs))
print(mypatterns)

myfunc <- function(i) {
  message("Working on chr ", mychrs[i], "...")
  print(mypatterns[i, ])
  tt <- getAlleicDepthAndGenomicLocationMatrixFromVCFfile(fls = mypatterns[i,])
  tt2 <- AD2NumericMatrix(tt, prefixIdxs = 1:maxPrefixIdx)
  saveRDS(object = tt2, sprintf("%s/alleles.ADs.v5.chr%s.rds", outpath, mychrs[i]))
}
require(snow)

if(numOfCluster > 2) {
  cl <- makeCluster(as.numeric(numOfCluster),type = "SOCK")
  #clusterEvalQ(cl, library(GenomicRanges))
  ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
  clusterExport(cl,c("genotypes","mypatterns", "mychrs", "outpath", "maxPrefixIdx", ex))
  parLapply(cl, 1:nrow(mypatterns), myfunc)
  stopCluster(cl)
} else {
  lapply(1:nrow(mypatterns), myfunc)
}




cat("Task completed.\n")
  