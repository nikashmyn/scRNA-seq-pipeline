#generateRun_STARcmds.R
#Rscript ~ejacob/WORK/RNA/Workflows/generateRun_STARcmds.R /czlab/etai/ExperimentsData/Stamatis/run.v2/fastq/ /czlab/etai/ExperimentsData/Stamatis/run.v2/star/ /path/to/scripts/ > star.cmds
args <- commandArgs(trailingOnly = TRUE)

baseDir <- args[1]
outdir <- args[2]
scriptsdir <- args[3]

if(length(args) > 3) {
  genomeDir <- args[4]
  refDir <- args[5]
}

getRun_STARcmds <- function(baseDir = "/singlecellcenter/etai/ExperimentsData/Stamatis/run.v2/fastq/",
                            outdir = "/singlecellcenter/etai/ExperimentsData/Stamatis/run.v2/star/",
                            genomeDir="/singlecellcenter/etai/ReferenceData/STAR/GRCh38_Gencode.v25.ERCC92",
                            refDir="/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.gtf") {
 
  
  fastqFiles <- list.files(path = baseDir, recursive = F, 
                           pattern = ".*\\.R1\\..*fastq.gz$|.*\\.R1\\..*fastq$|.*_R1.*fastq$|.*_R1.*fastq.gz$|.*_1\\.fastq.gz$|.*\\.1\\.fastq.gz$", full.names = T)
  isPair <- function(fname) {
    if(!file.exists(fname)) {
      stop("Error R1 not found!")
    }
    R2 <- gsub(pattern = "\\.R1\\.", replacement = ".R2.", x = fname)
    if(R2 == fname) {
      R2 <- gsub(pattern = "_1\\.fastq", replacement = "_2.fastq", x = fname)
      if(R2 == fname) {
        R2 <- gsub(pattern = "\\.1\\.fastq", replacement = ".2.fastq", x = fname)
        if(R2 == fname) {
          R2 <- gsub(pattern = "_R1_001.fastq", replacement = "_R2_001.fastq", x = fname)
          if(R2 == fname) {
            message("Cannot create file..")
            stop("ERROR - no file found!")
          }
        }
      }
    }
    if(file.exists(R2)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  printCMD <- function(fname) {
    outprefix <- sprintf("%s/%s", outdir, basename(fname))
    if(isPair(fname)) {
      return(sprintf("python2 %s/run_STAR.py %s %s %s %s --isPairedEnd 1 1> %s.run_STAR.stdouterr 2>&1\n", 
                     scriptsdir, genomeDir, refDir, fname, outprefix, outprefix))
    } else {
      return(sprintf("python2 %s/run_STAR.py %s %s %s %s --isPairedEnd 0 1> %s.run_STAR.stdouterr 2>&1\n", 
                     scriptsdir, genomeDir, refDir, fname, outprefix, outprefix))
    }
  }
  cmds <- sapply(fastqFiles, printCMD)
  return(cmds)
}

if(length(args) > 3) {
  cat(getRun_STARcmds(baseDir = baseDir, outdir = outdir, genomeDir = genomeDir, refDir = refDir))
} else {
  #use default references
  cat(getRun_STARcmds(baseDir = baseDir, outdir = outdir))  
}

  
  
  

