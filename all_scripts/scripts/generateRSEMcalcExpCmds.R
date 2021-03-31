args <- commandArgs(trailingOnly = TRUE)
#Rscript ~/WORK/RNA/Workflows/generateRSEMcalcExpCmds.R /czlab/etai/ExperimentsData/Stamatis/run.v2/fastq/ /czlab/etai/ExperimentsData/Stamatis/run.v2/rsem/ /path/to/scriptsdir/ > run.v2/rsem.cmds
baseDir <- args[1]
rsemoutdir <- args[2]
scriptsdir <- args[3]
addRefGenome <- F
#print("fasdfadadsfads")
#print(args)
#print(length(args))
if(length(args) == 6) {
  addRefGenome <- T
  gtffile <- args[4] 
  rsemrefpath <- args[5] 
  refgenome <- args[6]
}


#quit("no")

#NOTE: "./" was always used as no arugment was provided in last line of code "cat(gen..."  check behavior if no scripts dir argument is specified 
#  	 	I would then need to add an if statement for the scriptsdir or not in the last line of code. 
#Mind that this function works only on .R2. and .R1. file format
generateRSEMcmds <- function(baseDir = "/singlecellcenter/etai/ExperimentsData/Stamatis/run.v2/fastq",
                             rsemoutdir = "/singlecellcenter/etai/ExperimentsData/Stamatis/run.v2/rsem/",
                             scriptsdir = "./") {
  fastqFiles <- list.files(path = baseDir, recursive = F, 
                           pattern = ".*\\.R1\\..*fastq.gz$|.*\\.R1\\..*fastq$|.*_R1.*fastq$|.*_R1.*fastq.gz$|.*_1\\.fastq.gz$|.*\\.1\\.fastq.gz$", full.names = T)
  
  isPair <- function(fname) {
    if(!file.exists(fname)) {
      stop("rsem Error R1 not found!")
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
    outprefix <- sprintf("%s/%s", rsemoutdir, basename(fname))
    if(isPair(fname)) {
      R2 <- gsub(pattern = "\\.R1\\.", replacement = ".R2.", x = fname)
      if(R2 == fname) {
        message("Cannot produce R2 file for paired end read files: ", R2, "\n", fname, "\nTrying another format..\n")
        R2 <- gsub(pattern = "_R1_", replacement = "_R2_", x = fname)
        if(R2 == fname) {
          message("Cannot produce R2 file for paired end read files - again..")  
          R2 <- gsub(pattern = "_1\\.fastq", replacement = "_2.fastq", x = fname)
          if(R2 == fname) {
            R2 <- gsub(pattern = "\\.1\\.fastq", replacement = ".2.fastq", x = fname)
            if(R2 == fname) {
              message("Cannot produce R2 file for paired end read files - quiting.")  
              return(NA)
            }
          }
        }
      }
# I REMOVED THE / BETWEEN THE "%s" AND "runRSEMcalcExp.sh" 
      if(addRefGenome) {
        return(sprintf("bash %s/runRSEMcalcExp.sh 2 %s %s %s %s %s %s 1>%s.runRSEMcalcExp.stdouterr 2>&1\n", 
                       scriptsdir, fname, R2, outprefix, gtffile, rsemrefpath, refgenome, outprefix))    
      } else {
        return(sprintf("bash %s/runRSEMcalcExp.sh 2 %s %s %s 1>%s.runRSEMcalcExp.stdouterr 2>&1\n", 
                       scriptsdir, fname, R2, outprefix, outprefix))    
      }
    } else {
      if(addRefGenome) {
        return(sprintf("bash %s/runRSEMcalcExp.sh 1 %s %s %s %s %s 1>%s.runRSEMcalcExp.stdouterr 2>&1\n", 
                       scriptsdir, fname, outprefix, gtffile, rsemrefpath, refgenome, outprefix))
      } else {
        return(sprintf("bash %s/runRSEMcalcExp.sh 1 %s %s 1>%s.runRSEMcalcExp.stdouterr 2>&1\n", 
                       scriptsdir, fname, outprefix, outprefix))    
      }
    }
  }
  cmds <- sapply(fastqFiles, printCMD)
  return(cmds)
}

cat(generateRSEMcmds(baseDir = baseDir, rsemoutdir = rsemoutdir, scriptsdir = scriptsdir))
