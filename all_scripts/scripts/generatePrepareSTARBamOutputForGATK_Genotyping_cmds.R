#generateRun_STARcmds.R
#Rscript ~ejacob/WORK/RNA/Workflows/generateRun_STARcmds.R /czlab/etai/ExperimentsData/Stamatis/run.v2/fastq/ /czlab/etai/ExperimentsData/Stamatis/run.v2/star/ > star.cmds
args <- commandArgs(trailingOnly = TRUE)

baseDir <- args[1] # /czlab/etai/ExperimentsData/PRJEB20121/star
outdir <- args[2] # /czlab/etai/ExperimentsData/PRJEB20121/gatk
pattern <- args[3] # *.fastq.gz.Aligned.sortedByCoord.out.bam$
scriptsdir <- args[4]

stopifnot(length(args) >= 2)

if(length(args) == 2)
  pattern <- "*.R1.fastq.gz.Aligned.sortedByCoord.out.bam$"

#bash ~/WORK/RNA/Workflows/PrepareSTARBamOutputForGATK_Genotyping.sh 
#/czlab/etai/ExperimentsData/Stamatis/run.v2/star/170216_A9.R1.fastq.gz.Aligned.sortedByCoord.out.bam 
#/czlab/etai/ExperimentsData/Stamatis/run.v2/gatk


getPrepare_GATKcmds <- function(baseDir = "/singlecellcenter/etai/ExperimentsData/Stamatis/run.v2/star/",
                                outdir = "/singlecellcenter/etai/ExperimentsData/Stamatis/run.v2/gatk/",
                                pattern = "*.R1.fastq.gz.Aligned.sortedByCoord.out.bam$",
                                scriptsdir = "./") {
 
  
  bamFiles <- list.files(path = baseDir, recursive = F, pattern = pattern, full.names = T)
  
  printCMD <- function(fname) {
    cmd <- sprintf("bash %s/PrepareSTARBamOutputForGATK_Genotyping.sh %s %s\n",
                   scriptsdir, fname, outdir)
    return(cmd)
  }
  cmds <- sapply(bamFiles, printCMD)
  return(cmds)
}

cat(getPrepare_GATKcmds(baseDir = baseDir, outdir = outdir, pattern = pattern, scriptsdir = scriptsdir))
  
  
  

