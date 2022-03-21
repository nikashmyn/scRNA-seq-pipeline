##########################
### Script Explanation ###
##########################

#-------------------------------------------------------
#This script was written by Nikos Mynhier 

########################################
### Input Desired script Directories ###
########################################

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/rerun_6_9_2021/data", "/pellmanlab/nikos/Stam_Etai_Data", "SIS1025a",  "SIS1025b", "SIS1025c", "SIS1025d", "SIS1025e", "SIS1025f_Lane1", "SIS1025f_Lane2", "SIS1025g_Lane1", "SIS1025g_Lane2", "SIS1025misc", "SIS1025targ")
print(args)
scriptsdir <- args[1]
wkpath <- dirpath <- args[2]
datadir <- args[3]
experiments <- args[4:length(args)]

controlSampleIDs <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs.rds", wkpath))


##################################
#### Read in ASE and TPM files ###
##################################
#
#TPM <- readRDS(file=sprintf("%s/aggregated_results/TPM.bygene.rds", dirpath))
#
#ASE <- readRDS(file=sprintf("%s/aggregated_results/ASE.bygene.rds", dirpath))
#
#########################
#### Calculate AS-TPM ###
#########################
#
#AS_TPM <- list()
#AS_TPM$A <- cbind(ASE$AF[,c(1:6)], (TPM[,-c(1:6)] * ASE$AF[,-c(1:6)]))
#AS_TPM$B <- cbind(ASE$AF[,c(1:6)], (TPM[,-c(1:6)] * (1-ASE$AF[,-c(1:6)])))
#
#saveRDS(AS_TPM, file=sprintf("%s/aggregated_results/AS-TPM.bygene.rds", dirpath))
#
#############################################################################################################
##End of Script


