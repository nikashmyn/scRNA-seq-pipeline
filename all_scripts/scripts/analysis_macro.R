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

#######################
### Source Packages ###
#######################

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(reshape)
require(data.table)
require(gtools)
require(readxl)
require(ggplot2)
require(matrixStats)
require(reshape2)
require(gridExtra)
require(dplyr)
require(ggthemes)
require(rlang)
require(stringr)

########################################
### Input Desired script Directories ###
########################################

#Collect and process ASE data
#WARNING: LONG RUN TIME
message("Running QC")
source( sprintf('%s/scripts/QC_calculations.R', scriptsdir) )

#Collect and process ASE data
#WARNING: LONG RUN TIME
message("Aggregating ASE")
source( sprintf('%s/scripts/ASE_calculations_v2.R', scriptsdir) )

#Collect and process mRNA transcript data
#WARNING: LONG RUN TIME
message("Aggregating TPM")
source( sprintf('%s/scripts/TPM_calculations.R', scriptsdir) )

#Calculation TPM ratios and adjust for cell specific effects
message("Normalizing data")
source( sprintf('%s/scripts/normalize_data.R', scriptsdir) )

#Aggregate by inverse variance at bin, arm, and chr levels
message("Binning data")
source( sprintf('%s/scripts/agg_by_inv_var.R', scriptsdir) )

#Curate reference distributions for each chr and arm 
message("Creating reference distributions")
source( sprintf('%s/scripts/create_ref_distr.R', scriptsdir) )

#Calculate p-values at chr and arm level
message("Calculating p-values")
source( sprintf('%s/scripts/pval_calculations_2.R', scriptsdir) )

#Automatically call abnormal CN families and global summary 
message("Analyzing results")
source( sprintf('%s/scripts/automated_analysis_3.R', scriptsdir) )

#Generate Visuals
message("Creating individual visuals")
source( sprintf('%s/plots/visual_generator_v2.R', scriptsdir) )

#Generate Summary Visuals
message("Creating summary visuals")
source( sprintf('%s/plots/summary_visual_generator.R', scriptsdir) )













