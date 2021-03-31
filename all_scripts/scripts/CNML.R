
# Positive set - total del, monosomy, trisomy, two disomy + one mono and trisomy,
# Negative set - native disomy, disomy from another control swapped with native disomy, disomy from the same pair samples (swapped and swapping) but different chrs.
# Check PR within chr
# Check PR between chrs - that would be like between different cell types.
# Check PR with mix sizes of deletions or same sizes.
# 
# 
# Cross validation between pairs of control and swappers.
# Features - window of 11 genes

#TODO:
# show that the reg model captures non linear correlations also
# show examples of emrichments in co-expression that are not achieved using pearson correlation
# 0. do feature selection - DONE
# 0. check performance of feature selected 
# 0. implement on regression model also - DONE - see if there is an improvement
# 1. need to do regression for only 0-1-2-3CN model. maybe it would work better
# 2. add geneSpecific feature of distance between genes and GC and gene length - interspace and exonGC have a minor effect - DONE.
# 7. validate set by genotypes
# 8. run on tirosh data - maybe need to zscore before train - DONE - some improvments but maybe need to train on the same data set before. 
# 9. generate population with the same alteration but different sample background that will show applicability - how small changes are sporadic
#    also will enable to say to what extend we can detect a sub clone with a specific alteration?
# 10. run using deep learning : https://tensorflow.rstudio.com/keras/articles/examples/      
# 12. Implement train using lgbm as part of the main function generateAndEstimatePredictionMultiClassModel in buildCNMLlearningAndTestSet or divide the generateAndEstimatePredictionMultiClassModel to two functions: build_training_Set and train
# 13. refine predict_aneuploidy.R for tri's and mono's such that allele level data is more usfult (e.g. tri -> A > 1.15B & A > 0.8 & B > 0.8) - DONE (no specific improvement but better to use)
# 5. add global features as % of genes detected, % of zeros from chr in sample - ? not sure might effect generalization - maybe if we train on each dataset before prediction then it should be added.
# 0. add option to do 4C or not since not sure if it works - DONE 
#History
########
#run bst_predict.regression and bst.predict on all dataset and analyze GSEA and pathview of tpm and lgsm - which method to use for DE?
# 1. add doRatioBeforeScaleAndZscore in buildLearningOrTestSet - DONE - has some improvements to scale only
# 2. Win ftrs size of 21 or 25 works worse than 31 - DONE
# 3. do CN insert size of variable length - DONE (no specific improvment with insert size 31 of 600) - works good when inseret size is 100 for all 5-31, ftrs importance looks good 5q works a little bit less good but mono/tri/di works better
# 2. win size of 41 works less good than 31. Try 25. - DONE
# 3. add MA to ftrs - works worse DONE
# 4. SSM features
# 5. Try learning with higher gene mean expression then 1 to see if it performs better on tirosh's data
# 2. *2gain - trisomy + monosomy, disomy + disomy (4CN) - DONE - works approx like the 0-3 CNs maybe less good.
# 3. 3gain - trisomy + disomy (5CN) - since 4CN is not performing better, no plan to apply
# 4. 4gain - trisomy and trisomy (6CN) - since 4CN is not performing better, no plan to apply
# 5. *do DE based on regression model instead of multi classification problem - DONE - need to do regression for only 0-1-2-3CN model. maybe it would work better
# 6. change CV gene specific features to mean tpm - DONE - no improvemnt, even worse
#1. add myTCV and mean for gene features - DONE
#1. implement cross validation (one sample out) - done 20/80 - takes too long to run
#2*. 2loss - add vector for doCalibration for each sample and add zero sample vector - DONE
#4. join chrs - DONE - evaluating performance
#6*. run on targMN 5q and other cases of aneuploidies - DONE
#8. mix chrs in learning - DONE
#9. add capping to learning set - DONE
#10. try implement with lightgbm: https://www.kaggle.com/nschneider/gbm-vs-xgboost-vs-lightgbm - DONE
#11. round features to 3 digits after point to avoid over fitting - DONE

#see:
#http://fastml.com/what-is-better-gradient-boosted-trees-or-random-forest/


subsampleSampleSetsByField <- function(cs, field = "classes", method = c("median", "min"), minFactor = 1.0) {
  require(data.table)
  if(method[1] == "median") {
    q <- median(table(cs[[field]]))*minFactor
  } else if(method[1] == "min") {
    q <- mean(table(cs[[field]]))*minFactor
  }
  
  return(cs[,.SD[sample(.N, min(q,.N))],by = field])
}
#training set construcution:
# F (no change) | T (point of change)
# -----------------------------------
# Zero: 0C      | 0 to 1, 0 to 2, 0 to 3, 0 to 4
# Mono: 1C      | 1 to 2, 1 to 3, 1 to 4
# diploid: 2C   | 2 to 3, 2 to 4
# Trisomy: 3C   | 3 to 4
# Tetrasomy: 4C | -
#------------------------------------
# In total we have 14 sub groups of training set we need to construct
generateLearningAndTestSamples.pointOfChange <- function(adt, minDetectedGenes = 4000, chrToExclude = c("chrX", "chrY", "chrM"), excludeIDs, trainProb = 0.8,
                                           totalLossCount = 20, do4C = T, fuzzyRegion = 2, 
                                           file.1C = "/singlecellcenter/etai/data/WORK/secondaryanalysis/methods_paper_results/data/list_of_predicted_monosomies3.csv",
                                           file.2C = "/singlecellcenter/etai/data/WORK/secondaryanalysis/methods_paper_results/data/list_of_predicted_disomies3.csv",
                                           file.3C = "/singlecellcenter/etai/data/WORK/secondaryanalysis/methods_paper_results/data/list_of_predicted_trisomies3.csv") {
  require(data.table)
  
  ##################
  #Data preparation:
  ##################
  samp_pool <- names(which(colSums(adt[,-c(1:4)]>0) > minDetectedGenes))
  samp_pool <- samp_pool[!samp_pool %in% excludeIDs]
  
  c1 = data.table(read.csv(file.1C, stringsAsFactors = F), stringsAsFactors = F)
  c1 <- c1[(!seqnames %in% chrToExclude) & (id %in% samp_pool)]
  c2 = data.table(read.csv(file.2C, stringsAsFactors = F), stringsAsFactors = F)
  c2 <- c2[(!seqnames %in% chrToExclude) & (id %in% samp_pool)]
  c3 = data.table(read.csv(file.3C, stringsAsFactors = F), stringsAsFactors = F)
  c3 <- c3[(!seqnames %in% chrToExclude) & (id %in% samp_pool)]
  
  
  
  c1$train <- F
  c1$train[sample(1:nrow(c1), replace = F, size = round(nrow(c1)*trainProb))] <- T
  
  c2$train <- F
  c2$train[sample(1:nrow(c2), replace = F, size = round(nrow(c2)*trainProb))] <- T
  
  c3$train <- F
  c3$train[sample(1:nrow(c3), replace = F, size = round(nrow(c3)*trainProb))] <- T
  
  ##################################
  #FALSE group (no point of change)
  ##################################
  #This gorup has typeOfEventToGenerate = "amplitude"
  #1C
  nativeIDs <- c1$id[unlist(lapply(unique(as.character(c1$seqnames)), function(chr) sample(which(c1$seqnames == chr), size = sum(c1$seqnames == chr), replace = F)))]
  pseudoIDs <- c1$id[unlist(lapply(unique(as.character(c1$seqnames)), function(chr) sample(which(c1$seqnames == chr), size = sum(c1$seqnames == chr), replace = F)))]
  c1a <- cbind(c1, nativeIDs = nativeIDs, pseudoIDs = pseudoIDs)
  setnames(c1a, old = "id", "altIDs")
  c1a$altIDs2 <- NA
  c1a$nativeIDs2 <- NA
  c1a$pseudoIDs2 <- NA
  c1a$classes <- "1C"
  
  #2C
  nativeIDs <- c2$id[unlist(lapply(unique(as.character(c2$seqnames)), function(chr) sample(which(c2$seqnames == chr), size = sum(c2$seqnames == chr), replace = F)))]
  pseudoIDs <- c2$id[unlist(lapply(unique(as.character(c2$seqnames)), function(chr) sample(which(c2$seqnames == chr), size = sum(c2$seqnames == chr), replace = F)))]
  c2a <- cbind(c2, nativeIDs = nativeIDs, pseudoIDs = pseudoIDs)
  setnames(c2a, old = "id", "altIDs")
  c2a$altIDs2 <- NA
  c2a$nativeIDs2 <- NA
  c2a$pseudoIDs2 <- NA
  c2a$classes <- "2C"
  
  #3C
  nativeIDs <- c3$id[unlist(lapply(unique(as.character(c3$seqnames)), function(chr) sample(which(c3$seqnames == chr), size = sum(c3$seqnames == chr), replace = F)))]
  pseudoIDs <- c3$id[unlist(lapply(unique(as.character(c3$seqnames)), function(chr) sample(which(c3$seqnames == chr), size = sum(c3$seqnames == chr), replace = F)))]
  c3a <- cbind(c3, nativeIDs = nativeIDs, pseudoIDs = pseudoIDs)
  setnames(c3a, old = "id", "altIDs")
  c3a$altIDs2 <- NA
  c3a$nativeIDs2 <- NA
  c3a$pseudoIDs2 <- NA
  c3a$classes <- "3C"
  #4C
  if(do4C) {
    #We 
    c4.1 <- merge(c1[,-1], c3[,-1], by = c("seqnames", "train"), allow.cartesian = T)
    print(dim(c4.1))
    print((c4.1[1:3,]))
    c4.2 <- merge(c2[,-1], c2[,-1], by = c("seqnames", "train"), allow.cartesian = T)
    print(dim(c4.2))
    print((c4.2[1:3,]))
    
    c4size <- round(mean(c(nrow(c1), nrow(c3)))/2)
    c4 <- rbindlist(list(c4.1[sample(x = nrow(c4.1), size = c4size),], c4.2[sample(x = nrow(c4.2), size = c4size),]))
    setnames(c4, old = c("id.x", "id.y"), new =  c("altIDs", "altIDs2"))
    c4n <- rbindlist(list(c4.1[sample(x = nrow(c4.1), size = c4size),], c4.2[sample(x = nrow(c4.2), size = c4size),]))
    c4$nativeIDs <- c4n$id.x
    c4$nativeIDs2 <- c4n$id.y
    c4p <- rbindlist(list(c4.1[sample(x = nrow(c4.1), size = c4size),], c4.2[sample(x = nrow(c4.2), size = c4size),]))
    c4$pseudoIDs <- c4p$id.x
    c4$pseudoIDs2 <- c4p$id.y
    c4$classes <- "4C"
    c4$X <- 1:nrow(c4)
    setcolorder(c4, names(c3a))
    
    cs <- rbindlist(list(c1a, c2a, c3a, c4))
  } else {
    cs <- rbindlist(list(c1a, c2a, c3a))
  }
  #data preparation for 4 copies:
  if(do4C)
    c4 <- data.table(X = 1:nrow(c4), seqnames = c4$seqnames, id = c4$altIDs, train = c4$train, altIDs2 = c4$altIDs2)
  
  #0C
  c0 <- cs[sample(nrow(cs), size = totalLossCount),]
  c0$altIDs2 <- NA
  #c4$pseudoIDs <- set up by cs sampling
  c0$nativeIDs2 <- NA
  c0$pseudoIDs2 <- NA
  c0$classes <- "0C"
  c0$nativeIDs <- "FullLoss"
  c0$altIDs <- "FullLoss"  #paste("FullLoss", 1:totalLossCount, sep = "")
  cs <- rbind(cs,c0)
  
  cs$calibrate = F 
  cs$typeOfEventToGenerate <- "amplitude"
  cs$fuzzyRegion <- 0
  
  ########################################
  #TRUE group (includes a point of change)
  ########################################
  #This group has typeOfEventToGenerate = "pointOfChange"
  #nat/alt should be a data.table with:    X seqnames         id train
  get_my_change_of_point_sample_pair <- function(nat, alt, myclass = "1to3", calibrate = F) {
    #returns the header:
    #X seqnames    altIDs train nativeIDs pseudoIDs altIDs2 nativeIDs2 pseudoIDs2 classes calibrate
    g_mix <- merge(nat[,-1], alt[,-1], by = c("seqnames", "train"), allow.cartesian = T)
    print(dim(g_mix))
    print((g_mix[1:3,]))
    setnames(g_mix, old = c("id.x", "id.y"), new =  c("nativeIDs", "altIDs"))
    if(length(alt$altIDs2) == 0) {
      g_mix$altIDs2 <- NA
    } 
    g_mix$nativeIDs2 <- NA
    g_mix$pseudoIDs2 <- NA
    #pseudoID has the same ploidy as the nativeID one:
    g_mix$pseudoIDs <- g_mix$nativeIDs[unlist(lapply(unique(as.character(g_mix$seqnames)), function(chr) sample(which(g_mix$seqnames == chr), size = sum(g_mix$seqnames == chr), replace = F)))]
    g_mix$classes <- myclass
    g_mix$typeOfEventToGenerate <- "pointOfChange"
    g_mix$fuzzyRegion <- fuzzyRegion
    g_mix$calibrate <- calibrate
    g_mix$X <- 1:nrow(g_mix)
    return(g_mix)
  }
  
  
  # F (no change) | T (point of change)
  # -----------------------------------
  # Zero: 0C      | 0 to 1, 0 to 2, 0 to 3, 0 to 4
  # Mono: 1C      | 1 to 2, 1 to 3, 1 to 4
  # diploid: 2C   | 2 to 3, 2 to 4
  # Trisomy: 3C   | 3 to 4
  # Tetrasomy: 4C | -
  tcs <- NULL
  
  #zeros to
  #--------
  #0to1
  tbl <- rbind(data.table(X = 1, seqnames = unique(c1$seqnames), id = "FullLoss", train = T), 
               data.table(X = 1, seqnames = unique(c1$seqnames), id = "FullLoss", train = F))
  tc <- get_my_change_of_point_sample_pair(nat = c1, alt = tbl, myclass = "0to1")
  tcs <- rbind(tcs, tc)
  #0to2
  tbl <- rbind(data.table(X = 1, seqnames = unique(c2$seqnames), id = "FullLoss", train = T), 
               data.table(X = 1, seqnames = unique(c2$seqnames), id = "FullLoss", train = F))
  tc <- get_my_change_of_point_sample_pair(nat = c2, alt = tbl, myclass = "0to2")
  tcs <- rbind(tcs, tc)
  #0to3
  tbl <- rbind(data.table(X = 1, seqnames = unique(c3$seqnames), id = "FullLoss", train = T), 
               data.table(X = 1, seqnames = unique(c3$seqnames), id = "FullLoss", train = F))
  tc <- get_my_change_of_point_sample_pair(nat = c3, alt = tbl, myclass = "0to3")
  tcs <- rbind(tcs, tc)
  #0to4
  if(do4C) {
    tbl <- rbind(data.table(X = 1, seqnames = unique(c4$seqnames), id = "FullLoss", train = T), 
               data.table(X = 1, seqnames = unique(c4$seqnames), id = "FullLoss", train = F))
    tc <- get_my_change_of_point_sample_pair(nat = tbl, alt = c4, myclass = "0to4")
    tcs <- rbind(tcs, tc)
  }
  
  #ones to
  #-------
  #1to2
  tc <- get_my_change_of_point_sample_pair(nat = c1, alt = c2, myclass = "1to2")
  tcs <- rbind(tcs, tc)
  #1to3
  tc <- get_my_change_of_point_sample_pair(nat = c1, alt = c3, myclass = "1to3")
  tcs <- rbind(tcs, tc)
  #1to4
  if(do4C) {
    tc <- get_my_change_of_point_sample_pair(nat = c1, alt = c4, myclass = "1to4")
    tcs <- rbind(tcs, tc)
  }
  #twos to
  #-------
  #2to3
  tc <- get_my_change_of_point_sample_pair(nat = c2, alt = c3, myclass = "2to3")
  tcs <- rbind(tcs, tc)
  #2to4
  if(do4C) {
    tc <- get_my_change_of_point_sample_pair(nat = c2, alt = c4, myclass = "2to4")
    tcs <- rbind(tcs, tc)
  }
  #threes to
  #---------
  #3to4
  if(do4C) {
    tc <- get_my_change_of_point_sample_pair(nat = c3, alt = c4, myclass = "3to4")
    tcs <- rbind(tcs, tc)
  }
  #completing..
  setcolorder(tcs, names(cs))
  cs.all <- rbind(cs, tcs)

  return(cs.all)
}




generateLearningAndTestSamples <- function(adt, minDetectedGenes = 4000, chrToExclude = c("chrX", "chrY", "chrM"), excludeIDs, trainProb = 0.66,
                                           totalLossCount = 20, do4C = T,
                                           file.1C = "/singlecellcenter/etai/data/WORK/secondaryanalysis/methods_paper_results/data/list_of_predicted_monosomies3.csv",
                                           file.2C = "/singlecellcenter/etai/data/WORK/secondaryanalysis/methods_paper_results/data/list_of_predicted_disomies3.csv",
                                           file.3C = "/singlecellcenter/etai/data/WORK/secondaryanalysis/methods_paper_results/data/list_of_predicted_trisomies3.csv") {
  require(data.table)
  
  samp_pool <- names(which(colSums(adt[,-c(1:4)]>0) > minDetectedGenes))
  samp_pool <- samp_pool[!samp_pool %in% excludeIDs]
  
  c1 = data.table(read.csv(file.1C, stringsAsFactors = F), stringsAsFactors = F)
  c1 <- c1[(!seqnames %in% chrToExclude) & (id %in% samp_pool)]
  c2 = data.table(read.csv(file.2C, stringsAsFactors = F), stringsAsFactors = F)
  c2 <- c2[(!seqnames %in% chrToExclude) & (id %in% samp_pool)]
  c3 = data.table(read.csv(file.3C, stringsAsFactors = F), stringsAsFactors = F)
  c3 <- c3[(!seqnames %in% chrToExclude) & (id %in% samp_pool)]
  
  
  
  c1$train <- F
  c1$train[sample(1:nrow(c1), replace = F, size = round(nrow(c1)*trainProb))] <- T
  
  c2$train <- F
  c2$train[sample(1:nrow(c2), replace = F, size = round(nrow(c2)*trainProb))] <- T
  
  c3$train <- F
  c3$train[sample(1:nrow(c3), replace = F, size = round(nrow(c3)*trainProb))] <- T
  
  #1C
  nativeIDs <- c2$id[unlist(lapply(unique(as.character(c1$seqnames)), function(chr) sample(which(c2$seqnames == chr), size = sum(c1$seqnames == chr), replace = F)))]
  pseudoIDs <- c2$id[unlist(lapply(unique(as.character(c1$seqnames)), function(chr) sample(which(c2$seqnames == chr), size = sum(c1$seqnames == chr), replace = F)))]
  c1a <- cbind(c1, nativeIDs = nativeIDs, pseudoIDs = pseudoIDs)
  setnames(c1a, old = "id", "altIDs")
  c1a$altIDs2 <- NA
  c1a$classes <- "1C"
  
  #3C
  nativeIDs <- c2$id[unlist(lapply(unique(as.character(c3$seqnames)), function(chr) sample(which(c2$seqnames == chr), size = sum(c3$seqnames == chr), replace = F)))]
  pseudoIDs <- c2$id[unlist(lapply(unique(as.character(c3$seqnames)), function(chr) sample(which(c2$seqnames == chr), size = sum(c3$seqnames == chr), replace = F)))]
  c3a <- cbind(c3, nativeIDs = nativeIDs, pseudoIDs = pseudoIDs)
  setnames(c3a, old = "id", "altIDs")
  c3a$altIDs2 <- NA
  c3a$classes <- "3C"
  
  #4C
  if(do4C) {
    c4.1 <- merge(c1[,-1], c3[,-1], by = c("seqnames", "train"), allow.cartesian = T)
    print(dim(c4.1))
    print((c4.1[1:3,]))
    c4.2 <- merge(c2[,-1], c2[,-1], by = c("seqnames", "train"), allow.cartesian = T)
    print(dim(c4.2))
    print((c4.2[1:3,]))
    
    c4size <- round(mean(c(nrow(c1), nrow(c3)))/2)
    c4 <- rbindlist(list(c4.1[sample(x = nrow(c4.1), size = c4size),], c4.2[sample(x = nrow(c4.2), size = c4size),]))
    setnames(c4, old = c("id.x", "id.y"), new =  c("altIDs", "altIDs2"))
    c4$nativeIDs <- c2$id[unlist(lapply(unique(as.character(c4$seqnames)), function(chr) sample(which(c2$seqnames == chr), size = sum(c4$seqnames == chr), replace = F)))]
    c4$pseudoIDs <- c2$id[unlist(lapply(unique(as.character(c4$seqnames)), function(chr) sample(which(c2$seqnames == chr), size = sum(c4$seqnames == chr), replace = F)))]
    c4$classes <- "4C"
    c4$X <- 1:nrow(c4)
    setcolorder(c4, names(c3a))
  
    cs <- rbindlist(list(c1a, c3a, c4))
  } else {
    cs <- rbindlist(list(c1a, c3a))
  }
  c0 <- cs[sample(nrow(cs), size = totalLossCount),]
  c0$altIDs2 <- NA
  c0$classes <- "0C"
  c0$altIDs <- "FullLoss"  #paste("FullLoss", 1:totalLossCount, sep = "")

  cs <- rbind(cs,c0)
  
  cs$calibrate = F  
  return(cs)
}



generateTrainingAndTestSet.pointOfChange <- function(adt, controlSampleIDs = controlSampleIDs, 
                                                     cpus = 21, output_prefix_fname = NA, 
                                                     doZscore = F, doScale = T, doRatioBeforeScaleAndZscore = F,
                                                     pseudoInsertAsNegativeControl = T, geneSpecificFtrs = F, useMA = F, useSSM = F, useGlobalFtrs = F, 
                                                     insertSize = 25, extraRegion = 31, ftrsWindowSize = 15, 
                                                     numOfAlteredInserts = 30, numOfPseudoInserts = 10, 
                                                     chrs = c(rep("chr1", 12)), 
                                                     nativeIDs = c("161005c_A5", "170222_B2b", "161228_A1", "170208_A10", "170208_A3", "170208_D1",
                                                                   "161005c_A5", "170222_B2b", "161228_A1", "170208_A10", "0613_B3", "161023_C4"), 
                                                     nativeIDs2 = NULL,
                                                     altIDs = c("171205_3D", "0705_B12", "161003_A2", "180228_1A", "180228_1B", "0627_C8",
                                                                "0627_B1", "160627_B1", "170202_A1", "171205_2G", "180301_5D", "171205_3B"), 
                                                     altIDs2 = NULL,
                                                     pseudoIDs = c("170425_A12", "170425_A7", "170320_B11", "170222_B1b", "170208_C7b", "170222_B1b",
                                                                   "170425_A12", "170425_A7", "170320_B11", "170222_B1b", "170120_A3", "170208_B9"),
                                                     pseudoIDs2 = NULL,
                                                     classes = c(rep("3C", 6),
                                                                 rep("1C", 6)),
                                                     typeOfEventToGenerate = NULL,
                                                     fuzzyRegion = NULL,
                                                     calibrateAlt = NULL, calibratePseudo = NULL,
                                                     nativeClass = "2C",
                                                     trainNativeIdxs = c(3:6,9:12), trainAltIdxs = c(3:6,9:12), trainPseudoIdxs = c(3:6,9:12),
                                                     amplitudeShouldBeNegative = T) {
  
  require(data.table)
  require(randomForest)
  require(doMC)
  require(ROCR)
  require(zoo)
  require(matrixStats)
  
  if(is.null(nativeIDs2)) {
    nativeIDs2 <- c(rep(NA, length(nativeIDs)))
  }
  if(is.null(altIDs2)) {
    altIDs2 <- c(rep(NA, length(altIDs)))
  }
  if(is.null(pseudoIDs2)) {
    pseudoIDs2 <- c(rep(NA, length(pseudoIDs)))
  }
  if(is.null(calibrateAlt)) {
    calibrateAlt <- c(rep(F, length(altIDs)))
  }
  if(is.null(calibratePseudo)) {
    calibratePseudo <- c(rep(F, length(pseudoIDs)))
  }
  
  if(is.null(typeOfEventToGenerate)) {
    typeOfEventToGenerate <- c(rep("amplitute", length(nativeIDs)))
  }
  if(is.null(fuzzyRegion)) {
    fuzzyRegion <- c(rep(0, length(pseudoIDs)))
  }
  
  
  #########################
  #Generating training set:
  #########################
  message("Generating learning set..")
  if(cpus > 1) {
    require(snow)
    cl <- snow::makeCluster(cpus, type = "SOCK")
    clusterEvalQ(cl, library(data.table))
    clusterEvalQ(cl, library(matrixStats))
    clusterEvalQ(cl, library(zoo))
    clusterEvalQ(cl, source('~/WORK/secondaryanalysis/methods_paper_results/R/CNML.R'))
    snow::clusterExport(cl = cl, list = c(names(environment())), envir = environment())
    myLSs <- snow::parLapply(cl, 1:length(trainNativeIdxs), function(i) buildLearningOrTestSet(adt = adt, controlSampleIDs = controlSampleIDs, useGlobalFtrs = useGlobalFtrs, 
                                                                                               native_id = nativeIDs[trainNativeIdxs[i]], native_id2 = nativeIDs2[trainNativeIdxs[i]],  
                                                                                               alt_id = altIDs[trainAltIdxs[i]], alt_id2 = altIDs2[trainAltIdxs[i]], 
                                                                                               pseudo_id = pseudoIDs[trainPseudoIdxs[i]], pseudo_id2 = pseudoIDs2[trainPseudoIdxs[i]], 
                                                                                               typeOfEventToGenerate = typeOfEventToGenerate[trainNativeIdxs[i]],
                                                                                               calibrateAlt = calibrateAlt[trainAltIdxs[i]],
                                                                                               calibratePseudo = calibratePseudo[trainPseudoIdxs[i]],
                                                                                               chr = chrs[trainAltIdxs[i]], label = classes[trainAltIdxs[i]], 
                                                                                               fuzzyRegion = fuzzyRegion[trainAltIdxs[i]],
                                                                                               doRatioBeforeScaleAndZscore = doRatioBeforeScaleAndZscore, geneSpecificFtrs = geneSpecificFtrs,
                                                                                               insertSize = insertSize, extraRegion = extraRegion, ftrsWindowSize = ftrsWindowSize,
                                                                                               numOfAlteredInserts = numOfAlteredInserts, numOfPseudoInserts = numOfPseudoInserts, 
                                                                                               doScale = doScale, useMA = useMA, useSSM = useSSM, 
                                                                                               doZscore = doZscore, amplitudeShouldBeNegative = amplitudeShouldBeNegative))
    
  } else {
    myLSs <- lapply(1:length(trainNativeIdxs), function(i) buildLearningOrTestSet(adt = adt, controlSampleIDs = controlSampleIDs, useGlobalFtrs = useGlobalFtrs,
                                                                                  native_id = nativeIDs[trainNativeIdxs[i]], native_id2 = nativeIDs2[trainNativeIdxs[i]],  
                                                                                  alt_id = altIDs[trainAltIdxs[i]], alt_id2 = altIDs2[trainAltIdxs[i]], 
                                                                                  pseudo_id = pseudoIDs[trainPseudoIdxs[i]], pseudo_id2 = pseudoIDs2[trainPseudoIdxs[i]],  
                                                                                  typeOfEventToGenerate = typeOfEventToGenerate[trainNativeIdxs[i]],
                                                                                  calibrateAlt = calibrateAlt[trainAltIdxs[i]],
                                                                                  calibratePseudo = calibratePseudo[trainPseudoIdxs[i]],
                                                                                  chr = chrs[trainAltIdxs[i]], label = classes[trainAltIdxs[i]], 
                                                                                  fuzzyRegion = fuzzyRegion[trainAltIdxs[i]],
                                                                                  doRatioBeforeScaleAndZscore = doRatioBeforeScaleAndZscore, geneSpecificFtrs = geneSpecificFtrs,
                                                                                  insertSize = insertSize, extraRegion = extraRegion, ftrsWindowSize = ftrsWindowSize,
                                                                                  numOfAlteredInserts = numOfAlteredInserts, numOfPseudoInserts = numOfPseudoInserts, 
                                                                                  doScale = doScale, useMA = useMA, useSSM = useSSM,
                                                                                  doZscore = doZscore, amplitudeShouldBeNegative = amplitudeShouldBeNegative))
  }
  
  ALS <- list(real = rbindlist(lapply(myLSs, function(x) x$alteredInsert)),
              pseudo = rbindlist(lapply(myLSs, function(x) data.table(x$pseudoInsert))))
  
  rm(myLSs)
  
  message("Learning set size = ", nrow(ALS$real))
  print(table(ALS$real$label))
  
  if(!is.na(output_prefix_fname)) {
    saveRDS(ALS, file = sprintf("%s.ALS.rds", output_prefix_fname))
  }
  ########################
  #generating testing set:
  ########################
  message("Generating test set..")
  testNativeIdxs <- (1:length(nativeIDs))[-trainNativeIdxs]
  testAltIdxs <- (1:length(nativeIDs))[-trainAltIdxs]
  testPseudoIdxs <- (1:length(nativeIDs))[-trainPseudoIdxs]
  if(cpus > 1) {
    snow::clusterExport(cl = cl, list = c(names(environment())), envir = environment())
    myTSs <- snow::parLapply(cl, 1:length(testNativeIdxs), function(i) buildLearningOrTestSet(adt = adt, controlSampleIDs = controlSampleIDs, useGlobalFtrs = useGlobalFtrs,
                                                                                              native_id = nativeIDs[testNativeIdxs[i]], native_id2 = nativeIDs2[testNativeIdxs[i]],  
                                                                                              alt_id = altIDs[testAltIdxs[i]], alt_id2 = altIDs2[testAltIdxs[i]], 
                                                                                              pseudo_id = pseudoIDs[testPseudoIdxs[i]], pseudo_id2 = pseudoIDs2[testPseudoIdxs[i]], 
                                                                                              typeOfEventToGenerate = typeOfEventToGenerate[testNativeIdxs[i]],
                                                                                              calibrateAlt = calibrateAlt[testAltIdxs[i]],
                                                                                              calibratePseudo = calibratePseudo[testPseudoIdxs[i]],
                                                                                              chr = chrs[testAltIdxs[i]], label = classes[testNativeIdxs[i]], 
                                                                                              fuzzyRegion = fuzzyRegion[testAltIdxs[i]],
                                                                                              doRatioBeforeScaleAndZscore = doRatioBeforeScaleAndZscore, geneSpecificFtrs = geneSpecificFtrs,
                                                                                              insertSize = insertSize, extraRegion = extraRegion, ftrsWindowSize = ftrsWindowSize,
                                                                                              numOfAlteredInserts = numOfAlteredInserts, numOfPseudoInserts = numOfPseudoInserts, 
                                                                                              doScale = doScale, useMA = useMA, useSSM = useSSM,
                                                                                              doZscore = doZscore, amplitudeShouldBeNegative = amplitudeShouldBeNegative))
    
    stopCluster(cl)
  } else {
    myTSs <- lapply(1:length(testNativeIdxs), function(i) buildLearningOrTestSet(adt = adt, controlSampleIDs = controlSampleIDs, useGlobalFtrs = useGlobalFtrs,
                                                                                 native_id = nativeIDs[testNativeIdxs[i]], native_id2 = nativeIDs2[testNativeIdxs[i]],  
                                                                                 alt_id = altIDs[testAltIdxs[i]], alt_id2 = altIDs2[testAltIdxs[i]], 
                                                                                 pseudo_id = pseudoIDs[testPseudoIdxs[i]], pseudo_id2 = pseudoIDs2[testPseudoIdxs[i]], 
                                                                                 typeOfEventToGenerate = typeOfEventToGenerate[testNativeIdxs[i]],
                                                                                 calibrateAlt = calibrateAlt[testAltIdxs[i]],
                                                                                 calibratePseudo = calibratePseudo[testPseudoIdxs[i]],
                                                                                 chr = chrs[testAltIdxs[i]], label = classes[testNativeIdxs[i]], 
                                                                                 fuzzyRegion = fuzzyRegion[testAltIdxs[i]],
                                                                                 doRatioBeforeScaleAndZscore = doRatioBeforeScaleAndZscore, geneSpecificFtrs = geneSpecificFtrs,
                                                                                 insertSize = insertSize, extraRegion = extraRegion, ftrsWindowSize = ftrsWindowSize,
                                                                                 numOfAlteredInserts = numOfAlteredInserts, numOfPseudoInserts = numOfPseudoInserts, 
                                                                                 doScale = doScale, useMA = useMA, useSSM = useSSM,
                                                                                 doZscore = doZscore, amplitudeShouldBeNegative = amplitudeShouldBeNegative))
  }
  ATS <- list(real = rbindlist(lapply(myTSs, function(x) x$alteredInsert)),
              pseudo = rbindlist(lapply(myTSs, function(x) data.table(x$pseudoInsert))))
  
  message("Learning set size = ", nrow(ATS$real))
  print(table(ATS$real$label))
  
  if(!is.na(output_prefix_fname)) {
    saveRDS(ATS1, file = sprintf("%s.ATS1.rds", output_prefix_fname))
  }
  
  return(list(train = ALS, test = ATS))
}

getExpRatioforML <- function(mydt, doZscore, mySds, doRatioBeforeScaleAndZscore, min_sd = 1, mySds2, myMeans2, 
                             normalizationFactor, doScale = T, zscore.th = 3, ratio.th = 5) {
  
  mydt$ratio <- log2((mydt$y + 1)/normalizationFactor)
  
  if(doZscore && !is.null(mySds)) {
    if(doRatioBeforeScaleAndZscore) {
      mydt$y <- t(scale(t(mydt$ratio), center = myMeans2, scale = mySds2))
    } else {  
      mySds <- ifelse(mySds < min_sd, min_sd, mySds)
      mydt$y <- t(scale(t(mydt$y), center = normalizationFactor, scale = mySds))
    }
  } else {
    mydt$y <- log2((mydt$y + 1)/normalizationFactor)
  }
  if(doScale) {
    mydt$y <- mydt$y - median(mydt$y)
  }
  if(doZscore) {
    th <- zscore.th
  } else {
    th <- ratio.th
  }
  q01 <- max(-th, unname(quantile(mydt$y, 0.01)))
  q99 <- min(th, unname(quantile(mydt$y, 0.99)))
  mydt[y > q99]$y <- q99
  mydt[y < q01]$y <- q01
  
  mydt <- mydt[order(seqnames, start)]
  return(mydt)
}


#the same as "fromRawTPMtoMLInput" but mean on log space and not mean before log transform
fromRawTPMtoMLInput2 <- function(rsemtpm, geneRanges, controlSampleIDs, max_pos_col = 6, plusOne = 1, min_sd = 0.0001, 
                                 zerosAsNA = F, normBySd = F,
                                maxExp = 1200, quantileOfExp = 0.99, minDetectionLevel = 5, minNumOfSamplesToDetect = 3, 
                                geneInterspaceLogBase = 8, interspaceCoef = 5, doNormalizeByControls = F) {
  require(data.table)
  require(matrixStats)
  
  wi <- which(rowSums(rsemtpm[, controlSampleIDs] >= minDetectionLevel) >= minNumOfSamplesToDetect)
  dt <- data.table(rsemtpm[wi,], keep.rownames = "id")
  if(length(names(which(table(colnames(dt))>1))) > 0)
    dt <- dt[, -which(colnames(dt) %in% names(which(table(colnames(dt))>1))), with=F]
  setkey(geneRanges, id)
  setkey(dt, id)
  dt <- merge(geneRanges, dt)
  dt <- dt[order(seqnames, start, end)]
  dt <- dt[seqnames != "chrM" & seqnames != "chrY"]
  
  M <- copy(dt[, controlSampleIDs, with = F])
  if(zerosAsNA) 
    M[M==0] <- NA  
  
  message("NAs = ", sum(is.na(M)))
 
  q <- min(maxExp, unname(quantile(as.matrix(M), quantileOfExp, na.rm=T)))
  M[M > q] <- q
  message("Quantile q = ", q)
  dt2go <- copy(dt[, -c(1:max_pos_col), with=F])
  dt2go[dt2go>q] <- q
  
  #log before mean:
  dt2go <- log2(dt2go + plusOne)

  M <- log2(M + plusOne)
  message("NAs = ", sum(is.na(M)))
  #stats if doRatioBeforeScaleAndZscore == F
  means <- rowMeans(M, na.rm = T)
  sds <- rowSds(as.matrix(M), na.rm = T)
  
  #stats if doRatioBeforeScaleAndZscore == T
  M2 <- sweep(M, 1, means, "-")
  if(normBySd)
    M2 <- sweep(M2, 1, sds, "/") #* sqrt(mean(sds^2, na.rm=T))
  
  means2 <- rowMeans(M2, na.rm = T)
  sds2 <- rowSds(as.matrix(M2), na.rm = T)
  sds2 <- ifelse(sds2 < min_sd, min_sd, sds2)
  
  #gene interspace:
  interspace <- unlist(lapply(unique(dt$seqnames), 
                              function(chr) c(dt[seqnames == chr]$start[-1], dplyr::last(dt[seqnames == chr]$end)+1) - dt[seqnames == chr]$end + geneInterspaceLogBase))
  interspace[interspace <= 0] <- 1
  interspace <- round(log(x = interspace, base = geneInterspaceLogBase), digits = 2)
  interspace <- exp(-interspaceCoef*interspace/max(interspace))
  
  #calculating final dt:
  if(doNormalizeByControls) {
    dt <- cbind(dt[, 1:max_pos_col, with = F], sweep(dt2go, 1, means, "-") )
    if(normBySd)
      dt <- cbind(dt[, 1:max_pos_col, with = F], sweep(dt2go, 1, sds, "/")) #* sqrt(mean(sds^2, na.rm=T)) )
  } else {
    dt <- cbind(dt[, 1:max_pos_col, with = F], dt2go)
  }
  
  res <- list(adt = dt, means = means, means2 = means2, sds = sds, sds2 = sds2, interspace = interspace, controlSampleIDs = controlSampleIDs)
  return(res)
}


#this function is equivalent to "prepareTotalExpData" in runRNAseqPipelineFor_BF9VP.utils.R that is used in preparing the trainig set in trainAndTestLGBM_for_POC_CN1.R
#it outputs organized exp data for ML and statistics for getExpRatioforML
fromRawTPMtoMLInput <- function(rsemtpm, geneRanges, controlSampleIDs, max_pos_col = 6, plusOne = 1, min_sd = 0.0001,
                                maxExp = 1200, quantileOfExp = 0.99, minDetectionLevel = 5, minNumOfSamplesToDetect = 3, 
                                geneInterspaceLogBase = 8, interspaceCoef = 5, doNormalizeByControls = F) {
  require(data.table)
  require(matrixStats)
  
  wi <- which(rowSums(rsemtpm[, controlSampleIDs] >= minDetectionLevel) >= minNumOfSamplesToDetect)
  dt <- data.table(rsemtpm[wi,], keep.rownames = "id")
  if(length(names(which(table(colnames(dt))>1))) > 0)
    dt <- dt[, -which(colnames(dt) %in% names(which(table(colnames(dt))>1))), with=F]
  setkey(geneRanges, id)
  setkey(dt, id)
  dt <- merge(geneRanges, dt)
  dt <- dt[order(seqnames, start, end)]
  dt <- dt[seqnames != "chrM" & seqnames != "chrY"]
  
  M <- copy(dt[, controlSampleIDs, with = F] + plusOne)
  q <- min(maxExp, unname(quantile(as.matrix(M), quantileOfExp)))
  M[M > q] <- q
  dt2go <- copy(dt[, -c(1:max_pos_col), with=F])
  dt2go[dt2go>q] <- q
  
  #stats if doRatioBeforeScaleAndZscore == F
  means <- rowMeans(M)
  sds <- rowSds(as.matrix(M), na.rm = T)
  
  #stats if doRatioBeforeScaleAndZscore == T
  M2 <- log2(sweep(M, 1, rowMeans(as.matrix(M)), "/"))
  means2 <- rowMeans(M2, na.rm = T)
  sds2 <- rowSds(as.matrix(M2), na.rm = T)
  sds2 <- ifelse(sds2 < min_sd, min_sd, sds2)
  
  #gene interspace:
  interspace <- unlist(lapply(unique(dt$seqnames), 
                              function(chr) c(dt[seqnames == chr]$start[-1], dplyr::last(dt[seqnames == chr]$end)+1) - dt[seqnames == chr]$end + geneInterspaceLogBase))
  interspace[interspace <= 0] <- 1
  interspace <- round(log(x = interspace, base = geneInterspaceLogBase), digits = 2)
  interspace <- exp(-interspaceCoef*interspace/max(interspace))
  
  #calculating final dt:
  if(doNormalizeByControls) {
    dt <- cbind(dt[, 1:max_pos_col, with = F], log2(sweep(dt2go + 1, 1, means, "/") ))
  } else {
    dt <- cbind(dt[, 1:max_pos_col, with = F], dt2go)
  }
  
  res <- list(adt = dt, means = means, means2 = means2, sds = sds, sds2 = sds2, interspace = interspace, controlSampleIDs = controlSampleIDs)
  return(res)
}


getMySegmentationLiklihood <- function(mypreds, chr = "chr5") {
  tt <- mypreds[seqnames == chr,]
  rtt <- rle(tt$prediction)
  
  df <- data.frame(from = c(0, cumsum(rtt$lengths)[-length(rtt$lengths)])+1, 
                   to=c(0, cumsum(rtt$lengths)[-length(rtt$lengths)]) + 1 + rtt$lengths - 1, 
                   width = rtt$lengths)
  
  ll <- unlist(lapply(1:nrow(df), function(i) rep(sum(log(1 - tt[df$from[i]:df$to[i],]$vote)), df$width[i])))
  return(ll)
  
}

#TODO: 
#1. write lgbm_predict - DONE
#2.and run it on all samples
#2. run ssm, mv-ssm, ma and imf on all samples (including new ones)
#3. produce reports for stamatis
#4. organize source code in one folder including: ML, TDA (time domain analysis), normalization and pre processing, reports and figures

#input:
#rsemtpm should contain at least the sample under analysis and the controlSampleIDs
#example:
#mypreds <- runMyCheck2(id = "170512_B8", rsemtpm = rsemtpm, )
predict_cn_for_sample <- function(id, rsemtpm, geneRanges, trainingRunID ="id17", controlSampleIDs = controlSampleIDs2, GC_length_rds_fname = NULL,
                                  chrsToExcludeFromNormalization = c("chrX", "chrM", "chr10"), cpus = 1, 
                                  min_sd = 0.0001, doModelOfSelectedFeatures = T, geneInterspaceLogBase = 8, max_pos_col = 4,
                                  pathToTrainingOutput = "/singlecellcenter/etai/tmp/CNML/pointOfChange/test2/", MLtype = "multiclass") {
  require(data.table)
  require(ggplot2)
  
  config_list <- readRDS(sprintf("%s/%s.config_list.rds", pathToTrainingOutput, trainingRunID))
  ftrsWindowSize <- config_list$ftrsWindowSize
  training_insertSize <- config_list$insertSize
  
  maxExp = config_list$maxExp
  quantileOfExp = config_list$quantileOfExp
  minDetectionLevel = config_list$minDetectionLevel
  minNumOfSamplesToDetect = config_list$minNumOfSamplesToDetect
  useGeneSpecificFeatures <- config_list$geneSpecificFtrs
  useSSM <- config_list$useSSM
  useMA <- config_list$useMA
  
  mlinput <- fromRawTPMtoMLInput(rsemtpm = rsemtpm, geneRanges = geneRanges, controlSampleIDs, max_pos_col = max_pos_col, plusOne = 1, min_sd = min_sd,
                                 maxExp = maxExp, quantileOfExp = quantileOfExp, minDetectionLevel = minDetectionLevel, minNumOfSamplesToDetect = minNumOfSamplesToDetect, 
                                 geneInterspaceLogBase = geneInterspaceLogBase, doNormalizeByControls = F)
  
  adt <- mlinput$adt
  if(useGeneSpecificFeatures) {
    geneSpecific <- generateGeneSpecificFtrs(adt = adt, myMeans = mlinput$means, mySds = mlinput$sds, base = 2, digits = 2, GC_length_rds_fname = GC_length_rds_fname)
    myTCV <- geneSpecific$TCV
  } else {
    geneSpecific = NA
  }
  
  if(MLtype == "multiclass") {
    typeOfPrediction <- "classification"
  } else {
    typeOfPrediction <- "regression"
  }
  
  
  if(doModelOfSelectedFeatures) {
    bst.fname <- sprintf("%s/%s.bst.ftrs_selected_by_0.0050_th.%s.rds", pathToTrainingOutput, trainingRunID, MLtype)
    selected.ftrs <- readRDS(sprintf("%s/%s.names_of_selected_features_by_0.0050_th.%s.rds", pathToTrainingOutput, trainingRunID, MLtype))
    
    mydt <- data.table(adt[,1:max_pos_col], y = adt[[id]]) #clear C1
    mypreds <- lgbm_predict(mydt = mydt, bst.fname = bst.fname, chr = chr, cpus = cpus,
                            useMA = F, useSSM = useSSM, ftrsWindowSize = ftrsWindowSize, doZscore = T, featureNames = selected.ftrs,
                            doRatioBeforeScaleAndZscore = T, 
                            myMeans2 = mlinput$means2, mySds2 = mlinput$sds2,
                            normalizationFactor = mlinput$means, mySds = mlinput$sds, 
                            geneSpecific = geneSpecific, 
                            training_insertSize = training_insertSize, typeOfPrediction = typeOfPrediction)
    
  } else { 
    
    bst.fname <- sprintf("%s/%s.bst.all_ftrs.%s.rds", pathToTrainingOutput, trainingRunID, MLtype)
    mydt <- data.table(adt[,1:max_pos_col], y = adt[[id]]) #clear C1
    
    mypreds <- lgbm_predict(mydt = mydt, bst.fname = bst.fname, chr = chr,
                            useMA = F, useSSM = useSSM, ftrsWindowSize = ftrsWindowSize, doZscore = T,  
                            doRatioBeforeScaleAndZscore = T, 
                            myMeans2 = mlinput$means2, mySds2 = mlinput$sds2,
                            normalizationFactor = mlinput$means, mySds = mlinput$sds,
                            geneSpecific = geneSpecific, 
                            training_insertSize = training_insertSize, typeOfPrediction = typeOfPrediction)
  }  
  return(mypreds)
}

#bst.fname overrides bst if not null
lgbm_predict <- function(mydt, bst, bst.fname = NULL, normalizationFactor, mySds = NULL, chr = "chr5", useMA = F, useSSM = F,
                         doPlot = T, cpus = 22, shouldNormalizeExpData = T, max_pos_col = 4,
                         featureNames = NULL, doZscore = F,  imf_level = 6, chrsToExcludeFromNormalization = c("chrX", "chrM", "chr10"),
                         ftrsWindowSize = 31, training_insertSize = 35, global = NA, 
                         geneSpecific = NA, plotAllGenome = F, MAWinSize = 100, allFeaturesPrecisionNumOfDigits = 3,
                         doRatioBeforeScaleAndZscore = T, mySds2 = NULL, myMeans2 = NULL, min_sd = 0.0001, MA_winSize = 50,
                         typeOfPrediction = "classification") {
  
  require(lightgbm)
  require(data.table)
  require(zoo)
  
  if(!is.null(bst.fname))
    bst <- lgb.load(filename = bst.fname)
  
  tbl <- table(mydt$seqnames)
  tbl <- tbl[tbl!=0]
  minSize <- min(tbl)
  if(MAWinSize > minSize) {
    message("MA win size is greater than smallest chr num of genes")
    print(tbl)
    MAWinSize <- minSize - 1
  }
  
  if(shouldNormalizeExpData)
    mydt <- getExpRatioforML(mydt = mydt, doZscore = doZscore, mySds = mySds, 
                             doRatioBeforeScaleAndZscore = doRatioBeforeScaleAndZscore, 
                             min_sd = min_sd, mySds2 = mySds2, myMeans2 = myMeans2, 
                             normalizationFactor = normalizationFactor)
  chrs <- unique(mydt$seqnames)
  
  
  if(useMA) {
    ma <- mydt[, lapply(.SD, rollmean, k = MA_winSize, fill = "extend", align = "center"),
               .SDcols = names(mydt)[-c(1:max_pos_col)], by = seqnames]
    #cell normalization:
    dtm <- ma[, lapply(.SD, median),
              .SDcols = c("y"), by = "seqnames"]
    ma$y <- ma$y - median(dtm[!(seqnames %in% chrsToExcludeFromNormalization)]$y)
    
  }
  message("Running local SSM..")
  if(useSSM) {
    #generate ssm for all chromosomes in the cell for normalization:
    ssm.all <- mydt[, lapply(.SD, getMyGSSM), .SDcols = names(mydt)[-c(1:4)], by = seqnames]
    dtm <- ssm.all[, lapply(.SD, median),
                   .SDcols = c("y"), by = "seqnames"]
    ssmNorm <- median(dtm[!(seqnames %in% chrsToExcludeFromNormalization)]$y)
    if(is.null(cpus) || is.na(cpus) || cpus < 2) {
      localSSM <- lapply(chrs, 
                         function(chr) rollapply(data = mydt[seqnames==chr]$y, width=training_insertSize, FUN = getMyGSSM, fill = "extend", align = "center"))
    } else {
      require(snow)
      cl <- snow::makeCluster(cpus, type = "SOCK")
      clusterEvalQ(cl, library(data.table))
      clusterEvalQ(cl, library(zoo))
      clusterEvalQ(cl, source('~/WORK/secondaryanalysis/methods_paper_results/R/CNML.R'))
      snow::clusterExport(cl = cl, list = c("getMyGSSM", "mydt", "training_insertSize"), envir = environment())
      localSSM <- snow::parLapply(cl, chrs, 
                                  function(chr) rollapply(data = mydt[seqnames==chr]$y, width=training_insertSize, FUN = getMyGSSM, fill = "extend", align = "center"))
      stopCluster(cl)
      
    }
    localSSM <- lapply(localSSM, function(x) x - ssmNorm)
    names(localSSM) <- as.character(unique(mydt$seqnames))
  } 
  message("Local SSM calculations completed.")
  
  res <- lapply(chrs, function(chr) {
    if(!is.na(geneSpecific[1])) {
      gs <- list()
      for(gn in names(geneSpecific)) {
        gs[[gn]] <- geneSpecific[[gn]][mydt$seqnames == chr]
      }
    } else {
      gs <- NA
    }
    
    Ys <- list()
    inputYasMatrix <- NULL
    Ys$tpm <- mydt[seqnames == chr]$y
    if(useMA)
      Ys$MA <- ma[seqnames == chr]$y
    if(useSSM) {
      Ys$SSM <- localSSM[[chr]]
      inputYasMatrix <- list(SSM = T)
    }
    
    
    if(is.na(global[1])) {
      tt <- generateFeaturesForExpressionVector(Ys = Ys, inputYasMatrix = inputYasMatrix, ftrsWindowSize = ftrsWindowSize, 
                                                geneSpecific = gs)
    } else {
      tt <- generateFeaturesForExpressionVector(Ys = Ys, inputYasMatrix = inputYasMatrix, ftrsWindowSize = ftrsWindowSize, 
                                                global = mydt[seqnames == chr]$y, 
                                                geneSpecific = gs)
    }
    return(tt)
    
  })
  names(res) <- chrs
  #return(res)
  res2 <- list()
  idxs <- lapply(chrs, function(chr) {
    idx <- res[[chr]]$idxs
    return(idx)
    #tt <- generateFeaturesForExpressionVector(Y = mydt[seqnames == chr]$y, global = mydt[seqnames == chr]$y, ftrsWindowSize = ftrsWindowSize)
    #return(as.numeric(tt[[1]]))
  })
  names(idxs) <- chrs
  
  for(n in names(res))
    res[[n]] <- round(res[[n]][,-1, with=F], digits = allFeaturesPrecisionNumOfDigits)
  
  message("Performing predictions..")
  #classification
  if(typeOfPrediction == "classification") {
    maxcs <- lapply(chrs, function(chr) {
      dm <- data.matrix(res[[chr]])
      if(!is.null(featureNames)) {
        dm <- dm[, featureNames]
        colnames(dm) <- featureNames
      }
      preds <- predict(bst, dm, reshape=T)
      y <- max.col(preds)
      p <- sapply(1:nrow(preds), function(i) preds[i, y[i]])
      return(data.frame(y = y - 1, p = p))
    })
    names(maxcs) <- chrs
    
    
    #return(maxcs)
    myres <- lapply(chrs, function(chr) {
                data.table(id = mydt[seqnames == chr][idxs[[chr]]]$id,
                           seqnames = mydt[seqnames == chr]$seqnames[idxs[[chr]]], 
                           start = mydt[seqnames == chr]$start[idxs[[chr]]],
                           end = mydt[seqnames == chr]$end[idxs[[chr]]], 
                           tpm = mydt[seqnames == chr][idxs[[chr]]]$y,
                           prediction = maxcs[[chr]]$y,
                           vote = maxcs[[chr]]$p) })
              
    
    mypreds <- rbindlist(myres)
    mypreds$ll <- unlist(lapply(unique(mypreds$seqnames), function(chr) getMySegmentationLiklihood(mypreds = mypreds, chr = chr)))
  } else {
    maxcs <- lapply(chrs, function(chr) {
      dm <- data.matrix(res[[chr]])
      if(!is.null(featureNames)) {
        dm <- dm[, featureNames]
        colnames(dm) <- featureNames
      }
      preds <- predict(bst, dm, reshape=T)
      y <- preds
      return(data.frame(y = y))
    })
    names(maxcs) <- chrs
    
    
    #return(maxcs)
    myres <- lapply(chrs, function(chr) {
      data.table(id = mydt[seqnames == chr][idxs[[chr]]]$id,
                 seqnames = mydt[seqnames == chr]$seqnames[idxs[[chr]]], 
                 start = mydt[seqnames == chr]$start[idxs[[chr]]],
                 end = mydt[seqnames == chr]$end[idxs[[chr]]], 
                 tpm = mydt[seqnames == chr][idxs[[chr]]]$y,
                 prediction = maxcs[[chr]]$y) })
    
    
    mypreds <- rbindlist(myres)
  }
  
  return(mypreds)
}  


generateAndEstimatePredictionMultiClassModel <- function(adt, controlSampleIDs = controlSampleIDs, splitTrainingSet = 0, cpus = 21,
                                                         output_prefix_fname = NA, noLearning = F, geneSpecificFtrs = F,
                                                         doScale = T, pseudoInsertAsNegativeControl = T, runInParallel = T, 
                                                         insertSize = 25, extraRegion = 31, ftrsWindowSize = 15, 
                                                         numOfAlteredInserts = 30, numOfPseudoInserts = 10, useMA = F, 
                                                         chrs = c(rep("chr1", 12)), useGlobalFtrs = F, doRatioBeforeScaleAndZscore = F,
                                                         nativeIDs = c("161005c_A5", "170222_B2b", "161228_A1", "170208_A10", "170208_A3", "170208_D1",
                                                                       "161005c_A5", "170222_B2b", "161228_A1", "170208_A10", "0613_B3", "161023_C4"), 
                                                         nativeIDs2 = NULL,
                                                         altIDs = c("171205_3D", "0705_B12", "161003_A2", "180228_1A", "180228_1B", "0627_C8",
                                                                    "0627_B1", "160627_B1", "170202_A1", "171205_2G", "180301_5D", "171205_3B"), 
                                                         altIDs2 = NULL,
                                                         pseudoIDs = c("170425_A12", "170425_A7", "170320_B11", "170222_B1b", "170208_C7b", "170222_B1b",
                                                                       "170425_A12", "170425_A7", "170320_B11", "170222_B1b", "170120_A3", "170208_B9"),
                                                         pseudoIDs2 = NULL,
                                                         classes = c(rep("3C", 6),
                                                                     rep("1C", 6)),
                                                         doZscore = F, do4C = T,
                                                         calibrateAlt = c(rep(TRUE, 12)), calibratePseudo = c(rep(TRUE, 12)),
                                                         nativeClass = "2C",
                                                         trainNativeIdxs = c(3:6,9:12), trainAltIdxs = c(3:6,9:12), trainPseudoIdxs = c(3:6,9:12),
                                                         ntree = 300, ncores = 18, amplitudeShouldBeNegative = T) {
  
  require(data.table)
  require(randomForest)
  require(doMC)
  require(ROCR)
  require(zoo)
  require(matrixStats)
  
  if(is.null(nativeIDs2)) {
    nativeIDs2 <- c(rep(NA, length(nativeIDs)))
  }
  if(is.null(altIDs2)) {
    altIDs2 <- c(rep(NA, length(altIDs)))
  }
  if(is.null(pseudoIDs2)) {
    pseudoIDs2 <- c(rep(NA, length(pseudoIDs)))
  }
  
  #########################
  #Generating training set:
  #########################
  message("Generating learning set..")
  if(cpus > 1) {
    require(snow)
    cl <- snow::makeCluster(cpus, type = "SOCK")
    clusterEvalQ(cl, library(data.table))
    clusterEvalQ(cl, library(matrixStats))
    clusterEvalQ(cl, library(zoo))
    clusterEvalQ(cl, source('~/WORK/secondaryanalysis/methods_paper_results/R/CNML.R'))
    snow::clusterExport(cl = cl, list = c(names(environment())), envir = environment())
    myLSs <- snow::parLapply(cl, 1:length(trainNativeIdxs), function(i) buildLearningOrTestSet(adt = adt, controlSampleIDs = controlSampleIDs, useGlobalFtrs = useGlobalFtrs,
                                                                                  native_id = nativeIDs[trainNativeIdxs[i]], native_id2 = nativeIDs2[trainNativeIdxs[i]],  
                                                                                  alt_id = altIDs[trainAltIdxs[i]], alt_id2 = altIDs2[trainAltIdxs[i]], 
                                                                                  pseudo_id = pseudoIDs[trainPseudoIdxs[i]], pseudo_id2 = pseudoIDs2[trainPseudoIdxs[i]],  
                                                                                  doRatioBeforeScaleAndZscore = doRatioBeforeScaleAndZscore,
                                                                                  calibrateAlt = calibrateAlt[trainAltIdxs[i]],
                                                                                  calibratePseudo = calibratePseudo[trainPseudoIdxs[i]],
                                                                                  chr = chrs[trainAltIdxs[i]], label = classes[trainAltIdxs[i]], geneSpecificFtrs = geneSpecificFtrs,
                                                                                  insertSize = insertSize, extraRegion = extraRegion, ftrsWindowSize = ftrsWindowSize,
                                                                                  numOfAlteredInserts = numOfAlteredInserts, numOfPseudoInserts = numOfPseudoInserts, 
                                                                                  doScale = doScale, useMA = useMA, doZscore = doZscore, amplitudeShouldBeNegative = amplitudeShouldBeNegative))
    
  } else {
    myLSs <- lapply(1:length(trainNativeIdxs), function(i) buildLearningOrTestSet(adt = adt, controlSampleIDs = controlSampleIDs, useGlobalFtrs = useGlobalFtrs,
                                                                                  native_id = nativeIDs[trainNativeIdxs[i]], native_id2 = nativeIDs2[trainNativeIdxs[i]],  
                                                                                  alt_id = altIDs[trainAltIdxs[i]], alt_id2 = altIDs2[trainAltIdxs[i]], 
                                                                                  pseudo_id = pseudoIDs[trainPseudoIdxs[i]], pseudo_id2 = pseudoIDs2[trainPseudoIdxs[i]],  
                                                                                  doRatioBeforeScaleAndZscore = doRatioBeforeScaleAndZscore,
                                                                                  calibrateAlt = calibrateAlt[trainAltIdxs[i]],
                                                                                  calibratePseudo = calibratePseudo[trainPseudoIdxs[i]],
                                                                                  chr = chrs[trainAltIdxs[i]], label = classes[trainAltIdxs[i]], geneSpecificFtrs = geneSpecificFtrs,
                                                                                  insertSize = insertSize, extraRegion = extraRegion, ftrsWindowSize = ftrsWindowSize,
                                                                                  numOfAlteredInserts = numOfAlteredInserts, numOfPseudoInserts = numOfPseudoInserts, 
                                                                                  doScale = doScale, useMA = useMA, doZscore = doZscore, amplitudeShouldBeNegative = amplitudeShouldBeNegative))
  }
  if(pseudoInsertAsNegativeControl) {
    ALS <- rbindlist(lapply(myLSs, function(myLS) {
      myLS$pseudoInsert$label <- "0"
      learningSet <- rbindlist(myLS)
      learningSet$label[learningSet$label == "0"] <- nativeClass
      learningSet$label <- factor(learningSet$label, levels = sort(as.character(unique(learningSet$label))))
      return(learningSet)
    }))
  } else {
    message("Not yet supported.")
    return(NA)
  }
  
  rm(myLSs)
  
  message("Learning set size = ", nrow(ALS))
  print(table(ALS$label))
  
  if(!is.na(output_prefix_fname)) {
    saveRDS(ALS, file = sprintf("%s.ALS.rds", output_prefix_fname))
  }
  ########################
  #generating testing set:
  ########################
  message("Generating test set..")
  testNativeIdxs <- (1:length(nativeIDs))[-trainNativeIdxs]
  testAltIdxs <- (1:length(nativeIDs))[-trainAltIdxs]
  testPseudoIdxs <- (1:length(nativeIDs))[-trainPseudoIdxs]
  if(cpus > 1) {
    snow::clusterExport(cl = cl, list = c(names(environment())), envir = environment())
    myTSs <- snow::parLapply(cl, 1:length(testNativeIdxs), function(i) buildLearningOrTestSet(adt = adt, controlSampleIDs = controlSampleIDs,
                                                                                 native_id = nativeIDs[testNativeIdxs[i]], native_id2 = nativeIDs2[testNativeIdxs[i]],  
                                                                                 alt_id = altIDs[testAltIdxs[i]], alt_id2 = altIDs2[testAltIdxs[i]], 
                                                                                 pseudo_id = pseudoIDs[testPseudoIdxs[i]], pseudo_id2 = pseudoIDs2[testPseudoIdxs[i]], 
                                                                                 doRatioBeforeScaleAndZscore = doRatioBeforeScaleAndZscore,
                                                                                 calibrateAlt = calibrateAlt[testAltIdxs[i]],
                                                                                 calibratePseudo = calibratePseudo[testPseudoIdxs[i]],
                                                                                 chr = chrs[testAltIdxs[i]], label = classes[testNativeIdxs[i]], geneSpecificFtrs = geneSpecificFtrs,
                                                                                 insertSize = insertSize, extraRegion = extraRegion, ftrsWindowSize = ftrsWindowSize,
                                                                                 numOfAlteredInserts = numOfAlteredInserts, numOfPseudoInserts = numOfPseudoInserts, 
                                                                                 doScale = doScale, useMA = useMA, doZscore = doZscore, amplitudeShouldBeNegative = amplitudeShouldBeNegative))
    
    stopCluster(cl)
  } else {
    myTSs <- lapply(1:length(testNativeIdxs), function(i) buildLearningOrTestSet(adt = adt, controlSampleIDs = controlSampleIDs,
                                                                                 native_id = nativeIDs[testNativeIdxs[i]], native_id2 = nativeIDs2[testNativeIdxs[i]],  
                                                                                 alt_id = altIDs[testAltIdxs[i]], alt_id2 = altIDs2[testAltIdxs[i]], 
                                                                                 pseudo_id = pseudoIDs[testPseudoIdxs[i]], pseudo_id2 = pseudoIDs2[testPseudoIdxs[i]], 
                                                                                 doRatioBeforeScaleAndZscore = doRatioBeforeScaleAndZscore,
                                                                                 calibrateAlt = calibrateAlt[testAltIdxs[i]],
                                                                                 calibratePseudo = calibratePseudo[testPseudoIdxs[i]],
                                                                                 chr = chrs[testAltIdxs[i]], label = classes[testNativeIdxs[i]], geneSpecificFtrs = geneSpecificFtrs,
                                                                                 insertSize = insertSize, extraRegion = extraRegion, ftrsWindowSize = ftrsWindowSize,
                                                                                 numOfAlteredInserts = numOfAlteredInserts, numOfPseudoInserts = numOfPseudoInserts, 
                                                                                 doScale = doScale, useMA = useMA, doZscore = doZscore, amplitudeShouldBeNegative = amplitudeShouldBeNegative))
  }
  if(pseudoInsertAsNegativeControl) {
    ATS1 <- rbindlist(lapply(myTSs, function(myTS) {
      myTS$pseudoInsert$label <- nativeClass
      testSet <- rbindlist(myTS)
      testSet$label[testSet$label == "0"] <- nativeClass
      testSet$label <- factor(testSet$label, levels = sort(as.character(unique(testSet$label)))) 
      return(testSet)
    }))
  } else {
    message("Not yet supported.")
    return(NA)
  }
  
  if(!is.na(output_prefix_fname)) {
    saveRDS(ATS1, file = sprintf("%s.ATS1.rds", output_prefix_fname))
  }
  
  ATS2 <- rbindlist(lapply(myTSs, function(myTS) {
    testSet <- myTS$alteredInsert
    testSet$label[testSet$label == "0"] <- nativeClass
    testSet$label <- factor(testSet$label, levels = sort(as.character(unique(testSet$label)))) 
    return(testSet)
  }))
  
  rm(myTSs)
  message("Testing set 1 size = ", nrow(ATS1))
  print(table(ATS1$label))
  message("Testing set 2 size = ", nrow(ATS2))
  print(table(ATS2$label))
  
  if(!is.na(output_prefix_fname)) {
    saveRDS(ATS2, file = sprintf("%s.ATS2.rds", output_prefix_fname))
  }
  
  if(noLearning == T) {
    message("Completed (no learning requested).")
    return(T)
  }
  
  #return(as.list(environment()))
  
  
  
  ######################################################
  #training using regular Random Forest - Can be slow!!:
  #####################################################
  #TODO: exclude this code for learning from this function to an outside function like lgbm
  message("Training..")
  if(runInParallel) {
    if(splitTrainingSet > 0) {
      Nsplits <- sample(1:splitTrainingSet, nrow(ALS), replace = T)
      print(table(Nsplits))
      registerDoMC()
      rfs <- list()
      message("Spliting votes to ", splitTrainingSet, " datasets..")
      for(i in 1:splitTrainingSet) {
        message("Doing split: ", i, " size = ", sum(Nsplits == i))
        system.time(rfs[[i]] <- foreach(ntree=rep(ntree, ncores), .combine=combine, .multicombine=TRUE,
                                  .packages='randomForest') %dopar% {
                                    randomForest(label ~ . , 
                                                 data = ALS[Nsplits == i, ][sample(sum(Nsplits == i), size = min(table(Nsplits))), ], 
                                                 proximity = F, importance = F, norm.votes = F)
                                  })
        if(!is.na(output_prefix_fname)) {
          saveRDS(rfs[[i]], file = sprintf("%s.rf%d.rds", output_prefix_fname, i))
        }
        gc()
      }
      message("Combining rfs..")
      rf <- do.call(combine, rfs)
      saveRDS(rf, file = sprintf("%s.rf.all.rds", output_prefix_fname))
    } else {
      registerDoMC()
      system.time(rf <- foreach(ntree=rep(ntree, ncores), .combine=combine, .multicombine=TRUE,
                                .packages='randomForest') %dopar% {
                                  randomForest(label ~ . , data = ALS , proximity = F, importance = F)
                                })
    }
    
  } else {
    system.time(rf <- randomForest(label ~ . , ntree = ntree, data = ALS , ntreeproximity = F, importance = T))
  }
  
  #######################
  #Performacne evaluation
  #######################
  message("Performance evaluation..")
  
  
  #On learning set:
  if(splitTrainingSet == 0)
    evaluateMultiClassPerformance(votes = rf$votes, trueLabels = ALS$label, main = "Learning set")
  
  message("Predicting on test set 1..")
  #On test set 1 (with psuedoAlts):
  
  if(splitTrainingSet > 0) {
    Nsplits <- sample(1:splitTrainingSet, nrow(ATS1), replace = T)
    print(table(Nsplits))
    preds1 <- list()
    for(i in 1:splitTrainingSet) {
      message("Doing split: ", i, " size = ", sum(Nsplits == i))
      preds1[[i]] <- data.table(predict(rf, ATS1[Nsplits == i, ], type = "prob"))
      preds1[[i]]$label <- ATS1[Nsplits == i, ]$label
    }
    pred1 <- rbindlist(preds1)
    saveRDS(pred1, file = sprintf("%s.pred1.rds", output_prefix_fname))
  } else {
    pred1 <- data.table(predict(rf, ATS1, type = "prob"))
    pred1$label <- ATS1$label
  }
  evaluateMultiClassPerformance(votes = pred1[ ,-which(names(pred1) %in% "label"), with = F ], trueLabels = pred1$label, main = "Test set 1 (including pseudo in test)")
  
  #On test set 2 (with psuedoAlts):
  message("Predicting on test set 2..")
  
  
  if(splitTrainingSet > 0) {
    Nsplits <- sample(1:splitTrainingSet, nrow(ATS2), replace = T)
    print(table(Nsplits))
    preds2 <- list()
    for(i in 1:splitTrainingSet) {
      message("Doing split: ", i, " size = ", sum(Nsplits == i))
      preds2[[i]] <- data.table(predict(rf, ATS2[Nsplits == i, ], type = "prob"))
      preds2[[i]]$label <- ATS2[Nsplits == i, ]$label
    }
    pred2 <- rbindlist(preds2)
    saveRDS(pred2, file = sprintf("%s.pred2.rds", output_prefix_fname))
  } else {
    pred2 <- data.table(predict(rf, ATS2, type = "prob"))
    pred2$label <- ATS2$label
  }
  evaluateMultiClassPerformance(votes = pred2[ ,-which(names(pred2) %in% "label"), with = F ], trueLabels = pred2$label, main = "Test set 2 (no pseudo in test)")
  
  #pred2 <- predict(rf, ATS2, type = "prob")
  #evaluateMultiClassPerformance(votes = pred2, trueLabels = ATS2$label, main = "Test set 2 (no pseudo in test)")
  
  message("Completed.")
  return(as.list(environment()))
  
}

evaluateMultiClassPerformance <- function(votes, trueLabels, ROCplot = F, main = "") {
  require(ROCR)
  votes <- data.table(votes)
  classes <- unique(colnames(votes))
  trueLabels <- as.character(trueLabels)
  i <- 0
  rocs <- list()
  prs <- list()
  
  aucs.roc <- NULL
  aucs.pr <- NULL
  for(label in classes) {
    i <- i + 1
    myLabels <- as.factor((trueLabels == label) + 0)
    mypreds <- data.table(preds = votes[, label, with=F], labels = myLabels)
    
    pred <- prediction(mypreds$preds, labels = mypreds$label)
    AUC <- performance(pred,"auc")@y.values[[1]]
    aucs.roc <- c(aucs.roc, AUC)
    aucs.pr <- c(aucs.pr, pr.curve(scores.class0 = as.vector(as.matrix(mypreds[,1])), weights.class0 = as.numeric(mypreds$labels)-1)$auc.integral)
    
    prs[[label]] <- pr.curve(scores.class0 = as.vector(as.matrix(mypreds[,1])), weights.class0 = as.numeric(mypreds$labels)-1, curve = T)
    
    if(ROCplot) {
      perf_ROC <- performance(pred,"tpr","fpr") #plot the actual ROC curve
      plot(perf_ROC, main = sprintf("ROC %s", main), col=rainbow(length(classes))[i], lwd=3, add = ifelse(i == 1, F, T), 
           cex = 2, cex.lab = 1.4, cex.axis = 1.6)
    #text(0.5,0.5,paste("ROC AUC = ",format(AUC, digits=5, scientific=FALSE)))
    } else {
      perf_PR <- performance(pred,"prec","rec") #plot the actual ROC curve
      plot(perf_PR, main = sprintf("PR %s", main), col=rainbow(length(classes))[i], lwd=3, add = ifelse(i == 1, F, T), ylim = c(0,1), xlim = c(0,1), 
           cex = 2, cex.lab = 1.4, cex.axis = 1.6)
    
      #text(0.5,0.5,paste("ROC AUC = ",format(AUC, digits=5, scientific=FALSE)))
    }
  }
  
  allaucs <- list(aucs.roc = aucs.roc, aucs.pr = aucs.pr, pr = prs)
  if(ROCplot) {
    plot_aucs <- aucs.roc
    legend("bottomleft", inset=.05, title=NULL, cex = 1.5,
           legend = sprintf("%s - AUC = %.2f", classes, plot_aucs),
           fill=rainbow(length(classes)), horiz=F)
  } else {
    plot_aucs <- aucs.pr
    legend("bottomleft", inset=.05, title=NULL, cex = 1.5,
           legend = sprintf("%s - AUC = %.2f", classes, plot_aucs),
           fill=rainbow(length(classes)), horiz=F)
  }
  
  
  return(allaucs)
}


getFeatureImprotanceTbl <- function(mybst = "/singlecellcenter/etai/tmp/CNML/extraRegion32_rnd4-31_insert_31_rnd3-30_win31_gs_sc_zs_br2.0_v14/extraRegion32_rnd4-31_insert_31_rnd3-30_win31_gs_sc_zs_br2.0_v14.bst.all_ftrs.rds") {
  require(data.table)
  require(lightgbm)
  require(dplyr)
  
  bst <- lgb.load(filename = mybst)
  
  # multilogloss of final iteration on test set
  paste("# Rounds:", bst$current_iter())
  paste("Multilogloss of best model:", bst$record_evals$test$multi_logloss$eval %>% unlist() %>% tail(1))
  
  #Calculate variable importance. (note: takes awhile since single-threaded)
  df_imp <- tbl_df(lgb.importance(bst, percentage = TRUE))
  
  #df_imp$Feature <- factor(df_imp$Feature, levels=rev(df_imp$Feature))
  return(data.table(df_imp))
}
train_using_lightbgm <- function(trainingSet.fname = "/singlecellcenter/etai/tmp/CNML/learn_extraRegion32_insert_31_win31_gs_zs_v5.ALS.rds",
                                 testSet.fname = "/singlecellcenter/etai/tmp/CNML/learn_extraRegion32_insert_31_win31_gs_zs_v5.ATS1.rds",
                                 bst.fname = "/singlecellcenter/etai/tmp/CNML/learn_extraRegion32_insert_31_win31_gs_zs_v5.bst.all_ftrs.rds",
                                 report.fname = "/singlecellcenter/etai/tmp/CNML/learn_extraRegion32_insert_31_win31_gs_zs_v5.trainTestReport.pdf",
                                 selected.features.names = NA, allFeaturesPrecisionNumOfDigits = 3,
                                 fraction_of_total_trainSet = 1.0, fraction_of_total_testSet = 1.0, balance_ratio = NA,
                                 bagging_fraction = 1.0, bagging_freq = 0, nrounds = 1200,
                                 num_threads = 24, learning_rate = 0.1, early_stopping_rounds = 12) {
  require(data.table)
  require(lightgbm)
  require(dplyr)
  require(gtools)
  require(ggplot2)
  
  source('~ejacob/WORK/secondaryanalysis/methods_paper_results/R/lgbm.multiclass_example.Rstart.R')
  
  #load full region in CNV:
  message("Having ", length(trainingSet.fname), " training files to read..")
  if(length(trainingSet.fname) > 1) {
    df.train <- rbindlist(lapply(trainingSet.fname, function(f) readRDS(file = f)))
  } else {
    df.train <- readRDS(file = trainingSet.fname)
  }
  print(dim(df.train))
  
  all.levels <- mixedsort(as.character(unique(df.train$label)))
  df.train <- data.frame(round(df.train[,-ncol(df.train), with=F], digits = allFeaturesPrecisionNumOfDigits), label = factor(df.train$label, all.levels))
  
  #load full region in CNV:
  message("Having ", length(testSet.fname), " test files to read..")
  if(length(testSet.fname) > 1) {
    df.test <- rbindlist(lapply(testSet.fname, function(f) readRDS(file = f)))
  } else {
    df.test <- readRDS(file = testSet.fname)
  }
  print(dim(df.test))
  
  df.test <- data.frame(round(df.test[,-ncol(df.test), with=F], digits = allFeaturesPrecisionNumOfDigits), label = factor(df.test$label, all.levels))
  
  #adjust labels to lgbm foramt:
  df.train$label <- as.numeric(factor(df.train$label, levels = all.levels)) - 1
  df.test$label <- as.numeric(factor(df.test$label, levels = all.levels)) - 1
  
  #The data must be randomized for lightgbm to give unbiased scores. 
  df.train <- df.train %>% sample_frac(size = fraction_of_total_trainSet)
  df.test <- df.test %>% sample_frac(size = fraction_of_total_testSet)
  
  #TODO: implement a function to adjust balancing of classes outside the training 
  #Adjust classes for balancing:
  table(df.train$label)
  
  if(!is.na(balance_ratio)) {
    sampSize <- balance_ratio*round(sum(df.train$label != "2" & df.train$label != "0")/(length(unique(df.train$label - 2)))) #round(sum(df.train$label == "1" | df.train$label == "3")/2)
    if(sampSize < sum(df.train$label == "2")) {
      message("Sample size ", sampSize, " (including balance factor - ", balance_ratio, ") of 1C and 3C is smaller than 2C.")
      sampIdxs <- sample(x = which(df.train$label == "2"),
                         size = sum(df.train$label == "2") - min(sum(df.train$label == "2"), sampSize) + 1,
                         replace = F)
      df.train <- df.train[-sampIdxs, ]
      table(df.train$label)
    } else {
      message("Sample size ", sampSize, " (including balance factor - ", balance_ratio, ") of 1C and 3C is greater than 2C.")
    }
  }
  if(is.na(selected.features.names[1])) {
    selected.features.names <- colnames(df.train)[-ncol(df.train)]
    message("Running using all features.")
  } else {
    stopifnot(sum(!selected.features.names %in% colnames(df.train)) == 0)
    message("Selected features:")
    print(selected.features.names)
  }
  
  ####################################
  # (1)  - all or selected features:
  ####################################
  dtrain <- lgb.Dataset(data = data.matrix(df.train[, selected.features.names]), 
                        colnames = selected.features.names,
                        categorical_feature = NULL,
                        label = df.train$label, free_raw_data=T)
  dtest <- lgb.Dataset.create.valid(dataset = dtrain, data = data.matrix(df.test[, selected.features.names]), label = df.test$label)
  
  
  ####################################
  # (2) - specific features selected:
  ####################################
  
  # selected_features <- as.character(df_imp[df_imp$Gain > 0.005,]$Feature)
  # 
  # dtrain <- lgb.Dataset(data = data.matrix(df.train[, selected_features]), 
  #                       colnames = selected_features,
  #                       categorical_feature = NULL,
  #                       label = df.train$label, free_raw_data=T)
  # dtest <- lgb.Dataset.create.valid(dataset = dtrain, data = data.matrix(df.test[, selected_features]), label = df.test$label)
  # 
  
  
  valids <- list(test=dtest)
  num_classes <- length(unique(all.levels))
  
  #Training:
  params <- list(objective = "multiclass", metric = "multi_logloss", 
                 bagging_fraction = bagging_fraction, bagging_freq = bagging_freq)
  
  message("Begin training...")
  # training output printed. (verbose = 0 + record = T)
  bst <- lgb.train(params = params, data = dtrain, nrounds = nrounds, valids = valids,
                   num_threads = num_threads, learning_rate = learning_rate,
                   num_class = num_classes,
                   verbose = 1,
                   record = T,
                   early_stopping_rounds = early_stopping_rounds)
  
  lgb.save(booster = bst, filename = bst.fname)
  message("Completed training...")
  
  # multilogloss of final iteration on test set
  message(paste("# Rounds:", bst$current_iter()))
  message(paste("Multilogloss of best model:", bst$record_evals$test$multi_logloss$eval %>% unlist() %>% tail(1)))
  
  #Calculate variable importance. (note: takes awhile since single-threaded)
  df_imp <- tbl_df(lgb.importance(bst, percentage = TRUE))
  
  #Visualizations
  ###############
  pdf(file = report.fname)
  df_imp$Feature <- factor(df_imp$Feature, levels=rev(df_imp$Feature))
  p <- ggplot(df_imp, aes(x=Feature, y=Gain)) +
    geom_bar(stat="identity", fill="#34495e", alpha=0.9) +
    geom_text(aes(label=sprintf("%0.1f%%", Gain*100)), color="#34495e", hjust=-0.25, size=2.5) +
    #fte_theme() +
    coord_flip() +
    scale_y_continuous(limits = c(0, 0.9), labels=percent) +
    theme(plot.title=element_text(hjust=0.5), axis.title.y=element_blank()) +
    labs(title="Feature Importance", y="% of Total Gain in LightGBM Model")
  
  print(p)
  print(df_imp)
  
  message("Begin testing...")
  preds_matrix <- predict(bst, data.matrix(df.test[, selected.features.names]), reshape=T)
  preds_cor <- cor(preds_matrix)
  
  
  #get prediction labels (max prob)
  maxc <- max.col(preds_matrix)
  results <- do.call(rbind, lapply(1:nrow(preds_matrix), function(i) cbind(maxc[i] - 1, preds_matrix[i, maxc[i]])))
  
  df_results <- data.frame(results, label_act = df.test$label) %>%
    tbl_df() %>%
    transmute(label_pred = X1, prod_pred = X2, label_act)
  
  df_results %>% arrange(desc(prod_pred)) %>% head(20)
  
  require(caret)
  new.levels <- factor(sort(unique(df_results$label_act)), levels = as.numeric(as.factor(all.levels))-1)
  cm <- confusionMatrix(factor(df_results$label_pred, levels = new.levels), factor(df_results$label_act, levels = new.levels))
  print(data.frame(cm$overall))
  
  
  require(PRROC)
  #all aginast all PR
  pr <- pr.curve(scores.class0 = df_results$prod_pred, weights.class0 = (df_results$label_pred == df_results$label_act), curve = T)
  plot(pr)
  
  pm <- data.table(preds_matrix)
  setnames(pm, old = names(pm), new = as.character(new.levels))
  evaluateMultiClassPerformance(votes = pm, trueLabels = df.test$label)
  
  dev.off()
  message("All completed.")
}




train_using_lightbgm.regression <- function(trainingSet.fname = "/singlecellcenter/etai/tmp/CNML/learn_extraRegion32_insert_31_win31_gs_zs_v5.ALS.rds",
                                            testSet.fname = "/singlecellcenter/etai/tmp/CNML/learn_extraRegion32_insert_31_win31_gs_zs_v5.ATS1.rds",
                                            bst.fname = "/singlecellcenter/etai/tmp/CNML/learn_extraRegion32_insert_31_win31_gs_zs_v5.bst.all_ftrs.reg.rds",
                                            report.fname = "/singlecellcenter/etai/tmp/CNML/learn_extraRegion32_insert_31_win31_gs_zs_v5.trainTestReport.reg.pdf",
                                            selected.features.names = NA, allFeaturesPrecisionNumOfDigits = 3,
                                            fraction_of_total_trainSet = 1.0, fraction_of_total_testSet = 1.0, balance_ratio = NA,
                                            bagging_fraction = 1.0, bagging_freq = 0, nrounds = 1200,
                                            num_threads = 24, learning_rate = 0.1, early_stopping_rounds = 12) {
  require(data.table)
  require(lightgbm)
  require(dplyr)
  require(gtools)
  require(ggplot2)
  
  source('~ejacob/WORK/secondaryanalysis/methods_paper_results/R/lgbm.multiclass_example.Rstart.R')
  
  #load full region in CNV:
  message("Having ", length(trainingSet.fname), " training files to read..")
  if(length(trainingSet.fname) > 1) {
    df.train <- rbindlist(lapply(trainingSet.fname, function(f) readRDS(file = f)))
  } else {
    df.train <- readRDS(file = trainingSet.fname)
  }
  print(dim(df.train))
  
  all.levels <- mixedsort(as.character(unique(df.train$label)))
  df.train <- data.frame(round(df.train[,-ncol(df.train), with=F], digits = allFeaturesPrecisionNumOfDigits), label = factor(df.train$label, all.levels))
  
  #load full region in CNV:
  message("Having ", length(testSet.fname), " test files to read..")
  if(length(testSet.fname) > 1) {
    df.test <- rbindlist(lapply(testSet.fname, function(f) readRDS(file = f)))
  } else {
    df.test <- readRDS(file = testSet.fname)
  }
  print(dim(df.test))
  
  df.test <- data.frame(round(df.test[,-ncol(df.test), with=F], digits = allFeaturesPrecisionNumOfDigits), label = factor(df.test$label, all.levels))
  
  #adjust labels to lgbm foramt:
  df.train$label <- as.numeric(factor(df.train$label, levels = all.levels)) - 1
  df.test$label <- as.numeric(factor(df.test$label, levels = all.levels)) - 1
  
  #The data must be randomized for lightgbm to give unbiased scores. 
  df.train <- df.train %>% sample_frac(size = fraction_of_total_trainSet)
  df.test <- df.test %>% sample_frac(size = fraction_of_total_testSet)
  
  
  #Adjust classes for balancing:
  table(df.train$label)
  if(!is.na(balance_ratio)) {
    sampSize <- balance_ratio*round(sum(df.train$label != "2" & df.train$label != "0")/(length(unique(df.train$label - 2)))) #round(sum(df.train$label == "1" | df.train$label == "3")/2)
    if(sampSize < sum(df.train$label == "2")) {
      message("Sample size ", sampSize, " (including balance factor - ", balance_ratio, ") of 1C and 3C is smaller than 2C.")
      sampIdxs <- sample(x = which(df.train$label == "2"), 
                         size = sum(df.train$label == "2") - min(sum(df.train$label == "2"), sampSize) + 1, 
                         replace = F)
      df.train <- df.train[-sampIdxs, ]
      table(df.train$label)
    } else {
      message("Sample size ", sampSize, " (including balance factor - ", balance_ratio, ") of 1C and 3C is greater than 2C.")
    }
  }
  if(is.na(selected.features.names[1])) {
    selected.features.names <- colnames(df.train)[-ncol(df.train)] 
    message("Running using all features.")
  } else {
    stopifnot(sum(!selected.features.names %in% colnames(df.train)) == 0)
    message("Selected features:")
    print(selected.features.names)
  }
  
  
  ####################################
  # (1)  - all or selected features:
  ####################################
  dtrain <- lgb.Dataset(data = data.matrix(df.train[, selected.features.names]), 
                        colnames = selected.features.names,
                        categorical_feature = NULL,
                        label = df.train$label, free_raw_data=T)
  dtest <- lgb.Dataset.create.valid(dataset = dtrain, data = data.matrix(df.test[, selected.features.names]), label = df.test$label)
  
  # dtrain <- lgb.Dataset(data = data.matrix(df.train[, -ncol(df.train)]), 
  #                       colnames = colnames(df.train)[-ncol(df.train)],
  #                       categorical_feature = NULL,
  #                       label = df.train$label, free_raw_data=T)
  # dtest <- lgb.Dataset.create.valid(dataset = dtrain, data = data.matrix(df.test[, -ncol(df.test)]), label = df.test$label)
  # 
  valids <- list(test=dtest)
  
  
  #Training:
  params <- list(objective = "regression", metric = "l1", 
                 bagging_fraction = bagging_fraction, bagging_freq = bagging_freq)
  
  message("Begin training...")
  # training output printed. (verbose = 0 + record = T)
  bst <- lgb.train(params = params, data = dtrain, nrounds = nrounds, valids = valids,
                   num_threads = num_threads, learning_rate = learning_rate,
                   verbose = 1,
                   record = T,
                   early_stopping_rounds = early_stopping_rounds)
  
  lgb.save(booster = bst, filename = bst.fname)
  message("Completed training...")
  
  # multilogloss of final iteration on test set
  message(paste("# Rounds:", bst$current_iter()))
  message(paste("Multilogloss of best model:", bst$record_evals$test$l1$eval %>% unlist() %>% tail(1)))
  
  #Calculate variable importance. (note: takes awhile since single-threaded)
  df_imp <- tbl_df(lgb.importance(bst, percentage = TRUE))
  
  #Visualizations
  ###############
  df_imp$Feature <- factor(df_imp$Feature, levels=rev(df_imp$Feature))
  print(df_imp)
  
  pdf(file = report.fname)
  
  p <- ggplot(df_imp, aes(x=Feature, y=Gain)) +
    geom_bar(stat="identity", fill="#34495e", alpha=0.9) +
    geom_text(aes(label=sprintf("%0.1f%%", Gain*100)), color="#34495e", hjust=-0.25, size=2.5) +
    #fte_theme() +
    coord_flip() +
    scale_y_continuous(limits = c(0, 0.9), labels=percent) +
    theme(plot.title=element_text(hjust=0.5), axis.title.y=element_blank()) +
    labs(title="Feature Importance", y="% of Total Gain in LightGBM Model")
  
  print(p)
  
  message("Begin testing...")
  #preds_matrix <- predict(bst, data.matrix(df.test[, -ncol(df.test)]), reshape=T)
  preds_matrix <- predict(bst, data.matrix(df.test[, selected.features.names]), reshape=T)
  
  
  #get prediction labels (max prob)
  results <- data.frame(label_pred = preds_matrix, label_act = df.test$label)
  message(paste("# Corr test preds:", cor(results)[1,2]))
  p <- ggplot(results, aes(x=as.factor(label_act), y=label_pred)) + geom_boxplot(outlier.colour="darkgray", outlier.shape=16, outlier.size=0.5, notch=T) + xlab("Assigned CN") + ylab("Predicted Expression") + 
        ggtitle(sprintf("Cor = %.2f", round(cor(results)[1,2], digits = 2)))
  
  print(p)
  
  
  dev.off()
  message("All completed.")
}

evaluate_bst_performance <- function(bst.fname = "/singlecellcenter/etai/tmp/CNML/pointOfChange/test2/id17_noSSM.bst.ftrs_selected_by_0.0050_th.multiclass.rds", 
                                     testSet.fname = "/singlecellcenter/etai/tmp/CNML/pointOfChange/test2/id17_noSSM.testing_set.resampled.rds",
                                     selected_ftrs_fname = "/singlecellcenter/etai/tmp/CNML/pointOfChange/test2/id17_noSSM.names_of_selected_features_by_0.0050_th.multiclass.rds", 
                                     selected.features.names = NA, loadgraphicalParams = F) {
  if(loadgraphicalParams)
    source('~ejacob/WORK/secondaryanalysis/methods_paper_results/R/lgbm.multiclass_example.Rstart.R')
  require(lightgbm)
  require(dplyr)
  
  if(!is.na(selected_ftrs_fname)) {
    selected.features.names <- readRDS(selected_ftrs_fname)
  }
  bst <- lgb.load(filename = bst.fname)
  #load full region in CNV:
  message("Having ", length(testSet.fname), " test files to read..")
  if(length(testSet.fname) > 1) {
    df.test <- rbindlist(lapply(testSet.fname, function(f) readRDS(file = f)))
  } else {
    df.test <- readRDS(file = testSet.fname)
  }
  print(dim(df.test))
  
  all.levels <- mixedsort(as.character(unique(df.test$label)))
  df.test <- data.frame(round(df.test[,-ncol(df.test), with=F], digits = 3), label = factor(df.test$label, all.levels))
  
  #adjust labels to lgbm foramt:
  df.test$label <- as.numeric(factor(df.test$label, levels = all.levels)) - 1
  
  df.test <- df.test %>% sample_frac(size = 1.0)
  
  
  print(table(df.test$label))
  
  #Calculate variable importance. (note: takes awhile since single-threaded)
  df_imp <- tbl_df(lgb.importance(bst, percentage = TRUE))
  print(df_imp)
  #Visualizations
  ###############
  
  df_imp$Feature <- factor(df_imp$Feature, levels=rev(df_imp$Feature))
  p <- ggplot(df_imp[1:10,], aes(x=Feature, y=Gain)) +
    geom_bar(stat="identity", fill="#34495e", alpha=0.9) +
    geom_text(aes(label=sprintf("%0.1f%%", Gain*100)), color="#34495e", hjust=-0.25, size=8) +
    #fte_theme() +
    coord_flip() +
    scale_y_continuous(limits = c(0, 0.9), labels=percent) +
    theme(plot.title=element_text(hjust=0.5, size=20), 
          axis.title.x = element_text(size=15), axis.title.y = element_blank(), 
          axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) +
    labs(y="% of Total Gain in LightGBM Model")
  
  print(p)
  print(df_imp)
  
  if(is.na(selected.features.names[1])) {
    selected.features.names <- colnames(df.test)[-ncol(df.test)] 
    message("Running using all features.")
  } else {
    stopifnot(sum(!selected.features.names %in% colnames(df.test)) == 0)
    message("Selected features:")
    print(selected.features.names)
  }
  
  
  message("Begin testing...")
  preds_matrix <- predict(bst, data.matrix(df.test[, selected.features.names]), reshape=T)
  preds_cor <- cor(preds_matrix)
  
  
  #get prediction labels (max prob)
  maxc <- max.col(preds_matrix)
  results <- do.call(rbind, lapply(1:nrow(preds_matrix), function(i) cbind(maxc[i] - 1, preds_matrix[i, maxc[i]])))
  
  df_results <- data.frame(results, label_act = df.test$label) %>%
    tbl_df() %>%
    transmute(label_pred = X1, prod_pred = X2, label_act)
  
  df_results %>% arrange(desc(prod_pred)) %>% head(20)
  
  require(caret)
  new.levels <- factor(sort(unique(df_results$label_act)), levels = as.numeric(as.factor(all.levels))-1)
  cm <- confusionMatrix(factor(df_results$label_pred, levels = new.levels), factor(df_results$label_act, levels = new.levels))
  print(data.frame(cm$overall))
  
  
  require(PRROC)
  #all aginast all PR
  pr <- pr.curve(scores.class0 = df_results$prod_pred, weights.class0 = (df_results$label_pred == df_results$label_act), curve = T)
  plot(pr)
  
  pm <- data.table(preds_matrix)
  setnames(pm, old = names(pm), new = as.character(new.levels))
  allaucs <- evaluateMultiClassPerformance(votes = pm, trueLabels = df.test$label)
  
  return(allaucs)
}
#normalization of expression to ratio of ref IDs is not by default - 
#this function gets the adt un-normalized and uses the function getExpRatioforML to normalize after artificial substitution of the relevant segment
buildLearningOrTestSet <- function(adt, label = 1, controlSampleIDs, geneSpecificFtrs = F, doZscore = F,
                                   native_id = "161005c_A5", native_id2 = NULL, 
                                   alt_id =  "0705_B12", alt_id2 = NULL, 
                                   pseudo_id = "170425_A1", pseudo_id2 = NULL,
                                   calibrateAlt = T, calibratePseudo = T, doRatioBeforeScaleAndZscore = F, amplitudeShouldBeNegative = T,
                                   chr = "chr1", doScale = F, useGlobalFtrs = F, min_sd = 0.0001,
                                   insertSize = 20, extraRegion = 31, ftrsWindowSize = 11, 
                                   QCplotMe = F, typeOfEventToGenerate = c("amplitude", "pointOfChange"), fuzzyRegion = 0,
                                   numOfAlteredInserts = 300, numOfPseudoInserts = 150, 
                                   useMA = F, useSSM = F,
                                   getSimulatedSample = F, issueRndForSimu = 0) {
  
  message(chr)
  
  message("Alt: ", alt_id, ", native: ", native_id, ", pseudo: ", pseudo_id)
  M <- copy(adt[, controlSampleIDs, with = F] + 1)
  q <- min(1500, unname(quantile(as.matrix(M), 0.99)))
  M[M > q] <- q
  
  myMeans <- rowMeans(M, na.rm = T)
  message("myMeans has ", sum(myMeans == 0), " zeros, ", sum(is.na(myMeans)), " NAs and ", sum(is.nan(myMeans)), " not a number")
  mySds <- rowSds(as.matrix(M), na.rm = T)
  mySds <- ifelse(mySds < min_sd, min_sd, mySds)
  message("mySds has ", sum(mySds == 0), " zeros, ", sum(is.na(mySds)), " NAs and ", sum(is.nan(mySds)), " not a number")
  
  if(doRatioBeforeScaleAndZscore) {
    M2 <- log2(sweep(M, 1, rowMeans(as.matrix(M)), "/"))
    myMeans2 <- rowMeans(M2, na.rm = T)
    message("myMeans2 has ", sum(myMeans2 == 0), " zeros, ", sum(is.na(myMeans2)), " NAs and ", sum(is.nan(myMeans2)), " not a number")
    mySds2 <- rowSds(as.matrix(M2), na.rm = T)
    mySds2 <- ifelse(mySds2 < min_sd, min_sd, mySds2)
    message("mySds2 has ", sum(mySds2 == 0), " zeros, ", sum(is.na(mySds2)), " NAs and ", sum(is.nan(mySds2)), " not a number")
  }
  
  geneSpecific <- NA
  if(geneSpecificFtrs) {
    geneSpecific <- generateGeneSpecificFtrs(adt, myMeans, mySds, chr, base = 2, digits = 2)
  }
  
  
  #Alt assignment:
  if(calibrateAlt) {
    if(is.null(alt_id2) || is.na(alt_id2)) {
      alterInsertDT <- calibrating_pair_of_samples(adt = adt, insert_id = alt_id, native_id = native_id, chr = chr, doPlot = F, adtColumns = 1:4)
    } else {
      message("Doing +alt_id2")
      alterInsertDT1 <- calibrating_pair_of_samples(adt = adt, insert_id = alt_id, native_id = native_id, chr = chr, doPlot = F, adtColumns = 1:4)
      alterInsertDT2 <- calibrating_pair_of_samples(adt = adt, insert_id = alt_id2, native_id = native_id, chr = chr, doPlot = F, adtColumns = 1:4)
      alterInsertDT <- copy(alterInsertDT1)
      alterInsertDT$y <- alterInsertDT$y + alterInsertDT2$y
    }
  } else {
    if(is.null(alt_id2) || is.na(alt_id2)) {
      alterInsertDT <- data.table(adt[, 1:4, with = F], y = adt[[alt_id]])
    } else {
      message("Doing +alt_id2")
      alterInsertDT <- data.table(adt[, 1:4, with = F], y = adt[[alt_id]] + adt[[alt_id2]])
    }
  }
  
 
  alterInsertDT <- getExpRatioforML(mydt = alterInsertDT, doZscore = doZscore, mySds = mySds, doRatioBeforeScaleAndZscore = doRatioBeforeScaleAndZscore, 
                                    min_sd = min_sd, mySds2 = mySds2, myMeans2 = myMeans2, normalizationFactor = myMeans, doScale = F)
  
  
  #Pseudo:
  if(calibratePseudo) {
    if(is.null(pseudo_id2) || is.na(pseudo_id2)) {
      pseudoInsertDT <- calibrating_pair_of_samples(adt = adt, insert_id = pseudo_id, native_id = native_id, chr = chr, doPlot = F, adtColumns = 1:4)
    } else {
      message("Doing +pseudo_id2")
      pseudoInsertDT1 <- calibrating_pair_of_samples(adt = adt, insert_id = pseudo_id, native_id = native_id, chr = chr, doPlot = F, adtColumns = 1:4)
      pseudoInsertDT2 <- calibrating_pair_of_samples(adt = adt, insert_id = pseudo_id2, native_id = native_id, chr = chr, doPlot = F, adtColumns = 1:4)
      pseudoInsertDT <- copy(pseudoInsertDT1)
      pseudoInsertDT$y <- pseudoInsertDT$y + pseudoInsertDT2$y
    }
  } else {
    if(is.null(pseudo_id2) || is.na(pseudo_id2)) {
      pseudoInsertDT <- data.table(adt[, 1:4, with = F], y = adt[[pseudo_id]])
    } else {
      message("Doing +pseudo_id2")
      pseudoInsertDT <- data.table(adt[, 1:4, with = F], y = adt[[pseudo_id]] + adt[[pseudo_id2]])
    }
  }
  
  if(doScale) {
    if(sum(abs(adt[[alt_id]])) < 1) { #in case the vector is whole genome loss (all zeros)
      message(alt_id, " is zero - using ", pseudo_id, " for median.")
      tmpy <- getExpRatioforML(mydt = pseudoInsertDT, doZscore = doZscore, mySds = mySds, doRatioBeforeScaleAndZscore = doRatioBeforeScaleAndZscore, 
                               min_sd = min_sd, mySds2 = mySds2, myMeans2 = myMeans2, normalizationFactor = myMeans, doScale = F)
      
      alterInsertDT$y <- alterInsertDT$y - median(tmpy$y)
    } else {
      alterInsertDT$y <- alterInsertDT$y - median(alterInsertDT$y)
    }
  }
  
  pseudoInsertDT <- getExpRatioforML(mydt = pseudoInsertDT, doZscore = doZscore, mySds = mySds, doRatioBeforeScaleAndZscore = doRatioBeforeScaleAndZscore, 
                                     min_sd = min_sd, mySds2 = mySds2, myMeans2 = myMeans2, normalizationFactor = myMeans, doScale = doScale)
  
  
  
  #Native:
  nativeDT <- data.table(adt[, 1:4, with = F], y = adt[[native_id]])
  
  
  if(is.null(native_id2) || is.na(native_id2)) {
    #(no need to calibrate since everything else is calibrated to native)
    nativeDT <- data.table(adt[, 1:4, with = F], y = adt[[native_id]])
  } else {
    message("Doing +native_id2")
    if(calibratePseudo) {
      #(need to calibrate only native_id2)
      nativeDT2 <- calibrating_pair_of_samples(adt = adt, insert_id = native_id2, native_id = native_id, chr = chr, doPlot = F, adtColumns = 1:4)
      nativeDT$y <- nativeDT$y + nativeDT2$y
    } else {
      nativeDT <- data.table(adt[, 1:4, with = F], y = adt[[native_id]] + adt[[native_id2]])
    }
    
  }
  
  nativeDT <- getExpRatioforML(mydt = nativeDT, doZscore = doZscore, mySds = mySds, doRatioBeforeScaleAndZscore = doRatioBeforeScaleAndZscore, 
                               min_sd = min_sd, mySds2 = mySds2, myMeans2 = myMeans2, normalizationFactor = myMeans, doScale = doScale)
  
  useAltSampleForSSMnorm = T
  if(alt_id == "FullLoss")
    useAltSampleForSSMnorm = F
    
  mySet <- generate_learning_entries(nativeDT = nativeDT, alterInsertDT = alterInsertDT, pseudoInsertDT = pseudoInsertDT, useGlobalFtrs = useGlobalFtrs, QCplotMe = QCplotMe,
                                     chr = chr, insertSize = insertSize, extraRegion = extraRegion, ftrsWindowSize = ftrsWindowSize, geneSpecific = geneSpecific,
                                     numOfAlteredInserts = numOfAlteredInserts, numOfPseudoInserts = numOfPseudoInserts, 
                                     useMA = useMA, useSSM = useSSM, useAltSampleForSSMnorm = useAltSampleForSSMnorm,
                                     typeOfEventToGenerate = typeOfEventToGenerate, fuzzyRegion = fuzzyRegion, amplitudeShouldBeNegative = amplitudeShouldBeNegative,
                                     getSimulatedSample = getSimulatedSample, issueRndForSimu = issueRndForSimu)
  
  if(getSimulatedSample)
    return(mySet)
  
  if(length(mySet$pseudoInsert) >= 1 && !is.na(mySet$pseudoInsert)) {
    #mySet$pseudoInsert$label[mySet$pseudoInsert$label == 1] <- label
    mySet$pseudoInsert$label <- paste(label, mySet$pseudoInsert$label, sep = "_")
  }
  
  if(length(mySet$alteredInsert) >= 1 && !is.na(mySet$alteredInsert)) {
    #mySet$alteredInsert$label[mySet$alteredInsert$label == 1] <- label
    mySet$alteredInsert$label <- paste(label, mySet$alteredInsert$label, sep = "_")
  }
  return(mySet)
  
} 

#mySet is a data.table with at least one column of name "label"
resample_set <- function(mySet, minimalCountLabels = c("1to2_1", "2to3_1"), countFactor = 1.0) {
  tbl <- table(mySet$label)
  message("Before resampling:")
  print(tbl)  
  if(sum(!minimalCountLabels %in% names(tbl)) > 0) {
    stop("ERROR - minimal count labels are not in the set")
  }
  
  minCount <- round(min(tbl[minimalCountLabels]))*countFactor
  
  reSampledSet <- rbind(mySet[label %in% minimalCountLabels], mySet[!label %in% minimalCountLabels,.SD[sample(.N, min(minCount,.N))],by = "label"])
  message("Following resampling:")
  print(table(reSampledSet$label))
  return(reSampledSet)
}


QC.generate_learning_entries <- function() {
  posAndFilterFile <- "/singlecellcenter/etai/ExperimentsData/Stamatis/May072018/MVSSM/GSSM.1/data/dt.1.rsemtpm.rds"
  rawTPMfile <- "/singlecellcenter/etai/ExperimentsData/Stamatis/Stamatis_list_v12_180404.rsemtpm.rds"
  min_sd <- 0.0001
  rsemtpm <- readRDS(rawTPMfile)
  dt <- readRDS(posAndFilterFile)[, 1:4]
  dt <- dt[order(seqnames, start)]
  adt <- merge(dt, data.table(rsemtpm[dt$id, ], keep.rownames = "id"), by = "id")
  adt <- adt[order(seqnames, start)]
  adt <- adt[seqnames != "chrY" & seqnames != "chrM"]
  
  
  file.1C = "/singlecellcenter/etai/data/WORK/secondaryanalysis/methods_paper_results/data/list_of_predicted_monosomies4.csv"
  file.2C = "/singlecellcenter/etai/data/WORK/secondaryanalysis/methods_paper_results/data/list_of_predicted_disomies4.csv"
  file.3C = "/singlecellcenter/etai/data/WORK/secondaryanalysis/methods_paper_results/data/list_of_predicted_trisomies4.csv"
  
  test.1C <- fread(file.1C, header = T)
  test.2C <- fread(file.2C, header = T)
  test.3C <- fread(file.3C, header = T)
  
  #training set construcution:
  # F (no change) | T (point of change)
  # -----------------------------------
  # Zero          | 0 to 1, 0 to 2, 0 to 3, 0 to 4
  # Mono          | 1 to 2, 1 to 3, 1 to 4
  # diploid       | 2 to 3, 2 to 4
  # Trisomy       | 3 to 4
  #------------------------------------
  # In total we have 14 sub groups of training set we need to construct
  
  cs <- generateLearningAndTestSamples.pointOfChange(adt = adt, minDetectedGenes = 4000, excludeIDs = excludeIDs, trainProb = 0.8, totalLossCount = 40,
                                                     chrToExclude = c("chrX", "chrY", "chrM", "chr10", "chr8", "chr21"), do4C = T,
                                                     file.1C = "/singlecellcenter/etai/data/WORK/secondaryanalysis/methods_paper_results/data/list_of_predicted_monosomies4.csv",
                                                     file.2C = "/singlecellcenter/etai/data/WORK/secondaryanalysis/methods_paper_results/data/list_of_predicted_disomies4.csv",
                                                     file.3C = "/singlecellcenter/etai/data/WORK/secondaryanalysis/methods_paper_results/data/list_of_predicted_trisomies4.csv")
                  
  #mono1
  #As in the ML:
  #positive learning set:
  par(mfrow=c(2,1))
  tmp <- buildLearningOrTestSet(adt = adt, label = 1, controlSampleIDs = controlSampleIDs, geneSpecificFtrs = F, doZscore = T,
                                 native_id = test.3C[seqnames=="chr1"]$id[2], alt_id =  test.1C[seqnames=="chr1"]$id[3], alt_id2 = NULL, pseudo_id = "170425_A1", 
                                 calibrateAlt = F, calibratePseudo = F, doRatioBeforeScaleAndZscore = T,
                                 chr = "chr1", doScale = T, useGlobalFtrs = F, min_sd = 0.0001,
                                 insertSize = 42, extraRegion = 0, ftrsWindowSize = 31, 
                                 QCplotMe = T, typeOfEventToGenerate = c("amplitude", "pointOfChange")[2], fuzzyRegion = 2,
                                 numOfAlteredInserts = 1, numOfPseudoInserts = 0, useMA = F)
  
  #negative learrning set
  tmp <- buildLearningOrTestSet(adt = adt, label = 1, controlSampleIDs = controlSampleIDs, geneSpecificFtrs = F, doZscore = T,
                                native_id = test.2C[seqnames=="chr1"]$id[6], alt_id =  test.1C[seqnames=="chr1"]$id[2], alt_id2 = NULL, pseudo_id = "170425_A1", 
                                calibrateAlt = F, calibratePseudo = F, doRatioBeforeScaleAndZscore = T,
                                chr = "chr1", doScale = T, useGlobalFtrs = F, min_sd = 0.0001,
                                insertSize = 30, extraRegion = 0, ftrsWindowSize = 27, 
                                QCplotMe = T, typeOfEventToGenerate = c("amplitude", "pointOfChange")[1], fuzzyRegion = 0,
                                numOfAlteredInserts = 1, numOfPseudoInserts = 0, useMA = F)
  
  #with SSM:
  
  tmp <- buildLearningOrTestSet(adt = adt, label = 1, controlSampleIDs = controlSampleIDs, geneSpecificFtrs = F, doZscore = T,
                                native_id = test.2C[seqnames=="chr1"]$id[2], alt_id =  test.1C[seqnames=="chr1"]$id[3], alt_id2 = NULL, pseudo_id = "170425_A1", 
                                calibrateAlt = F, calibratePseudo = F, doRatioBeforeScaleAndZscore = T,
                                chr = "chr1", doScale = T, useGlobalFtrs = F, min_sd = 0.0001,
                                insertSize = 32, extraRegion = 0, ftrsWindowSize = 31, 
                                QCplotMe = T, typeOfEventToGenerate = c("amplitude", "pointOfChange")[2], fuzzyRegion = 2,
                                numOfAlteredInserts = 1, numOfPseudoInserts = 0, useMA = T, useSSM = T)
  
  tmp <- buildLearningOrTestSet(adt = adt, label = 1, controlSampleIDs = controlSampleIDs, geneSpecificFtrs = F, doZscore = T,
                                native_id = test.2C[seqnames=="chr1"]$id[2], alt_id =  test.3C[seqnames=="chr1"]$id[3], alt_id2 = NULL, pseudo_id = "170425_A1", 
                                calibrateAlt = F, calibratePseudo = F, doRatioBeforeScaleAndZscore = T,
                                chr = "chr1", doScale = T, useGlobalFtrs = F, min_sd = 0.0001,
                                insertSize = 32, extraRegion = 0, ftrsWindowSize = 31, 
                                QCplotMe = T, typeOfEventToGenerate = c("amplitude", "pointOfChange")[1], fuzzyRegion = 2,
                                numOfAlteredInserts = 1, numOfPseudoInserts = 0, useMA = T, useSSM = T)
  
  
  
  
}

generateGeneSpecificFtrs <- function(adt, myMeans, mySds, chr = NULL, base = 2, digits = 2, GC_length_rds_fname = NULL) {
  if(is.null(GC_length_rds_fname)) {
    glen <- fread("/singlecellcenter/etai/ReferenceData/Gencode/GC_lengths.tsv") #produced by: GTF2LengthGC.R in annotations dir - the sum of all exons' length
    setnames(glen, old = "V1", new = "id")
  } else {
    message("Reading GC and length from input file..")
    glen <- readRDS(GC_length_rds_fname)
    names(glen)[1] <- "id"
  }
  
  message("Gene specific features: ", nrow(glen), " from tpm data: ", nrow(adt))
  setkey(glen, "id")
  glen <- glen[adt$id]
  message("Following reduction - Gene specific features: ", nrow(glen), " from tpm data: ", nrow(adt))
  myTCV <- log2((mySds + 0.1)/(myMeans))
  myTCV[which(myTCV < -5)] <- -5
  myTCV[which(myTCV > 5)] <- 5
  
  myTCV <- round(scale(myTCV, center = T, scale = T), digits = digits)
  geneSpecific <- myTCV #myMeans => v7
  #TODO: maybe use the new probability based 1-exp transformation function
  interspace <- unlist(lapply(unique(adt$seqnames), function(chr) c(adt[seqnames == chr]$start[-1], data.table::last(adt[seqnames == chr]$end)+1) - adt[seqnames == chr]$end + 1.5))
  interspace[interspace <= 0] <- 1
  interspace <- round(log(x = interspace, base = base), digits = digits)
  
  geneSpecific <- list(TCV = myTCV, interspace = interspace, exonsLenSum = glen$Length, exonsGC = glen$GC)
  return(geneSpecific)
}

#Input:
#ftrsWindowSize - odd integer
generate_learning_entries <- function(nativeDT, alterInsertDT, pseudoInsertDT, useGlobalFtrs = F, geneSpecific = NA, amplitudeShouldBeNegative = T,
                                      chr, insertSize = 21, extraRegion = 21, ftrsWindowSize = 11, chrsToExcludeFromNormalization = c("chrX", "chrM", "chr10"),
                                      typeOfEventToGenerate = c("amplitude", "pointOfChange"), fuzzyRegion = 0, QCplotMe = F, getSimulatedSample = F, issueRndForSimu = 0,
                                      numOfAlteredInserts = 5, numOfPseudoInserts = 5, useMA = F, useSSM = F, useAltSampleForSSMnorm = T, MA_winSize = 50) {

  require(data.table)
  
  
  
  getMyLearningSet <- function(mydt, altdt, N) {
    
    rng <- range(which(mydt$seqnames == chr))
    
    mySizeLimit <- max(insertSize, ftrsWindowSize)
    z <- ceiling(((mySizeLimit - 1)/2))
    wi <- which(mydt$seqnames == chr)
    wi <- wi[-c(1:(z), (length(wi)-z + 1):length(wi))]
    
    L <- length(wi)
    if(N > L) {
      message("Sampling cannot be bigger than chromosome length.\nAdjusting to chromosome size*0.8 (", L, ")")
      N <- round(L*0.8)
    }
    
    
    
    if(extraRegion > 0) {
      rnd <- sample(x = wi, size = N, replace = F)  
    } else {
      rnd <- sample(x = wi, size = N, replace = F)
    }
    
    if(issueRndForSimu > 0 & getSimulatedSample == T) {
      rnd <- issueRndForSimu
      rnd <- data.frame(from = rnd - round((insertSize - 1)/2), 
                        to = rnd + round((insertSize - 1)/2),
                        midPoint = rnd)  
    } else {
      rnd <- data.frame(from = rnd - round((insertSize - 1)/2), 
                        to = rnd + round((insertSize - 1)/2),
                        midPoint = rnd)
    }
    #calculating norm factor outside the random loop to decrease computation time.
    #the disadvantage is that we do not include the simulated abberation in the normalizaton but since it is small (e.g. 31 genes) and the
    #norm factor a median of the whole chrs then it should be neglegable
    if(useSSM) {
      
      #IMF:
      #ssm <- mydt2[, lapply(.SD, runMyIterativeMedianFilter, levels = #), .SDcols = names(mydt2)[-c(1:4)], by = seqnames]
      #SSM global:
      if(useAltSampleForSSMnorm) {
        #in case 3/4 of the alt sample is almost unexpressed in comparisson to the host sample we take the host sample anyway:
        if(sum(altdt$y < quantile(mydt$y, 0.50))/nrow(mydt) > 0.75) {
          ssm <- mydt[, lapply(.SD, getMyGSSM), .SDcols = names(mydt)[-c(1:4)], by = seqnames]
        } else {
          ssm <- altdt[, lapply(.SD, getMyGSSM), .SDcols = names(altdt)[-c(1:4)], by = seqnames]  
        }
      } else {
        ssm <- mydt[, lapply(.SD, getMyGSSM), .SDcols = names(mydt)[-c(1:4)], by = seqnames]
      }
      #cell normalization:
      dtm <- ssm[, lapply(.SD, median),
                 .SDcols = c("y"), by = "seqnames"]
      ssmNorm <- median(dtm[!(seqnames %in% chrsToExcludeFromNormalization)]$y)
    }
    
    
    
    getMyLSBatch <- function(mypos) {
      from <- mypos[1]
      to <- mypos[2]
      midPointPos <- mypos[3]
      
      grab <- max(rng[1], (from - extraRegion)):min(rng[2], (to + extraRegion))
      
      
      C <- rep(0, length(mydt$y))
      
      flip = 1
      if(typeOfEventToGenerate[1] == "amplitude") {
        #class assignment:
        if(amplitudeShouldBeNegative) {
          C[from:to] <- 0
        } else {
          C[from:to] <- 1
        }
        #values assignments:
        Y <- mydt$y 
        Y[from:to] <- altdt$y[from:to]
      } else if(typeOfEventToGenerate[1] == "pointOfChange") {
        #class assignment:
        C[(midPointPos - fuzzyRegion):(midPointPos + fuzzyRegion)] <- 1
        #values assignments:
        Y <- mydt$y 
        if(sample(2)[1] == 1) { #flip left
          Y[from:midPointPos] <- altdt$y[from:midPointPos]
        } else { #flip right
          Y[midPointPos:to] <- altdt$y[midPointPos:to]
          flip = 2
        }
      } else {
        stop("No such event ", typeOfEventToGenerate[1], " quiting..")
      }
      
      if(getSimulatedSample) {
        simu.tmp <- copy(mydt)
        simu.tmp$y <- Y
        simu.tmp$ratio <- Y
        simuSample <- list(mydt = simu.tmp, pos = list(from = mydt[from, 1:4], to = mydt[to, 1:4], midPointPos = mydt[midPointPos, 1:4]))
      }
      
      C <- C[grab]
      Yalt <- Y
      Y <- Y[grab]
      
      if(QCplotMe) {
        plot(Y, lwd = 3, col = "blue", main = sprintf("Window Size = %d", ftrsWindowSize), cex.lab=1.5, cex = 1.5, 
             ylab = "Relative Score Ratio", xlab = "Gene index in Window")
        abline(v=range(which(C == 1)), col = "green", lwd = 5, lty = "dashed")
        if(flip == 1) {
          #abline(v=min(which(C == 1)) - 3, col = "red", lwd = 3, lty = "dashed")
          rect(0, min(Y)-1, min(which(C == 1)), max(Y)+1, density = NULL, angle = 45,
               col = alpha(colour = "red", alpha = 0.3), border = NA)
        } else {
          #abline(v=max(which(C == 1)) + 3, col = "red", lwd = 3, lty = "dashed")
          rect(max(which(C == 1)), min(Y)-1, length(Y), max(Y)+1, density = NULL, angle = 45,
               col = alpha(colour = "red", alpha = 0.3), border = NA)
        }
      }
      Ys <- list()
      Ys$tpm <- Y
      if(useMA) {
        mydt2 <- copy(mydt)
        if(typeOfEventToGenerate[1] == "amplitude") {
          mydt2$y[from:to] <- altdt$y[from:to] 
        } else if(typeOfEventToGenerate[1] == "pointOfChange") {
          mydt2$y <- Yalt
        } else {
          stop("No such event ", typeOfEventToGenerate[1], " quiting..")
        }
        ma <- mydt2[, lapply(.SD, rollmean, k = MA_winSize, fill = "extend", align = "center"),
                    .SDcols = names(mydt2)[-c(1:4)], by = seqnames]
        #cell normalization:
        dtm <- ma[, lapply(.SD, median),
                  .SDcols = c("y"), by = "seqnames"]
        ma <- ma$y - median(dtm[!(seqnames %in%chrsToExcludeFromNormalization)]$y)
        
        #ma <- ma$y - median(ma$y, na.rm = T) #this is the old normalization 
        ma <- ma[grab]
        Ys$MA <- ma
        if(QCplotMe) {
          #lines(ma, lwd = 3, col = "magenta")
          # lines(rollmean(Ys$tpm - median(dtm[!(seqnames %in%chrsToExcludeFromNormalization)]$y), k=ftrsWindowSize, fill = "extend", align = "center"), 
          #       lwd = 3, col = "magenta")
        }
        #Y <- ma
      } 
      
      if(useSSM) {
        # mydt2 <- copy(mydt)
        # if(typeOfEventToGenerate[1] == "amplitude") {
        #   mydt2$y[from:to] <- altdt$y[from:to] 
        # } else if(typeOfEventToGenerate[1] == "pointOfChange") {
        #   if(sample(2)[1] == 1) { #flip left
        #     mydt2$y[from:midPointPos] <- altdt$y[from:midPointPos]
        #   } else { #flip right
        #     mydt2$y[midPointPos:to] <- altdt$y[midPointPos:to]
        #   }
        # } else {
        #   stop("No such event ", typeOfEventToGenerate[1], " quiting..")
        # }
        # 
        
        #IMF:
        #ssm <- mydt2[, lapply(.SD, runMyIterativeMedianFilter, levels = #), .SDcols = names(mydt2)[-c(1:4)], by = seqnames]
        #SSM global:
        #ssm <- mydt2[, lapply(.SD, getMyGSSM), .SDcols = names(mydt2)[-c(1:4)], by = seqnames]
        #SSM local:
        ssm_local <- getMyGSSM(mv1 = Ys$tpm)
        
        #cell normalization:
        #dtm <- ssm[, lapply(.SD, median),
        #            .SDcols = c("y"), by = "seqnames"]
        #ssm <- ssm$y - median(dtm[!(seqnames %in% chrsToExcludeFromNormalization)]$y)
        
      
        Ys$SSM <- ssm_local - ssmNorm
        #ssm <- ssm_local - median(dtm[!(seqnames %in% chrsToExcludeFromNormalization)]$y)#ssm[grab]
        #Ys$SSM <- ssm
        if(QCplotMe) {
          lines(Ys$SSM, lwd = 6, col = "darkgreen")
          abline(h=0, lwd=3, lty="dashed", col="gray")
          # legend("bottomleft", fill = c("blue", "magenta", "darkgreen", "green"),
          #        inset = 0.02, legend = c("Gene Relative TPM", "MA Score", "SSM Score", "Point of Change"), cex=1)
          legend("bottomleft", fill = c("blue","darkgreen", "green"),
                 inset = 0.02, legend = c("Gene Relative TPM", "SSM Score", "Point of Change"), cex=1)
        }
        
        #Y <- ma
      } 
      if(getSimulatedSample) {
        message("Returning simu sample")
        return(simuSample)
      }
      gs <- NA
      if(!is.na(geneSpecific[1])) {
        gs <- lapply(names(geneSpecific), function(x) geneSpecific[[x]][grab])
        names(gs) <- names(geneSpecific)
      }
      
      midPoint <- (ftrsWindowSize - 1)/2 + 1
      #ftrs <- do.call(rbind, wapply(Y = Y, width = ftrsWindowSize, by = 1))
      if(useGlobalFtrs) {
        ftrs <- generateFeaturesForExpressionVector(Ys = Ys, global = Y, ftrsWindowSize = ftrsWindowSize, geneSpecific = gs, isLogged = T)[,-1]
      } else {
        ftrs <- generateFeaturesForExpressionVector(Ys = Ys, ftrsWindowSize = ftrsWindowSize, geneSpecific = gs, isLogged = T)[,-1]
      }
      classes <- do.call(rbind, wapply(Y = C, width = ftrsWindowSize, by = 1))[,midPoint]
      
      myLSB <- data.table(ftrs, label = classes)
      
      return(myLSB)
    }
    if(getSimulatedSample) {
      print(as.vector(as.matrix(rnd[1, ])))
      return(getMyLSBatch(as.vector(as.matrix(rnd[1, ]))))
    }
    
    res <- rbindlist(apply(rnd, 1, getMyLSBatch))
    return(res)
  }
  
  
  #Positive set:
  if(numOfAlteredInserts > 0) {
    mydt <- copy(nativeDT)
    altdt <- copy(alterInsertDT)
    if(getSimulatedSample)
      return(getMyLearningSet(mydt = mydt, altdt = altdt, N = numOfAlteredInserts))
    
    alteredInsertLS <- getMyLearningSet(mydt = mydt, altdt = altdt, N = numOfAlteredInserts)
  } else {
    alteredInsertLS <- NA
  }
  
  
  #Pseudo positive set:
  if(numOfPseudoInserts > 0) {
    mydt <- copy(nativeDT)
    altdt <- copy(pseudoInsertDT)
    pseudoInsertLS <- getMyLearningSet(mydt = mydt, altdt = altdt, N = numOfPseudoInserts)
  } else {
    pseudoInsertLS <- NA
  }
  #Negative set:
  #? should we use a negative set of just the nativeDT?
  #? should we use a negative set of the alteredDT of another diploidic chr?
  
  return(list(pseudoInsert = pseudoInsertLS, alteredInsert = alteredInsertLS))
  
}

runMyIterativeMedianFilter <- function(x, t = NA, min_win = 15, levels = 3, plotMe = F) {
  require(signal)
  m <- 0
  if(plotMe)
    plot(t, x)
  xx <- list()
  for(i in 1:max(levels)) {
    if(plotMe)
      message(min_win*i + m)
    x <- signal::filter(filt = MedianFilter(min_win*i + m), x = x)
    if(i %in% levels)
      xx[[sprintf("level%d", i)]] <- x
    m <- 1
  }
  if(plotMe) {
    lines(t, x, lwd = 3, col="red")
    abline(v = 65600208, lty = "dashed", col = "green", lwd=2)
  }
  
  return(xx[[1]])
}


getMyGSSM <- function(mv1) {
  require(KFAS)
  Zt <- matrix(c(1, 0), 1, 2)
  Ht <- matrix(NA)
  Tt <- matrix(c(1, 0, 1, 1), 2, 2)
  Rt <- matrix(c(1, 0), 2, 1)
  Qt <- matrix(NA)
  a1 <- matrix(c(1, 0), 2, 1)
  P1 <- matrix(0, 2, 2)
  P1inf <- diag(2)
  
  if(sum(mv1, na.rm = T) == 0)
    mv1[sample(length(mv1), 2)] <- 0.00001
  
  model_gaussian <- SSModel(mv1 ~ -1 +
                              SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf), 
                            H = Ht)
  fit_gaussian <- fitSSM(model_gaussian, inits = c(0, 0), method = "BFGS")
  out_gaussian <- KFS(fit_gaussian$model)
  return(as.numeric(coef(out_gaussian)[,1]))
}


calibrating_pair_of_samples <- function(adt, insert_id, native_id, chr, doPlot = T, adtColumns = 1:5) {
  print(native_id)
  print(insert_id)
  if(doPlot) {
    real_insert <- log2(adt[seqnames!=chr][[insert_id]] + 1)
    diploid <- log2(adt[seqnames!=chr][[native_id]] + 1)
    fz <- sum(diploid == 0)
    rz <- sum(real_insert == 0)
    
    par(mfrow=c(1,1))
    plot(diploid, real_insert, cex.axis = 1.5, cex.lab = 1.2, lwd=3,
         main = sprintf("Calibration - Control (|zeros| = %d) - %s\naberration (|zeros| = %d) - %s", fz, native_id, rz, insert_id), 
         xlab ="Diploid sample (log2(TPM + 1))", 
         ylab ="Aberrated sample (log2(TPM + 1))", col = alpha("lightblue", alpha = 0.85))
    abline(0,1, col="red", lty = "dashed", lwd = 2)
    abline(lm(real_insert ~ diploid), col = "blue", lwd = 4)
    
    ids <- union(which(diploid == 0), which(real_insert == 0))
    diploid <- diploid[-ids]
    real_insert <- real_insert[-ids]
    points(diploid, real_insert, col = alpha("lightgreen", alpha = 0.5), lwd=1)
    abline(0,1, col="red", lty = "dashed", lwd = 2)
    fit <- lm(real_insert ~ diploid)
    abline(fit, col = "darkgreen", lwd = 4)
    prd <- predict(fit, newdata=data.frame(real_insert = real_insert),interval = c("confidence"), level = 0.90,type="response")
    lines(diploid, prd[,2],col = alpha(colour = "darkgreen", alpha = 0.6), lwd = 3, lty=2)
    lines(diploid, prd[,3],col = alpha(colour = "darkgreen", alpha = 0.6), lwd = 3, lty=2)
  }
  #fitting
  real_insert <- log2(adt[seqnames!=chr][[insert_id]] + 1)
  diploid <- log2(adt[seqnames!=chr][[native_id]] + 1)
  ids <- union(which(diploid == 0), which(real_insert == 0))
  diploid <- diploid[-ids]
  real_insert <- real_insert[-ids]
  fit <- lm(real_insert ~ diploid)
  
  if(doPlot) {
    real_insert.test <- 1/fit$coefficients[2]*(-fit$coefficients[1] + real_insert)
    abline(lm(real_insert.test ~ diploid), col = alpha(colour = "pink", alpha = 0.7), lwd = 7)
  }
  # Actual calibration
  real_insert.all <- log2(adt[[insert_id]] + 1)
  diploid.all <- log2(adt[[native_id]] + 1)
  real_insert.all2 <- 1/fit$coefficients[2]*(-fit$coefficients[1] + real_insert.all)
  real_insert.all2[real_insert.all2 < 0] <- 0
  
  if(doPlot) {
    ids <- union(which(diploid.all == 0), which(real_insert.all2 == 0))
    diploid <- diploid.all[-ids]
    real_insert2 <- real_insert.all2[-ids]
    abline(lm(real_insert2 ~ diploid), col = "magenta", lwd = 4)
    legend("topleft", inset=.05, title=NULL,
           c("Diploid ~ Aberrated (zeros excluded)", "Diploid ~ Aberrated (including zeros)", 
             "Diploid ~ Aberrated (Following Calibration) Excluding zeros", "x = y", "Calibration test"), 
           fill=c("darkgreen", "blue", "magenta", "red", alpha(colour = "pink", alpha = 0.7)), horiz=F)
  }
  
  y <- 2^real_insert.all2 - 1
  tt <- data.table(adt[, adtColumns, with = F], y = ifelse(test = y > 0, yes = y, no = 0))
  
  return(tt)
}




#input of Ys should be TPM and any other- no log!
#Both Ys and geneSpecific are lists of vectors
#e.g. inputYasMatrix <- list(SSM=T)
generateFeaturesForExpressionVector <- function(Ys, global = NA, geneSpecific = NA, ftrsWindowSize, by = 1, inputYasMatrix = NULL, isLogged = T) {
  #check if window size is even:
  if((ftrsWindowSize %% 2) == 0) {
    stop("ftrsWindowSize is even.")
  }
  for(n in names(Ys)) {
    Y <- Ys[[n]]
    
    if(is.null(inputYasMatrix[[n]])) {
      if(ftrsWindowSize > length(Y)) {
        stop("ftrsWindowSize of ", n, " is larger than vector length.")
      }
    } else {
      if(ftrsWindowSize > ncol(Y)) {
        stop("ftrsWindowSize of ", n, " is larger than vector length.")
      }
    }
  }
  midPoint <- (ftrsWindowSize - 1)/2 + 1
  allftrs <- list()
  for(n in names(Ys)) {
    Y <- Ys[[n]]
    message("Working on ", n, " ftrs..")
    if(!is.null(inputYasMatrix[[n]])) {
      message("Using ", n, " ftrs as matrix.")
      message("Number of NAs: ", sum(is.na(Y)), " out of: ", length(Y))
      exwin <- round((ncol(Y) - ftrsWindowSize)/2)
      idxs <- do.call(rbind, wapply(Y = 1:nrow(Y), width = ftrsWindowSize, by = by))[, midPoint]
      ftrs <- Y[idxs, (exwin + 1):(ftrsWindowSize + exwin) ]
      
    } else {
      message("Number of NAs: ", sum(is.na(Y)))
      ftrs <- do.call(rbind, wapply(Y = Y, width = ftrsWindowSize, by = by))
      idxs <- do.call(rbind, wapply(Y = 1:length(Y), width = ftrsWindowSize, by = by))[, midPoint]
    }
    print(dim(ftrs))
    rownames(ftrs) <- idxs
    colnames(ftrs) <- c(sprintf("N%d", (midPoint-1):1), "M0", sprintf("P%d", 1:(midPoint-1))) #paste(n, c(sprintf("N%d", (midPoint-1):1), "0", sprintf("P%d", 1:(midPoint-1))), sep = "_")
    
    Mean <- rowMeans(as.matrix(ftrs), na.rm = T)
    Var <- rowVars(as.matrix(ftrs), na.rm = T)
    #print(ftrs)
    quant50 <- matrix(apply(ftrs, 1, quantile, 0.5, na.rm = T))
    quant25 <- matrix(apply(ftrs, 1, quantile, 0.25, na.rm = T))
    quant75 <- matrix(apply(ftrs, 1, quantile, 0.75, na.rm = T))
    if(isLogged) {
      midRatio <- matrix(apply(ftrs, 1, function(v) log2((sum(2^v[1:(midPoint - 1)]) + 0.1)/(sum(2^v[(midPoint + 1):ncol(ftrs)]) + 0.1))))
    } else {
      midRatio <- matrix(apply(ftrs, 1, function(v) log2((sum(v[1:(midPoint - 1)]) + 0.1)/(sum(v[(midPoint + 1):ncol(ftrs)]) + 0.1))))
    }
    ftrs <- data.frame(ftrs, Mean = Mean, Var = Var, quant50 = quant50, quant25 = quant25, midRatio = midRatio)
    #colnames(ftrs)[(ncol(ftrs) - 6 + 1):ncol(ftrs)] <- paste(n, colnames(ftrs)[(ncol(ftrs) - 6 + 1):ncol(ftrs)], sep = "_")
    allftrs[[n]] <- ftrs
  }
  mftrs <- do.call(cbind, allftrs)
  
  if(!is.na(geneSpecific[1])) {
    #gftrs <- do.call(rbind, wapply(Y = geneSpecific, width = ftrsWindowSize, by = by))
    #colnames(gftrs) <- c(sprintf("N%d", (midPoint-1):1), "M0", sprintf("P%d", 1:(midPoint-1)))
    allgftrs <- list()
    for(n in names(geneSpecific)) {
      Y <- geneSpecific[[n]]
      message("Number of NAs: ", sum(is.na(Y)), " out of: ", length(Y))
      ftrs <- do.call(rbind, wapply(Y = Y, width = ftrsWindowSize, by = by))
      colnames(ftrs) <- c(sprintf("N%d", (midPoint-1):1), "M0", sprintf("P%d", 1:(midPoint-1)))
      
      Mean <- rowMeans(as.matrix(ftrs), na.rm = T)
      Var <- rowVars(as.matrix(ftrs), na.rm = T)
      #print(ftrs)
      quant50 <- matrix(apply(ftrs, 1, quantile, 0.5, na.rm = T))
      quant25 <- matrix(apply(ftrs, 1, quantile, 0.25, na.rm = T))
      quant75 <- matrix(apply(ftrs, 1, quantile, 0.75, na.rm = T))
      
      ftrs <- data.frame(ftrs, Mean = Mean, Var = Var, quant50 = quant50, quant25 = quant25)
      allgftrs[[n]] <- ftrs
    }
    gftrs <- do.call(cbind, allgftrs)
    
  }
  #return(ftrs)
  if(is.na(global[1])) {
    
    if(!is.na(geneSpecific[1])) {
      ftrs <- data.table(check.names = F,
                         idxs = as.numeric(rownames(ftrs)),
                         mftrs, 
                         geneSpecific = gftrs)
    } else {
      ftrs <- data.table(check.names = F,
                         idxs = as.numeric(rownames(ftrs)),
                         mftrs)
    }
  } else {
    if(!is.na(geneSpecific[1])) {
      ftrs <- data.table(check.names = F,
                         idxs = as.numeric(rownames(ftrs)), 
                         mftrs, 
                         geneSpecific = gftrs,
                         global.Mean = mean(global), 
                         global.quant25 = quantile(global, 0.25, na.rm = T), 
                         global.quant50 = quantile(global, 0.5, na.rm = T),  
                         global.quant75 = quantile(global, 0.75, na.rm = T))
    } else {
      ftrs <- data.table(check.names = F,
                         idxs = as.numeric(rownames(ftrs)),
                         mftrs, 
                         global.Mean = mean(global), 
                         global.quant25 = quantile(global, 0.25, na.rm = T), 
                         global.quant50 = quantile(global, 0.5, na.rm = T),  
                         global.quant75 = quantile(global, 0.75, na.rm = T))
    }
  }
  return(ftrs)
}

generateFeaturesForExpressionVector.old <- function(Ys, global = NA, geneSpecific = NA, ftrsWindowSize, by = 1, isLogged = T) {
  #check if window size is even:
  if((ftrsWindowSize %% 2) == 0) {
    stop("ftrsWindowSize is even.")
  }
  for(n in names(Ys)) {
    Y <- Ys[[n]]
    if(ftrsWindowSize > length(Y)) {
      stop("ftrsWindowSize of ", n, " is larger than vector length.")
    }
  }
  midPoint <- (ftrsWindowSize - 1)/2 + 1
  allftrs <- list()
  for(n in names(Ys)) {
    Y <- Ys[[n]]
    ftrs <- do.call(rbind, wapply(Y = Y, width = ftrsWindowSize, by = by))
    
    idxs <- do.call(rbind, wapply(Y = 1:length(Y), width = ftrsWindowSize, by = by))[, midPoint]
    rownames(ftrs) <- idxs
    colnames(ftrs) <- c(sprintf("N%d", (midPoint-1):1), "M0", sprintf("P%d", 1:(midPoint-1))) #paste(n, c(sprintf("N%d", (midPoint-1):1), "0", sprintf("P%d", 1:(midPoint-1))), sep = "_")
    
    Mean <- rowMeans(as.matrix(ftrs), na.rm = T)
    Var <- rowVars(as.matrix(ftrs), na.rm = T)
    #print(ftrs)
    quant50 <- matrix(apply(ftrs, 1, quantile, 0.5))#, na.rm = T))
    quant25 <- matrix(apply(ftrs, 1, quantile, 0.25))#, na.rm = T))
    quant75 <- matrix(apply(ftrs, 1, quantile, 0.75))#, na.rm = T))
    if(isLogged) {
      midRatio <- matrix(apply(ftrs, 1, function(v) log2((sum(2^v[1:(midPoint - 1)]) + 0.1)/(sum(2^v[(midPoint + 1):ncol(ftrs)]) + 0.1))))
    } else {
      midRatio <- matrix(apply(ftrs, 1, function(v) log2((sum(v[1:(midPoint - 1)]) + 0.1)/(sum(v[(midPoint + 1):ncol(ftrs)]) + 0.1))))
    }
    ftrs <- data.frame(ftrs, Mean = Mean, Var = Var, quant50 = quant50, quant25 = quant25, midRatio = midRatio)
    #colnames(ftrs)[(ncol(ftrs) - 6 + 1):ncol(ftrs)] <- paste(n, colnames(ftrs)[(ncol(ftrs) - 6 + 1):ncol(ftrs)], sep = "_")
    allftrs[[n]] <- ftrs
  }
  mftrs <- do.call(cbind, allftrs)
  
  if(!is.na(geneSpecific[1])) {
    #gftrs <- do.call(rbind, wapply(Y = geneSpecific, width = ftrsWindowSize, by = by))
    #colnames(gftrs) <- c(sprintf("N%d", (midPoint-1):1), "M0", sprintf("P%d", 1:(midPoint-1)))
    allgftrs <- list()
    for(n in names(geneSpecific)) {
      Y <- geneSpecific[[n]]
      ftrs <- do.call(rbind, wapply(Y = Y, width = ftrsWindowSize, by = by))
      colnames(ftrs) <- c(sprintf("N%d", (midPoint-1):1), "M0", sprintf("P%d", 1:(midPoint-1)))
      
      Mean <- rowMeans(as.matrix(ftrs), na.rm = T)
      Var <- rowVars(as.matrix(ftrs), na.rm = T)
      #print(ftrs)
      quant50 <- matrix(apply(ftrs, 1, quantile, 0.5))#, na.rm = T))
      quant25 <- matrix(apply(ftrs, 1, quantile, 0.25))#, na.rm = T))
      quant75 <- matrix(apply(ftrs, 1, quantile, 0.75))#, na.rm = T))
      
      ftrs <- data.frame(ftrs, Mean = Mean, Var = Var, quant50 = quant50, quant25 = quant25)
      allgftrs[[n]] <- ftrs
    }
    gftrs <- do.call(cbind, allgftrs)
    
  }
  #return(ftrs)
  if(is.na(global[1])) {
   
    if(!is.na(geneSpecific[1])) {
      ftrs <- data.table(check.names = F,
                         idxs = as.numeric(rownames(ftrs)),
                         mftrs, 
                         geneSpecific = gftrs)
    } else {
      ftrs <- data.table(check.names = F,
                         idxs = as.numeric(rownames(ftrs)),
                         mftrs)
    }
  } else {
    if(!is.na(geneSpecific[1])) {
      ftrs <- data.table(check.names = F,
                         idxs = as.numeric(rownames(ftrs)), 
                         mftrs, 
                         geneSpecific = gftrs,
                         global.Mean = mean(global), 
                         global.quant25 = quantile(global, 0.25), 
                         global.quant50 = quantile(global, 0.5),  
                         global.quant75 = quantile(global, 0.75))
    } else {
      ftrs <- data.table(check.names = F,
                         idxs = as.numeric(rownames(ftrs)),
                         mftrs, 
                         global.Mean = mean(global), 
                         global.quant25 = quantile(global, 0.25), 
                         global.quant50 = quantile(global, 0.5),  
                         global.quant75 = quantile(global, 0.75))
    }
  }
  return(ftrs)
}


#taken from: https://www.r-bloggers.com/wapply-a-faster-but-less-functional-rollapply-for-vector-setups/
wapply <- function(Y, width, by = NULL)
{
  if (is.null(by)) by <- width
  
  lenX <- length(Y)
  SEQ1 <- seq(1, lenX - width + 1, by = by)
  SEQ2 <- lapply(SEQ1, function(x) Y[x:(x + width - 1)])
  return(SEQ2)
}


