##############
# PARAMETERS #
##############
dstdir <- sprintf("%s/ML_data", dirpath)

normalize_expression <- T
NUM_OF_FEATURES <- 100 #50

base_model_fname <- sprintf("%s/NN/base_model_v1.WIN%d.tf", dirpath, NUM_OF_FEATURES)
base_model_avg_fname <- sprintf("%s/NN/base_model_v1.AVG%d.tf", dirpath, NUM_OF_FEATURES)
olr_model_fname <- sprintf("%s/NN/olr_model_v1.AVG%d.rds", dirpath, NUM_OF_FEATURES)

olr_model <- readRDS(olr_model_fname)
#olr_model <- load_model_tf(base_model_fname)

arms <- copy(ganno) # <- readRDS(file = sprintf("%s/aggregated_results/ganno.rds", dirpath))

#########################################################
# Normalization of expression before feature generation #
#########################################################
#This is the chosen expression normalization. 
#For the other methods studied see the script: explore_models_and_data_preprocessing.R
if(normalize_expression == T) {
  m <- as.matrix(adt.default[,-c(1:4)])
  agg <- adt.default[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), by="seqnames", .SDcols = -c(1:4)] 
  require(matrixStats)
  m <- sweep(x = m, MARGIN = 2, STATS = colMedians(as.matrix(agg[,-1])), FUN = "-")
  adt <- cbind(adt.default[, 1:4], m)
}

generate_features <- function(sample_id = "170512_B8", adt, winSize = NUM_OF_FEATURES) {
  myadt <- copy(adt[, c(colnames(adt)[1:4], sample_id), with=F])
  ftrs <- lapply(unique(myadt$seqnames), function(chr) 
                 data.table(rollapply(data = myadt[seqnames == chr][[sample_id]], partial = T, #added partial
                                      width = winSize, FUN = identity, fill = NA, align = "center")))
  ftrs <- cbind(myadt[,1:4], rbindlist(ftrs))
  return(ftrs)
}

#Using the OLR model:
model <- olr_model

#generate_features_and_run_integer_cn_predictions_for_multiple_samples <- function(adt, winSize = NUM_OF_FEATURES, FUN = "mean") {
#  ftrs <- adt[, lapply(.SD, rollapply, width = winSize, FUN = FUN, fill = NA, align = "center"),
#              .SDcols = names(adt)[-c(1:4)], by = seqnames]
#
#  tmp <- ftrs[, lapply(.SD, function(x) as.numeric(predict(object = model, data.table(x = x)))), .SDcols = names(ftrs)[-1]]
#  preds <- cbind(adt[,1:4], tmp)
#
#  return(preds)
#}

#prediction of CN and CN probabilities:
#--------------------------------------
generate_features_and_run_integer_cn_predictions_for_multiple_samples <- function(adt, winSize = NUM_OF_FEATURES, FUN = "mean") {
  ftrs <- adt[, lapply(.SD, rollapply, width = winSize, FUN = FUN, fill = NA, partial = T, align = "center"),  #added partial
              .SDcols = names(adt)[-c(1:4)], by = seqnames]
  
  tmp <- ftrs[, lapply(.SD, function(x) as.numeric(predict(object = model, data.table(x = x)))), .SDcols = names(ftrs)[-1]]
  preds <- cbind(adt[,1:4], tmp)
  
  return(preds)
}

ourcns <- generate_features_and_run_integer_cn_predictions_for_multiple_samples(adt = adt)

#changed FUN = FUN -> "mean" and width = winSize -> "mean"
ourftrs <- adt[, lapply(.SD, rollapply, width = NUM_OF_FEATURES, FUN = "mean", fill = NA, partial = T, align = "center"), #added partial
               .SDcols = names(adt)[-c(1:4)], by = seqnames]

ourprobs <- lapply(samples_to_use, function(i) predict(object = model, data.table(x = ourftrs[[i]]), type = "p"))
names(ourprobs) <- samples_to_use

ourpreds <- list(rowdata = adt[, 1:4], preds = ourprobs, cns = ourcns)
saveRDS(ourpreds, file = sprintf("%s/ML_data/preds.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES))  #"/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2/latest_model/preds.probs_olr_model.WIN50.rds")

plot(max.col(ourpreds$preds[[22]]))

preads <- readRDS(file = sprintf("%s/ML_data/preds.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES))
ourpreds <- readRDS(file = sprintf("%s/ML_data/preds.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES))

#--------------------------------------------
#calculate statistics for p-val calculations:
#--------------------------------------------

setkeyv(arms, c("seqnames", "start", "end"))
setkeyv(ourpreds$cns, c("seqnames", "start", "end"))
tmp <- foverlaps(arms, ourpreds$cns)
tmp <- tmp[, c("id", "seqnames", "start", "end", "arm", samples_to_use), with=F]

#Integer copy number stats:
#--------------------------
stats.cn.arm <- tmp[, lapply(.SD, function(x) mean(x, na.rm = T)), .SDcols = colnames(tmp)[-c(1:5)], by = c("arm")]
stats.cn.chr <- tmp[, lapply(.SD, function(x) mean(x, na.rm = T)), .SDcols = colnames(tmp)[-c(1:5)], by = c("seqnames")]

#Intermediate copy number stats:
#-------------------------------
#This part gives you the fraction of intermediate predictions (not clearly one CN state)
#out of all the predictions in a bin (bin = chr or arm)
maxprob <- do.call(cbind, lapply(ourpreds$preds, function(x) apply(x, 1, max)))
setkeyv(ourpreds$rowdata, c("seqnames", "start", "end", "id"))
setkeyv(arms, c("seqnames", "start", "end", "id"))
# I changed this to merge instead foverlap
tmp <- merge(ourpreds$rowdata, arms)
stopifnot(sum(tmp$id != ourpreds$rowdata$id) == 0)
tmp <- tmp[, c("id", "seqnames", "start", "end", "arm"), with=F]
maxprob <- cbind(tmp, maxprob)

#arm
stats.intermediate.arm <- maxprob[, lapply(.SD, function(x) mean(x < 0.67, na.rm = T)), .SDcols = colnames(maxprob)[-c(1:5)], by = c("arm")]
si <- melt(stats.intermediate.arm, id.vars="arm")
setnames(si, old = c("variable", "value"), new = c("sample_id", "newpred"))
setcolorder(x = si, neworder = c("sample_id", "arm", "newpred"))
setnames(si, old = "newpred", new = "intermediate_frac")
saveRDS(si, file = sprintf("%s/ML_data/interstat.armlevel.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES))

#chr
stats.intermediate.chr <- maxprob[, lapply(.SD, function(x) mean(x < 0.67, na.rm = T)), .SDcols = colnames(maxprob)[-c(1:5)], by = c("seqnames")]
si <- melt(stats.intermediate.chr, id.vars="seqnames")
setnames(si, old = c("variable", "value"), new = c("sample_id", "newpred"))
setcolorder(x = si, neworder = c("sample_id", "seqnames", "newpred"))
setnames(si, old = "newpred", new = "intermediate_frac")
saveRDS(si, file = sprintf("%s/ML_data/interstat.chrlevel.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES))

interstat.chr <- readRDS(file = sprintf("%s/ML_data/interstat.chrlevel.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES))
interstat.arm <- readRDS(file = sprintf("%s/ML_data/interstat.armlevel.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES))

#This script has the calculate_fraction... function used below
source( sprintf('%s/plots/plot_cn_and_intermediate_prediction_latest_model.R', scriptsdir) )

#fracstat: (the fraction of loci that are classified as each CN state)
tmp <- lapply(controlIDs, function(myid) {
  message(myid)
  rbindlist(lapply(unique(interstat.chr$seqnames), function(chr, myid) {
    data.table(sample_id = myid, 
               seqnames = chr, 
               t(as.vector(calculate_fraction_of_cns_for_chr(preads$cns[id == myid & seqnames == chr])))) }, myid)) 
})
tmp2 <- rbindlist(tmp)
setnames(tmp2, old = names(tmp2)[-c(1:2)], new = c("1","1.5", "2", "2.5", "3"))
require(entropy)
tmp2$entropy <- apply(tmp2[, -c(1:2)], 1, function(x) entropy(as.numeric(x)))
fracstat <- tmp2
saveRDS(fracstat, sprintf("%s/ML_data/fracstat.chrlevel.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES))
#fracstat <- readRDS(file = sprintf("%s/ML_data/fracstat.chrlevel.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES))

controlIDs <- Reduce(intersect, list(controlIDs, colnames(stats.cn.chr), colnames(stats.intermediate.arm)))

saveRDS(controlIDs, file = sprintf("%s/aggregated_results/reduced_controlIDs.rds", dirpath))

#---------------------
#End interstat section
#---------------------
