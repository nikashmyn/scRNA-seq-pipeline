
require(data.table)

#loading data used in this module:
#---------------------------------
#dbftrs <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/tmp/db.copy_number_preds_from_SSM_5_12_0.05_with_ftrs.zscored.winSize50.rds")
dbftrs <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/tmp/db.copy_number_preds_from_SSM_5_12_0.10b_with_ftrs.zscored.winSize50.rds")
#dbftrs <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/tmp/db.copy_number_preds_from_SSM_5_12_0.25_with_ftrs.zscored.winSize50.rds")
dball <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/tmp/db.all_samples_ftrs.zscored.winSize50.rds")
CN <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/CN.ssm.arms_5_12.rds")

adt <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/adt.ftrs.winSize50.raw.rds")
controlSampleIDs <- readRDS("/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/controlSampleIDs")
col_anno <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/col_anno.Apr122018.rds")
rsemtpm <- readRDS("/singlecellcenter/etai/tmp/CNML/pointOfChange/run_id17/data/Stamatis_list_v14_181025.rsemtpm.rds")

controlIDs <- col_anno[experimental_label %in% c("control_LCM_notreatment", "control_gen2_notreatment", "noco_gen2_noMN_lcm", "noco_gen1_noMN")]$WTA.plate
triIDs <- col_anno[experimental_label %in% c("Tri12_LCM", "Tri12_facs")]$WTA.plate[col_anno[experimental_label %in% c("Tri12_LCM", "Tri12_facs")]$WTA.plate %in% names(which(colSums(rsemtpm>0)>6000))]


source('~/WORK/Papers/MLpaper/R/utilities/performance_evaluation_of_CN_prediction.R')

meta_columns = c("sample_id", "seqnames", "contig", "size", "whole_chr_copy", "experimental_label", "id")
whole_chr_copy = T
class_column = "type"
experimental_label_to_exclude = c("bunch_of_grapes", "reinc_in_use")
chrs_to_exclude = names(which(table(adt$MEAN$seqnames)<=400)) #sprintf("chr%d", 16:22)
use_cluster = F
types = c(1, 2, 3)
precision = 0.05

#choise 1 for vars
all_vars <- names(dbftrs)[9:ncol(dbftrs)]


manual_selected_vars6 <- c("geneSpecific.MEDIAN.MEDIAN.GC",
                           "geneSpecific.MEDIAN.MEDIAN.cv.nz",
                           "geneSpecific.MEDIAN.MEDIAN.Length",
                           "geneSpecific.MEAN.MEAN.interspace",
                           "TE.all.SD",                     
                           "TE.onlyExpressed.MAXSTRETCH",      
                           "Frac.nonzeros.MAXSTRETCH",         
                           "TE.all.Q95",                       
                           "TE.onlyExpressed.Q95",             
                           "Frac.nonzeros.SD",                 
                           "TE.onlyExpressed.Q25",             
                           "TE.onlyExpressed.Q75",             
                           "TE.all.Q05",                       
                           "Frac.nonzeros.MEAN",               
                           "TE.onlyExpressed.MEDIAN",         
                           "TE.onlyExpressed.MEAN",           
                           "TE.all.Q25",                      
                           "TE.all.Q75",                      
                           "TE.all.MEDIAN",                   
                           "TE.all.MEAN")                     

#decide which vars to take:
variable_names <- manual_selected_vars6



#taking disomies only from controls
mat <- na.omit(dbftrs[whole_chr_copy == T & type %in% types &
                        ((sample_id %in% controlSampleIDs & type == 2) | type != 2) &
                        !(experimental_label %in% experimental_label_to_exclude) &
                        experimental_label != "Reinc_in_use" &
                        !(seqnames %in% chrs_to_exclude),
                      c("sample_id", "seqnames", "size", class_column, variable_names), with = F])

fitmat <- mat

dim(fitmat)
table(fitmat$type)


#light gbm model:
#----------------

source('~/WORK/Papers/MLpaper/R/utilities/train_using_lightgbm.R')


#dividing to train and test set:
test_data <- NULL
train_data <- NULL

for(i in unique(mat$type)) {
  sample_names <- unique(mat[type == i]$sample_id)
  N <- length(sample_names)
  idxs <- sample(x = N, size = round(N*0.3))
  test_data <- rbind(test_data, mat[type == i & sample_id %in% sample_names[idxs]])
  train_data <- rbind(train_data, mat[type == i & sample_id %in% sample_names[-idxs]])
}

table(train_data$type)
table(test_data$type)

use_double_size = T

#equal class size
if(use_double_size == F) {
  test_data <- test_data[,.SD[sample(.N, max(table(test_data$type)), replace = ifelse(max(table(test_data$type)) > .N, T, F))],by = type]
  train_data <- train_data[,.SD[sample(.N, max(table(train_data$type)), replace = ifelse(max(table(train_data$type)) > .N, T, F))],by = type]
} else {
#double size for 2 copies
  test_data <- test_data[,.SD[sample(.N, ifelse(.N == max(table(test_data$type)), .N, max(table(test_data$type))/2), 
                                     replace = ifelse(max(table(test_data$type))/2 > .N, T, F))],by = type]
  train_data <- train_data[,.SD[sample(.N, ifelse(.N == max(table(train_data$type)), .N, max(table(train_data$type))/2), 
                                       replace = ifelse(max(table(train_data$type))/2 > .N, T, F))],by = type]
}
table(train_data$type)
table(test_data$type)


fit.lm <- lm(formula = type ~ ., data = train_data[,-c(2:4)])
boxplot(fitted(fit.lm)~train_data$type, notch=F)

#train - classification
train_using_lightgbm(df.train = train_data[, -c(1:4)], train_labels = train_data$type,
                     df.tes = test_data[, -c(1:4)], test_labels = test_data$type, 
                     model_type = "classification", 
                     bst.fname = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/tmp/bst.basic.vars6.SSM_0.10b.classification.rds",
                     report.fname = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/tmp/bst.basic.vars6.SSM_0.10b.classification.report.pdf",)

#predict:

#all data:
all_data <- na.omit(dball[, c(names(dball)[1:6], variable_names), with=F])
all_meta <- all_data[, 1:6]
all_data <- as.matrix(all_data[, -c(1:6)])
dim(all_data)

source('~/WORK/Papers/MLpaper/R/utilities/predict_using_lightgbm.R')

predictions <- predict_using_lightgbm(new_data = all_data, 
                                      bst.fname = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/tmp/bst.basic.vars6.SSM_0.10b.classification.rds")
colnames(predictions) <- c("c1", "c2", "c3")
lgbmpreds5 <- cbind(all_meta, predictions)
saveRDS(lgbmpreds5, file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/tmp/lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")


#train - regression
train_using_lightgbm(df.train = train_data[, -c(1:4)], train_labels = train_data$type,
                     df.tes = test_data[, -c(1:4)], test_labels = test_data$type, 
                     model_type = "regression", 
                     bst.fname = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/tmp/bst.basic.vars6.SSM_0.10b.regression.l12.rds",
                     report.fname = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/tmp/bst.basic.vars6.SSM_0.10b.regression.report.l12.pdf",)


#predict:

#all data:
all_data <- na.omit(dball[, c(names(dball)[1:6], variable_names), with=F])
all_meta <- all_data[, 1:6]
all_data <- as.matrix(all_data[, -c(1:6)])
dim(all_data)

source('~/WORK/Papers/MLpaper/R/utilities/predict_using_lightgbm.R')

#regression mode
predictions <- predict_using_lightgbm(new_data = all_data, 
                                      bst.fname = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/tmp/bst.basic.vars6.SSM_0.10b.regression.rds")


lgbmpreds6 <- cbind(all_meta, pred = predictions + 1)
saveRDS(lgbmpreds6, file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/tmp/lgbmpreds.bst.basic.vars6.SSM_0.10b.regression.rds")




