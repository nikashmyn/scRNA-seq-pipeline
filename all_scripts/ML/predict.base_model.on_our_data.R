#predict.base_model.R

##############
# PARAMETERS #
##############
normalize_expression <- T
NUM_OF_FEATURES <- 50 

base_model_fname <- sprintf("%s/NN/base_model_v1.WIN%d.tf", dirpath, NUM_OF_FEATURES)
base_model_avg_fname <- sprintf("%s/NN/base_model_v1.AVG%d.tf", dirpath, NUM_OF_FEATURES)
olr_model_fname <- sprintf("%s/NN/olr_model_v1.AVG%d.rds", dirpath, NUM_OF_FEATURES)

#loading base model and data for predictions:
base_model <- load_model_tf(filepath = base_model_fname)
base_avg_model <- load_model_tf(filepath = base_model_avg_fname)
olr_model <- readRDS(olr_model_fname)

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
                                          data.table(rollapply(data = myadt[seqnames == chr][[sample_id]], 
                                                               width = winSize, FUN = identity, fill = NA, align = "center")))
  ftrs <- cbind(myadt[,1:4], rbindlist(ftrs))
  return(ftrs)
}

model <- base_model
#some validation on known cases:
ftrs <- generate_features(adt = adt, sample_id = "170512_B8", winSize = NUM_OF_FEATURES)
mypreds <- predict(object = model, x = as.matrix(ftrs[seqnames == "chr5"][, -c(1:4)]))
plot(max.col(mypreds))

## OLR
model <- olr_model
chr <- "chr5"
ftrs <- generate_features(adt = adt, sample_id = "170512_B4", winSize = NUM_OF_FEATURES)
mypreds <- predict(object = olr_model, data.table(x=as.numeric(rowMeans(as.matrix(ftrs[seqnames == chr][,-c(1:4)])))), type="p")
table(max.col(mypreds))
plot(max.col(mypreds))

## NN
model <- base_avg_model
mypreds <- predict(object = base_avg_model, rowMeans(as.matrix(ftrs[seqnames == chr][,-c(1:4)])), type="p")
table(max.col(mypreds))
#!!!!! Should this be mypreds? it was mypreds1.

plot(max.col(mypreds), ylim=c(0.5,3.5))
#points(mypreds[,1]+1, col="blue", cex = 0.5)
#points(mypreds[,2]+1, col="green", cex = 0.5)
points(rowMeans(as.matrix(ftrs[seqnames == chr][,-c(1:4)]))-min(rowMeans(as.matrix(ftrs[seqnames == chr][,-c(1:4)])), na.rm = T)+1, col="seagreen")
#mypreds <- predict(object = model, x = as.matrix(ftrs[seqnames == chr][, -c(1:4)]))
points(max.col(mypreds), col="purple")
points(mypreds[,1]+1, col="blue", cex = 0.5)
points(mypreds[,2]+1, col="green", cex = 0.5)
points(mypreds[,3]+1, col="red", cex = 0.5)
abline(h=1.67)
abline(h=1.33)

# I added this model change line. Before it was the NN model.
model <- base_model
# end change
ftrs <- generate_features(adt = adt, sample_id = controlSampleIDs[2])
mypreds <- predict(object = model, x = as.matrix(ftrs[seqnames == "chr2"][, -c(1:4)]))
plot(max.col(mypreds))

ftrs <- generate_features(adt = adt, sample_id = "170126_A3")
mypreds <- predict(object = model, x = as.matrix(ftrs[seqnames == "chr1"][, -c(1:4)]))
plot(max.col(mypreds))


