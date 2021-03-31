#data_preparation.base_model.R

################################################
# PARAMETERS INITIZALIZATION / CONFIGURATION   #
################################################

semi_supervised_method_to_use <- "SSM"  #"SSM" #MA100, IMF

#monosomies and total loss:
#Restricted thresholds: (0.8 accuracy)
LOSS_TH <- 0.6
NORMAL_HIGH_TH <- 1.1
NORMAL_LOW_TH <- 0.9
GAIN_TH <- 1.4
#Premissive ths: (0.75 accuracy)
# LOSS_TH <- 0.7
# NORMAL_HIGH_TH <- 1.2
# NORMAL_LOW_TH <- 0.8
# GAIN_TH <- 1.3
#Hard ths: (0.75 accuracy)
# LOSS_TH <- 0.55
# NORMAL_HIGH_TH <- 1.05
# NORMAL_LOW_TH <- 0.95
# GAIN_TH <- 1.45


#NUM_OF_FEATURES <- 10 #accuracy of 0.67
#NUM_OF_FEATURES <- 25 #accuracy of 0.78
#NUM_OF_FEATURES <- 31 #accuracy of 0.83
NUM_OF_FEATURES <- 50 #accuracy of 0.88
#NUM_OF_FEATURES <- 100 #accuracy of 0.963
number_of_entries_per_instance_for_training <- 150 #an instance is a chromosome in a sample that is used to be sampled for training entries
downsample_training_frac <- 0.7

#We decide training factor based on the accuracy in the predictions of each copy number group:
#for example, when all training factors are 1, the accuracy of OLR is 0.91 for loss, 0.8 for normal and 0.74 for gain.
#when we increase the learning factor of the gain group to 1.22 we increase the prediction accuracy to 0.778
gain_training_factor <- 2.0
loss_training_factor <- 1.0
normal_training_factor <- 1.3


######################################
# BUILDING THE TRAINING AND TEST SET #
######################################

if(semi_supervised_method_to_use == "SSM") {
  #Let's use large cnv preds using SSM: 0.85
  samples_to_use <- intersect(samples_to_use, colnames(arms$SSM.TE))
  #arms <- arms$SSM.TE[!arm %in% arms_to_exclude, c("arm", samples_to_use), with=F]
} else if(semi_supervised_method_to_use == "MA100") {
  #Let's try MA for large cnv predcitions and see how better or worse we perform: 0.86
  samples_to_use <- intersect(samples_to_use, colnames(arms$MA100))
  #arms <- arms$MA100[!arm %in% arms_to_exclude, c("arm", samples_to_use), with=F]
} else if(semi_supervised_method_to_use == "IMF") {
  #Let's try IMF for large cnv predcitions and see how better or worse we perform: 0.80
  samples_to_use <- intersect(samples_to_use, colnames(arms$IMF))
  #arms <- arms$IMF[!arm %in% arms_to_exclude, c("arm", samples_to_use), with=F]
} else {
  stop("ERROR - no method to use")
}


#################################################
# CREATION OF A GOLDEN STANDARD SHARED TEST SET #
#################################################

#Let's create a golden standard test set that will be used to test SSM and MA
SSM <- arms$SSM.TE[!arm %in% arms_to_exclude, c("arm", samples_to_use), with=F]
MA <- arms$MA100[!arm %in% arms_to_exclude, c("arm", samples_to_use), with=F]
A <- arms$raw.A[!arm %in% arms_to_exclude, c("arm", samples_to_use), with=F]
B <- arms$raw.B[!arm %in% arms_to_exclude, c("arm", samples_to_use), with=F]

stopifnot(sum(colnames(MA) != colnames(SSM)) == 0 &
            sum(colnames(MA) != colnames(A)) == 0 &
            sum(colnames(MA) != colnames(B)) == 0 &
            sum(rownames(MA) != rownames(B)) == 0 &
            sum(rownames(MA) != rownames(A)) == 0 &
            sum(rownames(MA) != rownames(SSM)) == 0)



#monosomies and total loss:
which.loss <- which((A < 0.1 | B < 0.1) & SSM < 0.6 & MA < 0.6, arr.ind = TRUE)
#excluding the cases of arm columns 
if(sum(which.loss[,2] == 1) > 0) {
  message("There are: ", sum(which.loss[,2] == 1), " non relevant cases.")
  which.loss <- which.loss[-which(which.loss[,2] == 1),]
}
which.loss <- data.table(arm = SSM$arm[which.loss[,1]], sample_id = colnames(SSM)[which.loss[,2]])
dim(which.loss)
head(which.loss)
sort(table(which.loss$arm))
#normal:
which.normal <- which((A > 0.4 &  A < 0.6 & B > 0.4 & B < 0.6) & (SSM < 1.1 & SSM > 0.9) & (MA < 1.1 & MA > 0.9), arr.ind = TRUE)
if(sum(which.normal[,2] == 1) > 0) {
  message("There are: ", sum(which.normal[,2] == 1), " non relevant cases.")
  which.normal <- which.normal[-which(which.normal[,2] == 1),]
}
which.normal <- data.table(arm = SSM$arm[which.normal[,1]], sample_id = colnames(SSM)[which.normal[,2]])
dim(which.normal)
head(which.normal)
sort(table(which.normal$arm))
#Trisomies and additional gains:
which.gain <- which((A > 0.65 | B > 0.65) & SSM > 1.4 & MA > 1.4, arr.ind = TRUE)
if(sum(which.gain[,2] == 1) > 0) {
  message("There are: ", sum(which.gain[,2] == 1), " non relevant cases.")
  which.gain <- which.gain[-which(which.gain[,2] == 1),]
}
which.gain <- data.table(arm = SSM$arm[which.gain[,1]], sample_id = colnames(SSM)[which.gain[,2]])
dim(which.gain)
head(which.gain)
sort(table(which.gain$arm))
#let's keep only cases with 1 or 2 aneuploides at most to avoid biases:
too_many_aneuploidies <- names(which(table(c(which.loss$sample_id, which.gain$sample_id))>2))
length(too_many_aneuploidies)
golden.which.gain <- copy(which.gain[!sample_id %in% too_many_aneuploidies, ])
golden.which.loss <- copy(which.loss[!sample_id %in% too_many_aneuploidies, ])
golden.which.normal <- copy(which.normal[!sample_id %in% too_many_aneuploidies, ])
message("Total of: ", nrow(golden.which.gain), " gains, ", nrow(golden.which.loss), " losses, and ", nrow(golden.which.normal), " normals.")



##################################################
# CREATION OF A TRAINING AND TEST SET PER METHOD #
##################################################

if(semi_supervised_method_to_use == "SSM") {
  #BY SSM: for DNN with size = 25 we get 0.785 accuracy
  arms <- arms$SSM.TE[!arm %in% arms_to_exclude, c("arm", samples_to_use), with=F]
} else if(semi_supervised_method_to_use == "MA100") {
  #BY MA100: for NN with size = 25 we get 0.782 accuracy
  arms <- arms$MA100[!arm %in% arms_to_exclude, c("arm", samples_to_use), with=F]
} else if(semi_supervised_method_to_use == "IMF") {
  arms <- arms$IMF[!arm %in% arms_to_exclude, c("arm", samples_to_use), with=F]
}
dim(arms)

  
which.loss <- which(arms < LOSS_TH, arr.ind = TRUE)
if(sum(which.loss[,2] == 1) > 0)
  which.loss <- which.loss[-which(which.loss[,2] == 1),]
which.loss <- data.table(arm = arms$arm[which.loss[,1]], sample_id = colnames(arms)[which.loss[,2]])
dim(which.loss)
#normal:
which.normal <- which(arms < NORMAL_HIGH_TH & arms > NORMAL_LOW_TH, arr.ind = TRUE)
if(sum(which.normal[,2] == 1) > 0)
  which.normal <- which.normal[-which(which.normal[,2] == 1),]
which.normal <- data.table(arm = arms$arm[which.normal[,1]], sample_id = colnames(arms)[which.normal[,2]])
dim(which.normal)
#Trisomies and additional gains:
which.gain <- which(arms > GAIN_TH, arr.ind = TRUE)
if(sum(which.gain[,2] == 1) > 0)
  which.gain <- which.gain[-which(which.gain[,2] == 1),]
which.gain <- data.table(arm = arms$arm[which.gain[,1]], sample_id = colnames(arms)[which.gain[,2]])
dim(which.gain)
#let's keep only cases with 1 or 2 aneuploides at most to avoid biases:
too_many_aneuploidies <- names(which(table(c(which.loss$sample_id, which.gain$sample_id))>2))
which.gain <- which.gain[!sample_id %in% too_many_aneuploidies, ]
which.loss <- which.loss[!sample_id %in% too_many_aneuploidies, ]
which.normal <- which.normal[!sample_id %in% too_many_aneuploidies, ]

#There are more losses in the golden standard set because there are less aneuploidies in total 
#Since there are less aneuploidies, lower number of samples are excluded due to "too_many_aneuoploides"
message("GOLDEN - Total of: ", nrow(golden.which.gain), " gains, ", nrow(golden.which.loss), " losses, and ", nrow(golden.which.normal), " normals.")
message("METHOD SPECIFIC - Total of: ", nrow(which.gain), " gains, ", nrow(which.loss), " losses, and ", nrow(which.normal), " normals.")

#let's divide the samples to training and test sets where the test set includes only golden standard samples in order to be able to compare
#between the explicit methods used for labeling:
golden.samples.all <- unique(c(golden.which.gain$sample_id, golden.which.loss$sample_id, golden.which.normal$sample_id))
length(golden.samples.all)
samples.all <- unique(c(which.gain$sample_id, which.loss$sample_id, which.normal$sample_id, golden.samples.all))

length(samples.all)
samples.train <- unique(c(sample(unique(which.gain$sample_id), size = length(unique(which.gain$sample_id))/2),
                          sample(unique(which.loss$sample_id), size = length(unique(which.loss$sample_id))/2),
                          sample(unique(which.normal$sample_id), size = length(unique(which.normal$sample_id))/2)))
length(samples.train)
#for the test sample group we want it to include samples that are not in the training set and only in the golden standard set
samples.test.bymethod <- samples.all[!samples.all %in% samples.train]
message("When not restricting the set to be from the golden standard only, we get ", length(samples.test.bymethod), " samples.")
samples.test <- samples.all[(!samples.all %in% samples.train) & (samples.all %in% golden.samples.all)]
message("When we restrict the set to be from the golden standard, we get ", length(samples.test), " samples.")
message("The training set includes ", length(samples.train), " samples.")

#QC feature:
QC <- round(scale((colSums(rsemtpm[, samples.all] > 0)), center = T, scale = T), digits = 0)[,1]

#########################################################
# Normalization of expression before feature generation #
#########################################################
#This is the chosen expression normalization. For the other methods studied see the script: explore_models_and_data_preprocessing.R
m <- as.matrix(adt.default[,-c(1:4)])
agg <- adt.default[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), by="seqnames", .SDcols = -c(1:4)] 
require(matrixStats)
m <- sweep(x = m, MARGIN = 2, STATS = colMedians(as.matrix(agg[,-1])), FUN = "-")
adt <- cbind(adt.default[, 1:4], m)


####################################
# GENERATING FEATURE DATA          #
####################################
#we will generate feature data based on the number of training samples for each group such that it is balanced.
#the strategy is to sample more entries for training from the smaller group and less entries from the larger groups
#loss and gain have less entries than normal:
max_grp_size <- max(nrow(which.gain), nrow(which.loss))
gain_N_entries <- min(500, (max_grp_size/nrow(which.gain))*number_of_entries_per_instance_for_training)*gain_training_factor
loss_N_entries <- min(500, (max_grp_size/nrow(which.loss))*number_of_entries_per_instance_for_training)*loss_training_factor
normal_N_entries <- max(2, (max_grp_size/nrow(which.normal))*(loss_N_entries + gain_N_entries)*2)*normal_training_factor
message("After balancing strategy we will sample, ", gain_N_entries, " from the gain instances, ", loss_N_entries, " from the loss instances and ",
        normal_N_entries, " from the normal instances.")

#let's generate the features data for training and testing:
get_random_seq <- function(sample_id, arm, size = 20, N = 10, dt, ganno, label = "gain", addQC = F, verbose = F) {
  mydt <- dt[which(ganno$arm == arm), sample_id, with = F][[1]]
  if(verbose)
    message("Working on: ", arm, ", ", sample_id, " with ", length(mydt), " genes")
  if((length(mydt) - size + 1) < N) {
    newN <- (length(mydt) - size + 1)
    if(verbose)
      message(" Downgrading to ", newN, " instead of ", N)
    if(newN < 1) {
      message("segement length is too long to be sampled from this arm.")
      return(NA)
    }
    N <- newN
  }
  rnds <- sample(x = 1:(length(mydt) - size + 1), size = N)
  rnd_segs <- sapply(rnds, function(x) mydt[x:(x + size - 1)]) 
  #print(dim(rnd_segs))
  if(addQC) {
    res <- data.table(sample_id = sample_id, arm = arm, start_idx = rnds, end_idx = rnds + size - 1, t(rnd_segs), QC = QC[sample_id], label = label)
  } else {
    res <- data.table(sample_id = sample_id, arm = arm, start_idx = rnds, end_idx = rnds + size - 1, t(rnd_segs), label = label)
  }
  return(res)
}
################
# Training set #
################

SIZE <- NUM_OF_FEATURES
number_of_entries_per_sample <- number_of_entries_per_instance_for_training

#gain
ll <- apply(which.gain[sample_id %in% samples.train, ], 1, function(xx) { 
  get_random_seq(sample_id = xx[2], arm = xx[1], dt = adt, ganno = ganno, size = SIZE, N = gain_N_entries, label = "gain")
})
train_set <- rbindlist(ll[!is.na(ll)])
#loss
ll <- apply(which.loss[sample_id %in% samples.train, ], 1, function(xx) { 
  get_random_seq(sample_id = xx[2], arm = xx[1], dt = adt, ganno = ganno, size = SIZE, N = loss_N_entries, label = "loss")
})
train_set <- rbindlist(list(train_set, rbindlist(ll[!is.na(ll)])))
table(train_set$label)
#normal
normal_sampled <- sample(intersect(which.normal$sample_id, samples.train), length(unique(train_set$sample_id)))
ll <- apply(which.normal[sample_id %in% normal_sampled, ], 1, function(xx) { 
  get_random_seq(sample_id = xx[2], arm = xx[1], dt = adt, ganno = ganno, size = SIZE, N = normal_N_entries, label = "normal")
})
#since each sample has multiple arms with normal expression we shuffle the normal entries of the training set to sample a variety of sample_ids and not only a couple:
train_set <- rbindlist(list(train_set, head(sample_frac(rbindlist(ll[!is.na(ll)]), size = 1), n = nrow(train_set))))
table(train_set$label)
#let's shuffle the dataset before training:
train_set <- sample_frac(train_set, size = 1)

################
# Test set     #
################
#gain
ll <- apply(which.gain[sample_id %in% samples.test, ], 1, function(xx) { 
  get_random_seq(sample_id = xx[2], arm = xx[1], dt = adt, ganno = ganno, size = SIZE, N = gain_N_entries, label = "gain")
})

test_set <- rbindlist(ll[!is.na(ll)])
#loss
ll <- apply(which.loss[sample_id %in% samples.test, ], 1, function(xx) { 
  get_random_seq(sample_id = xx[2], arm = xx[1], dt = adt, ganno = ganno, size = SIZE, N = loss_N_entries, label = "loss")
})
test_set <- rbindlist(list(test_set, rbindlist(ll[!is.na(ll)])))
table(test_set$label)
#normal
normal_sampled <- sample(intersect(which.normal$sample_id, samples.test), length(unique(test_set$sample_id)))
ll <- apply(which.normal[sample_id %in% normal_sampled, ], 1, function(xx) { 
  get_random_seq(sample_id = xx[2], arm = xx[1], dt = adt, ganno = ganno, size = SIZE, N = normal_N_entries, label = "normal")
})

#since each sample has multiple arms with normal expression we shuffle the normal entries of the testing set to sample a variety of sample_ids and not only a couple:
test_set <- rbindlist(list(test_set, head(sample_n(rbindlist(ll[!is.na(ll)]), size = nrow(rbindlist(ll[!is.na(ll)]))), n = nrow(test_set))))
table(test_set$label)

#let's shuffle the dataset before training:
test_set <- as.data.table(sample_frac(test_set, size = 1))

#QC - let's make sure we don't have the same sample ids in train and test:
stopifnot(length(intersect(unique(test_set$sample_id), unique(train_set$sample_id))) == 0)
levels <- unique(c(train_set$label, test_set$label))
levels <- c("loss", "normal", "gain")

#########################################
##### Define training and test sets #####
#########################################

x_train <- train_set[, -c(1:4, ncol(train_set)), with=F]
y_train <- train_set$label
x_train <- as.matrix(x_train) #dense_layer input
dim(x_train)
y_train <- to_categorical(as.numeric(factor(y_train, levels = levels)) - 1)
dim(y_train)

x_test <- test_set[, -c(1:4, ncol(test_set)), with=F]
y_test <- test_set$label
x_test <- as.matrix(x_test) #dense_layer input
dim(x_test)
y_test <- to_categorical(as.numeric(factor(y_test, levels = levels)) - 1)
dim(y_test)

