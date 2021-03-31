#training.base_model.R

#this script should be run after executing: data_preparation.base_model.R

##############
# PARAMETERS #
##############

base_model_fname <- sprintf("%s/NN/base_model_v1.WIN%d.tf", dirpath, NUM_OF_FEATURES)
base_model_avg_fname <- sprintf("%s/NN/base_model_v1.AVG%d.tf", dirpath, NUM_OF_FEATURES)
olr_model_fname <- sprintf("%s/NN/olr_model_v1.AVG%d.rds", dirpath, NUM_OF_FEATURES)

require(data.table)
require(readxl)
require(dplyr)


##############################
# Building a model - Training #
###############################

# Initialize model
#these model reaches 0.74 (with feature size = 25) and 0.84 (with feature size = 50) - consider increasing number of layer for a larger training set
model <- keras_model_sequential()
early_stopping <- callback_early_stopping(monitor ='loss', min_delta=0.0001)
#LINEAR REGRESSION MODEL+sigmoid: 0.887 (size = 50), 0.78 for size = 25
#this is the best performing and simplest model that was used out of several models explored. 
#See the script: explore_models_and_data_preprocessing.R

model %>% 
  layer_dense(units = 1, input_shape = dim(x_train)[2], activation='linear') %>%
  layer_dense(units = 3, activation = 'softmax') %>%
  #layer_dense(units = 3, activation = 'linear') %>% #doesn't do any learning - that proofs that it is a non linear relashionship
  summary(model)

model %>% compile(
  loss = 'categorical_crossentropy',
  #optimizer = optimizer_adadelta(lr = 0.2), #all optimizers have about the same performance of about 0.85 for size = 50
  optimizer = optimizer_adam(lr = 0.0005),
  #optimizer = optimizer_rmsprop(lr = 0.0005, rho = 0.95),
  metrics = c('accuracy')
)

# set global default to show metrics
options(keras.view_metrics = TRUE)

# Train model over four epochs
cat('Train...\n')
history <- model %>% fit(
  x_train, y_train,
  batch_size = 32,
  epochs = 64,
  callbacks = early_stopping,
  validation_data = list(x_test, y_test)
)
#Saving the base model (for predictions or transfer learning):
save_model_tf(object = model, filepath = base_model_fname)

#BUILDING MODEL BASED ON AVERAGE OF THE VALUES IN THE WINDOW:
x_train_avg <- rowMeans(x_train)
x_test_avg <- rowMeans(x_test)
# Initialize model
#these model reaches 0.74 (with feature size = 25) and 0.84 (with feature size = 50) - consider increasing number of layer for a larger training set
model_avg <- keras_model_sequential()
early_stopping <- callback_early_stopping(monitor ='loss', min_delta=0.0001)
#LINEAR REGRESSION MODEL+sigmoid: 0.887 (size = 50), 0.78 for size = 25
#this is the best performing and simplest model that was used out of several models explored. 
#See the script: explore_models_and_data_preprocessing.R

model_avg %>% 
  layer_dense(units = 1, input_shape = 1, activation='linear') %>%
  layer_dense(units = 3, activation = 'softmax') %>%
  #layer_dense(units = 3, activation = 'linear') %>% #doesn't do any learning - that proofs that it is a non linear relashionship
  summary(model_avg)

model_avg %>% compile(
  loss = 'categorical_crossentropy',
  #optimizer = optimizer_adadelta(lr = 0.2), #all optimizers have about the same performance of about 0.85 for size = 50
  optimizer = optimizer_adam(lr = 0.0005),
  #optimizer = optimizer_rmsprop(lr = 0.0005, rho = 0.95),
  metrics = c('accuracy')
)

# set global default to show metrics
options(keras.view_metrics = TRUE)

# Train model over four epochs
cat('Train...\n')
history <- model_avg %>% fit(
  x_train_avg, y_train,
  batch_size = 32,
  epochs = 64,
  callbacks = early_stopping,
  validation_data = list(x_test_avg, y_test)
)

save_model_tf(object = model_avg, filepath = base_model_avg_fname)

require(MASS)
require(data.table)
dt <- data.table(y = rowMeans(x_train), x = factor(max.col(y_train)))
boxplot(y ~ x, data = dt)
require(ggplot2)
require(rms)
require(matrixStats)

p <- ggplot(dt, aes(x = x, y = y, fill = x)) +
     geom_violin(trim = F) + geom_boxplot(width = 0.1) + theme_bw() + #theme_minimal() + theme_classic() +
     xlab(label = "Arm level scale copy number") + 
     ylab(label = "small scale (31 genes) copy number") +
     theme(axis.text=element_text(size = 12),
           axis.title=element_text(size = 14,face = "bold"))
  
p + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

# I switched this to orm instead of polr
olr_model <- polr( y ~ ., data = data.table(x_train, y = factor(max.col(y_train))) )
#summary(olr_model)
preds <- predict(object = olr_model, newdata = data.table(x_test))
tbl <- table(preds, factor(max.col(y_test)))
print(tbl)
message("Accuracy of OLR model with all x variables as input: ", sum(diag(tbl/sum(tbl)))) #accuracy

olr_model <- polr(y ~ x.mean + x.sd, data = data.table(x.mean = rowMeans(x_train), x.sd = rowSds(x_train), y = factor(max.col(y_train))))
#summary(olr_model)
preds <- predict(object = olr_model, newdata = data.table(x.mean = rowMeans(x_test), x.sd = rowSds(x_test)))
tbl <- table(preds, factor(max.col(y_test)))
print(tbl)
message("Accuracy of OLR model with sd and mean as input: ", sum(diag(tbl/sum(tbl)))) #accuracy

olr_model <- polr(y ~ x, data = data.table(x = rowMeans(x_train), y = factor(max.col(y_train))))
#summary(olr_model)

preds <- predict(object = olr_model, newdata = data.table(x = rowMeans(x_test)))
tbl <- table(preds, factor(max.col(y_test)))
print(tbl)
message("Accuracy of OLR model with mean: ", sum(diag(tbl/sum(tbl)))) #accuracy
print(diag(tbl)/rowSums(tbl))

preds <- max.col(predict(model_avg, x = rowMeans(x_test)))
tbl <- table(preds, factor(max.col(y_test)))
print(tbl)
message("NN 1 layer softmax accuracy: ", sum(diag(tbl/sum(tbl)))) #accuracy
print(diag(tbl)/rowSums(tbl))

saveRDS(olr_model, file = olr_model_fname)



