
# Based on the following kaggle competition:
# https://www.kaggle.com/c/lish-moa/data?select=train_targets_scored.csv
#
# Summary: Predict probability of response for each MoA (Mechanism of Action) of
# different drug samples given features such as: 772 genes expresions, cell 
# viability (% alive after tx) in 100 types of cells, treatment (tx or control)
# duration and dose of drug.
# Training data has 23814 rows
# Training targets are binary for each mechanism of action (207 MoAs)
# 

# simple mean

# mahalanobis distance

# k nearest neighbour

# logistic regression

## Targeting?

## Read data and preprocessing ####
rm(list=ls()) #CAUTION: clears working directory
setwd("~/Berkeley Classes/2021Spring/Big data A public health perspective/Project 4")
source("P4_project_functions.R")

wd <- getwd()
setwd("data")

features <- read.csv("train_features.csv")
#drugs    <- read.csv("train_drug.csv") #non-informative and not in testing set
targets  <- read.csv("train_targets_scored.csv")

# test <- read.csv("test_features.csv")
# all(names(features)==names(test_features))
# rm(test_features)
setwd(wd)




# all rows match in table, can erase ids in training / testing
if(all(features$sig_id==targets$sig_id)){
  targets <- targets[,-1]
  features <- features[,-1]  
}

# keep originals
orig_features <- features
orig_targets <- targets

# change low and high dose for -1 and 1
features$cp_dose <- 1-2*as.numeric(features$cp_dose=="D1")

# remove experiments used for control, no control had a MoA activated
ids_trt <- which(features$cp_type=="trt_cp")
id_rm <- which(colnames(features)=="cp_type") # MoA in control == 0 always, no variability
features <- features[ids_trt,-id_rm]
targets <- targets[ids_trt,]

## Fold structure and loss function to crossvalidate algorithms ####
#

library(origami)
# set as strata ID whether drug has none or any MoA (balanced sets)
v <- 2
folds <- make_folds(n = nrow(targets),fold_fun = folds_vfold, V=v,
                    strata_ids = rowSums(targets)>0)


print(timestamp())
# euclidean distance
risks_euc <- matrix(NA,5,v)
for(i in c(1:v)){
  
  train_targets  <- targets[folds[[i]]$training_set,]
  train_features <- features[folds[[i]]$training_set,]
  
  test_targets  <- targets[folds[[i]]$validation_set,]
  test_features <- features[folds[[i]]$validation_set,]
  
  # Regular clusters in "raw data"
  clusters_moa    <- clusters_characterization(targets = train_targets,
                                               features =   train_features,
                                               activeMOA=TRUE)
  
  clusters_no_moa <- clusters_characterization(targets = train_targets,
                                               features = train_features,
                                               activeMOA=FALSE)
  
  clusters <- list(MOA=clusters_moa,NoMOA=clusters_no_moa)
  
  # euclidian distance with number of classes
  preds1 <- predict_distance_based(train_ft = NULL,
                                   train_tgts = NULL,
                                   test_ftrs =  test_features,
                                   dist_type ="euclidean",
                                   trained_clusters = clusters,
                                   weights_fun = "nclass",pow = 2)$preds
  
  # euclidean distance with inverse
  preds2 <- predict_distance_based(train_ft = NULL,
                                   train_tgts = NULL,
                                   test_ftrs =  test_features,
                                   dist_type ="euclidean",
                                   trained_clusters = clusters,
                                   weights_fun = "inverse",pow = 1)$preds
  # euclidean distance with inverse
  preds3 <- predict_distance_based(train_ft = NULL,
                                   train_tgts = NULL,
                                   test_ftrs =  test_features,
                                   dist_type ="euclidean",
                                   trained_clusters = clusters,
                                   weights_fun = "inverse",pow = 5)$preds
  
  # euclidean distance with softmax
  preds4 <- predict_distance_based(train_ft = NULL,
                                   train_tgts = NULL,
                                   test_ftrs =  test_features,
                                   dist_type ="euclidean",
                                   trained_clusters = clusters,
                                   weights_fun = "softmax",pow = 1)$preds
  
  preds5 <- predict_distance_based(train_ft = NULL,
                                   train_tgts = NULL,
                                   test_ftrs =  test_features,
                                   dist_type ="euclidean",
                                   trained_clusters = clusters,
                                   weights_fun = "softmax",pow = 5)$preds
  
  risks_euc[1,i] <- log_loss_score(preds1,test_targets)
  risks_euc[2,i] <- log_loss_score(preds2,test_targets)
  risks_euc[3,i] <- log_loss_score(preds3,test_targets)
  risks_euc[4,i] <- log_loss_score(preds4,test_targets)
  risks_euc[5,i] <- log_loss_score(preds5,test_targets)
}

print(timestamp())
# mahalanobis distance
risks_maha <- matrix(NA,5,v)
for(i in c(1:v)){
  
  train_targets  <- targets[folds[[i]]$training_set,]
  train_features <- features[folds[[i]]$training_set,]
  
  test_targets  <- targets[folds[[i]]$validation_set,]
  test_features <- features[folds[[i]]$validation_set,]
  
  # Regular clusters in "raw data"
  clusters_moa    <- clusters_characterization(targets = train_targets,
                                               features =   train_features,
                                               activeMOA=TRUE)
  
  clusters_no_moa <- clusters_characterization(targets = train_targets,
                                               features = train_features,
                                               activeMOA=FALSE)
  
  clusters <- list(MOA=clusters_moa,NoMOA=clusters_no_moa)
  
  # euclidian distance with number of classes
  preds1 <- predict_distance_based(train_ft = NULL,
                                   train_tgts = NULL,
                                   test_ftrs =  test_features,
                                   dist_type ="mahalanobis",
                                   trained_clusters = clusters,
                                   weights_fun = "nclass",pow = 2)$preds
  
  # euclidean distance with inverse
  preds2 <- predict_distance_based(train_ft = NULL,
                                   train_tgts = NULL,
                                   test_ftrs =  test_features,
                                   dist_type ="mahalanobis",
                                   trained_clusters = clusters,
                                   weights_fun = "inverse",pow = 1)$preds
  # euclidean distance with inverse
  preds3 <- predict_distance_based(train_ft = NULL,
                                   train_tgts = NULL,
                                   test_ftrs =  test_features,
                                   dist_type ="mahalanobis",
                                   trained_clusters = clusters,
                                   weights_fun = "inverse",pow = 5)$preds
  
  # euclidean distance with softmax
  preds4 <- predict_distance_based(train_ft = NULL,
                                   train_tgts = NULL,
                                   test_ftrs =  test_features,
                                   dist_type ="mahalanobis",
                                   trained_clusters = clusters,
                                   weights_fun = "softmax",pow = 1)$preds
  
  preds5 <- predict_distance_based(train_ft = NULL,
                                   train_tgts = NULL,
                                   test_ftrs =  test_features,
                                   dist_type ="mahalanobis",
                                   trained_clusters = clusters,
                                   weights_fun = "softmax",pow = 5)$preds
  
  risks_maha[1,i] <- log_loss_score(preds1,test_targets)
  risks_maha[2,i] <- log_loss_score(preds2,test_targets)
  risks_maha[3,i] <- log_loss_score(preds3,test_targets)
  risks_maha[4,i] <- log_loss_score(preds4,test_targets)
  risks_maha[5,i] <- log_loss_score(preds5,test_targets)
}

print(timestamp())




## With PCA transform ####
# get pca transformed features
# "training"
# transform_object <- pca_features(features,comps = 55)
# train_fts_pca <- transform_object$newCords
# 
# # Transform new data (for testing)
# test_fts_pca <- pca_new_coords(features,transform_object)
# 
# pca_clusters_moa    <- clusters_characterization(targets = targets,
#                                                  features =   train_fts_pca,
#                                                  activeMOA=TRUE)
# 
# pca_clusters_no_moa <- clusters_characterization(targets = targets,
#                                                  features =   train_fts_pca,
#                                                  activeMOA=FALSE)
# 
# pca_clusters <- list(MOA=pca_clusters_moa,NoMOA=pca_clusters_no_moa)
# 
# 
# preds_pca <- predict_distance_based(train_ft = train_fts_pca,train_tgts = targets,
#                                  test_ftrs =  test_fts_pca,dist_type ="mahalanobis",
#                                  trained_clusters = pca_clusters,
#                                  weights_fun = "nclass",pow = 2)
# 
# loss_nclass_maha <- log_loss_score(preds_pca$preds,targets)
# 
# p_moa <- distance_mahalanobis(test_fts_pca,pca_clusters[[1]])$probability
# p_no_moa <- distance_mahalanobis(test_fts_pca,pca_clusters[[2]])$probability
# 
# preds_p <- predict_f_weights(p_moa,p_no_moa,
#                              pca_clusters_moa$nobs,pca_clusters_no_moa$nobs)
# 
# loss_nclass_p <- log_loss_score(preds_p,targets)







