### Function to run the Super Learner on data set
library(origami)
library(sl3)
library(SuperLearner)
library(forecast)

prediction_sl_origin <- function(dataset,outname,varnames){
  
  ## Define cross validation fold structure apt for time series
  #
  # It will be made to have 10 folds, 
  # 20% of data not seen in training
  val_size <- floor(nrow(dataset)*0.15)
  gap_size <- floor(nrow(dataset)*0.05)
  #
  #begin training with half of data and make 10 folds
  windw_size <- floor(nrow(dataset)*0.5) # 50% of data
  batch_size <- floor(nrow(dataset)*0.03) #30% remaining of data in 10 folds
  #
  # Rolling origin that grows the training set in time and keeps a constant gap 
  # and a constant validation set
  sl_folds_origin <- make_folds(n=nrow(dataset),
                                fold_fun = folds_rolling_origin,
                                first_window = windw_size,
                                validation_size = val_size,
                                gap = gap_size,
                                batch = batch_size)
  
  ## Define the machine learning task 
  sl_task_origin <- make_sl3_Task(data = dataset, 
                                  covariates = varnames,
                                  outcome = outname,
                                  outcome_type = "continuous",
                                  folds=sl_folds_origin)

  # initialize learners: ARIMA, splines(with sin basis expansion), SVM, xgboosting
  lrn_arima <- make_learner(Lrnr_arima)
  lrn_splne <- Lrnr_pkg_SuperLearner$new("SL.polymars")
  lrn_svm_s <- make_learner(Lrnr_svm, kernel = "sigmoid")
  lrn_svm_p <- make_learner(Lrnr_svm, kernel = "polynomial",degree = 3)
  
  lrn_nnt <- make_learner(Lrnr_nnet)
  
  # gradient boosting, regression tree with various eta, depths and rounds
  grid_params <- list(eta = seq(0.2,0.4,by=0.1), max_depth = seq(5,7,by=1))
  grid <- expand.grid(grid_params,KEEP.OUT.ATTRS = FALSE)
  xgb_learners <- apply(grid, MARGIN = 1, function(tuning_params) {
    do.call(Lrnr_xgboost$new, as.list(tuning_params))
  })

  
  SL.library = c(lrn_arima,lrn_splne,lrn_svm_p,lrn_svm_s,xgb_learners,lrn_nnt)
  
  super_learner <- make_learner(Lrnr_sl,learners=SL.library)
  
  sl_trained <- super_learner$train(sl_task_origin)
  
  result <- list("SuperLearner" = sl_trained,
                 "ML_Task" = sl_task_origin,
                 "Preds" = sl_trained$predict(sl_task_origin))
  
  return(result)
  
}

prediction_sl_window <- function(dataset,outname,varnames){
  
  ## Define cross validation fold structure apt for time series
  #
  # It will be made to have 10 folds, 
  # 20% of data not seen in training
  val_size <- floor(nrow(dataset)*0.15)
  gap_size <- floor(nrow(dataset)*0.05)
  #
  #begin training with half of data and make 10 folds
  windw_size <- floor(nrow(dataset)*0.5) # 50% of data
  batch_size <- floor(nrow(dataset)*0.03) #30% remaining of data in 10 folds
  #
  # Rolling origin that grows the training set in time and keeps a constant gap 
  # and a constant validation set
  sl_folds_origin <- make_folds(n=nrow(dataset),
                                fold_fun = folds_rolling_window,
                                window_size = windw_size,
                                validation_size = val_size,
                                gap = gap_size,
                                batch = batch_size)
  
  ## Define the machine learning task 
  sl_task_origin <- make_sl3_Task(data = dataset, 
                                  covariates = varnames,
                                  outcome = outname,
                                  outcome_type = "continuous",
                                  folds=sl_folds_origin)
  
  # initialize learners: ARIMA, splines(with sin basis expansion), SVM, xgboosting
  lrn_arima <- make_learner(Lrnr_arima)
  lrn_splne <- Lrnr_pkg_SuperLearner$new("SL.polymars")
  lrn_svm_s <- make_learner(Lrnr_svm, kernel = "sigmoid")
  lrn_svm_p <- make_learner(Lrnr_svm, kernel = "polynomial",degree = 3)
  
  lrn_nnt <- make_learner(Lrnr_nnet)
  
  # gradient boosting, regression tree with various eta, depths and rounds
  grid_params <- list(eta = seq(0.2,0.4,by=0.1), max_depth = seq(5,7,by=1))
  grid <- expand.grid(grid_params,KEEP.OUT.ATTRS = FALSE)
  xgb_learners <- apply(grid, MARGIN = 1, function(tuning_params) {
    do.call(Lrnr_xgboost$new, as.list(tuning_params))
  })
  
  
  SL.library = c(lrn_arima,lrn_splne,lrn_svm_p,lrn_svm_s,xgb_learners,lrn_nnt)
  
  super_learner <- make_learner(Lrnr_sl,learners=SL.library)
  
  sl_trained <- super_learner$train(sl_task_origin)
  
  result <- list("SuperLearner" = sl_trained,
                 "ML_Task" = sl_task_origin,
                 "Preds" = sl_trained$predict(sl_task_origin))
  
  return(result)
  
}



cross_validate_SL <- function(dataset,SL,outname,varnames){
  
  # make fold structure
  val_size <- floor(nrow(dataset)*0.15)
  gap_size <- floor(nrow(dataset)*0.05)
  #
  #begin training with half of data and make 10 folds
  windw_size <- floor(nrow(dataset)*0.5) # 50% of data
  batch_size <- floor(nrow(dataset)*0.1) #30% remaining of data in 4 folds
  #
  # Rolling origin that grows the training set in time and keeps a constant gap 
  # and a constant validation set
  folds_cv <- make_folds(n=nrow(dataset),
                                fold_fun = folds_rolling_window,
                                window_size = windw_size,
                                validation_size = val_size,
                                gap = gap_size,
                                batch = batch_size)
  

  
  # make task of cross validation
  cv_task <- make_sl3_Task(data=dataset, covariates=varnames, outcome=outname,
                           outcome_type = "continuous",
                           folds =folds_cv)
  
  CV_sl <- CV_lrnr_sl(lrnr_sl = SL,cv_task,loss_fun = loss_squared_error)
  
  return(CV_sl)
}