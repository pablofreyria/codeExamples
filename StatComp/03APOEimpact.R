
#required data on file: dataPax, dataNoId, optim_eta,optim_cp



### stratify by # of e4 apoe allele and remove column
id_apoe <- which(colnames(dataNoId)%in% c("NACCNE4S","NACCAPOE")) # NACCAPOE is the type of mutation
dataNoApoe <- dataNoId[,-id_apoe]
dataApoe0 <- dataNoId[dataNoId$NACCNE4S==0,-id_apoe]
dataApoe1 <- dataNoId[dataNoId$NACCNE4S==1,-id_apoe]
dataApoe2 <- dataNoId[dataNoId$NACCNE4S==2,-id_apoe]

true_covs <- which(colnames(dataNoApoe)!="NACCALZD")
### model in all the population and mean sq error

train_xgmodel <- function(dataF,covs,eta=optim_eta,thres=0.5){
  # change of format to match xgboost
  data_all <- data.matrix(dataF)
  xgboost_mod <- xgboost(data= data_all[,true_covs],
                         params=list(eta=optim_eta,gamma=0.001),
                         label = data_all[,"NACCALZD"],
                         nrounds=150, # with verbose showed good convergence at 50
                         verbose=0)
  
  #predictions
  preds_all <- predict(xgboost_mod,newdata = data_all[,true_covs],type = "prob")
  
  mse <- mean((data_all[,"NACCALZD"] - preds_all)^2)
  pos_AD <- mean(data_all[,"NACCALZD"])
  pos_preds <- mean(preds_all>=thres)
  prevalence <- data.frame(Original = pos_AD,Predicted=pos_preds)
  
  n <- nrow(data_all)
  missclass <- mean(as.numeric(preds_all>=thres) !=  data_all[,"NACCALZD"])
  fp <- mean(preds_all>=thres &  data_all[,"NACCALZD"]==0)
  fn <- mean(preds_all<thres &  data_all[,"NACCALZD"]==1)
  err_report <- data.frame(MissClas = missclass,FP = fp,FN=fn)
  
  results <- list(model = xgboost_mod,MSE=mse,N=n,Report=err_report,Prev = prevalence)
  return(results)
}

pop_all <- train_xgmodel(dataNoApoe,true_covs)
pop_0e4 <- train_xgmodel(dataApoe0,true_covs)
pop_1e4 <- train_xgmodel(dataApoe1,true_covs)
pop_2e4 <- train_xgmodel(dataApoe2,true_covs)

### test model against the other population
test_model_prevalence <- function(xgmodel,test_data,covs=true_covs,thres=0.5){
  
  test_form <- data.matrix(test_data)
  #predictions
  preds_all <- predict(xgmodel,newdata = test_form[,true_covs],type = "prob")
  
  pred_prev <- mean(preds_all>=thres)
  obs_prev <- mean(test_form[,"NACCALZD"])
  
  res <- data.frame(Obs=obs_prev,Pred=pred_prev)
  return(res)
}
 
test_model_prevalence(pop_0e4$model,dataApoe1) 
