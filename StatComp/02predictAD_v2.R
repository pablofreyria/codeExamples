library(rpart)
library(rpart.plot)
library(e1071)
library(xgboost)
library(glmnet)
library(origami)
library(ggplot2)

#v1 original analysis
#v2 changed to factor CVHATT, CVCHF, CBSTROKE,DIABETES,B12DEF, DEP2YRS, DEPOTHR, ALCOHOL,ABUSEOTHR

# remove scores and years of events that are unlikely (heart attack, parkinsons, stroke) 
# drop hispanic and keep race
# quitsmoke is left as 888 because around 50% have never smoked

#rm_vars <- c(scores2,"HATTYEAR","PDYR","NACCSTYR") 
rm_vars <- c("HATTYEAR","PDYR","NACCSTYR")
vars2 <- lapply(vars, function(x) setdiff(x,rm_vars))

dataPax <- finalData[,unlist(vars2)]

# distribution of variable types
# numer <- c("BIRTHYR","NACCAGE","SMOKYRS","QUITSMOK","EDUC","BMI")
# categor <- c("RACE","NACCCOGF","COGMODE","MEMORY","ORIENT","PERSCARE","SATIS","BORED")
# ordinal <- c("COGMODE","MEMORY","ORIENT","PERSCARE","SATIS","BORED")
# binar <- c("SEX","DEMENTED","NACCALZD","CVHATT","CVCHF","CBSTROKE","PD","DIABETES",
#            "B12DEF","DEP2YRS","DEPOTHR","ALCOHOL","ABUSOTHR","ANYMEDS")

#CVHATT, CVCHF, CBSTROKE,DIABETES,B12DEF, DEP2YRS, DEPOTHR, ALCOHOL,ABUSEOTHR

dataPax[,scores2][dataPax[,scores2]==-4] <- NA
#imput BMI as the median
dataPax$NACCBMI[dataPax$NACCBMI==888.8] <- median(dataPax$NACCBMI[dataPax$NACCBMI!=888.8])
dataPax$RACE <- as.factor(dataPax$RACE)
dataPax$NACCCOGF <- as.factor(dataPax$NACCCOGF)
dataPax$SATIS <- as.factor(dataPax$SATIS)
dataPax$BORED <- as.factor(dataPax$BORED)
dataPax$NACCAPOE <- as.factor(dataPax$NACCAPOE)
dataPax$PD <- as.factor(dataPax$PD)

# set the QUITSMOK as the interaction between SMOKYRS and QUITSMOK
dataPax$QUITSMOK <-dataPax$QUITSMOK*dataPax$SMOKYRS

# V2 change
dataPax$CVHATT <- as.factor(dataPax$CVHATT)
dataPax$CVCHF <- as.factor(dataPax$CVCHF)
dataPax$CBSTROKE <- as.factor(dataPax$CBSTROKE)
dataPax$DIABETES <- as.factor(dataPax$DIABETES)
dataPax$B12DEF <- as.factor(dataPax$B12DEF)
dataPax$DEP2YRS <- as.factor(dataPax$DEP2YRS)
dataPax$DEPOTHR <- as.factor(dataPax$DEPOTHR)
dataPax$ALCOHOL <- as.factor(dataPax$ALCOHOL)
dataPax$ABUSOTHR <- as.factor(dataPax$ABUSOTHR)


dataNoId <- dataPax[,-1]

### Cross validation procedure to tune parameters ####

# select data for procedure
dataCV <- dataNoId

# main variables
V <- 10
nparams <- 20
cp_seq <- seq(0.00001,0.015,length = nparams)
eta_seq <- seq(0.2, 0.6, length = nparams)
n_lambda <- 80

# shuffling data and creating folds
nCV <- nrow(dataCV)
dataCV <- dataCV[sample(1:nCV,nCV,replace = F),]
cv_folds <- make_folds(n=nCV,fold_fun = folds_vfold,V=V,strata_ids = dataCV$NACCALZD)

loss_cp <- matrix(NA,V,nparams) 
loss_eta <-  matrix(NA,V,nparams)

loss_lambda <- vector("list",V)
lambda_seqs <- vector("list",V)
n_complete <- numeric(V)
true_covs <- which(colnames(dataCV)!="NACCALZD")
# try parameters
for (i in 1:V) {
  y_obj <- dataCV$NACCALZD[cv_folds[[i]]$validation_set]
  train_data <- dataCV[cv_folds[[i]]$training_set,]
  valid_data <- dataCV[cv_folds[[i]]$validation_set,]
  # coercing to double the data frame
  data_tr <- data.matrix(train_data)
  data_va <- data.matrix(valid_data)
  
  data_tr_compl <- data_tr[complete.cases(data_tr),]
  data_va_compl <- data_va[complete.cases(data_va),]
  aux_glmnet <- glmnet(data_tr_compl[,true_covs],data_tr_compl[,"NACCALZD"], 
                       family = "binomial",
                       alpha=1,
                       nlambda = n_lambda)
  y_pred_glmnet <- predict(aux_glmnet,newx = data_va_compl[,true_covs],type = "response")
  n_complete[i] <- sum(complete.cases(data_va[,true_covs]))
  
  lambda_seqs[[i]] <- aux_glmnet$lambda
  loss_lambda[[i]] <- colSums(apply(y_pred_glmnet, 2, function(x) x - y_obj[complete.cases(data_va)])^2)

  for (j in 1:nparams) {
    # tree
    aux_tree <- rpart(NACCALZD~.,
                      data = train_data,
                      method = "class",
                      cp=cp_seq[j])
    y_pred_tre <- predict(aux_tree,
                      newdata=valid_data,
                      type="prob")[,"1"]

    loss_cp[i,j] <- sum((y_obj - y_pred_tre)^2)

    ## xgboost
    aux_xgbost <- xgboost(data= data_tr[,true_covs],
                          params=list(eta=eta_seq[j],gamma=0.001),
                          label = data_tr[,"NACCALZD"],
                          nrounds = 150,
                          verbose=0)
    y_pred_xg <- predict(aux_xgbost,newdata = data_va[,true_covs])

    loss_eta[i,j] <- sum((y_obj - y_pred_xg)^2)
    print(sprintf("Param_id:%i, Fold=%i ",j,i))
  }
}

#
## retrieve optimal parameters
optim_cp  <- cp_seq[which.min(colSums(loss_cp))]
optim_eta <- eta_seq[which.min(colSums(loss_eta))] 

optim_lam <-  unlist(lambda_seqs)[which.min(unlist(loss_lambda))]
    
### Make plots of the MSE ###
#
# change format
dat_cp  <- data.frame(param = cp_seq,MSE = colSums(loss_cp)/nCV)
dat_eta <- data.frame(param = eta_seq, MSE = colSums(loss_eta)/nCV)

# glmnet : mse by fold using linear interpol to address different lambda seqs
# and averaging the folds
MSE_fold <- sapply(1:V,function(x) loss_lambda[[x]]/n_complete[x])
funs <- sapply(1:V, function(x) approxfun(lambda_seqs[[x]],MSE_fold[[x]],rule = 2))
all_lams <- unlist(lambda_seqs)
loss_lambda_matrix <- t(sapply(funs, function(x) x(all_lams)))
MSE_lambda <- colMeans(loss_lambda_matrix)

dat_lam <-  data.frame(param = all_lams, MSE = MSE_lambda)

## make plots

plot_cp <- ggplot(data = dat_cp, aes(x=param,y=MSE)) + geom_line() + 
  ggtitle("Classification tree") + xlab("cp-param") + ylab("CV-MSE") +  
  theme(text=element_text(size=18))

plot_eta <- ggplot(data = dat_eta, aes(x=param,y=MSE)) + geom_line() + 
  ggtitle("XGboosting") + xlab("eta-param") + ylab("CV-MSE") +
  theme(text=element_text(size=18))

plot_lambda <- ggplot(data = dat_lam, aes(x=param,y=MSE)) + geom_line() + 
  ggtitle("LASSO logistic regression") + xlab("lambda-param") + ylab("CV-MSE") +
  theme(text=element_text(size=18))

## print to system
png(filename="TreeParamTune.png", width = 480, height = 480)
plot(plot_cp)
dev.off()

png(filename="XgboostParamTune.png", width = 480, height = 480)
plot(plot_eta)
dev.off()

png(filename="LassoParamTune.png", width = 480, height = 480)
plot(plot_lambda)
dev.off()
 

### ROC curve with best model ###
# best model is xgboost and tree, no point in bringing glmnet
dataROC <- dataNoId

#train xgboost model
data_tr <- data.matrix(dataROC)
xgboost_opt <- xgboost(data= data_tr[,true_covs],
                       params=list(eta=optim_eta,gamma=0.001),
                       label = data_tr[,"NACCALZD"],
                       nrounds=80, # with verbose showed good convergence at 50
                       verbose=0)
xgboost_preds <- predict(xgboost_opt,newdata = data_tr[,true_covs],type = "prob")

mean((data_tr[,"NACCALZD"] - xgboost_preds)^2)

# train class tree
tree_opt <- rpart(NACCALZD~.,
                  data = dataROC,
                  method = "class",
                  cp=optim_cp)
tree_preds <- predict(tree_opt,
                      newdata=dataROC,
                      type="prob")[,"1"]

mean((dataROC$NACCALZD - tree_preds)^2)

## ROC
#xgboost
seq_thres <- seq(1,0,length=20)

compute_tp_fp <- function(preds,targets,thres){
  all_preds_pos <- sapply(thres, function(x) as.numeric(preds>=x)) 
  total_pos <- colSums(all_preds_pos)
  
  fp <- apply(all_preds_pos, 2, function(x) sum(targets==0 & x==1))
  tp <- apply(all_preds_pos, 2, function(x) sum(targets==1 & x==1))
  fn <- apply(all_preds_pos, 2, function(x) sum(targets==1 & x==0))
  tn <- apply(all_preds_pos, 2, function(x) sum(targets==0 & x==0))
  
  fpr <- fp/(fp+tn) # false positives among real negatives
  tpr <- tp/(tp+fn) # trie positives among real positives
  return(data.frame(fpr,tpr))
}

roc_xgboost <- compute_tp_fp(xgboost_preds,dataROC$NACCALZD,seq_thres)
roc_tree <- compute_tp_fp(tree_preds,dataROC$NACCALZD,seq_thres)
  
  
## plot xaxis=fp; yaxis=tp
roc_xgboost$Algorithm = rep("xgboost",length(seq_thres))
roc_tree$Algorithm = rep("Tree",length(seq_thres))
plot_data <- rbind.data.frame(roc_xgboost,roc_tree)

roc_plot <- ggplot(data=plot_data,aes(x=fpr,y=tpr,col=Algorithm)) + geom_line() + 
  ggtitle("ROC curves") + xlab("False positive rate") + ylab("True positive rate") +
  theme(text=element_text(size=15))

## print to system
png(filename="ROCcurves.png", width = 480, height = 480)
plot(roc_plot)
dev.off()



### clean all variables ###
keep <- c("dataPax","dataNoId","vars2","optim_cp","optim_eta")

rm(list=setdiff(ls(),keep))

