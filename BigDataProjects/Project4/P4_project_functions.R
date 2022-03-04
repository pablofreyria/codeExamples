
### PCA, projections, classes and pseudoinverse functions ####

# wrapper function to do pca and keep only relevant parts
pca_features <- function(features,comps=148,return_coords = TRUE){
  
  aux <- prcomp(features,rank. = comps)
  
  center <- aux$center
  scale <- aux$scale
  rot_matrix <- aux$rotation
  
  if (return_coords) {
    result <- list(center=center,scale=scale,rotation=rot_matrix,
                   newCords = aux$x)
  }else{
    result <- list(center=center,scale=scale,rotMatrix=rot_matrix)
  }
  
  return(result)
}

# return new coordinates from a pca projection
pca_new_coords <- function(new_features,pca_obj){
  
  new_coords <- as.matrix(scale(new_features,
                                center = pca_obj$center,
                                scale = pca_obj$scale)
  )%*%pca_obj$rotation
  
  return(new_coords)
}

pseudoinverse_variance <- function(covariance,tol=1e-5){
  
  # tol serves to ignore all eigenvalues below the threshold
  
  svd_comps <- svd(covariance)
  n <- nrow(covariance)
  m <- ncol(covariance)
  
  comp_tol <- which(svd_comps$d>=tol)
  n_comp_tol <- length(comp_tol)
  
  U_prime <- svd_comps$u[,comp_tol]
  d_prime <- diag(1/svd_comps$d[comp_tol],n_comp_tol,n_comp_tol)
  V_prime <- svd_comps$v[,comp_tol]
  
  #reconstructed_inv <- svd_comps$v%*%d_prime%*%svd_comps$u
  
  reconstructed_inv <- V_prime%*%d_prime%*%t(U_prime)

  
  return(reconstructed_inv)
}


clusters_characterization <- function(targets,features,activeMOA=TRUE,getVariance=TRUE){
  # activeMOA is for each class cluster observations with moa =1 vs moa = 0 
  
  n_classes <- ncol(targets)
  n_features <- ncol(features)
  
  class_mean <- matrix(NA,nrow = n_features,ncol = n_classes)
  if(getVariance){
    class_vars <- vector("list",length = n_classes)  
  }
  
  n <- numeric(n_classes)
  for (i in c(34:n_classes)) {
    ids_moa_active <- which(targets[,i]==as.numeric(activeMOA))
    n[i] <- length(ids_moa_active)  
    if (n[i]==0) {
      class_mean[,i] <- colMeans(features)
      
    } else if(n[i]>1) {
      class_mean[,i] <- colMeans(features[ids_moa_active,])
    } else {
      class_mean[,i] <-t(features[ids_moa_active,])
    }
      
    if(getVariance){
      if(length(ids_moa_active)>1){
        covariance <- cov(features[ids_moa_active,])
        if (det(covariance)<1e-10) {
          # compute the inverse with svd decomposition 
          class_vars[[i]] <- pseudoinverse_variance(covariance,tol = 1e-5)  
        }else{
          class_vars[[i]] <- solve(covariance)
        }
      }else{
        class_vars[[i]] <- 1 # no variance for single observation
      }
    }
    
  }
  
  if (getVariance) {
    result <- list(centers = class_mean, inv_variances = class_vars, nobs = n)
  }else{
    result <- list(centers = class_mean, nobs = n)
  }
  
  return(result)
}




### Distance functions ####

### Distance to weights ###

distance_euclidean <- function(features,centers){
  
  # outcome variable: distance to each point
  distances <- matrix(NA,nrow =nrow(features) ,
                      ncol = ncol(centers)
  )
  
  for (i in c(1:ncol(centers))) {
    distances[,i] <- t(colSums((t(features)-centers[,i])^2))
  }
  
  return(distances)
  
}

distance_mahalanobis <- function(features,clusters,computeP=TRUE,df=14){
# df corresponds to the meaningfull components in the psuedoinverse of covariance

  #outcome variable: distance to each point
  distances <- matrix(NA,nrow =nrow(features) ,
                      ncol = ncol(clusters$centers)
  )
  
  for (i in c(1:ncol(clusters$centers))) {
    
    if(clusters$nobs[i]>1){
      distances[,i] <- mahalanobis(as.matrix(features),clusters$centers[,i],
                                   clusters$inv_variances[[i]],
                                   inverted = TRUE)
    }else{
      print(sprintf("Warning: Only 1 observation in cluster %i, euclidian distance used",i))
      distances[,i] <- t(colSums((t(features)-clusters$centers[,i])^2))
    }
    
  }
  
  if (computeP) {
    p <- ncol(features)
    probs <- pchisq(distances,df,lower.tail = FALSE)
    result <- list(distance = distances,probability = probs)
  }else{
    result <- distances
  }
  
  return(result)
}


# Using powers of inverse distance as probability
distance_weights_inverse<- function(distances,power=1){
  
  weights <- (1/distances^power)
  weights[is.na(weights)] <- Inf #distance=0 (it is the stereotypical class)
  
  return(weights)
}

# Using softmax function (standarized exponential scale)
distance_weights_softmax <- function(distances,power=1){
  
  weights <- exp(-power*distances)
  
  return(weights)
  
}


distance_weights_nclass <- function(distances,n_clases){
  
  weights <- t(n_clases/t(distances^2))
  
  return(weights)
}



### Predict functions ####

predict_f_mean <- function(train_ft=NULL,train_tgts,test_ftrs){
  
  p.mean <- colMeans(train_tgts)
  
  n <- dim(test_ftrs)[1]
  m <- dim(train_tgts)[2]
  test_preds <- t(matrix(p.mean,nrow = m,ncol=n))
  
  return(test_preds)
}

predict_distance_based <- function(train_ft=features,train_tgts=targets,
                                   test_ftrs=test, dist_type="euclidean",
                                   trained_clusters = NULL,
                                   weights_fun = c("inverse","nclass","softmax"),
                                   pow=2)
  {
  if(is.null(trained_clusters)){
    # first get the clusters
    clusters_moa    <- clusters_characterization(targets = train_tgts,
                                                 features =   train_ft,
                                                 activeMOA=TRUE)
    
    clusters_no_moa <- clusters_characterization(targets = train_tgts,
                                                 features =   train_ft,
                                                 activeMOA=FALSE)
    trained_clusters <- list(MOA=clusters_no_moa,NoMOA=trained_clusters)
    
  }else{
    clusters_moa <- trained_clusters[[1]]
    clusters_no_moa <- trained_clusters[[2]]
  }
  
  # compute the distance
  if (dist_type=="euclidean") {
    dist_moa <- distance_euclidean(test_ftrs,clusters_moa$centers)
    dist_no_moa <- distance_euclidean(test_ftrs,clusters_no_moa$centers)
  }else{
    dist_moa <- distance_mahalanobis(test_ftrs,clusters_moa,computeP = FALSE)
    dist_no_moa <- distance_mahalanobis(test_ftrs,clusters_no_moa,computeP = FALSE)
  }
  
  # compute weights (input directly translated to probability)
  if (weights_fun=="inverse") {
    weights1 <- distance_weights_inverse(dist_moa,power = pow)
    weights2 <- distance_weights_inverse(dist_no_moa,power = pow)
    
  }else if (weights_fun=="softmax") {
    weights1 <- distance_weights_softmax(dist_moa,power = pow)
    weights2 <- distance_weights_softmax(dist_no_moa,power = pow)
  }else if(weights_fun=="nclass"){
    weights1 <- distance_weights_nclass(dist_moa,
                                        n_clases = clusters_moa$nobs)
    weights2 <- distance_weights_nclass(dist_no_moa,
                                        n_clases = clusters_no_moa$nobs)
  }
  
  predicted_p <- predict_f_weights(weights1,weights2,
                                   clusters_moa$nobs,clusters_no_moa$nobs)
  
  res <- list(preds = predicted_p,clusters=trained_clusters)
  return(res)
}

# use weights to predict probability of belonging to class 1
predict_f_weights <- function(weights1,weights2,n1,n2){
  
  ids_inf <- is.infinite(weights1)
  
  p_class1 <- t(t(weights1)*n1/(t(weights1)*n1 + t(weights2)*n2))
  
  id_indef <- which(is.na(p_class1))
  p_class1[id_indef] <- 0.5
  
  p_class1[ids_inf] <- 1
  
  return(p_class1)
}


### Error and cross validation functions ####

#loss function
log_loss_score <- function(pred,observed,bounds=1e-8){
  
  if (!all(dim(pred)==dim(observed))) {
    stop("Arguments have different dimensions")
  }
  n <- nrow(pred)
  m <- ncol(pred)
  
  # get as vector
  pred1 <- matrix(as.matrix(pred),nrow = m*n,ncol = 1)
  obs1  <- matrix(as.matrix(observed),nrow = m*n,ncol = 1)
  
  loss <- numeric(m*n)
  
  # loss when predict should be 0 and 1 respectively
  loss[obs1==0] <- -log(1 - pmin(pred1[obs1==0],1-bounds))
  loss[obs1==1] <- -log(pmax(pred1[obs1==1],bounds))
  
  # adjust for prediction exact matches
  ids_same <- which(pred1==obs1) 
  loss[ids_same] <- 0
  
  risk_est <- mean(loss)
  return(risk_est)
}

# cross validated risk for predicting with mean
cv_risk_f_mean <- function(targets,features,folds){
  
  risk <- numeric(length(folds))
  
  for (i in c(1:length(folds))){
    train_targets  <- targets[folds[[i]]$training_set,]
    train_features <- features[folds[[i]]$training_set,]
    
    test_targets  <- targets[folds[[i]]$validation_set,]
    test_features <- features[folds[[i]]$validation_set,]
    
    preds <- predict_f_mean(train_tgts = train_targets,test_ftrs = test_features)
    
    risk[i] <- log_loss_score(preds,test_targets)
    
  }
  
  return(risk)
}
