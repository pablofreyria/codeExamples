
# probability with random guessing
p_ini <- colMeans(targets[,-1]) 
plot(sort(p_ini))
summary(p_ini) #median is ~1 in 1000
# log base seems like a much better choice
logp_ini <- log(p_ini)
plot(sort(logp_ini))
summary(logp_ini)

# see if the classification is: exhaustive and exclusive
summary(rowSums(targets[,-1]))
sum(rowSums(targets[,-1])>1)/nrow(targets)
sum(rowSums(targets[,-1])==0)/nrow(targets)
# it is neither: having one MoA does not exclude from having other and 
# we can have drugs without any MoA. Only 8% have more than 1 and 40% have none


# check that there is no MoA activated with control vehicle
sum(rowSwums(orig_targets[orig_features$cp_type=="trt_cp",])==0)
# and if all drugs with treatment activate a MoA
sum(colSums(orig_targets[orig_features$cp_type=="ctl_vehicle",])==0)
## first rule if(type=control)->P=0, else P=E[p|x,cp_type=1]

## analysis on covariance matrix and its eigen values, singularity issues

# center and scale features
feat_center <- colMeans(features)
feat_var    <- apply(features,2,sd)

feat_scales <- scale(features,
                     center = feat_center,
                     scale = feat_var)

n_feat <- ncol(features)
# analyze eigen values to identify how many meaningful directions
eigen_vals <- eigen(cov(feat_scales))$values
#eigen values quickly drop to zero
plot(eigen_vals)

# see were relationship is flat and validate that eigen is above a threshold
id_delta <- which((eigen_vals[1:(n_feat-1)]-eigen_vals[2:n_feat])<1e-3) # second value is already too close
min(id_delta)
sum(eigen_vals[1:min(id_delta)])/sum(eigen_vals)
# main 148 components are going to be used
# preserves 76% of the standardized features matrix (1/6 of vars) and
# next component is within 1e-3 so it seems that curve is to flat already

#
