### MCMC sampling method applied to a sequential screening method####

# Box GEP, Tiao GC. Bayesian Inference in Statistical Analysis. John Wiley and Sons: New York, 1992.
# DOI : 10.1002/9781118033197 

# y1<-numeric()
# mat<-numeric()
load("Ejercicio_Cap_Jones_Completo.RData")
rm(list=setdiff(ls(),c("y1","mat")))
MCMC_mio <-function(y1,mat,ids){
#Parametros y datos
burnin <- 1000
step <- 10
epsilon <- 0.00
k <- length(ids)
n <- 21000
n_exp <- nrow(mat)
sigma <- numeric(n)
beta <- matrix(0,nrow=k,ncol=n)
X <- mat[,ids]
gam <- 1
# distribucion a priori (Ideas de Box-Meyer)
sigma[1] <- exp(runif(1))^2
#sigma[1] <-1
beta[,1] <- rnorm(k,0,sqrt(sigma[1]))
p0_sigma <- 0
p0_beta <- log(dnorm(beta[,1],0,gam*sqrt(sigma[1])))

# verosimilitud y a posteriori
# res_sq <- as.numeric(t(y1-X%*%beta[,1])%*%(y1-X%*%beta[,1]))
L <- Log_lik_normal(y1,mat,ids,beta[,1],sqrt(sigma[1]))
p_sigma <- p0_sigma + L
p_beta <- p0_beta + L
L <- Log_lik_normal(y1,mat,ids,beta[,1],sqrt(sigma[1]))
p_beta <- p0_beta + L


#iteracion
for(i in 2:n){
  # SIGMA: A priori
  aux_beta <- beta[,(i-1)]
  aux_sigma <- exp(runif(1))^2
  p0_sigma <- 0
  #aux_sigma<-1
  #p0_sigma <- log(1/aux_sigma)
  
  #SIGMA: A posteriori y verosimilitud
  res_sq <- as.numeric(t(y1-X%*%aux_beta)%*%(y1-X%*%aux_beta))
  L <-  Log_lik_normal(y1,mat,ids,aux_beta,sqrt(aux_sigma))
  p1_sigma <- L+p0_sigma  
  
  if(runif(1) < exp(p1_sigma-p_sigma)){
    sigma[i] <- aux_sigma
    p_sigma <- p1_sigma
  }else{
    sigma[i] <- sigma[i-1]
  }
  #BETA, ITERANDO
  for(j in c(1:k)){
    aux_beta[j] <- rnorm(1,beta[j,(i-1)],sqrt(sigma[i]))
    p0_beta <- log(dnorm(aux_beta[j],beta[j,(i-1)],sqrt(sigma[i])))
    
    res_sq <- as.numeric(t(y1-X%*%aux_beta)%*%(y1-X%*%aux_beta))
    L <-  Log_lik_normal(y1,mat,ids,aux_beta,sqrt(sigma[i]))
    p1_beta <- L + p0_beta
    
    if(runif(1)<exp(p1_beta-p_beta[j])){
      beta[j,i] <- aux_beta[j]
      p_beta[j] <- p1_beta
    }else{
      beta[j,i] <- beta[j,i-1]
      aux_beta[j] <- beta[j,i-1]
    }
  }
}

beta_fin <- beta[,-c(1:burnin)]
sigma_fin <- sigma[-c(1:burnin)]

beta_fin <- beta_fin[,seq(1,(n-burnin),step)]
sigma_fin <- sigma_fin[seq(1,(n-burnin),step)]
n_obs_fin <- length(beta_fin[1,])

beta_pos <- rowSums((beta_fin>epsilon)*1)
beta_neg <- rowSums((beta_fin<(-epsilon))*1)

p_beta_pos <- (beta_pos+1)/(n_obs_fin + 2)
p_beta_neg <- (beta_neg+1)/(n_obs_fin + 2)
odds_beta_pos <- (p_beta_pos)/(1-p_beta_pos)
odds_beta_neg <- (p_beta_neg+1)/(1-p_beta_neg)

odds <- list(Positive = odds_beta_pos,Negative=odds_beta_neg, Beta = beta_fin, sigma = sigma_fin,path=beta)
return(odds)
}

MLE_beta <-function(y1,mat,ids){
  X <- mat[,ids]
  beta_hat <- solve((t(X)%*%X),t(X)%*%y1)
  return(beta_hat)
}
Log_lik_normal <- function(y1,mat,ids,beta,sigma){
  X <- mat[,ids]
  n <- nrow(X)
  res_sq <- as.numeric(t(y1-X%*%beta)%*%(y1-X%*%beta))
  L<-log(sigma)*(-n)+(-1/(2*sigma^2)*(res_sq))
  return(L)
}

###1 paso todos los principales
ids1 <-c(1:7)
paso1 <- MCMC_mio(y1,mat,ids1)

rbind(colnames(mat)[ids1],paso1$Positive,paso1$Negative)

#Paso 2
# Significativos del pasado c(1,2,3,4,5)
  ids2 <- c(1,which(colnames(mat)=="A"),which(colnames(mat)=="B"),which(colnames(mat)=="C"),
            which(colnames(mat)=="D"),which(colnames(mat)=="E"))
  ids2 <- c(ids2,which(colnames(mat)=="AA"),which(colnames(mat)=="BB"),which(colnames(mat)=="CC"),
            which(colnames(mat)=="DD"),which(colnames(mat)=="EE"))

  paso2 <- MCMC_mio(y1,mat,ids2)
  rbind(colnames(mat)[ids2],paso2$Positive, paso2$Negative)
  
#Paso 2 con efectos de interaccion
  ids21 <- c(1,which(colnames(mat)=="A"),which(colnames(mat)=="B"),which(colnames(mat)=="C"),
            which(colnames(mat)=="D"))
  ids21 <- c(ids21,which(colnames(mat)=="AB"),which(colnames(mat)=="AC"),which(colnames(mat)=="AD"),
             which(colnames(mat)=="BC"),which(colnames(mat)=="BD"),which(colnames(mat)=="CD"))

  paso21 <- MCMC_mio(y1,mat,ids21)
  rbind(colnames(mat)[ids21],paso21$Positive, paso21$Negative)

#Paso 2 final
  ids22 <- c(1,which(colnames(mat)=="A"),which(colnames(mat)=="B"),which(colnames(mat)=="C"),
             which(colnames(mat)=="D"))
  ids22 <- c(ids22,which(colnames(mat)=="AA"))
  ids22 <- c(ids22,which(colnames(mat)=="BC"),which(colnames(mat)=="CD"))

  paso22 <- MCMC_mio(y1,mat,ids22)
  rbind(colnames(mat)[ids22],paso22$Positive, paso22$Negative)


 #paso 3 de segundo orden con todos los principales
 ids3 <- c(1:7)
 ids3 <- c(ids3,which(colnames(mat)=="AA"),which(colnames(mat)=="BC"),which(colnames(mat)=="EE"))
 paso3 <- MCMC_mio(y1,mat,ids3)

 rbind(colnames(mat)[ids3],paso3$Positive, paso3$Negative)

  #paso 4
  ids4 <- c(1:5)
  ids4 <- c(ids4,which(colnames(mat)=="AA"),which(colnames(mat)=="BC"))
  paso4 <- MCMC_mio(y1,mat,ids4)
  rbind(colnames(mat)[ids4],paso4$Positive, paso4$Negative)
  
# RESUMEN 
t1  <- rbind(colnames(mat)[ids1],paso1$Positive,paso1$Negative)
t2  <- rbind(colnames(mat)[ids2],paso2$Positive, paso2$Negative)
t21 <- rbind(colnames(mat)[ids21],paso21$Positive, paso21$Negative)
t22 <- rbind(colnames(mat)[ids22],paso22$Positive, paso22$Negative)
t3  <- rbind(colnames(mat)[ids3],paso3$Positive, paso3$Negative)
t4  <- rbind(colnames(mat)[ids4],paso4$Positive, paso4$Negative)

# latex
xtable(t(t1), digits = 2)
