# script to analyze design matrix by experiment

library('matrixStats')
# library('reshape2')
library(daewr)

matrizModeloCuadratico <- function(matriz_diseno){
  # construye a partir de una matriz de diseno la matriz del modelo
  # correspondiente a un modelo cuadratico completo
  
  A <- matriz_diseno
  r <- NROW(A)
  m <- NCOL(A)
  mSec <-m*(m+1)/2
  
  eff_main <- colnames(A)
  eff_sec <- character(length=mSec)
  B <- matrix(0,nrow=r,mSec)
  k <-1
  for(i in 1:m){
    B[,k]<-A[,i]*A[,i]
    eff_sec[k] <-paste(eff_main[i],eff_main[i],sep = "")
    k <- k+1
  }
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      B[,k]<-A[,i]*A[,j]
      eff_sec[k] <-paste(eff_main[i],eff_main[j],sep = "")
      k <- k+1
    }
  }
  colnames(B)<-eff_sec
  C <- cbind(A,B)
  A <- as.matrix(A)
  B <- as.matrix(B)
  C <- as.matrix(C)
  data = list(Main=A,SeconOrder=B,Model=C)
  return(data)
}




matrizDisenoDef<-function(m){
  #setwd("C:/Users/Pablo/Documents/TESIS2/R/DisenosDef")
  setwd("~/Documents/Tesis/Latex/R/DisenosDef")
  r <- 2*m+1
  
  name <- paste("JN",m,"F",r,"R",".csv",sep = "")
  
  A <-read.csv(name)
  mats <- matrizModeloCuadratico(A)
  
  return(mats)
}

disenoJones <- function(m){
  
  A <- DefScreen(m)
  mats <- matrizModeloCuadratico(A)
  
  return(mats)
}
