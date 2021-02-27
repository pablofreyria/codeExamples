# seleccion de modelos con jones
library("reshape2")  
library("ggplot2")
library("matrixStats")
library("Matrix")
library("AICcmodavg")
source('~/TESIS2/R/diseñoPorFactnum.R')

m <- 6
parMax <- 2*m

#modelo
mat<- matrizDisenoDef(m)
mat <- cbind(rep(1,NROW(mat$Model)),mat$Model)
colnames(mat)[1]<-"mu"
#mu,x1,x2,x3,x4, x11,x22,x33,x44, x12,x13,x14,x23,x24,x34 #
beta <- numeric(1+m + m*(m+1)/2)
beta[1] <- 20
beta[which(colnames(mat)=="A")] <- 4
beta[which(colnames(mat)=="B")] <- 3
beta[which(colnames(mat)=="C")] <- -2
beta[which(colnames(mat)=="D")] <- -1
beta[which(colnames(mat)=="BC")] <- 5
beta[which(colnames(mat)=="AA")] <- 6

y <- mat%*%beta
set.seed(1)
eps <- rnorm(NROW(mat),0,1)
y1 <- y+eps

## seleccion del modelo ####
id_x <- matrix(0,m,m+1)
id_x[1,] <- grep("A",colnames(mat))
id_x[2,] <- grep("B",colnames(mat))
id_x[3,] <- grep("C",colnames(mat))
id_x[4,] <- grep("D",colnames(mat))
id_x[5,] <- grep("E",colnames(mat))
id_x[6,] <- grep("F",colnames(mat))


# Build all possible models
gc()
est_models <- matrix(0,nrow=1,ncol=length(beta)-1)
notValid <- 0
for(i in c(1:(2^(length(beta)-2)))){
  conf <- as.integer(intToBits(i))[c(1:(length(beta)-1))]
  if (sum(conf)<=parMax-1){
    valid <- TRUE
    for(j in c(1:m)){
      main <- (conf[id_x[j,1]-1]==0)
      sec <- (sum(conf[id_x[j,-1]-1])>0)
      if(sec&&main){
        valid <- FALSE
      }
    }
    if(valid){
      est_models <- rbind(est_models,conf) 
    }else{
      notValid <- notValid
    }
  }else{
    notValid <- notValid + 1 
  }
}

# result of possible models
est_models <- est_models[c(1:103523),]

est_models <- cbind(rep(1,nrow(est_models)),est_models)
# all_models <- Matrix(0,nrow = 2^(length(beta)-1),ncol = length(beta)-1,sparse = TRUE)
# for(i in c(1:length(beta-1))){
#   ceros <- rep(0,2^(i-1))
#   unos  <- rep(1,2^(i-1))
#   bloque <- c(ceros,unos)
#   all_models[,length(beta)-i] <- rep(bloque,2^(length(beta)-i-1))
# }
# # quitar los modelos que tienen mas de 7 parametros
# all_models <- cbind(rep(1,nrow(all_models)),all_models)
# est_models <- all_models[which(rowSums(all_models)<=parMax),]
# 
# # quitar los que tengan interaccion o termino cuadratico pero no incluyan al lineal
# for(i in c(1:m)){
#   id_principal <- which(est_models[,id_x[i,1]]==0)
#   id_inter <- which(rowSums(est_models[,id_x[i,-1]])>0)
#   id_remove <- intersect(id_principal,id_inter)
#   est_models <- est_models[-id_remove,]
# }
# guardar el MSE para cada modelo, algunos no jalan (singular)
### Aca empieza el data ValidModelsJones6fJones ####
MSE <- numeric(nrow(est_models))
AIC <- numeric(nrow(est_models))
id_remove <- numeric()
for(i in c(1:nrow(est_models))){
  id_model <- as.logical(est_models[i,])
  X <- mat[,id_model]
  if(min(eigen(t(X)%*%X)$values)<1e-8){
    id_remove <- c(id_remove,i)
  }else{
    reg <- lm(y1 ~ X -1)
    y_hat <- reg$fitted.values
    AIC[i] <- AICc(reg)
    #AIC[i] <- AIC(reg)
    #AIC[i] <- AIC[i] + 2*sum(id_model)*(sum(id_model)+1)/(NROW(X)-sum(id_model)-1)
    MSE[i] <- sum((y1-y_hat)^2)/(2*m+1-sum(est_models[i,]))
  }
}
#MSE <- MSE[-id_remove]
#est_models <- est_models[-id_remove,]

fact_number <- rowSums(est_models)
# seleccionar modelo
MSE_min <- numeric(parMax)
Var_MSE_min <- numeric(parMax)
AIC_min <- numeric(parMax)
nr_facts <- numeric(parMax)
best_model_mse <- matrix(0,parMax,ncol = ncol(est_models))
best_model_aic <- matrix(0,parMax,ncol = ncol(est_models))
for(i in c(1:parMax)){
  MSE_min[i] <- min(MSE[which(fact_number==i)])
  Var_MSE_min[i] <- var(MSE[which(fact_number==i)])
  AIC_min[i] <- min(AIC[which(fact_number==i)])
  nr_facts[i] <- i
  best_model_mse[i,] <- est_models[which(MSE==MSE_min[i])[1],]
  best_model_aic[i,] <- est_models[which(AIC==AIC_min[i])[1],]
}
plot(nr_facts,AIC_min)

Var_MSE_min[1] <- 0 # solo hay modelo con solo la media
Var_MSE_min <- Var_MSE_min/max(Var_MSE_min)*max(MSE_min)
final_data <- as.data.frame(cbind(nr_facts,MSE_min,Var_MSE_min))
plot_data <- melt(final_data,id.vars=1)

fig<-ggplot(plot_data, aes(x=nr_facts,y=value, group = variable, colour = factor(variable)))+ 
  geom_point()+ geom_line()+ scale_x_continuous(breaks=seq(0,parMax,1))+
  scale_color_manual(labels = c("MSE mejor modelo","Varianza de MSE"),values=c("darkblue","black"))+
  labs(x="Número de parámetros",y="Error cuadrado medio ")+ theme_bw()+
  theme(legend.position = "bottom",axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+labs(colour = " ")

aic_data <- as.data.frame(cbind(nr_facts[c(1:9)], AIC_min[c(1:9)]))
names(aic_data) <- c("Parametros", "AICc")
fig2 <- ggplot(aic_data,aes(x=Parametros,y=AICc)) + theme_bw() + geom_line() + 
  labs(x="Número de parámetros",y=expression(AIC[c]) ) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))+
  scale_x_continuous(breaks=seq(0,9,1))
print(fig2)
  

# wd <- getwd()
# setwd("~/TESIS2/Latex/figure")
# png(filename = "ejJones_selec_modelo_AIC.png",width = 800,height = 480)
# print(fig2)
# dev.off()
# setwd(wd)
