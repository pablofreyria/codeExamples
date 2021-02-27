# Bayesian factor screening applied to example in Jones' paper


# leer los modelos simulados en el ejercicio anterior
#est_models <- matrix()
#mat <- matrix()
#y1 <- numeric()
#load("ValidModelsJones6f.RData")

# parametros del modelo de cribado
apriori <- 0.25
k <- 6
n_models <- NROW(est_models)
gam <- 3
n <-NROW(y1)
beta_0 <- sum(y1)/n
S_beta0 <- t(y1-beta_0)%*%(y1-beta_0)

# calcular la probabilidad apriori del modelo con base en el n?mero de factores
p_m <- numeric(length = n_models)
for(i in 1:n_models){
  f_i <- sum(est_models[i,c(1:k+1)])
  t_i <- sum(est_models[i,])-1
  
  gamma_mat <- matrix(0,nrow = t_i+1,ncol = t_i+1)
  diag(gamma_mat)<-1 
  diag(gamma_mat)[1]<- 0
  gamma_mat <- 1/gam^2*gamma_mat
  
  x_mat <- mat[,as.logical(est_models[i,])]
  beta_hat <- solve(gamma_mat + t(x_mat)%*%x_mat)%*%t(x_mat)%*%y1
  S_beta <- t(y1 -x_mat%*%beta_hat)%*%(y1 -x_mat%*%beta_hat)  
  
  p_m[i] <- (apriori/(1-apriori))^f_i*gam^(-t_i)*
            sqrt(n)/sqrt(det(gamma_mat + t(x_mat)%*%x_mat))*
            ((S_beta +t(beta_hat)%*%gamma_mat%*%beta_hat)/S_beta0)^(-(n-1)/2)
}

C <- sum(p_m)
p_m <- p_m/C

# suma de probabildad de modelos que incluyan al factor i
prob_facts <- numeric(length = k)
for(i in c(1:k)){
  ids <- which(est_models[,(i+1)]==1)
  prob_facts[i] <- sum(p_m[ids])
}


plot_bayes <- as.data.frame(cbind(c(1:k),prob_facts))
names(plot_bayes) <- c("factor","prob")

fig2<-ggplot(data = plot_bayes,aes(x=factor,y=prob)) + geom_col(fill="gray85",colour="black") + theme_bw()+
  labs(x="Factor",y="Probabilidad" ) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))+
  scale_x_continuous(breaks=seq(1,6,1),labels = c("A","B","C","D","E","F"))

print(fig2)

### identificar el modelo ####
# mayor probabilidad con los primeros num_fact factores
num_fact <- 4
fact_act <- rowSums(est_models[,c(1:k+1)])
p_m_max_fact <- max(p_m[which(fact_act==num_fact)])

best_model <- est_models[which(p_m==p_m_max_fact),]


# wd <- getwd()
# setwd("~/TESIS2/Latex/figure")
# png(filename = "graf_Bayes_Screen.png",width = 800,height = 480)
# print(fig2)
# dev.off()
# setwd(wd)
   

