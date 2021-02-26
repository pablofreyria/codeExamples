library("reshape2")
library("directlabels")
# Toy example of a complete model analysis with 2 factors and 3 levels each

#eps <- rnorm(confs,0,2)
load("DataEjemploCap2.RData")
eps <- dataEj[[1]]

efectos <- c(20,-4,0,4,-2,4,-2,-1,1,0,1,-3,2,0,2,-2)
facts <- c("u","t1","t2","t3","b1","b2","b3","tb11","tb12","tb13",
           "tb21","tb22","tb23","tb31","tb32","tb33")
names(efectos)<-facts
confs <- 3*3*4
nres  <- 1+1+2*3

X_inter <- matrix(0,nrow = confs,ncol = 9)
X_main <- matrix(0,nrow = confs,ncol = 1+3+3)
X_res <- matrix(0,nrow = nres,ncol = NCOL(X_main)+NCOL(X_inter))
# build design matrix ####
# main terms
X_main[,1]<-1
for(i in c(1:3)){
  X_main[c((12*(i-1)+1):(12*i)),i+1] <- 1
  for(j in c(1:3)){
    X_main[c((12*(i-1)+1+4*(j-1)):(12*(i-1)+4*j)),3+j+1]<-1
  }
}
# interaction terms
for(i in c(1:9)){
  X_inter[c((4*(i-1)+1):(4*i)),i]<-1
}

# Build model and "observed" response ####
X_model <- cbind(X_main,X_inter)
y <- X_model%*%efectos
y1 <- y+eps
y_reg <- rbind(y1,matrix(0,nres,1)) 
# promedio por config
y_bar_conf <- numeric(length = 9)
for(i in c(1:9)){
  y_bar_conf[i]<- mean(y_reg[c((4*i-3):(4*i))])
}

# Restricciones para calcular la solucion
for(i in c(1:5)){
  X_res[i,c((3*i-1):(3*i+1))]<-1
}
for(i in c(1:3)){
  X_res[i+5,7+i+3*c(0:2)] <- 1
}
# Solucion 
X_final <- rbind(X_model,X_res)

beta <- solve(t(X_final)%*%X_final)%*%t(X_final)%*%y_reg

y_hat <- X_model%*%beta
res <- y1-y_hat

# grafica cuantil-cuantil de los efectos ###
beta_s <- sort.int(beta[c(2:length(beta))],index.return = TRUE)
id_sort <- beta_s$ix
beta_s <- beta_s$x
dist_beta <- (c(1:length(beta_s))-0.5)/length(beta_s)
beta_teo <- qnorm(dist_beta)
m <- diff(quantile(beta_s,c(0.2,0.8)))/diff(qnorm(c(0.2,0.8)))
plot_data <- as.data.frame(cbind(beta_teo,beta_s))
names(plot_data) <- c("teo","emp")

### seccion ####
fig<-ggplot(plot_data,aes(x=teo,y=emp))+geom_point(shape=15)+  theme_bw()+
  geom_text(size=6,aes(label=facts[2:length(beta)][id_sort]),vjust=1)+ 
  geom_abline(slope = m,intercept = 0) + labs(x="Dist. teórica",y="Dist. empírica")+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=24,face="bold"))
wd <- getwd()
setwd("~/TESIS2/Latex/figure")
png(filename = "graf_qq_efectos.png",width = 800,height = 480)
print(fig)
dev.off()
setwd(wd)

# pruebas F 
y_bar <- mean(y1)
id_t <- c(2:4)
SS_t <- 0
for(i in id_t){
  id_fact <-which(X_model[,i]==1)
  SS_t<-SS_t+length(id_fact)*(mean(y1[id_fact])-y_bar)^2
}

id_b <- c(5:7)
SS_b <-0
for(i in id_b){
  id_fact <-which(X_model[,i]==1)
  SS_b<-SS_b+length(id_fact)*(mean(y1[id_fact])-y_bar)^2
}
SS_tb<-0
SS_e <-0
for(i in id_t){
  for(j in id_b){
    id.t <- which(X_model[,i]==1)
    id.b <- which(X_model[,j]==1)
    id_tb <- intersect(id.t,id.b)
    SS_tb<- SS_tb+length(id_tb)*(mean(y1[id_tb])-mean(y1[id.t])-mean(y1[id.b])+y_bar)^2
    SS_e <- SS_e + sum((y1[id_tb]-mean(y1[id_tb]))^2)
  }
}
SS_Total <- sum((y1-y_bar)^2)

Ftest <- as.data.frame(rbind(SS_t,SS_b,SS_tb,SS_e,SS_Total))
rownames(Ftest) <- c("Tau","Beta","TauBeta","Residual","Total")
df <- c(2,2,4,27,35)
MS <- Ftest[,1]/df
F_score <- MS/MS[4]
p_value <- pf(F_score, df,df[4],lower.tail = FALSE)
Ftest <- cbind(Ftest,df,MS,F_score,p_value)
colnames(Ftest) <- c("Suma de cuad","Grados Lib","Cuad medios",
                     "Estadistico","Significancia")

# graficas del ejemplo
# primera gráfica sobre la respuesta media por factor level
id_t <- c(2:4)
id_b <- c(5:7)
nivel_t <- c(0,1,3)
nivel_b <- c(1,3,5)
eje_t <- matrix(0,nrow = length(nivel_t),ncol = length(nivel_b))
contour_mat <- matrix(0,nrow = length(nivel_t)*length(nivel_b),ncol = 3)
k<-1
for(i in c(1:length(nivel_t))){
  for(j in c(1:length(nivel_b))){
    id.t <- which(X_model[,id_t[i]]==1)
    id.b <- which(X_model[,id_b[j]]==1)
    id_inter <- intersect(id.t,id.b)
    eje_t[i,j] <- mean(y1[id_inter])
    contour_mat[k,1]<-nivel_t[i]
    contour_mat[k,2]<-nivel_b[j]
    contour_mat[k,3]<-eje_t[i,j]
    k<-k+1
  }
}
cont_data <-as.data.frame(contour_mat)
names(cont_data)<- c("Tau","Beta","Value")
plot_data <- as.data.frame(cbind(nivel_t,eje_t))
names(plot_data) <- c("Tau","BETA1","BETA2","BETA3")
plot_data <- melt(plot_data,id.vars = 1)

# contour plot
fig<-ggplot(cont_data, aes(x=Tau, y=Beta, z = Value)) + theme_bw()+
  stat_contour(aes(colour = ..level..))+
  theme(axis.text=element_text(size=13),axis.title=element_text(size=24,face="bold"))
fig<-direct.label(fig,method="bottom.pieces")
wd <- getwd()
setwd("~/TESIS2/Latex/figure")
png(filename = "EJ2_curvas_nivel.png",width = 800,height = 480)
print(fig)
dev.off()
setwd(wd)

# relacion y detección de interacción
fig<-ggplot(plot_data,aes(x=Tau,y=value,group=variable,color=variable))+ geom_line() + 
  geom_point()+ scale_color_manual(values=c('black','#00BFC4','#F8766D') , 
                labels = c(expression( beta==1),expression( beta==3),expression( beta==5)))+
  labs(x="Tau",y=expression(bar(y)),colour =" ")+ theme_bw()+
  theme(legend.key = element_rect(size = 2),legend.key.size = unit(2, 'lines'),
        legend.text = element_text(size=24), legend.position = "bottom",
        axis.text=element_text(size=13),axis.title=element_text(size=24,face="bold"))

wd <- getwd()
setwd("~/TESIS2/Latex/figure")
png(filename = "EJ2_interaccion.png",width = 800,height = 480)
print(fig)
dev.off()
setwd(wd)

