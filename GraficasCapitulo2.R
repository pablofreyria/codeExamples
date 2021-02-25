# Graficas del ejemplo en el capitulo 2 
# library(tikzDevice)
# options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}",
# "\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"))
# 
 library(ggplot2)
 library(mvtnorm)
 library(MASS)

### funcion que grafique los residuos
plot_residuos <- function(dataset,xlab,ylab){
  # hace scatter plot de un data set de la forma [x,y]
  xname <- dataset[,1]
  yname <- dataset[,2]
  fig<-ggplot(data=dataset,aes(x=xname,y=yname)) + geom_point(size=4,shape=20)+
    theme_classic()+labs(x=xlab,y=ylab)+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold"))+
    scale_x_continuous(expand = c(0, 0), limits = c(min(xname)-1,max(xname)+1)) + 
    scale_y_continuous(expand = c(0, 0),limits = c(min(yname)-1,max(yname)+1))
  
  return(fig)  
}

# Mala distribución de los residuos ####
# Modelo 
n <- 50

x <- sample(seq(0,20,0.5),n,replace = TRUE)
x2 <- x^2
y <- 5+3*x+0.5*x^2

# Residuos buenos
eps <- rnorm(n,0,4)

y1 <- y + eps
reg <- stats::lm(y1 ~x+x2)
res <-reg$residuals
y_hat <- reg$fitted.values
plot_func <- as.data.frame(cbind(y_hat,res))

fig  <- plot_residuos(plot_func,expression(widehat(y)),
                      expression(y-widehat(y)))
wd <- getwd()
setwd("~/TESIS2/Latex/figure")
png(filename = "residuos_bien.png",width = 800,height = 480)
print(fig)
dev.off()
setwd(wd)

# Mala forma funcional
eps <- rnorm(n,0,5)

y1 <- y + eps
reg <- stats::lm(y1 ~x)
res <-reg$residuals
y_hat <- reg$fitted.values
plot_func <- as.data.frame(cbind(y_hat,res))

fig  <- plot_residuos(plot_func,expression(widehat(y)),
                      expression(y-widehat(y)))
wd <- getwd()
setwd("~/TESIS2/Latex/figure")
png(filename = "residuos_cuad.png",width = 800,height = 480)
print(fig)
dev.off()
setwd(wd)


# Heteroscedasticidad vs. yhat /trachtenberg
sigma <- matrix(0,nrow=n,ncol = n)
#diag(sigma) <- 0.5+0.5*floor(c(1:n)/10)
aux <-x
aux[x<6] <-1
diag(sigma) <- aux
eps <- mvrnorm(1,rep(0,n),sigma)
y1 <- y+eps
reg <- stats::lm(y1 ~ x+x2)
res <-reg$residuals
y_hat <- reg$fitted.values
plot_func <- as.data.frame(cbind(y_hat,res))

fig  <- plot_residuos(plot_func,expression(widehat(y)),
                      expression(y-widehat(y)))
print(fig)

wd <- getwd()
setwd("~/TESIS2/Latex/figure")
png(filename = "residuos_var.png",width = 800,height = 480)
print(fig)
dev.off()
setwd(wd)
print(fig)

# Autocorrelacion  pos /neg vs. orden 
sigma <- matrix(0,nrow=n,ncol = n)
diag(sigma) <- 1
id_all <- c(1:n^2)
id_diag <-seq(1,n^2,by = n+1)
ids <- setdiff(id_all,id_diag)
corr_vec <- -0.75
for(i in c(1:n)){
  for(j in c(i:n)){
    sigma[i,j] <- corr_vec^(j-i)
    sigma[j,i] <- sigma[i,j]
  }
}
eps <- mvrnorm(1,rep(0,n),sigma)
y1 <- y + eps
reg <- stats::lm(y1 ~ x+x2)
res <-reg$residuals
order <- c(1:n)
plot_func <- as.data.frame(cbind(order,res))

fig  <- plot_residuos(plot_func,"Orden", expression(e["i"]))
print(fig)
wd <- getwd()
setwd("~/TESIS2/Latex/figure")
png(filename = "residuos_orden_corrNeg.png",width = 800,height = 480)
print(fig)
dev.off()
setwd(wd)



# Cambio de l varianzaen el tiempo
sigma <- matrix(0,nrow=n,ncol = n)
diag(sigma) <- c(1:n)/2
eps <- mvrnorm(1,rep(0,n),sigma)
y1 <- y+eps
reg <- stats::lm(y1 ~ x+x2)
res <-reg$residuals
y_hat <- reg$fitted.values
plot_func1 <- as.data.frame(cbind(y_hat,res))
plot_func2 <- as.data.frame(cbind(c(1:n),res))

fig1  <- plot_residuos(plot_func1,expression(widehat(y)),
                      expression(y-widehat(y)))
fig2  <- plot_residuos(plot_func2,"Orden",expression(y-widehat(y)))
print(fig1)
print(fig2)

wd <- getwd()
setwd("~/TESIS2/Latex/figure")
png(filename = "residuos_ordenVar_y.png",width = 800,height = 480)
print(fig1)
dev.off()
setwd(wd)
print(fig)

# Res vs. res-1 (parte de autocorr))
sigma <- matrix(0,nrow=n,ncol = n)
diag(sigma) <- 1
id_all <- c(1:n^2)
id_diag <-seq(1,n^2,by = n+1)
ids <- setdiff(id_all,id_diag)
corr_vec <- 0.7
for(i in c(1:n)){
  for(j in c(i:n)){
    sigma[i,j] <- corr_vec^(j-i)
    sigma[j,i] <- sigma[i,j]
  }
}
eps <- mvrnorm(1,rep(0,n),sigma)
y1 <- 5+3*x + eps
reg <- stats::lm(y1 ~ x)
res <-reg$residuals
res_i <- res[c(2:n)]
res_minus <- res[c(1:(n-1))]

plot_func <- as.data.frame(cbind(res_minus,res_i))

fig  <- plot_residuos(plot_func,expression(e["i-1"]), expression(e["i"]))
print(fig)
wd <- getwd()
setwd("~/TESIS2/Latex/figure")
png(filename = "residuos_autocorr.png",width = 800,height = 480)
print(fig)
dev.off()
setwd(wd)


# Graficas normales de residuos (qqplot) ####
# Modelo 
n <- 50
x <- sample(seq(0,20,0.5),n,replace = TRUE)
y <- 2+4*x
eps <- rnorm(n,0,2)

y1 <- y+eps
reg <- lm(y1~x)
res <- reg$residuals
# distribucion empiríca con la correccion de yates 

res_sort <- sort(res)
dist_emp <- (c(1:n)-0.5)/n
res_teo <- qnorm(dist_emp)

# grafica normal
y_marks <- seq(1,n,6)
y_tik <- unname(res_teo[y_marks])
y_labels <- dist_emp[y_marks]
m_norm <- diff(res_teo[c(1,n)])/diff(res_sort[c(1,n)])
fig<-ggplot(as.data.frame(cbind(res_sort,dist_emp,res_teo)),aes(x=res_sort,y=res_teo))+
  geom_point()+geom_abline(slope =m_norm,intercept = 0)+ theme_bw()+
  scale_y_continuous(breaks=y_tik,labels = round(y_labels,2))+
  labs(x="Residuos observados",y="Distribución normal")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=22,face="bold"))
print(fig)
wd <- getwd()
setwd("~/TESIS2/Latex/figure")
png(filename = "graf_normal.png",width = 800,height = 480)
print(fig)
dev.off()
setwd(wd) 
print(fig)

# grafica cuantil cuantil
#para la linea del qqplot
vec_quant <- c(0.2,0.8)
coord_y <- quantile(reg$residuals,vec_quant)
coord_x <- qnorm(vec_quant)
m <- diff(coord_y)/diff(coord_x)
b <- coord_y[1]-m*coord_x[1]
plot_data <- as.data.frame(cbind(res_teo,res_sort))
names(plot_data) <- c("Teórica","Empírica")

fig<-ggplot(plot_data,aes(x=Teórica,y=Empírica))+geom_point()+ theme_bw()+
  geom_abline(slope = m,intercept = b) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=22,face="bold"))
print(fig)
wd <- getwd()
setwd("~/TESIS2/Latex/figure")
png(filename = "graf_qq.png",width = 800,height = 480)
print(fig)
dev.off()
setwd(wd) 
print(fig)


