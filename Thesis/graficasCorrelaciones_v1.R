library(ggplot2)
library(reshape2)
corrCuad <-function(m){
   y = 1/3-1/(m-1)
   return(y)
}
corrEff <- function(m){
  y = sqrt((2*m+1)/((3*m-3)*(m-2)))
  return(y)
}

m <- c(4:40)
m2<- seq(4,40)
cuad <- corrCuad(m)
inter <- corrEff(m)
aux1 <- matrix(0,length(m),1)
aux2 <- matrix(1,length(m2),1)
aux3 <- rbind(aux1,aux2)

data <- as.data.frame(cbind(c(m,m2),c(cuad,inter),aux3))
names(data)<-c("Factores","Correlacion","Tipo")


correlacionesCuad <- ggplot(data=data, aes(x=Factores,y=Correlacion, group=Tipo,colour=factor(Tipo),
                      shape=factor(Tipo)))+
  geom_line()+geom_point(aes(fill=factor(Tipo)),size=3) +
  guides(size=FALSE,fill=FALSE,linetype = guide_legend(ncol = 2,keywidth=4))+  
  scale_shape_manual(labels = c("Correlación entre efectos cuadráticos", 
                    "Correlación entre efecto cuadrático e interacción"),
                     values=c(16,17)) +
  scale_fill_manual(labels = c("Correlación entre efectos cuadráticos", 
                               "Correlación entre efecto cuadrático e interacción"),
                    values=c("darkblue","black"))+
  scale_color_manual(labels = c("Correlación entre efectos cuadráticos", 
                     "Correlación entre efecto cuadrático e interacción"),
                       values=c("darkblue","black"))+
  theme(legend.position = "bottom", legend.key = element_rect(size = 3),
        axis.title.y = element_text(margin = margin(0,0.6,0,0,unit="cm")),
        axis.title.x = element_text(margin = margin(0.5,0,0,0,unit="cm")),
        legend.key.size = unit(2, 'lines'),legend.text = element_text(size=18),
        axis.text=element_text(size=18), 
        axis.title=element_text(size=18,face="bold"))+
  coord_cartesian(xlim=c(3, 40), ylim=c(0, 0.8))+
  scale_x_continuous(breaks=seq(0,40,5))+
  scale_y_continuous(breaks = seq(0,1,0.1)) + 
  labs(x="Factores", y="Valor absoluto de la correlación",colour=" ",shape=" ",
       shape_size=" ")

wd <- getwd()
setwd("~/TESIS2/Latex/figure")
png(filename = "correlacionesCuad2.png",width = 900,height = 480)
print(correlacionesCuad)
dev.off()
setwd(wd)

