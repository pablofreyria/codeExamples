# Correlation in estimated effects from example in Jones' paper

library("reshape2")  
library("ggplot2")
source('~/TESIS2/R/diseñoPorFactnum.R')
  
  
  facts <- c(6:30)
  corrEstimadas <- matrix(0,nrow = length(facts),9)
  
  for(factNum in c(1:length(facts))){
    # leer matriz de diseño y de experimento
    m <- facts[factNum]
    mat<- matrizDisenoDef(m)  # Diseños dados por Jones
    #mat <- disenoJones(m)      # Diseños de R
    # encontrar todas las combinaciones mC3
    all_three <- matrix(0,nrow=choose(m,3),ncol = 3)
    k<-1
    for(i in c(1:(m-2))){
      for(j in c((i+1):(m-1))){
        for(l in c((j+1):m)){
          all_three[k,1] <- i
          all_three[k,2] <- j
          all_three[k,3] <- l
          k <- k+1
        }
      }
    }
    
    # inicializar las correlacioones en 0
    cor_qqss_max <- 0
    cor_qqss_avg <- 0
      cor_qqqs_max <- 0
      cor_qqqs_avg <- 0
    cor_qqst_max <- 0
    cor_qqst_avg <- 0
      cor_stuv_max <- 0
      cor_stuv_avg <- 0
    # iterar sobre todas las formas de elegir 3 parámetros de m
    for(i in c(1:choose(m,3))){
      # hacer la matriz de diseño y de experimento a partir de la de tres factores
      aux_dis_mat <- mat$Main[,all_three[i,]]
      aux_mod_mat <- matrizModeloCuadratico(aux_dis_mat)
      X1 <- cbind(matrix(1,nrow = NROW(aux_mod_mat$Model),ncol=1),aux_mod_mat$Model)
      colnames(X1) <- c("mu", colnames(aux_mod_mat$Model))
      X1 <- as.matrix(X1)
      # calcular matriz de correlaciones en valor absoluto
      norm_eq <- t(X1)%*%X1
      sol <- solve(norm_eq)
      sds <- sqrt(diag(sol))
      matvars <- sds%*%t(sds)
      matcors <- sol/matvars
      matcors <- abs(matcors)
      
      # calcular la correlacion por tipo, la matriz tiene la misma forma->ids fijos
        # entre cuadráticos
        id_qqss<-c(46,47,57)
        cor_qqss_max <-max(cor_qqss_max,matcors[id_qqss])
        cor_qqss_avg <-sum(cor_qqss_avg,matcors[id_qqss])
      
        # entre cuadráticos interaccion factor  comun
        id_qqqs <- c(75,85,76,96,87,97)
        cor_qqqs_max <-max(cor_qqqs_max,matcors[id_qqqs])
        cor_qqqs_avg <-sum(cor_qqqs_avg,matcors[id_qqqs])
      
        # entre cuadráticos interaccion factor no comun
        id_qqst <- c(95,86,77)
        cor_qqst_max <-max(cor_qqst_max,matcors[id_qqst])
        cor_qqst_avg <-sum(cor_qqst_avg,matcors[id_qqst])
        
        # entre interacciones
        id_stuv <- c(79,80,90)
        cor_stuv_max <-max(cor_stuv_max,matcors[id_stuv])
        cor_stuv_avg <-sum(cor_stuv_avg,matcors[id_stuv])
    }
    # poner los valores en la tabla por número de parámetro
      corrEstimadas[factNum,1] <- m
      corrEstimadas[factNum,2] <- cor_qqss_max
      corrEstimadas[factNum,3] <- cor_qqss_avg/(choose(m,3)*3)
      corrEstimadas[factNum,4] <- cor_qqqs_max
      corrEstimadas[factNum,5] <- cor_qqqs_avg/(choose(m,3)*6)
      corrEstimadas[factNum,6] <- cor_qqst_max
      corrEstimadas[factNum,7] <- cor_qqst_avg/(choose(m,3)*3)
      corrEstimadas[factNum,8] <- cor_stuv_max
      corrEstimadas[factNum,9] <- cor_stuv_avg/(choose(m,3)*3)
  } 
  corrEstimadas <- as.data.frame(corrEstimadas)
  names(corrEstimadas) <- c("m","r_qqss_max","r_qqss_avg","r_qqqs_max","r_qqqs_avg",
                            "r_qqst_max","r_qqst_avg","r_stuv_max","r_stuv_avg")

  corrEstimadas_avg <- corrEstimadas[,c(1,3,5,7,9)]
  corrEstimadas_max <- corrEstimadas[,c(1,2,4,6,8)]
 # plot_corr <- melt(corrEstimadas,id.vars = 1)
  plot_corr_avg <- melt(corrEstimadas_avg,id.vars = 1)
  plot_corr_max <- melt(corrEstimadas_max,id.vars = 1)

  fig1<-ggplot(plot_corr_max,aes(x=m,y=value,group=variable,color=variable))+
    geom_line() + geom_point()+
    # scale_color_manual(labels=names(corrEstimadas)
    #                    values = c('black','blue4','darkred',"#009E73",'orange'))+
    scale_fill_manual(labels=names(corrEstimadas))+
    # scale_shape_manual(labels=c('Cuadráticos','Cuad-inter común','Cuad-inter no común',
    #                             'Inter-inter común','Inter-inter no común'),
    #                    values = c('black','blue4','darkred','palegreen4','orangered4'))
    labs(x="Número de factores", y="Corr. Max estimadores",colour=" ",shape=' ',shape_size=" ")+
    scale_x_continuous(breaks=seq(6,facts[factNum],3))+
    scale_y_continuous(breaks = seq(0,1,0.1))+
    theme(legend.position = "bottom", legend.key = element_rect(size = 3),
          axis.title.y = element_text(margin = margin(0,0.6,0,0,unit="cm")),
          axis.title.x = element_text(margin = margin(0.5,0,0,0,unit="cm")),
          legend.key.size = unit(2, 'lines'),legend.text = element_text(size=11),
          axis.text=element_text(size=12),
          axis.title=element_text(size=13,face="bold"))

  # wd <- getwd()
  # setwd("~/TESIS2/Latex/figure")
  # png(filename = "Correlacion_Max_Estimadores_Jones_30.png",width = 800,height = 480)
  # print(fig1)
  # dev.off()
  # setwd(wd)