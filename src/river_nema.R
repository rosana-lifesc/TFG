library(iNEXT);       # Paquete para analisis de diversidad alfa por intrapolacion y extrapolacion (metodo de Chao et al. 2014)
library(vegan);       # Paquete para el analisis de comunidades
library(ggplot2);     # Paquete para elaboracion de graficas (necesario para graficas de iNEXT)
library(factoextra);  # Paquete para la visualizacion de analisis multivariantes (lo usaremos para los dendrogramas)
library(RColorBrewer)
library(plotly)
library(RVAideMemoire)





getwd()
#setwd("D:/R/rios")
getwd()






## Importamos los datos 

datos <- read.csv2("nema_river.csv")
datos$Tramo = factor(datos$Tramo, levels = c("1", "2", "3"))
datos$Cuenca = factor(datos$Cuenca, levels = c("1", "2", "3", "4", "5"))

cp_csv <- read.csv2("cp+feeding.csv")
CP <- subset(cp_csv, select = "CP")
feeding <- subset(cp_csv, select = "FEEDING")

datos1 <- datos[datos$Tramo != 2, ]


muestra <- rownames(datos1)              # Nombre de las comunidades
nMuestras <- length(muestra)                        # Numero de comunidades
nPlots <- as.vector(t(datos1[1,]))


bichos <- subset(datos1, select = c(9:108))
#muestra <- subset(datos1, select = c(1))

bichos_d <- subset(bichos[25:37,])
bichos_i <- subset(bichos[1:24,])

#ref_bichos <- cbind(muestra, bichos)





#alpha diversity metrics



N0    <- rowSums(bichos > 0)                    # Species richness
H     <- diversity(bichos)                      # Shannon diversity
N1    <- exp(H)                                 # Hill?s 1 
S     <- diversity(bichos, index = "simpson")   # Simpson diversity
N2    <- diversity(bichos, index = "inv")       # Hill?s 2
E21   <- N2/N1                                  # Hill?s Evenness index 


datos1$Richness <- N0

datos1$Shannon  <- H

datos1$Evenness <- E21




###

#plots 

###


str(datos1, list.len=1999)



library(ggpubr)


tiff(file="FIG . Riqueza.tiff", width=8, height=10, units="in", res=800, compression = "lzw")


richness.plot <- ggplot(datos1, aes(x=factor(Tramo), y=Richness, color=factor(Tramo) ))+ geom_boxplot(size=1.5)+
  facet_wrap(~Margen, ncol=2) 

richness.plot <- richness.plot + ggtitle("") + labs(x="", y ="Riqueza") + theme_light() +
  theme(panel.border = element_rect(size= 3, color="gray22")) + 
  theme(axis.text=element_text(size=20, face="bold"),  axis.title=element_text(size=27,face="bold"), 
  axis.ticks=element_line(size=1.5)) + theme( legend.position = "none",strip.text = element_text(size = 15) ) + 
  scale_y_continuous(limits = c(0,15)) + scale_x_discrete(name="", breaks=c(1, 3), labels=c("Cabecera", "Curso Bajo"))


richness.plot

dev.off()



tiff(file="FIG . Shannon", width=8, height=10, units="in", res=800, compression = "lzw")

shannon.plot <- ggplot(datos1, aes(x=factor(Tramo), y=Shannon, color=factor(Tramo) ))+ geom_boxplot(size=1.5)+
  facet_wrap(~Margen, ncol=2) 


shannon.plot <- shannon.plot + ggtitle("") + labs(x="", y ="Shannon") + theme_light() +theme(panel.border = element_rect(size= 3, color="gray22")) + 
  theme(axis.text=element_text(size=20, face="bold"),  axis.title=element_text(size=27,face="bold"), 
        axis.ticks=element_line(size=1.5)) + theme( legend.position = "none",strip.text = element_text(size = 15) ) + scale_y_continuous(limits = c(0,3)) +
        scale_x_discrete(name="", breaks=c(1, 3), labels=c("Cabecera", "Curso bajo"))


shannon.plot

dev.off()



#### ASSESSING EFFECTS! 



prueba0 <- perm.anova(Richness ~ Tramo  , data = datos1,nperm=999)





prueba <- lm(Richness ~ Tramo + Margen , datos1 )

summary(prueba)

anova(prueba)


prueba1 <- lm(Richness ~ Tramo + Margen + Tramo*Margen , datos1 )



prueba2 <- lm(Richness ~ Tramo + Margen + Tramo*Margen + River, datos1)

summary(prueba1)
summary(prueba2)


anova(prueba, prueba1, prueba2)










### DIVERSIYDAD BETA 

library(adespatial)

#1 We transformed matrix by using Hellinger (Please see Blanchet and )

str(bichos)

bichos.hel <- decostand(bichos, "hellinger")## It is including in the following code about beta diversity´s computation. 



#?beta.div# to know about outputs of functions.. 



bicho.var <- beta.div(bichos, "hellinger", nperm=10000,save.D=TRUE)# computing beta diversity as the overall variance.. Legendre and de Caceres 2013


names(bicho.var)# checking object composition 



bdtotal_bichos   <-  bicho.var$beta 



#SStotal    BDtotal 
#28.3570946  0.7876971

# I like a lot our results.. high total variance of community composition .. 0.78.. 

var.beta.bichos <- bicho.var$D # dissimilarity matrix


#How is the effect of agriculture on beta diversity? Here we are going on preliminary results (i.e. focusing on total variance of community composition)


#subseting 

View(datos1)

str(datos1)



datos_d <- subset(datos1[25:37,])
datos_i <- subset(datos1[1:24,])


bichos_d <- subset(bichos[25:37,])
bichos_i <- subset(bichos[1:24,])



### the same.. 


bichos.d.hel <- decostand(bichos_d, "hellinger")## It is including in the following code about beta diversity´s computation. 


bicho.d.var <- beta.div(bichos_d, "hellinger", nperm=10000,save.D=TRUE)# computing beta diversity as the overall variance.. Legendre and de Caceres 2013


names(bicho.d.var)# checking object composition 

#?beta.div# to know about outputs of functions.. 


bdtotal_bichos_d   <-  bicho.d.var$beta 

#SStotal    BDtotal 
#10.1261605  0.8438467 



var.beta.bichos_d <- bicho.d.var$D # dissimilarity matrix






bichos.i.hel <- decostand(bichos_i, "hellinger")## It is including in the following code about beta diversity´s computation. 


bicho.i.var <- beta.div(bichos_i, "hellinger", nperm=10000,save.D=TRUE)# computing beta diversity as the overall variance.. Legendre and de Caceres 2013


names(bicho.i.var)# checking object composition 

#?beta.div# to know about outputs of functions.. 


bdtotal_bichos_i   <-  bicho.i.var$beta 

#SStotal    BDtotal 
#16.0830063  0.6992611


  
  
  
var.beta.bichos_i <- bicho.i.var$D # dissimilarity matrix





#### INTERESANTES RESULTADOS! 






## Marti Anderson's PERMDISP2   Multivariate homogeneity of groups dispersions (variances)

print(bicho.i.var)
indicativo_derecha <- bicho.d.var$SCBD
indicativo_izquierda <- bicho.i.var$SCBD
indicativo_total <- bicho.var$SCBD
#write.csv2(indicativo_derecha, "indicativo_derecha.csv")
# write.csv2(indicativo_izquierda, "indicativo_izquierda.csv")
# write.csv2(indicativo_total, "indicativo_total.csv")

var.beta.bichos <- bicho.var$D # dissimilarity matrix


mod <- betadisper(var.beta.bichos, datos1$Margen, type="centroid")

plot(mod)

names(mod)



#aov.beta.margend <- aov( D ~ Tramo, data = bichos_d)


distances <- mod$distances



datos1$BD <- distances

print(datos1$BD)

#Anova for beta diversity depending on the side of the river (margen)
margend <- subset(datos1, Margen %in% c("Derecha"))
aov.beta.margend <- aov(BD ~ Tramo, data = margend)
summary(aov.beta.margend)
tab_model(aov.beta.margend)


margeni <- subset(datos1, Margen %in% c("Izquierda"))
aov.beta.margeni <- aov(BD ~ Tramo, data = margeni)
summary(aov.beta.margeni)
tab_model(aov.beta.margeni)


tiff(file="FIG . Beta.tiff", width=8, height=10, units="in", res=800, compression = "lzw")


bd.plot <- ggplot(datos1, aes(x=factor(Tramo), y=BD, color=factor(Tramo) )) + geom_boxplot(size=1.5)+
  facet_wrap(~Margen, ncol=2) 


bd.plot <- bd.plot + ggtitle("") + labs(x="", y ="Diversidad beta") + theme_light() +theme(panel.border = element_rect(size= 3, color="gray22")) + 
  theme(axis.text=element_text(size=20, face="bold"),  axis.title=element_text(size=27,face="bold"), 
        axis.ticks=element_line(size=3)) + theme( legend.position = "none",strip.text = element_text(size = 15) ) +
        scale_x_discrete(name="", breaks=c(1, 3), labels=c("Cabecera", "Curso bajo"))


bd.plot

dev.off()



## assessing effects! no echar cuenta 

anova(mod)
#tab_model(mod)
#Analysis of Variance Table


#Response: Distances
#Df  Sum Sq  Mean Sq F value Pr(>F)
#Groups     1 0.04966 0.049664  2.8072 0.1028
#Residuals 35 0.61921 0.017692               

  
permutest(mod, pairwise = TRUE, permutations = 999)


  
plot(mod)





#tramos
  
  
  
mod.t <- betadisper(var.beta.bichos, datos1$Tramo, type="centroid")


permutest(mod.t, pairwise = TRUE, permutations = 999)


plot(mod.t)



#cuenca
  
  
  
mod.c <- betadisper(var.beta.bichos, datos1$Cuenca, type="centroid")


permutest(mod.c, pairwise = TRUE, permutations = 999)



plot(mod.c)





####### GLOBAL PLOT



library(ggpubr)


tiff(file="FIG 1. NEMATODE.tiff", width=8, height=13, units="in", res=800, compression = "lzw")



figure <-  ggarrange(richness.plot, shannon.plot, bd.plot,  ncol = 1, nrow = 3, labels = c("A", "B", "C"), font.label = list(size = 30, color = "#0072B2"))


figure


dev.off()



pdf(file="FIG 1. NEMA.pdf", width=8, height=13)



figure <-  ggarrange(richness.plot, shannon.plot, bd.plot,  ncol = 1, nrow = 3, labels = c("A", "B", "C"), font.label = list(size = 30, color = "#0072B2"))


figure


dev.off()



pdf(file="ALFA_NEMA.pdf", width=8, height=13)



figure <-  ggarrange(richness.plot, shannon.plot,  ncol = 1, nrow = 2, labels = c("A", "B"), font.label = list(size = 30, color = "#0072B2"))


figure


dev.off()





pdf(file="BETANEMA.pdf", width=8, height=6)



figure <-  ggarrange(bd.plot,  ncol = 1, nrow = 1, font.label = list(size = 30, color = "black"))


figure


dev.off()
