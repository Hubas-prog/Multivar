#================================
#  Statistiques multivariées
#  Script C. Hubas (MNHN)
#  V2.0 - 2020
#  Script du cours
#================================

###########################
# Packages necessaires
###########################

library(MASS)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(shiny)
library(vegan)
library(ggplot2)
library(scatterplot3d)
library(pvclust)
library(ape)
library(vegan)

##########################
# Correlation 
##########################

data<-read.table("COR.txt",h=T) # charger le fichier COR.txt contenant vos données
plot(data$hauteur~data$opercule,
	xlab="Largeur de l'opercule (cm)",
	ylab="hauteur de la coquille (cm)") # pour faire un nuage de points
abline(h=mean(data$hauteur),lty=1,col="red")
abline(v=mean(data$opercule),lty=1,col="red")
cor.test(data$hauteur,data$opercule,method="pearson")

###########################
# Regression de Y en X
############################

h.o<-lm(data$hauteur~data$opercule) # pour calculer la droite de régression

# Dispersion résiduelle
plot(data$hauteur~data$opercule,
	xlab="Largeur de l'opercule (cm)",
	ylab="hauteur de la coquille (cm)",
	col="grey") # pour faire un nuage de points
abline(h.o,col="blue",lwd=2) # pour tracer la droite de régression h.o
yb<-predict(h.o)
segments(x0=data$opercule,x1=data$opercule,y0=data$hauteur,y1=yb,col="blue")

# Dispersion totale
plot(data$hauteur~data$opercule,
	xlab="Largeur de l'opercule (cm)",
	ylab="hauteur de la coquille (cm)",
	col="grey") # pour faire un nuage de points
abline(h=mean(data$hauteur),lty=1,col="purple")
segments(x0=data$opercule,x1=data$opercule,y0=data$hauteur,y1=mean(data$hauteur),col="purple")

# Dispersion expliquée
plot(data$hauteur~data$opercule,
	xlab="Largeur de l'opercule (cm)",
	ylab="hauteur de la coquille (cm)",
	col="grey") # pour faire un nuage de points
abline(h.o,col="blue",lwd=2) # pour tracer la droite de régression h.o
abline(h=mean(data$hauteur),lty=1,col="purple")
segments(x0=data$opercule,x1=data$opercule,y0=mean(data$hauteur),y1=yb,col="orange")

############################
# Regression de X en Y
############################

o.h<-lm(data$opercule~data$hauteur)
summary(o.h)
# Dispersion résiduelle
plot(data$opercule~data$hauteur,
	ylab="Largeur de l'opercule (cm)",
	xlab="hauteur de la coquille (cm)",
	col="grey") # pour faire un nuage de points
abline(o.h,col="green",lwd=2) # pour tracer la droite de régression h.o
ybp<-predict(o.h)
segments(x0=data$hauteur,x1=data$hauteur,y0=data$opercule,y1=ybp,col="green")


############################
# Les deux droites dans le même espace
############################

hauteur.scaled<-scale(data$hauteur)
opercule.scaled<-scale(data$opercule)
plot(hauteur.scaled~opercule.scaled,
	xlab="Largeur (centrée,réduite) de l'opercule (cm)",
	ylab="Hauteur (centrée,réduite) de la coquille (cm)",
	col="black")
abline(h=mean(hauteur.scaled),lty=1,col="red")
abline(v=mean(opercule.scaled),lty=1,col="red")
hsos<-lm(hauteur.scaled~opercule.scaled)
abline(hsos,col="blue",lwd=2)
oshs<-lm(opercule.scaled~hauteur.scaled)
oshs$coef
hauteur.bis=(opercule.scaled+oshs$coef[1])/oshs$coef[2]
points(hauteur.bis~opercule.scaled,type="l",col="green3",lwd=2)


############################
# Lien entre angle alpha & r
############################

# 1° × π/180 = 0,01745 rad

cos(0*pi/180)
round(cos(36.9*pi/180),1)
round(cos(66.4*pi/180),1)
round(cos(90*pi/180),1)
round(cos(113.6*pi/180),1)
round(cos(143.1*pi/180),1)
cos(180*pi/180)

############################
# Jouez avec les correlations
# Avec l'aide précieuse d'Elie Arnaud
# Ingénieur Données et Métadonnées, UMS Patrinat / PNDB
############################

ui <- fluidPage(
  # Side panel to chose n and r settings
  sidebarPanel(
    # n
    # numericInput("n", "Number of samples", max = 10000, min = 1, step = 1),
    sliderInput("n", "Number of samples", max = 10000, min = 1, step = 1, value = 10000),
    # r
    # numericInput("r", "Correlation", 0.23, min = 0, max = 1, step = 0.01)
    sliderInput("r", "Correlation", max = 1, min = -1, step = 0.01, value = 0.23),
  ),
  # Main panel to display plot
  mainPanel(
    plotOutput("plot")
  )
)

server <- function(input, output, session) {
  
  make.plot<-function(n,r){
    data.simul <- MASS::mvrnorm(
      n, 
      mu=c(0, 0), 
      Sigma=matrix(c(1, r, r, 1), nrow=2),
      empirical=TRUE
    )
    x<-data.simul[, 1]
    y<-data.simul[, 2]
    plot(y~x,xlim=c(-4,4),ylim=c(-4,4),col="grey")
    M1<-lm(y~x)
    M2<-lm(x~y)
    abline(M1,col="blue",lwd=3)
    points(
      (x+M2$coef[1])/M2$coef[2]~x,
      type="l",
      col="green3",
      lwd=3
    )
  }
  
  # Fournir le coefficient de corrélation r et le nombre d'échantillons n pour générer le graphique correspondant
  make.plot(n=10000,r=0.23)
  
  output$plot <- renderPlot(
    # make.plot(n=10000,r=0.23)
    make.plot(
      n=input$n,
      r=input$r
    )
  )
}

shinyApp(ui, server)


############################
# les composantes principales
############################

EIG<-eigen(cor(data))

plot(hauteur.scaled~opercule.scaled,
	xlab="Largeur (centrée,réduite) de l'opercule (cm)",
	ylab="Hauteur (centrée,réduite) de la coquille (cm)",
	col="black")
abline(h=mean(hauteur.scaled),lty=1,col="red")
abline(v=mean(opercule.scaled),lty=1,col="red")

arrows(x0=-1*EIG$vectors[1,1],
	y0=-1*EIG$vectors[2,1],
	x1=EIG$vectors[1,1],
	y1=EIG$vectors[2,1],lwd=2)

arrows(x0=-1*EIG$vectors[1,2],
	y0=-1*EIG$vectors[2,2],
	x1=EIG$vectors[1,2],
	y1=EIG$vectors[2,2],lwd=2)

############################
# Généralisation a k variables
############################

acp <- read.table("ACP.txt",,h=T)
corrplot(cor(acp[,1:13]), method="number")
pairs(acp[,1:13])
eigen(cor(acp[,1:13]))

# Vérification
ev <- eigen(cor(acp[,1:13]))$values
M <- cor(acp[,1:13])
P <- eigen(cor(acp[,1:13]))$vectors
round(solve(P)%*%as.matrix(M)%*%P,2) # diagonalistation P-1MP

# valeurs propre
plot(cumsum(ev)/13*100,type="b",ylab="valeurs propres cumulées (%)",xlab="axes factoriels", xaxp=c(1,13,12))
plot(ev,type="b",ylab="valeurs propres",xlab="axes factoriels", xaxp=c(1,13,12))

# règle de Kaiser-Guttman
abline(h=1,col="red")
text(x=10,y=1.2,"Seuil de Kaiser-Guttman",col="red")

# La règle de Karlis – Saporta - Spinaki
seuil.vp <- 1+2*sqrt((dim(acp[,1:13])[2]-1)/(dim(acp[,1:13])[1]-1))
abline(h=seuil.vp,col="darkred")
text(x=10,y=seuil.vp+0.2,"Seuil de Karlis–Saporta-Spinaki",col="darkred")

############################
# Correlations des variables
############################

COR <- t(t(P) * sqrt(ev)/apply(scale(acp[,1:13]),2,sd)) 
rownames(COR)<-colnames(acp[,1:13])
colnames(COR)<-paste("Dim.",1:13); COR

# La corrélation d’une variable Vi d'intérêt avec un axe j est donnée par la formule ci dessus. Avec P la matrice de vecteurs propre, ev les valeurs propres, et l’écart-type des variables (celle qui nous intéresse). A noter que dans le cas d’une ACP normée, l'écart-type vaut 1. A noter l'utilisation de la fonction t() pour tranposer la matrice de vecteurs propres et faciliter le calcul. Source : https://stats.stackexchange.com/questions/253718/correlation-between-an-original-variable-and-a-principal component

############################
# Cos2
############################

COS2 <- round(COR^2,5) ; COS2

############################
# Contribution des variables aux axes
############################

CONT <- P*P*100
rownames(CONT)<-colnames(acp[,1:13])
colnames(CONT)<-paste("Dim.",1:13); CONT
apply(CONT,1,sum)

############################
# Avec FactoMine R
############################

res.pca <- PCA(acp[,1:13], scale.unit=TRUE, ncp=5, graph=F)

# Variables
res.pca$var$cont
res.pca$var$cos2
res.pca$var$cor

# représentation graphique
fviz_pca_var(res.pca, col.var = "black")

# contribution des variables
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)

# identifier les variables > à contribution moyenne sur axe1
WAv1<-which(res.pca$var$cont[,1]>mean(res.pca$var$cont[,1])) ; WAv1
# identifier les variables > à contribution moyenne sur axe2
WAv2<-which(res.pca$var$cont[,2]>mean(res.pca$var$cont[,2])) ; WAv2

# calculer la contribution totale des variables dont la contibution individuelle est > à la cont. moyenne 
sum(res.pca$var$cont[,1][WAv1]) # axe1
sum(res.pca$var$cont[,2][WAv2]) # axe2

# Individus
res.pca$ind$cont
res.pca$ind$cos2

# représentation graphique
fviz_pca_ind(res.pca, col.var = "black")

# contribution des indiv
fviz_contrib(res.pca, choice = "ind", axes = 1, top = 20)
fviz_contrib(res.pca, choice = "ind", axes = 2, top = 20)

# identifier les indiv > à contribution moyenne sur axe1
WAi1<-which(res.pca$ind$cont[,1]>mean(res.pca$ind $cont[,1])) ; WAi1
# identifier les variables > à contribution moyenne sur axe2
WAi2<-which(res.pca$ind $cont[,2]>mean(res.pca$ind $cont[,2])) ; WAi2

# calculer la contribution totale des variables dont la contibution individuelle est > à la cont. moyenne 
sum(res.pca$ind$cont[,1][WAi1]) # axe1
sum(res.pca$ind$cont[,2][WAi2]) # axe2


############################
# Indices de dissimilarité
############################

afc<-read.table("AFC.txt",h=T)
names(afc)
dist.B<-vegdist(afc[,1:35],"bray")
dist.J<-vegdist(ifelse(afc[,1:35]>0,1,0),"jaccard")

# Représentation d'une matrice de dissimilarité
source("https://ichthyology.usm.edu/courses/multivariate/coldiss.R")
par(mfrow=c(1,2))
coldiss(dist.B)
coldiss(dist.J)

multi2 =metaMDS(afc[,1:35],"bray",k=2,trymax=500,plot=T) # pour faire une nMDS avec un indice de Bray-Curtis
summary(multi2)
multi2$stress # pour connaitre la valeur du "stress"
stressplot(multi2) # pour obtenir le diagrame de Shepard

# on représente l'ordination
Stations<-afc$sitesafc
ggplot(as.data.frame(multi2$points),aes(x=MDS1,y=MDS2,col=Stations))+
	geom_point(cex=4)+
	theme_bw()

# on peut représenter en 3 dimensions (package scatterplot3d obligatoire)
res.mds2<-metaMDS(afc[,1:35], distance="bray", k=3, trymax=50, plot=T) # on refait le calcul avec 3 dimensions

# on va faire ressortir les groupes 
c<-data.frame(cbind(res.mds2$points,afc$sitesafc))
c.s<-split(c,c[,4])
res<-scatterplot3d(res.mds2$points[,1],
	res.mds2$points[,2],
	res.mds2$points[,3],
	highlight.3d=T,
	type="h",
	main=c("stress=",res.mds2$stress))
res$points3d(c.s$"1"[,1],c.s$"1"[,2],c.s$"1"[,3],
	type="h",
	col="red",
	pch=16)
res$points3d(c.s$"2"[,1],c.s$"2"[,2],c.s$"2"[,3],
	type="h",
	col="green",
	pch=16)
res$points3d(c.s$"3"[,1],c.s$"3"[,2],c.s$"3"[,3],
	type="h",
	col="blue",
	pch=16)

# Arbre phylo
tree<-hclust(vegdist(afc[,1:35],"bray"),"ward.D2")
plot(tree)
plot(as.phylo(tree),
	type = "unrooted",
	cex =1,
	label.offset = 0,
	no.margin = TRUE,
	tip.color = rainbow(3)[cutree(tree,3)])