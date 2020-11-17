#================================
#  Statistiques multivariées
#  Script C. Hubas (MNHN)
#  V2.0 - 2020
#  Script du TD
#================================

###########################
# Packages necessaires 
###########################

library(corrplot)
library(PerformanceAnalytics)
library(ade4)
library(factoextra)
library(cowplot)

###########################
# aesthetics 
###########################

my.palette <- colorRampPalette(c("red3","orange","yellow","green3","royalblue"))

My_Theme <- theme(
   axis.text.x = element_text(angle = 45,hjust = 1,size = 12),
   axis.text.y = element_text(size = 12),
   axis.title = element_text(size = 12),
   legend.title=element_text(size = 12),
   legend.text=element_text(size = 12),
   strip.background =element_rect(fill="lightgrey"))

###########################
# ACP étude de cas 
###########################

acp <- read.table("ACP.txt",h=T) # on charge les fichiers ACP.txt
corrplot(cor(acp[,1:13]), method="number") # matrice de corrélations
pairs(acp[,1:13]) # petite pause paour jetter un oeil aux nuages de points
chart.Correlation(acp[,1:13], histogram=TRUE, pch=19) # petite pause 2 paour jetter un oeil aux nuages de points

boxplot(acp[,1:13]) # Est-ce que quelque chose vous choque ?
acp.scale<-scale(acp[,1:13]) # on centre et on réduit

par(mfrow=c(1,2), las=2)
boxplot(acp[,1:13])
boxplot(acp.scale) # Quelle est la principale différence entre les deux ?

res.acp<-dudi.pca(acp[,1:13])
3 # ici essayer d'abord 2 axes 

# règle de Kaiser-Guttman
abline(h=1,col="red")
text(x=10,y=1.2,"Seuil de Kaiser-Guttman",col="red")

# La règle de Karlis – Saporta - Spinaki
nb_col<-dim(acp[,1:13])[2]
nb_li<-dim(acp[,1:13])[1]
seuil.vp <- 1+2*sqrt((nb_col-1)/(nb_li-1))
abline(h=seuil.vp,col="darkred")
text(x=10,y=seuil.vp+0.2,"Seuil de Karlis–Saporta-Spinaki",col="darkred")

res.acp<-dudi.pca(acp[,1:13],scannf = F,nf=3) # arguments pour réaliser l'ACP en fixant le nbr d'axes

res.acp # prenez le temps de décortiquer cet objet
res.acp$eig # valeurs propres
get_eig(res.acp) # même chose avec factoextra

fviz_screeplot(res.acp, addlabels = TRUE, ylim = c(0, 50)) # pourcentgae d'explication porté par chaque axe

res.acp.var <- get_pca_var(res.acp) # obtenir les cos2/corel/contrib des variables
res.acp.ind <- get_pca_ind(res.acp) # obtenir les cos2/contrib des individus
res.acp.var$contrib # exemple : les contributions des variables aux axes

# Contributions des variables à l'axe 1
fviz_contrib(res.acp, choice = "var", axes = 1, top = 10)
# Contributions des variables to l'axe 2
fviz_contrib(res.acp, choice = "var", axes = 2, top = 10)


###########################
# Représenter les variables
###########################

fviz_pca_var(res.acp, col.var = "black") # représenter les variables
# variante
fviz_pca_var(res.acp, col.var="contrib", # colorier les variables en fonction de leur contrib
             gradient.cols = my.palette(3), # controler les couleurs
             repel = TRUE # éviter l'overlapping
)

###########################
# Représenter les individus
###########################

fviz_pca_ind(res.acp, col.ind = "cos2", # colorier les individus en fonction de leur contrib
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # une autre palette de couleur
             repel = TRUE # éviter l'overlapping
)

###########################
# Utiliser des variables supp
###########################

# sur la variable 14 (laissée de côté) le nom de sa station (A, B ou C) est noté
# sur la variable 15 (laissée de côté) le mois d'échantillonnage (1-12, i.e. janvier-décembre) est noté
acp[,14:15]
# petite retouche esthétique de la variable 15
acp$moisacp <- factor(
   c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov")[acp$moisacp],
   level=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"))
# ce sont des variables discrètes
# pour améliorer la lecture de l'ordination, nous pouvons representer ces groupes sur le plan factoriel
# chaque groupe sera repéré par son barycentre

fviz_pca_ind(res.acp,
             label = "none", # cacher les étiquettes
             habillage = acp$sitesacp, # couleur selon le facteur
             palette = c("#00AFBB", "#E7B800", "#FC4E07"), # controler les couleurs
             addEllipses = TRUE # ajouter des ellipses
)

fviz_pca_ind(res.acp,
             label = "none", # cacher les étiquettes
             habillage = acp$sitesacp, # couleur selon le facteur
             palette = c("#00AFBB", "#E7B800", "#FC4E07"), # controler les couleurs
             addEllipses = TRUE, # ajouter des ellipses
             ellipse.type = "convex" # controller le "type" des ellipses
)

# Une grande liberté dans la représentation

fviz_pca_biplot(res.acp,
                habillage = acp$sitesacp,
                label = "var",
                addEllipses = TRUE,
                ellipse.type = "convex"
)

fviz_pca_biplot(res.acp,
                habillage = acp$sitesacp,
                label = "var",
                addEllipses = TRUE,
                ellipse.type = "convex"
) +
theme_bw()

fviz_pca_biplot(res.acp,
                habillage = acp$sitesacp,
                label = "var",
                col.var="white",
                addEllipses = TRUE,
                ellipse.type = "convex"
) +
theme_dark()

###########################
# Compiler les figures
###########################

ind1 <- fviz_pca_ind(res.acp,
                     label = "none", # cacher les étiquettes
                     habillage = acp$sitesacp, # couleur selon le facteur
                     palette = c("#00AFBB", "#E7B800", "#FC4E07"), # controler les couleurs
                     addEllipses = TRUE, # ajouter des ellipses
                     ellipse.type = "convex" # controller le "type" des ellipses
)+
   My_Theme

ind2 <- fviz_pca_ind(res.acp,
                     label = "none", # cacher les étiquettes
                     habillage = as.factor(acp$moisacp), # couleur selon le facteur
                     palette = my.palette(11), # controler les couleurs
                     addEllipses = TRUE, # ajouter des ellipses
                     ellipse.type = "convex" # controller le "type" des ellipses
)+
   My_Theme


var <- fviz_pca_var(res.acp, col.var="contrib", # colorier les variables en fonction de leur contrib
                    gradient.cols = my.palette(3), # controler les couleurs
                    repel = TRUE # éviter l'overlapping
)+
   My_Theme

left <- plot_grid(ind1, ind2,labels=c("a)", "b)"),ncol=1)
plot_grid(left, var,labels=c("", "c)"),ncol=2)

# interpréter ..........................

###########################
# EXERCICE ACP 1
###########################

# Essayer d'explorer la fonction fviz_pca_var() et de jouer avec les arguments
# Essayer de trouver un moyen de représenter la variable NH4 sur le graphique
# quelle est la contribution de cette variable à l'axe 1 ?

# réponse : 
#fviz_pca_ind(res.acp,
#             geom = "point",
#             label="none",
#             axes = c(1, 2),
#             pointsize = acp$NH4
#)


###########################
# EXERCICE ACP 2
###########################

percent <- read.table("pourcentages.txt",h=T)

# Le jeu de données pourcentgae.txt présente les résulats d'une expérience
# des sacs à litière ont été disposés dans une mangrove
# certains sacs avec 100% d'algues (codé A00)
# certains sacs avec 100% de feuilles de palétuviers (codé F00)
# certains sacs avec 50% de chaque (codé M50)
# certains sacs avec 75% de feuilles (codé M75)
# les sacs ont été laissé entre 0 et 35 jours (codé T00 à T35)
# les codes sont dans la varaible "ID"
# Aidez le chercheur à y voir plus clair dans ses données

# Aide : il est possible d'extraire les facteurs de la variable "ID" avec la commande substr()
percent$treatment <- substr(percent$ID,1,3)
percent$time <- substr(percent$ID,4,6)

######################################
# NMDS
######################################

library(vegan)
?metaMDS

afc <- read.table("AFC.txt",h=T)

dist=vegdist(afc[,1:35],"bray")
list(dist)
afc$sitesafc
as.matrix(dist)

matrice=as.matrix(dist)
cols=substr(colnames(matrice),1,2)
rows=substr(rownames(matrice),1,2)
matrice[rows=="SA",cols=="SA"]

as.vector(dist)

   N <- attr(dist, "Size")
   irow <- as.vector(as.dist(row(matrix(nrow = N, ncol = N))))
   icol <- as.vector(as.dist(col(matrix(nrow = N, ncol = N))))
   group.row=afc$sitesafc[irow]
   group.col=afc$sitesafc[icol]
   new.factor=paste(group.row,group.col)
   site <- split(dist, new.factor)
names(site)
boxplot(site)

dist=vegdist(afc[,1:35],"euclidean")
list(dist)
afc$sitesafc

   N <- attr(dist, "Size")
   irow <- as.vector(as.dist(row(matrix(nrow = N, ncol = N))))
   icol <- as.vector(as.dist(col(matrix(nrow = N, ncol = N))))
   group.row=afc$sitesafc[irow]
   group.col=afc$sitesafc[icol]
   new.factor=paste(group.row,group.col)
   site <- split(dist, new.factor)
names(site)
boxplot(site)

?vegdist

######################################
# Tests Mentel
######################################
library(vegan)
data(varespec)
dim(varespec)
data(varechem)
dim(varechem)
veg.dist <- vegdist(varespec) # Bray-Curtis
env.dist <- vegdist(scale(varechem), "euclid")
mantel(veg.dist, env.dist)
mantel(veg.dist, env.dist, method="spear")
attr(veg.dist,"Size")
attr(env.dist,"Size")
mantel.rtest(veg.dist, env.dist, nrepet = 9999)
cor.test(as.vector(veg.dist),as.vector(env.dist))
geo=data.frame(LAT=sample(seq(from=-49.28217,to=-49.28159,by=0.00001),size=24),LONG=sample(seq(from=-16.59406,to=-16.59384,by=0.000001),size=24))
geo.dist <- vegdist(scale(geo), "euclid")
plot(geo)
mantel.partial(veg.dist, env.dist, geo.dist, method="pearson", permutations=999)
?sample


######################################
# Tests Mentel Partiel
######################################
comm<-read.table("/Users/cedrichubas/Documents/CED/2-Enseignement/SEP M2 : ED - multivari?/comm.txt",head=T)
soil<-read.table("/Users/cedrichubas/Documents/CED/2-Enseignement/SEP M2 : ED - multivari?/soil.txt",head=T,row.names=1)
geo<-read.table("/Users/cedrichubas/Documents/CED/2-Enseignement/SEP M2 : ED - multivari?/geo.txt",head=T,row.names=1)

library(vegan)
commdist<-vegdist(comm, method="bray")
soildist<-vegdist(soil, method="bray")
geodist<-vegdist(geo, method="euclidean")

### Now, we are going to use a partial mantel test to test if communities have more different species compositions as the soil composition in the sites get increasingly different. We will do that while removing any possible spatial autocorrelation. So, the third distance matrix in the analysis will be the spatial distance between sites.
mantel(commdist, soildist)
mantel.partial(commdist, soildist, geodist, method="pearson", permutations=999)


Ztheo=NULL
for (i in 1:1000) {
	new.veg.dist=sample(as.vector(veg.dist),replace=F)
	Ztheo[i]=sum(new.veg.dist*env.dist)
}

Zobs=sum(veg.dist*env.dist)
hist(Ztheo,xlim=c(800,870))
abline(v=Zobs)
nbr.sup=length(which(Ztheo>Zobs))
p.val=(1+nbr.sup)/1000

#######################################



res.mds=metaMDS(afc[,1:35], distance="bray", k=2, trymax=500, plot=T) # pour faire le calcul avec un indice de Bray-Curtis
res.mds
summary(res.mds)
res.mds$stress # pour connaitre la valeur du "stress"
stressplot(res.mds) # pour obtenir le diagrame de Shepard

# on repr?sente la matrice de similarit? 
par(mfrow=c(1,2))
plot(res.mds, display = c("site"), type="t", main=c("Stress=",res.mds$stress),cex=2)
s.class(res.mds$point, afc$sitesafc, add.plot=T,col=c("violetred","darkblue","lightseagreen"))
plot(res.mds, display = c("site"), type="t",main=c("Stress=",res.mds$stress))
s.class(res.mds$point, as.factor(afc$moisafc), add.plot=T,col=colors()[100:112])

# on peut repr?senter en 3 dimensions (package scatterplot3d obligatoire)
res.mds2=metaMDS(afc[,1:35], distance="bray", k=3, trymax=50, plot=T) # on refait le calcul avec 3 dimensions
library(scatterplot3d)

par(mfrow=c(1,2))
plot(res.mds, display = c("site"), type="t", main=c("Stress=",res.mds$stress))
s.class(res.mds$point, afc$sitesafc, add.plot=T,col=c("violetred","darkblue","lightseagreen"))
scatterplot3d(res.mds2$points[,1], res.mds2$points[,2],res.mds2$points[,3],highlight.3d=T,type="h", main=c("stress=",res.mds2$stress/100))

# on va faire ressortir les groupes 
c=data.frame(cbind(res.mds2$points,afc$sitesafc))
c.s=split(c,c[,4])
res=scatterplot3d(res.mds2$points[,1], res.mds2$points[,2],res.mds2$points[,3],highlight.3d=T,type="h", main=c("stress=",res.mds2$stress))
res$points3d(c.s$"1"[,1],c.s$"1"[,2],c.s$"1"[,3], type="h",col="violetred", pch=16)
res$points3d(c.s$"2"[,1],c.s$"2"[,2],c.s$"2"[,3], type="h",col="darkblue", pch=16)
res$points3d(c.s$"3"[,1],c.s$"3"[,2],c.s$"3"[,3], type="h",col="lightseagreen", pch=16)

# g?ce ? la fonction metaMDs on peut aussi repr?senter les esp?ces
plot(res.mds, display = c("site","species"), type="t", main=c("Stress=",res.mds$stress))
s.class(res.mds$point, sitesafc, add.plot=T,col=c("violetred","darkblue","lightseagreen"))


# EXERCICE avec AFC01
nmds01=read.table("/Users/cedrichubas/Documents/CED/2-Enseignement/SEP M2 : ED - multivari?/AFC01.txt", h=T)
attach(nmds01)
names(nmds01)
library(vegan)
multi.afc01=metaMDS(nmds01[,0:35],distance="jaccard",k=3, trymax=50, plot=T)
plot(multi.afc01, display = c("site","species"), type="t", main=c("Stress=", multi.afc01$stress))
s.class(res.mds$point, nmds01$sitesafc, add.plot=T,col=c("violetred","darkblue","lightseagreen"))


######################################
# PCoA Bonus !
######################################

multi=vegdist(afc[,1:35],"bray")
attach(afc)
library(vegan)
library(ade4)
multibis=cmdscale(multi,eig=T)
multibis$eig/sum(multibis$eig)*100
summary(multibis)
ordiplot(multibis)
s.class(multibis$points, afc$sitesafc,col=c("violetred","darkblue","lightseagreen"),add.plot=T)

par(mfrow=c(1,2))
