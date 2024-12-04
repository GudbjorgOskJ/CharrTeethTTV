###### Correlation between bone shape and teeth traits ##########
##Library packages 
library(geomorph)
library(shapes)
library(svd)
library(scatterplot3d)
library(rgl)
library(MASS)
library(ape)
library(vegan)
library(ggplot2)
library(ggfortify)
library(abind)
library(car)
library(patchwork)
library(emmeans)
library(lsmeans)


###################################################################
####### Correlation of Dentary shape to all teeth traits ##############

#Note this script uses shape data from a separate project, tps files containing shape data and csv files with individual data can be found here: https://github.com/GudbjorgOskJ/CharrBonesTVV 

###Reading data and general shape data preparation

##Reading data
den<-readland.tps("Dentary_LM.tps", negNA = TRUE, specID = c("imageID"))

#### Estimate missing landmarks
den1<-estimate.missing(den,method = "TPS")

#### Flipping of individuals where it is needed.
Need_flipping_den<-den1[,,40:41]
Flipped_den<-rotate.coords(Need_flipping_den, "flipX")
den1[,,40:41]<-Flipped_den

den1<-rotate.coords(den1, type = c("rotateC"))
den1<-rotate.coords(den1, type = c("flipX"))

#### Specifying which landmarks are sliding. LM = 10, 14 and 15 for Dentary
sliding <- matrix(c(9, 10, 11, 13, 14, 15, 14, 15, 1), ncol = 3, byrow = TRUE)

####Performing Generlized Procrustes analysis
gpa.lands<-gpagen(den1, curves = sliding)


### Adding character information
denCh<- read.csv("Dentary_Ch.csv", stringsAsFactors = TRUE, header=TRUE)


gdf.den<-geomorph.data.frame(gpa.lands, morph = denCh$morph, 
                             csize = gpa.lands$Csize,logcsize = log(gpa.lands$Csize),
                             ind = dimnames(gpa.lands$coords)[[3]], length = denCh$length_cm,
                             weight = denCh$weight_g, id = denCh$id_me, sex = denCh$sex, 
                             age = denCh$age, loglength = log(denCh$length_cm), teethd = denCh$teeth_den,
                             teethm =denCh$teeth_max, teethpr = denCh$teeth_pre, teethpa = denCh$teeth_pal,
                             teethv = denCh$teeth_vo, teethg = denCh$teeth_glos)


###############################################
### Correlation Analysis #############


#### Dentary landmark data ########

#### All morphs together
#Vomer
missingv=c(19,30,73,78,80,89,115,133,135,166,167,206,209)
#Glossohyal
missingg=c(17,26,167,178)

denxden<-two.b.pls(A1=gdf.den$coords, A2=gdf.den$teethd, iter =999, seed= NULL, print.progress =TRUE)
denxmax<-two.b.pls(A1=gdf.den$coords[,,-95], A2=gdf.den$teethm[-95], iter = 999, seed = NULL, print.progress = TRUE)
denxpre<-two.b.pls(A1=gdf.den$coords, A2=gdf.den$teethpr, iter = 999, seed = NULL, print.progress = TRUE)
denxpal<-two.b.pls(A1=gdf.den$coords[,,-139], A2=gdf.den$teethpa[-139], iter = 999, seed = NULL, print.progress = TRUE)
denxvom<-two.b.pls(A1=gdf.den$coords[,,-missingv], A2=gdf.den$teethv[-missingv], iter = 999, seed = NULL, print.progress = TRUE)
denxglos<-two.b.pls(A1=gdf.den$coords[,,-missingg], A2=gdf.den$teethg[-missingg], iter = 999, seed = NULL, print.progress = TRUE)

plot(denxden)
#Run for other models

compare.pls(denxden, denxmax, denxpre, denxpal, denxvom, denxglos)

### By morph
###Preparation
sb<- c(1:33, 90:101, 103, 104, 105, 116, 129, 136:167)
sbm<- c(1:33, 90:94, 96:101, 103, 104, 105, 116, 129, 136:167)
sbpal<- c(1:33, 90:101, 103, 104, 105, 116, 129, 136:138, 140:167)
sbv<- c(1:18, 20:29, 31:33, 90:101, 103, 104, 105, 116, 129, 136:165)
sbg<- c(1:16,18:25, 27:33, 90:101, 103, 104, 105, 116, 129, 136:166)
pl<-c(41, 52, 53, 54, 55, 60, 61, 62, 63, 64, 71:89, 102, 117, 118, 130:135, 195:216,237)
plv<-c(41, 52, 53, 54, 55, 60, 61, 62, 63, 64, 71:72, 74:77, 79, 81:88, 102, 117, 118, 130:132, 134, 195:205, 207:208, 210:216,237)
lb<-c(37:40, 43, 45, 46, 47, 51,58,59,67,106,109,114, 168:194, 217, 219, 221, 222, 225, 233, 234)
lbg<-c(37:40, 43, 45, 46, 47, 51, 58, 59, 67, 106, 109, 114, 168:177,179:194, 217, 219, 221, 222, 225, 233, 234)
pi<-c(34:36, 42, 44, 48, 49, 50, 56, 57, 65, 66, 68, 69, 70, 107, 108, 110, 111, 112, 113, 115, 119:128, 218, 220, 223, 224, 226:232, 235, 236, 238, 239)
piv<- c(34:36, 42, 44, 48, 49, 50, 56, 57, 65, 66, 68, 69, 70, 107, 108, 110, 111, 112, 113, 119:128, 218, 220, 223, 224, 226:232, 235, 236, 238, 239)


##denxden
denxdensb<-two.b.pls(A1=gdf.den$coords[,,sb], A2=gdf.den$teethd[sb], iter =999, seed= NULL, print.progress =TRUE)
denxdenpl<-two.b.pls(A1=gdf.den$coords[,,pl], A2=gdf.den$teethd[pl], iter =999, seed= NULL, print.progress =TRUE)
denxdenlb<-two.b.pls(A1=gdf.den$coords[,,lb], A2=gdf.den$teethd[lb], iter =999, seed= NULL, print.progress =TRUE)
denxdenpi<-two.b.pls(A1=gdf.den$coords[,,pi], A2=gdf.den$teethd[pi], iter =999, seed= NULL, print.progress =TRUE)

plot(denxdensb)
#Run for other models

compare.pls(denxdensb, denxdenpl, denxdenlb, denxdenpi)

##denxmax
denxmaxsb<-two.b.pls(A1=gdf.den$coords[,,sbm], A2=gdf.den$teethm[sbm], iter = 999, seed = NULL, print.progress = TRUE)
denxmaxpl<-two.b.pls(A1=gdf.den$coords[,,pl], A2=gdf.den$teethm[pl], iter = 999, seed = NULL, print.progress = TRUE)
denxmaxlb<-two.b.pls(A1=gdf.den$coords[,,lb], A2=gdf.den$teethm[lb], iter = 999, seed = NULL, print.progress = TRUE)
denxmaxpi<-two.b.pls(A1=gdf.den$coords[,,pi], A2=gdf.den$teethm[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(denxmaxsb)
#Run for other models

compare.pls(denxmaxsb, denxmaxpl, denxmaxlb, denxmaxpi)

##denxpre
denxpresb<-two.b.pls(A1=gdf.den$coords[,,sb], A2=gdf.den$teethpr[sb], iter = 999, seed = NULL, print.progress = TRUE)
denxprepl<-two.b.pls(A1=gdf.den$coords[,,pl], A2=gdf.den$teethpr[pl], iter = 999, seed = NULL, print.progress = TRUE)
denxprelb<-two.b.pls(A1=gdf.den$coords[,,lb], A2=gdf.den$teethpr[lb], iter = 999, seed = NULL, print.progress = TRUE)
denxprepi<-two.b.pls(A1=gdf.den$coords[,,pi], A2=gdf.den$teethpr[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(denxpresb)
#Run for other models

compare.pls(denxpresb, denxprepl, denxprelb, denxprepi)

##denxpal
denxpalsb<-two.b.pls(A1=gdf.den$coords[,,sbpal], A2=gdf.den$teethpa[sbpal], iter = 999, seed = NULL, print.progress = TRUE)
denxpalpl<-two.b.pls(A1=gdf.den$coords[,,pl], A2=gdf.den$teethpa[pl], iter = 999, seed = NULL, print.progress = TRUE)
denxpallb<-two.b.pls(A1=gdf.den$coords[,,lb], A2=gdf.den$teethpa[lb], iter = 999, seed = NULL, print.progress = TRUE)
denxpalpi<-two.b.pls(A1=gdf.den$coords[,,pi], A2=gdf.den$teethpa[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(denxpalsb)
#Run for other models

compare.pls(denxpalsb, denxpalpl, denxpallb, denxpalpi)

##denxvom
denxvomsb<-two.b.pls(A1=gdf.den$coords[,,sbv], A2=gdf.den$teethv[sbv], iter = 999, seed = NULL, print.progress = TRUE)
denxvompl<-two.b.pls(A1=gdf.den$coords[,,plv], A2=gdf.den$teethv[plv], iter = 999, seed = NULL, print.progress = TRUE)
denxvomlb<-two.b.pls(A1=gdf.den$coords[,,lb], A2=gdf.den$teethv[lb], iter = 999, seed = NULL, print.progress = TRUE)
denxvompi<-two.b.pls(A1=gdf.den$coords[,,piv], A2=gdf.den$teethv[piv], iter = 999, seed = NULL, print.progress = TRUE)

plot(denxvomsb)
#Run for other models

compare.pls(denxvomsb, denxvompl, denxvomlb, denxvompi)

##denxglos
denxglossb<-two.b.pls(A1=gdf.den$coords[,,sbg], A2=gdf.den$teethg[sbg], iter = 999, seed = NULL, print.progress = TRUE)
denxglospl<-two.b.pls(A1=gdf.den$coords[,,pl], A2=gdf.den$teethg[pl], iter = 999, seed = NULL, print.progress = TRUE)
denxgloslb<-two.b.pls(A1=gdf.den$coords[,,lbg], A2=gdf.den$teethg[lbg], iter = 999, seed = NULL, print.progress = TRUE)
denxglospi<-two.b.pls(A1=gdf.den$coords[,,pi], A2=gdf.den$teethg[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(denxglossb)
#Run for other models

compare.pls(denxglossb, denxglospl, denxgloslb, denxglospi)

#### Dentary dental palate length #####
#### All morphs together
lmks<-matrix(c(2,4,2,5,1,2,4,14), ncol = 2, byrow = T, dimnames = list(c("dental", "lingual","MS", "height"), c("start", "end")))
lineardists <- interlmkdist(gdf.den$coords, lmks)

DPdenxden<-two.b.pls(A1=log(lineardists[,1]), A2=gdf.den$teethd, iter =999, seed= NULL, print.progress =TRUE)
DPdenxmax<-two.b.pls(A1=log(lineardists[-95,1]), A2=gdf.den$teethm[-95], iter = 999, seed = NULL, print.progress = TRUE)
DPdenxpre<-two.b.pls(A1=log(lineardists[,1]), A2=gdf.den$teethpr, iter = 999, seed = NULL, print.progress = TRUE)
DPdenxpal<-two.b.pls(A1=log(lineardists[-139,1]), A2=gdf.den$teethpa[-139], iter = 999, seed = NULL, print.progress = TRUE)
DPdenxvom<-two.b.pls(A1=log(lineardists[-missingv,1]), A2=gdf.den$teethv[-missingv], iter = 999, seed = NULL, print.progress = TRUE)
DPdenxglos<-two.b.pls(A1=log(lineardists[-missingg,1]), A2=gdf.den$teethg[-missingg], iter = 999, seed = NULL, print.progress = TRUE)

plot(DPdenxden)
#Run for other models

compare.pls(DPdenxden, DPdenxmax, DPdenxpre,DPdenxpal,DPdenxvom,DPdenxglos)

## By morph

##DPdenxden
DPdenxdensb<-two.b.pls(A1=log(lineardists[sb,1]), A2=gdf.den$teethd[sb], iter =999, seed= NULL, print.progress =TRUE)
DPdenxdenpl<-two.b.pls(A1=log(lineardists[pl,1]), A2=gdf.den$teethd[pl], iter =999, seed= NULL, print.progress =TRUE)
DPdenxdenlb<-two.b.pls(A1=log(lineardists[lb,1]), A2=gdf.den$teethd[lb], iter =999, seed= NULL, print.progress =TRUE)
DPdenxdenpi<-two.b.pls(A1=log(lineardists[pi,1]), A2=gdf.den$teethd[pi], iter =999, seed= NULL, print.progress =TRUE)

plot(DPdenxdensb)
#Run for other models

compare.pls(DPdenxdensb, DPdenxdenpl, DPdenxdenlb, DPdenxdenpi)

##DPdenxmax
DPdenxmaxsb<-two.b.pls(A1=log(lineardists[sbm,1]), A2=gdf.den$teethm[sbm], iter = 999, seed = NULL, print.progress = TRUE)
DPdenxmaxpl<-two.b.pls(A1=log(lineardists[pl,1]), A2=gdf.den$teethm[pl], iter = 999, seed = NULL, print.progress = TRUE)
DPdenxmaxlb<-two.b.pls(A1=log(lineardists[lb,1]), A2=gdf.den$teethm[lb], iter = 999, seed = NULL, print.progress = TRUE)
DPdenxmaxpi<-two.b.pls(A1=log(lineardists[pi,1]), A2=gdf.den$teethm[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(DPdenxmaxsb)
#Run for other models

compare.pls(DPdenxmaxsb, DPdenxmaxpl, DPdenxmaxlb, DPdenxmaxpi)

##DPdenxpre
DPdenxpresb<-two.b.pls(A1=log(lineardists[sb,1]), A2=gdf.den$teethpr[sb], iter = 999, seed = NULL, print.progress = TRUE)
DPdenxprepl<-two.b.pls(A1=log(lineardists[pl,1]), A2=gdf.den$teethpr[pl], iter = 999, seed = NULL, print.progress = TRUE)
DPdenxprelb<-two.b.pls(A1=log(lineardists[lb,1]), A2=gdf.den$teethpr[lb], iter = 999, seed = NULL, print.progress = TRUE)
DPdenxprepi<-two.b.pls(A1=log(lineardists[pi,1]), A2=gdf.den$teethpr[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(DPdenxpresb)
#Run for other models

compare.pls(DPdenxpresb, DPdenxprepl, DPdenxprelb, DPdenxprepi)

##DPdenxpal
DPdenxpalsb<-two.b.pls(A1=log(lineardists[sbpal,1]), A2=gdf.den$teethpa[sbpal], iter = 999, seed = NULL, print.progress = TRUE)
DPdenxpalpl<-two.b.pls(A1=log(lineardists[pl,1]), A2=gdf.den$teethpa[pl], iter = 999, seed = NULL, print.progress = TRUE)
DPdenxpallb<-two.b.pls(A1=log(lineardists[lb,1]), A2=gdf.den$teethpa[lb], iter = 999, seed = NULL, print.progress = TRUE)
DPdenxpalpi<-two.b.pls(A1=log(lineardists[pi,1]), A2=gdf.den$teethpa[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(DPdenxpalsb)
#Run for other models

compare.pls(DPdenxpalsb, DPdenxpalpl, DPdenxpallb, DPdenxpalpi)

##DPdenxvom
DPdenxvomsb<-two.b.pls(A1=log(lineardists[sbv,1]), A2=gdf.den$teethv[sbv], iter = 999, seed = NULL, print.progress = TRUE)
DPdenxvompl<-two.b.pls(A1=log(lineardists[plv,1]), A2=gdf.den$teethv[plv], iter = 999, seed = NULL, print.progress = TRUE)
DPdenxvomlb<-two.b.pls(A1=log(lineardists[lb,1]), A2=gdf.den$teethv[lb], iter = 999, seed = NULL, print.progress = TRUE)
DPdenxvompi<-two.b.pls(A1=log(lineardists[piv,1]), A2=gdf.den$teethv[piv], iter = 999, seed = NULL, print.progress = TRUE)

plot(DPdenxvomsb)
#Run for other models

compare.pls(DPdenxvomsb, DPdenxvompl, DPdenxvomlb, DPdenxvompi)

##DPdenxglos
DPdenxglossb<-two.b.pls(A1=log(lineardists[sbg,1]), A2=gdf.den$teethg[sbg], iter = 999, seed = NULL, print.progress = TRUE)
DPdenxglospl<-two.b.pls(A1=log(lineardists[pl,1]), A2=gdf.den$teethg[pl], iter = 999, seed = NULL, print.progress = TRUE)
DPdenxgloslb<-two.b.pls(A1=log(lineardists[lbg,1]), A2=gdf.den$teethg[lbg], iter = 999, seed = NULL, print.progress = TRUE)
DPdenxglospi<-two.b.pls(A1=log(lineardists[pi,1]), A2=gdf.den$teethg[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(DPdenxglossb)
#Run for other models

compare.pls(DPdenxglossb, DPdenxglospl, DPdenxgloslb, DPdenxglospi)

### Dentary landmark data vs. dental palate length vs. fork length vs. bone C-size  ####
lenxden<-two.b.pls(A1=gdf.den$loglength, A2=gdf.den$teethd, iter =999, seed= NULL, print.progress =TRUE)
Csizexden<-two.b.pls(A1=gdf.den$logcsize, A2=gdf.den$teethd, iter =999, seed= NULL, print.progress =TRUE)

compare.pls(denxden, DPdenxden, lenxden, Csizexden)

####################################################################
#####################################################################
############### Maxilla morphology ################################

###Reading data and general data preparation

max<-readland.tps("Maxilla_LM.tps", negNA = TRUE, specID = c("imageID"))

max1<-estimate.missing(max,method = "TPS")

Need_flipping_max<-max1[,,131]
Flipped_max<-rotate.coords(Need_flipping_max, "flipY")
max1[,,131]<-Flipped_max

max1<-rotate.coords(max1, type = c("flipY"))

##### LM: 15, 16, 17
sliding1 <-matrix(c(14, 15, 16, 15, 16, 7, 7, 17, 5),ncol = 3,byrow = TRUE) # Landmarks 15, 16, 17
max.land1<-gpagen(max1, curves = sliding1)

maxCh<- read.csv("Maxilla_Ch.csv", stringsAsFactors = TRUE, header=TRUE)

gdf.max<-geomorph.data.frame(max.land1, morph = maxCh$morph, 
                             csize = max.land1$Csize,logcsize = log(max.land1$Csize),
                             ind = dimnames(max.land1$coords)[[3]], length = maxCh$length_cm,
                             weight = maxCh$weight_g, id = maxCh$id_me, sex = maxCh$sex, 
                             age = maxCh$age, loglength = log(maxCh$length_cm), 
                             teethd = maxCh$teeth_den, teethm =maxCh$teeth_max, 
                             teethpr = maxCh$teeth_pre, teethpa = maxCh$teeth_pal, 
                             teethv = maxCh$teeth_vo, teethg = maxCh$teeth_glos, 
                             teethang = maxCh$teeth_ang)

###############################################
### Correlation Analysis #############

### Tooth count ####

### Maxilla landmark data ####
#### All morphs together
#Vomer
missingv=c(19,30,74,79,81,90,115,133,135,166,167,206,209)
#Glossohyal
missingg=c(17,26,167,178)

maxxden<-two.b.pls(A1=gdf.max$coords[,,-49], A2=gdf.max$teethd[-49], iter =999, seed= NULL, print.progress =TRUE)
maxxmax<-two.b.pls(A1=gdf.max$coords, A2=gdf.max$teethm, iter = 999, seed = NULL, print.progress = TRUE)
maxxpre<-two.b.pls(A1=gdf.max$coords, A2=gdf.max$teethpr, iter = 999, seed = NULL, print.progress = TRUE)
maxxpal<-two.b.pls(A1=gdf.max$coords[,,-139], A2=gdf.max$teethpa[-139], iter = 999, seed = NULL, print.progress = TRUE)
maxxvom<-two.b.pls(A1=gdf.max$coords[,,-missingv], A2=gdf.max$teethv[-missingv], iter = 999, seed = NULL, print.progress = TRUE)
maxxglos<-two.b.pls(A1=gdf.max$coords[,,-missingg], A2=gdf.max$teethg[-missingg], iter = 999, seed = NULL, print.progress = TRUE)

plot(maxxden)
#Run for other models

compare.pls(maxxden,maxxmax,maxxpre,maxxpal, maxxvom,maxxglos)

## By morph
###Preparation
sb<- c(1:33, 91:101, 103, 104, 105, 116, 129, 136:167)
sbpal<- c(1:33, 91:101, 103, 104, 105, 116, 129, 136:138, 140:167)
sbv<- c(1:18, 20:29, 31:33, 91:101, 103, 104, 105, 116, 129, 136:165)
sbg<- c(1:16, 18:25, 27:33, 91:101, 103, 104, 105, 116, 129, 136:166)

pl<-c(41, 53, 54, 55,56, 61, 62, 63, 64, 65, 72:90, 102, 117, 118, 130:135, 195:216,237)
plv<- c(41, 53, 54, 55,56, 61, 62, 63, 64, 65, 72:73, 75:78, 80, 82:89, 102, 117, 118, 130:132, 134, 195:205, 207:208, 210:216,237)

lb<-c(37:40, 43, 45, 46, 47, 52, 59, 60, 68, 106, 109,114, 168:194, 217, 219, 221, 222, 225, 233, 234)
lbg<-c(37:40, 43, 45, 46, 47, 52, 59, 60, 68, 106, 109,114, 168:177, 179:194, 217, 219, 221, 222, 225, 233, 234)

pi<-c(34:36, 42, 44, 48, 49, 50, 51, 57, 58, 66, 67, 69, 70, 71, 107, 108, 110, 111, 112, 113, 115, 119:128, 218, 220, 223, 224, 226:232, 235, 236, 238, 239)
pid<-c(34:36, 42, 44, 48, 50, 51, 57, 58, 66, 67, 69, 70, 71, 107, 108, 110, 111, 112, 113, 115, 119:128, 218, 220, 223, 224, 226:232, 235, 236, 238, 239)
piv<-c(34:36, 42, 44, 48, 50, 51, 57, 58, 66, 67, 69, 70, 71, 107, 108, 110, 111, 112, 113, 119:128, 218, 220, 223, 224, 226:232, 235, 236, 238, 239)

##maxxden
maxxdensb<-two.b.pls(A1=gdf.max$coords[,,sb], A2=gdf.max$teethd[sb], iter =999, seed= NULL, print.progress =TRUE)
maxxdenpl<-two.b.pls(A1=gdf.max$coords[,,pl], A2=gdf.max$teethd[pl], iter =999, seed= NULL, print.progress =TRUE)
maxxdenlb<-two.b.pls(A1=gdf.max$coords[,,lb], A2=gdf.max$teethd[lb], iter =999, seed= NULL, print.progress =TRUE)
maxxdenpi<-two.b.pls(A1=gdf.max$coords[,,pid], A2=gdf.max$teethd[pid], iter =999, seed= NULL, print.progress =TRUE)

plot(maxxdensb)
#Run for other models

compare.pls(maxxdensb, maxxdenpl, maxxdenlb, maxxdenpi)

##maxxmax
maxxmaxsb<-two.b.pls(A1=gdf.max$coords[,,sb], A2=gdf.max$teethm[sb], iter = 999, seed = NULL, print.progress = TRUE)
maxxmaxpl<-two.b.pls(A1=gdf.max$coords[,,pl], A2=gdf.max$teethm[pl], iter = 999, seed = NULL, print.progress = TRUE)
maxxmaxlb<-two.b.pls(A1=gdf.max$coords[,,lb], A2=gdf.max$teethm[lb], iter = 999, seed = NULL, print.progress = TRUE)
maxxmaxpi<-two.b.pls(A1=gdf.max$coords[,,pi], A2=gdf.max$teethm[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(maxxmaxsb)
#Run for other models

compare.pls(maxxmaxsb, maxxmaxpl, maxxmaxlb, maxxmaxpi)

##maxxpre
maxxpresb<-two.b.pls(A1=gdf.max$coords[,,sb], A2=gdf.max$teethpr[sb], iter = 999, seed = NULL, print.progress = TRUE)
maxxprepl<-two.b.pls(A1=gdf.max$coords[,,pl], A2=gdf.max$teethpr[pl], iter = 999, seed = NULL, print.progress = TRUE)
maxxprelb<-two.b.pls(A1=gdf.max$coords[,,lb], A2=gdf.max$teethpr[lb], iter = 999, seed = NULL, print.progress = TRUE)
maxxprepi<-two.b.pls(A1=gdf.max$coords[,,pi], A2=gdf.max$teethpr[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(maxxpresb)
#Run for other models

compare.pls(maxxpresb, maxxprepl, maxxprelb, maxxprepi)

##maxxpal
maxxpalsb<-two.b.pls(A1=gdf.max$coords[,,sbpal], A2=gdf.max$teethpa[sbpal], iter = 999, seed = NULL, print.progress = TRUE)
maxxpalpl<-two.b.pls(A1=gdf.max$coords[,,pl], A2=gdf.max$teethpa[pl], iter = 999, seed = NULL, print.progress = TRUE)
maxxpallb<-two.b.pls(A1=gdf.max$coords[,,lb], A2=gdf.max$teethpa[lb], iter = 999, seed = NULL, print.progress = TRUE)
maxxpalpi<-two.b.pls(A1=gdf.max$coords[,,pi], A2=gdf.max$teethpa[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(maxxpalsb)
#Run for other models

compare.pls(maxxpalsb, maxxpalpl, maxxpallb, maxxpalpi)

##maxxvom
maxxvomsb<-two.b.pls(A1=gdf.max$coords[,,sbv], A2=gdf.max$teethv[sbv], iter = 999, seed = NULL, print.progress = TRUE)
maxxvompl<-two.b.pls(A1=gdf.max$coords[,,plv], A2=gdf.max$teethv[plv], iter = 999, seed = NULL, print.progress = TRUE)
maxxvomlb<-two.b.pls(A1=gdf.max$coords[,,lb], A2=gdf.max$teethv[lb], iter = 999, seed = NULL, print.progress = TRUE)
maxxvompi<-two.b.pls(A1=gdf.max$coords[,,piv], A2=gdf.max$teethv[piv], iter = 999, seed = NULL, print.progress = TRUE)

plot(maxxvomsb)
#Run for other models

compare.pls(maxxvomsb, maxxvompl, maxxvomlb, maxxvompi)

## maxxglos
maxxglossb<-two.b.pls(A1=gdf.max$coords[,,sbg], A2=gdf.max$teethg[sbg], iter = 999, seed = NULL, print.progress = TRUE)
maxxglospl<-two.b.pls(A1=gdf.max$coords[,,pl], A2=gdf.max$teethg[pl], iter = 999, seed = NULL, print.progress = TRUE)
maxxgloslb<-two.b.pls(A1=gdf.max$coords[,,lbg], A2=gdf.max$teethg[lbg], iter = 999, seed = NULL, print.progress = TRUE)
maxxglospi<-two.b.pls(A1=gdf.max$coords[,,pi], A2=gdf.max$teethg[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(maxxglossb)
#Run for other models

compare.pls(maxxglossb, maxxglospl, maxxgloslb, maxxglospi)

#### Maxilla dental palate length #####
#### All morphs together
lmks<-matrix(c(4,8,4,9,6,7), ncol = 2, byrow = T, dimnames = list(c("dental1", "dental2","width"), c("start", "end")))
lineardists <- interlmkdist(gdf.max$coords, lmks)

DPmaxxden<-two.b.pls(A1=log(lineardists[-49,1]), A2=gdf.max$teethd[-49], iter =999, seed= NULL, print.progress =TRUE)
DPmaxxmax<-two.b.pls(A1=log(lineardists[,1]), A2=gdf.max$teethm, iter = 999, seed = NULL, print.progress = TRUE)
DPmaxxpre<-two.b.pls(A1=log(lineardists[,1]), A2=gdf.max$teethpr, iter = 999, seed = NULL, print.progress = TRUE)
DPmaxxpal<-two.b.pls(A1=log(lineardists[-139,1]), A2=gdf.max$teethpa[-139], iter = 999, seed = NULL, print.progress = TRUE)
DPmaxxvom<-two.b.pls(A1=log(lineardists[-missingv,1]), A2=gdf.max$teethv[-missingv], iter = 999, seed = NULL, print.progress = TRUE)
DPmaxxglos<-two.b.pls(A1=log(lineardists[-missingg,1]), A2=gdf.max$teethg[-missingg], iter = 999, seed = NULL, print.progress = TRUE)


plot(DPmaxxden)
#Run for other models

compare.pls(DPmaxxden, DPmaxxmax, DPmaxxpre,DPmaxxpal,DPmaxxvom,DPmaxxglos)


## By morph
##maxxden
DPmaxxdensb<-two.b.pls(A1=log(lineardists[sb,1]), A2=gdf.max$teethd[sb], iter =999, seed= NULL, print.progress =TRUE)
DPmaxxdenpl<-two.b.pls(A1=log(lineardists[pl,1]), A2=gdf.max$teethd[pl], iter =999, seed= NULL, print.progress =TRUE)
DPmaxxdenlb<-two.b.pls(A1=log(lineardists[lb,1]), A2=gdf.max$teethd[lb], iter =999, seed= NULL, print.progress =TRUE)
DPmaxxdenpi<-two.b.pls(A1=log(lineardists[pid,1]), A2=gdf.max$teethd[pid], iter =999, seed= NULL, print.progress =TRUE)

plot(DPmaxxdensb)
#Run for other models

compare.pls(DPmaxxdensb, DPmaxxdenpl, DPmaxxdenlb, DPmaxxdenpi)

##maxxmax
DPmaxxmaxsb<-two.b.pls(A1=log(lineardists[sb,1]), A2=gdf.max$teethm[sb], iter = 999, seed = NULL, print.progress = TRUE)
DPmaxxmaxpl<-two.b.pls(A1=log(lineardists[pl,1]), A2=gdf.max$teethm[pl], iter = 999, seed = NULL, print.progress = TRUE)
DPmaxxmaxlb<-two.b.pls(A1=log(lineardists[lb,1]), A2=gdf.max$teethm[lb], iter = 999, seed = NULL, print.progress = TRUE)
DPmaxxmaxpi<-two.b.pls(A1=log(lineardists[pi,1]), A2=gdf.max$teethm[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(DPmaxxmaxsb)
#Run for other models

compare.pls(DPmaxxmaxsb, DPmaxxmaxpl, DPmaxxmaxlb, DPmaxxmaxpi)

##maxxpre
DPmaxxpresb<-two.b.pls(A1=log(lineardists[sb,1]), A2=gdf.max$teethpr[sb], iter = 999, seed = NULL, print.progress = TRUE)
DPmaxxprepl<-two.b.pls(A1=log(lineardists[pl,1]), A2=gdf.max$teethpr[pl], iter = 999, seed = NULL, print.progress = TRUE)
DPmaxxprelb<-two.b.pls(A1=log(lineardists[lb,1]), A2=gdf.max$teethpr[lb], iter = 999, seed = NULL, print.progress = TRUE)
DPmaxxprepi<-two.b.pls(A1=log(lineardists[pi,1]), A2=gdf.max$teethpr[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(DPmaxxpresb)
#Run for other models

compare.pls(DPmaxxpresb, DPmaxxprepl, DPmaxxprelb, DPmaxxprepi)

##maxxpal
DPmaxxpalsb<-two.b.pls(A1=log(lineardists[sbpal,1]), A2=gdf.max$teethpa[sbpal], iter = 999, seed = NULL, print.progress = TRUE)
DPmaxxpalpl<-two.b.pls(A1=log(lineardists[pl,1]), A2=gdf.max$teethpa[pl], iter = 999, seed = NULL, print.progress = TRUE)
DPmaxxpallb<-two.b.pls(A1=log(lineardists[lb,1]), A2=gdf.max$teethpa[lb], iter = 999, seed = NULL, print.progress = TRUE)
DPmaxxpalpi<-two.b.pls(A1=log(lineardists[pi,1]), A2=gdf.max$teethpa[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(DPmaxxpalsb)
#Run for other models

compare.pls(DPmaxxpalsb, DPmaxxpalpl, DPmaxxpallb, DPmaxxpalpi)

##maxxvom
DPmaxxvomsb<-two.b.pls(A1=log(lineardists[sbv,1]), A2=gdf.max$teethv[sbv], iter = 999, seed = NULL, print.progress = TRUE)
DPmaxxvompl<-two.b.pls(A1=log(lineardists[plv,1]), A2=gdf.max$teethv[plv], iter = 999, seed = NULL, print.progress = TRUE)
DPmaxxvomlb<-two.b.pls(A1=log(lineardists[lb,1]), A2=gdf.max$teethv[lb], iter = 999, seed = NULL, print.progress = TRUE)
DPmaxxvompi<-two.b.pls(A1=log(lineardists[piv,1]), A2=gdf.max$teethv[piv], iter = 999, seed = NULL, print.progress = TRUE)

plot(DPmaxxvomsb)
#Run for other models

compare.pls(DPmaxxvomsb, DPmaxxvompl, DPmaxxvomlb, DPmaxxvompi)

##maxxglos
DPmaxxglossb<-two.b.pls(A1=log(lineardists[sbg,1]), A2=gdf.max$teethg[sbg], iter = 999, seed = NULL, print.progress = TRUE)
DPmaxxglospl<-two.b.pls(A1=log(lineardists[pl,1]), A2=gdf.max$teethg[pl], iter = 999, seed = NULL, print.progress = TRUE)
DPmaxxgloslb<-two.b.pls(A1=log(lineardists[lbg,1]), A2=gdf.max$teethg[lbg], iter = 999, seed = NULL, print.progress = TRUE)
DPmaxxglospi<-two.b.pls(A1=log(lineardists[pi,1]), A2=gdf.max$teethg[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(DPmaxxglossb)
#Run for other models

compare.pls(DPmaxxglossb, DPmaxxglospl, DPmaxxgloslb, DPmaxxglospi)

### Maxilla landmark data vs. dental palate length vs. fork length vs. bone C-size  ####
lenxmax<-two.b.pls(A1=gdf.max$loglength, A2=gdf.max$teethm, iter = 999, seed = NULL, print.progress = TRUE)
Csizexmax<-two.b.pls(A1=gdf.max$logcsize, A2=gdf.max$teethm, iter =999, seed= NULL, print.progress =TRUE)

compare.pls(maxxmax,DPmaxxmax,lenxmax,Csizexmax)

### Tooth angle ####

missingang=c(9,18,26,28,30,31,34,35,41,68,71,79,81,87,102,104,105,106,114, 117,118,129,130,131,154,163,208,216,217)
maxxang<-two.b.pls(A1=gdf.max$coords[,,-missingang], A2=gdf.max$teethang[-missingang], iter = 999, seed = NULL, print.progress = TRUE)

DPmaxxang<-two.b.pls(A1=log(lineardists[-x,missingang]), A2=gdf.max$teethang[-missingang], iter = 999, seed = NULL, print.progress = TRUE)

lenxang<-two.b.pls(A1=gdf.max$loglength[-missingang], A2=gdf.max$teethang[-missingang], iter = 999, seed = NULL, print.progress = TRUE)
Csizexang<-two.b.pls(A1=gdf.max$logcsize[-missingang], A2=gdf.max$teethang[-missingang], iter =999, seed= NULL, print.progress =TRUE)

compare.pls(maxxang,DPmaxxang,lenxang,Csizexang)

####################################################################
#####################################################################
############### Premaxilla morphology ################################

###Reading data and general data preparation

prem<-readland.tps("Premaxilla_LM.tps", negNA = TRUE, specID = c("imageID"))

prem1<-estimate.missing(prem,method = "TPS")


x<-c(1,9,28,64,79,82,87,97,99,198,206,216,219)
Need_flipping_prem<-prem1[,,x]
Flipped_prem<-rotate.coords(Need_flipping_prem, "flipY")
prem1[,,x]<-Flipped_prem


prem1<-rotate.coords(prem1, type = c("rotateC"))

#### No sliding LM
pre.land<-gpagen(prem1)


preCh<- read.csv("Premaxilla_Ch.csv", stringsAsFactors = TRUE, header=TRUE) # Yes same data set as Articular angular


gdf.pre<-geomorph.data.frame(pre.land, morph = preCh$morph, csize = pre.land$Csize,
                             logcsize = log(pre.land$Csize),ind = dimnames(pre.land$coords)[[3]],
                             length = preCh$length_cm,weight = preCh$weight_g, id = preCh$id_me,
                             sex = preCh$sex, age = preCh$age, loglength = log(preCh$length_cm),
                             teethd = preCh$teeth_den, teethm =preCh$teeth_max,
                             teethpr = preCh$teeth_pre, teethpa = preCh$teeth_pal,
                             teethv = preCh$teeth_vo, teethg = preCh$teeth_glos)

###############################################
### Correlation Analysis #############

#### Premaxilla landmark data ########

#### All morphs together
#Vomer
missingv=c(19,30,74,79,81,90,116,134,136,167,168,207,210)
#Glossohyal
missingg=c(17,26,168,179)

prexden<-two.b.pls(A1=gdf.pre$coords[,,-49], A2=gdf.pre$teethd[-49], iter =999, seed= NULL, print.progress =TRUE)
prexmax<-two.b.pls(A1=gdf.pre$coords[,,-96], A2=gdf.pre$teethm[-96], iter = 999, seed = NULL, print.progress = TRUE)
prexpre<-two.b.pls(A1=gdf.pre$coords, A2=gdf.pre$teethpr, iter = 999, seed = NULL, print.progress = TRUE)
prexpal<-two.b.pls(A1=gdf.pre$coords[,,-140], A2=gdf.pre$teethpa[-140], iter = 999, seed = NULL, print.progress = TRUE)
prexvom<-two.b.pls(A1=gdf.pre$coords[,,-missingv], A2=gdf.pre$teethv[-missingv], iter = 999, seed = NULL, print.progress = TRUE)
prexglos<-two.b.pls(A1=gdf.pre$coords[,,-missingg], A2=gdf.pre$teethg[-missingg], iter = 999, seed = NULL, print.progress = TRUE)

plot(prexden)
#Run for other models

compare.pls(prexden,prexmax,prexpre,prexpal,prexvom,prexglos)

#### By morph
###Preparation
sb<- c(1:33, 91:102, 104, 105, 106, 117, 130, 137:168)
sbm<- c(1:33, 91:95, 97:102, 104, 105, 106, 117, 130, 137:168)
sbpal<- c(1:33, 91:102, 104, 105, 106, 117, 130, 137:139, 141:168)
sbv<- c(1:18, 20:29, 31:33, 91:102, 104, 105, 106, 117, 130, 137:166)
sbg<- c(1:16, 18:25, 27:33, 91:102, 104, 105, 106, 117, 130, 137:167)

pl<-c(41, 53, 54, 55, 56, 61, 62, 63, 64, 65, 72:90, 103, 118, 119, 131:136, 196:217,238)
plv<-c(41, 53, 54, 55, 56, 61, 62, 63, 64, 65, 72:73, 75:78, 80, 82:89, 103, 118, 119, 131:133, 135, 196:206, 208:209, 211:217,238)

lb<-c(37:40, 43, 45, 46, 47, 52,59,60,68,107,110,115, 169:195, 218, 220, 222, 223, 226, 234, 235)
lbg<-c(37:40, 43, 45, 46, 47, 52,59,60,68,107,110,115, 169:178, 180:195, 218, 220, 222, 223, 226, 234, 235)

pi<-c(34:36, 42, 44, 48, 49, 50, 51, 57, 58, 66, 67, 69, 70, 71, 108, 109, 111, 112, 113, 114, 116, 120:129, 219, 221, 224, 225, 227:233, 236, 237, 239, 240)
pid<-c(34:36, 42, 44, 48, 50, 51, 57, 58, 66, 67, 69, 70, 71, 108, 109, 111, 112, 113, 114, 116, 120:129, 219, 221, 224, 225, 227:233, 236, 237, 239, 240)
piv<-c(34:36, 42, 44, 48, 49, 50, 51, 57, 58, 66, 67, 69, 70, 71, 108, 109, 111, 112, 113, 114, 120:129, 219, 221, 224, 225, 227:233, 236, 237, 239, 240)


##prexden
prexdensb<-two.b.pls(A1=gdf.pre$coords[,,sb], A2=gdf.pre$teethd[sb], iter =999, seed= NULL, print.progress =TRUE)
prexdenpl<-two.b.pls(A1=gdf.pre$coords[,,pl], A2=gdf.pre$teethd[pl], iter =999, seed= NULL, print.progress =TRUE)
prexdenlb<-two.b.pls(A1=gdf.pre$coords[,,lb], A2=gdf.pre$teethd[lb], iter =999, seed= NULL, print.progress =TRUE)
prexdenpi<-two.b.pls(A1=gdf.pre$coords[,,pid], A2=gdf.pre$teethd[pid], iter =999, seed= NULL, print.progress =TRUE)

plot(prexdensb)
#Run for other models

compare.pls(prexdensb, prexdenpl, prexdenlb, prexdenpi)

##prexmax
prexmaxsb<-two.b.pls(A1=gdf.pre$coords[,,sbm], A2=gdf.pre$teethm[sbm], iter = 999, seed = NULL, print.progress = TRUE)
prexmaxpl<-two.b.pls(A1=gdf.pre$coords[,,pl], A2=gdf.pre$teethm[pl], iter = 999, seed = NULL, print.progress = TRUE)
prexmaxlb<-two.b.pls(A1=gdf.pre$coords[,,lb], A2=gdf.pre$teethm[lb], iter = 999, seed = NULL, print.progress = TRUE)
prexmaxpi<-two.b.pls(A1=gdf.pre$coords[,,pi], A2=gdf.pre$teethm[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(prexmaxsb)
#Run for other models

compare.pls(prexmaxsb, prexmaxpl, prexmaxlb, prexmaxpi)

##prexpre
prexpresb<-two.b.pls(A1=gdf.pre$coords[,,sb], A2=gdf.pre$teethpr[sb], iter = 999, seed = NULL, print.progress = TRUE)
prexprepl<-two.b.pls(A1=gdf.pre$coords[,,pl], A2=gdf.pre$teethpr[pl], iter = 999, seed = NULL, print.progress = TRUE)
prexprelb<-two.b.pls(A1=gdf.pre$coords[,,lb], A2=gdf.pre$teethpr[lb], iter = 999, seed = NULL, print.progress = TRUE)
prexprepi<-two.b.pls(A1=gdf.pre$coords[,,pi], A2=gdf.pre$teethpr[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(prexpresb)
#Run for other models

compare.pls(prexpresb, prexprepl, prexprelb, prexprepi)

##prexpal
prexpalsb<-two.b.pls(A1=gdf.pre$coords[,,sbpal], A2=gdf.pre$teethpa[sbpal], iter = 999, seed = NULL, print.progress = TRUE)
prexpalpl<-two.b.pls(A1=gdf.pre$coords[,,pl], A2=gdf.pre$teethpa[pl], iter = 999, seed = NULL, print.progress = TRUE)
prexpallb<-two.b.pls(A1=gdf.pre$coords[,,lb], A2=gdf.pre$teethpa[lb], iter = 999, seed = NULL, print.progress = TRUE)
prexpalpi<-two.b.pls(A1=gdf.pre$coords[,,pi], A2=gdf.pre$teethpa[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(prexpalsb)
#Run for other models

compare.pls(prexpalsb, prexpalpl, prexpallb, prexpalpi)

##prexom
prexvomsb<-two.b.pls(A1=gdf.pre$coords[,,sbv], A2=gdf.pre$teethv[sbv], iter = 999, seed = NULL, print.progress = TRUE)
prexvompl<-two.b.pls(A1=gdf.pre$coords[,,plv], A2=gdf.pre$teethv[plv], iter = 999, seed = NULL, print.progress = TRUE)
prexvomlb<-two.b.pls(A1=gdf.pre$coords[,,lb], A2=gdf.pre$teethv[lb], iter = 999, seed = NULL, print.progress = TRUE)
prexvompi<-two.b.pls(A1=gdf.pre$coords[,,piv], A2=gdf.pre$teethv[piv], iter = 999, seed = NULL, print.progress = TRUE)

plot(prexvomsb)
#Run for other models

compare.pls(prexvomsb, prexvompl, prexvomlb, prexvompi)

##prexglos
prexglossb<-two.b.pls(A1=gdf.pre$coords[,,sbg], A2=gdf.pre$teethg[sbg], iter = 999, seed = NULL, print.progress = TRUE)
prexglospl<-two.b.pls(A1=gdf.pre$coords[,,pl], A2=gdf.pre$teethg[pl], iter = 999, seed = NULL, print.progress = TRUE)
prexgloslb<-two.b.pls(A1=gdf.pre$coords[,,lbg], A2=gdf.pre$teethg[lbg], iter = 999, seed = NULL, print.progress = TRUE)
prexglospi<-two.b.pls(A1=gdf.pre$coords[,,pi], A2=gdf.pre$teethg[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(prexglossb)
#Run for other models

compare.pls(prexglossb, prexglospl, prexgloslb, prexglospi)

#### Premaxilla dental palate length #####
#### All morphs together
lmks<-matrix(c(1,10,8,3), ncol = 2, byrow = T, dimnames = list(c("dental","width"), c("start", "end")))
lineardists <- interlmkdist(gdf.pre$coords, lmks)

DPprexden<-two.b.pls(A1=log(lineardists[-49,1]), A2=gdf.pre$teethd[-49], iter =999, seed= NULL, print.progress =TRUE)
DPprexmax<-two.b.pls(A1=log(lineardists[-96,1]), A2=gdf.pre$teethm[-96], iter = 999, seed = NULL, print.progress = TRUE)
DPprexpre<-two.b.pls(A1=log(lineardists[,1]), A2=gdf.pre$teethpr, iter = 999, seed = NULL, print.progress = TRUE)
DPprexpal<-two.b.pls(A1=log(lineardists[-140,1]), A2=gdf.pre$teethpa[-140], iter = 999, seed = NULL, print.progress = TRUE)
DPprexvom<-two.b.pls(A1=log(lineardists[-missingv,1]), A2=gdf.pre$teethv[-missingv], iter = 999, seed = NULL, print.progress = TRUE)
DPprexglos<-two.b.pls(A1=log(lineardists[-missingg,1]), A2=gdf.pre$teethg[-missingg], iter = 999, seed = NULL, print.progress = TRUE)

plot(DPprexden)
#Run for other models

compare.pls(DPprexden, DPprexmax, DPprexpre, DPprexpal,DPprexvom,DPprexglos)

#### By morph

##prexden
DPprexdensb<-two.b.pls(A1=log(lineardists[sb,1]), A2=gdf.pre$teethd[sb], iter =999, seed= NULL, print.progress =TRUE)
DPprexdenpl<-two.b.pls(A1=log(lineardists[pl,1]), A2=gdf.pre$teethd[pl], iter =999, seed= NULL, print.progress =TRUE)
DPprexdenlb<-two.b.pls(A1=log(lineardists[lb,1]), A2=gdf.pre$teethd[lb], iter =999, seed= NULL, print.progress =TRUE)
DPprexdenpi<-two.b.pls(A1=log(lineardists[pid,1]), A2=gdf.pre$teethd[pid], iter =999, seed= NULL, print.progress =TRUE)

plot(DPprexdensb)
#Run for other models

compare.pls(DPprexdensb, DPprexdenpl, DPprexdenlb, DPprexdenpi)

##prexmax
DPprexmaxsb<-two.b.pls(A1=log(lineardists[sbm,1]), A2=gdf.pre$teethm[sbm], iter = 999, seed = NULL, print.progress = TRUE)
DPprexmaxpl<-two.b.pls(A1=log(lineardists[pl,1]), A2=gdf.pre$teethm[pl], iter = 999, seed = NULL, print.progress = TRUE)
DPprexmaxlb<-two.b.pls(A1=log(lineardists[lb,1]), A2=gdf.pre$teethm[lb], iter = 999, seed = NULL, print.progress = TRUE)
DPprexmaxpi<-two.b.pls(A1=log(lineardists[pi,1]), A2=gdf.pre$teethm[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(DPprexmaxsb)
#Run for other models

compare.pls(DPprexmaxsb, DPprexmaxpl, DPprexmaxlb, DPprexmaxpi)

##prexpre
DPprexpresb<-two.b.pls(A1=log(lineardists[sb,1]), A2=gdf.pre$teethpr[sb], iter = 999, seed = NULL, print.progress = TRUE)
DPprexprepl<-two.b.pls(A1=log(lineardists[pl,1]), A2=gdf.pre$teethpr[pl], iter = 999, seed = NULL, print.progress = TRUE)
DPprexprelb<-two.b.pls(A1=log(lineardists[lb,1]), A2=gdf.pre$teethpr[lb], iter = 999, seed = NULL, print.progress = TRUE)
DPprexprepi<-two.b.pls(A1=log(lineardists[pi,1]), A2=gdf.pre$teethpr[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(DPprexpresb)
#Run for other models

compare.pls(DPprexpresb, DPprexprepl, DPprexprelb, DPprexprepi)

##prexpal
DPprexpalsb<-two.b.pls(A1=log(lineardists[sbpal,1]), A2=gdf.pre$teethpa[sbpal], iter = 999, seed = NULL, print.progress = TRUE)
DPprexpalpl<-two.b.pls(A1=log(lineardists[pl,1]), A2=gdf.pre$teethpa[pl], iter = 999, seed = NULL, print.progress = TRUE)
DPprexpallb<-two.b.pls(A1=log(lineardists[lb,1]), A2=gdf.pre$teethpa[lb], iter = 999, seed = NULL, print.progress = TRUE)
DPprexpalpi<-two.b.pls(A1=log(lineardists[pi,1]), A2=gdf.pre$teethpa[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(DPprexpalsb)
#Run for other models

compare.pls(DPprexpalsb, DPprexpalpl, DPprexpallb, DPprexpalpi)

##prexvom
DPprexvomsb<-two.b.pls(A1=log(lineardists[sbv,1]), A2=gdf.pre$teethv[sbv], iter = 999, seed = NULL, print.progress = TRUE)
DPprexvompl<-two.b.pls(A1=log(lineardists[plv,1]), A2=gdf.pre$teethv[plv], iter = 999, seed = NULL, print.progress = TRUE)
DPprexvomlb<-two.b.pls(A1=log(lineardists[lb,1]), A2=gdf.pre$teethv[lb], iter = 999, seed = NULL, print.progress = TRUE)
DPprexvompi<-two.b.pls(A1=log(lineardists[piv,1]), A2=gdf.pre$teethv[piv], iter = 999, seed = NULL, print.progress = TRUE)

plot(DPprexvomsb)
#Run for other models

compare.pls(DPprexvomsb, DPprexvompl, DPprexvomlb, DPprexvompi)

#prexglos
DPprexglossb<-two.b.pls(A1=log(lineardists[sbg,1]), A2=gdf.pre$teethg[sbg], iter = 999, seed = NULL, print.progress = TRUE)
DPprexglospl<-two.b.pls(A1=log(lineardists[pl,1]), A2=gdf.pre$teethg[pl], iter = 999, seed = NULL, print.progress = TRUE)
DPprexgloslb<-two.b.pls(A1=log(lineardists[lbg,1]), A2=gdf.pre$teethg[lbg], iter = 999, seed = NULL, print.progress = TRUE)
DPprexglospi<-two.b.pls(A1=log(lineardists[pi,1]), A2=gdf.pre$teethg[pi], iter = 999, seed = NULL, print.progress = TRUE)

plot(DPprexglossb)
#Run for other models

compare.pls(DPprexglossb, DPprexglospl, DPprexgloslb, DPprexglospi)

### Premaxilla landmark data vs. dental palate length vs. fork length vs. bone C-size  ####
lenxpre<-two.b.pls(A1=gdf.pre$loglength, A2=gdf.pre$teethpr, iter = 999, seed = NULL, print.progress = TRUE)
Csizexpre<-two.b.pls(A1=gdf.pre$logcsize, A2=gdf.pre$teethpr, iter =999, seed= NULL, print.progress =TRUE)

compare.pls(prexpre,DPprexpre,lenxpre,Csizexpre)

##################################################
#### Comparing between bones #### 
compare.pls(denxden,maxxden,prexden)
compare.pls(denxmax, maxxmax, prexmax)
compare.pls(denxpre,maxxpre,prexpre)

compare.pls(DPdenxden,DPmaxxden,DPprexden)
compare.pls(DPdenxmax, DPmaxxmax, DPprexmax)
compare.pls(DPdenxpre,DPmaxxpre,DPprexpre)
