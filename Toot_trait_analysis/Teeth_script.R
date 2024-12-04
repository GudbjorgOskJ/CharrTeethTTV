
##Library packages 
library(ggplot2)
library(readxl)
library(tidyverse)
library(patchwork)
library(AER)
library(MASS)
library(ggcorrplot)
library(GGally)
library(Hmisc)
library(Rmisc)
library(knitr)
library(moments)
library(lme4)
library(glmmTMB)
library(emmeans)
library(MCMCglmm)
library(dplyr)
library("broom")


###################################################################
################# Teeth count analysis ############################

###Reading data and general data preparation

#Data in wide format
datW<-read.csv("Teeth_countW.csv", na.strings = "NA",sep = ",", stringsAsFactors = TRUE)

#Data in long format
datL<-read.csv("Teeth_countL.csv", na.strings = "NA", sep = ",",stringsAsFactors = TRUE)

head(datW)
head(datL)
dim(datW)
dim(datL)
summary(datL)

#More on datL later on for now, wide format is used.

#############################################################
### Repeatability of March and May data count #########

### Overall correlation
#### dentary
#Mean
cor.test(datW$march_dentary_mean, datW$may_dentary_mean, method=c("pearson"))
#Right
cor.test(datW$march_dentary_rt, datW$may_dentary_rt, method=c("pearson"))
#Left
cor.test(datW$march_dentary_lt, datW$may_dentary_lt, method=c("pearson"))

##### By morph
#LB
#Mean
cor.test(datW[datW$morph == "LB",]$march_dentary_mean, datW[datW$morph == "LB",]$may_dentary_mean, method=c("pearson"))
#Right
cor.test(datW[datW$morph == "LB",]$march_dentary_rt, datW[datW$morph == "LB",]$may_dentary_rt, method=c("pearson"))
#Left
cor.test(datW[datW$morph == "LB",]$march_dentary_lt, datW[datW$morph == "LB",]$may_dentary_lt, method=c("pearson"))
#PI
#Mean
cor.test(datW[datW$morph == "PI",]$march_dentary_mean, datW[datW$morph == "PI",]$may_dentary_mean, method=c("pearson"))
#Right
cor.test(datW[datW$morph == "PI",]$march_dentary_rt, datW[datW$morph == "PI",]$may_dentary_rt, method=c("pearson"))
#Left
cor.test(datW[datW$morph == "PI",]$march_dentary_lt, datW[datW$morph == "PI",]$may_dentary_lt, method=c("pearson"))
#PL
#Mean
cor.test(datW[datW$morph == "PL",]$march_dentary_mean, datW[datW$morph == "PL",]$may_dentary_mean, method=c("pearson"))
#Right
cor.test(datW[datW$morph == "PL",]$march_dentary_rt, datW[datW$morph == "PL",]$may_dentary_rt, method=c("pearson"))
#Left
cor.test(datW[datW$morph == "PL",]$march_dentary_lt, datW[datW$morph == "PL",]$may_dentary_lt, method=c("pearson"))
#SB
#Mean
cor.test(datW[datW$morph == "SB",]$march_dentary_mean, datW[datW$morph == "SB",]$may_dentary_mean, method=c("pearson"))
#Right
cor.test(datW[datW$morph == "SB",]$march_dentary_rt, datW[datW$morph == "SB",]$may_dentary_rt, method=c("pearson"))
#Left
cor.test(datW[datW$morph == "SB",]$march_dentary_lt, datW[datW$morph == "SB",]$may_dentary_lt, method=c("pearson"))


#### maxilla
#Mean
cor.test(datW$march_maxilla_mean, datW$may_maxilla_mean, method=c("pearson"))
#Right
cor.test(datW$march_maxilla_rt, datW$may_maxilla_rt, method=c("pearson"))
#Left
cor.test(datW$march_maxilla_lt, datW$may_maxilla_lt, method=c("pearson"))
##### By morph
#LB
#Mean
cor.test(datW[datW$morph == "LB",]$march_maxilla_mean, datW[datW$morph == "LB",]$may_maxilla_mean, method=c("pearson"))
#Right
cor.test(datW[datW$morph == "LB",]$march_maxilla_rt, datW[datW$morph == "LB",]$may_maxilla_rt, method=c("pearson"))
#Left
cor.test(datW[datW$morph == "LB",]$march_maxilla_lt, datW[datW$morph == "LB",]$may_maxilla_lt, method=c("pearson"))
#PI
#Mean
cor.test(datW[datW$morph == "PI",]$march_maxilla_mean, datW[datW$morph == "PI",]$may_maxilla_mean, method=c("pearson"))
#Right
cor.test(datW[datW$morph == "PI",]$march_maxilla_rt, datW[datW$morph == "PI",]$may_maxilla_rt, method=c("pearson"))
#Left
cor.test(datW[datW$morph == "PI",]$march_maxilla_lt, datW[datW$morph == "PI",]$may_maxilla_lt, method=c("pearson"))
#PL
#Mean
cor.test(datW[datW$morph == "PL",]$march_maxilla_mean, datW[datW$morph == "PL",]$may_maxilla_mean, method=c("pearson"))
#Right
cor.test(datW[datW$morph == "PL",]$march_maxilla_rt, datW[datW$morph == "PL",]$may_maxilla_rt, method=c("pearson"))
#Left
cor.test(datW[datW$morph == "PL",]$march_maxilla_lt, datW[datW$morph == "PL",]$may_maxilla_lt, method=c("pearson"))
#SB
#Mean
cor.test(datW[datW$morph == "SB",]$march_maxilla_mean, datW[datW$morph == "SB",]$may_maxilla_mean, method=c("pearson"))
#Right
cor.test(datW[datW$morph == "SB",]$march_maxilla_rt, datW[datW$morph == "SB",]$may_maxilla_rt, method=c("pearson"))
#Left
cor.test(datW[datW$morph == "SB",]$march_maxilla_lt, datW[datW$morph == "SB",]$may_maxilla_lt, method=c("pearson"))

#### premaxilla
#Mean
cor.test(datW$march_premaxilla_mean, datW$may_premaxilla_mean, method=c("pearson"))
#Right
cor.test(datW$march_premaxilla_rt, datW$may_premaxilla_rt, method=c("pearson"))
#Left
cor.test(datW$march_premaxilla_lt, datW$may_premaxilla_lt, method=c("pearson"))

##### By morph
#LB
#Mean
cor.test(datW[datW$morph == "LB",]$march_premaxilla_mean, datW[datW$morph == "LB",]$may_premaxilla_mean, method=c("pearson"))
#Right
cor.test(datW[datW$morph == "LB",]$march_premaxilla_rt, datW[datW$morph == "LB",]$may_premaxilla_rt, method=c("pearson"))
#Left
cor.test(datW[datW$morph == "LB",]$march_premaxilla_lt, datW[datW$morph == "LB",]$may_premaxilla_lt, method=c("pearson"))
#PI
#Mean
cor.test(datW[datW$morph == "PI",]$march_premaxilla_mean, datW[datW$morph == "PI",]$may_premaxilla_mean, method=c("pearson"))
#Right
cor.test(datW[datW$morph == "PI",]$march_premaxilla_rt, datW[datW$morph == "PI",]$may_premaxilla_rt, method=c("pearson"))
#Left
cor.test(datW[datW$morph == "PI",]$march_premaxilla_lt, datW[datW$morph == "PI",]$may_premaxilla_lt, method=c("pearson"))
#PL
#Mean
cor.test(datW[datW$morph == "PL",]$march_premaxilla_mean, datW[datW$morph == "PL",]$may_premaxilla_mean, method=c("pearson"))
#Right
cor.test(datW[datW$morph == "PL",]$march_premaxilla_rt, datW[datW$morph == "PL",]$may_premaxilla_rt, method=c("pearson"))
#Left
cor.test(datW[datW$morph == "PL",]$march_premaxilla_lt, datW[datW$morph == "PL",]$may_premaxilla_lt, method=c("pearson"))
#SB
#Mean
cor.test(datW[datW$morph == "SB",]$march_premaxilla_mean, datW[datW$morph == "SB",]$may_premaxilla_mean, method=c("pearson"))
#Right
cor.test(datW[datW$morph == "SB",]$march_premaxilla_rt, datW[datW$morph == "SB",]$may_premaxilla_rt, method=c("pearson"))
#Left
cor.test(datW[datW$morph == "SB",]$march_premaxilla_lt, datW[datW$morph == "SB",]$may_premaxilla_lt, method=c("pearson"))

#### palatine
#Mean
cor.test(datW$march_palatine_mean, datW$may_palatine_mean, method=c("pearson"))
#Right
cor.test(datW$march_palatine_rt, datW$may_palatine_rt, method=c("pearson"))
#Left
cor.test(datW$march_palatine_lt, datW$may_palatine_lt, method=c("pearson"))
##### By morph
#LB
#Mean
cor.test(datW[datW$morph == "LB",]$march_palatine_mean, datW[datW$morph == "LB",]$may_palatine_mean, method=c("pearson"))
#Right
cor.test(datW[datW$morph == "LB",]$march_palatine_rt, datW[datW$morph == "LB",]$may_palatine_rt, method=c("pearson"))
#Left
cor.test(datW[datW$morph == "LB",]$march_palatine_lt, datW[datW$morph == "LB",]$may_palatine_lt, method=c("pearson"))
#PI
#Mean
cor.test(datW[datW$morph == "PI",]$march_palatine_mean, datW[datW$morph == "PI",]$may_palatine_mean, method=c("pearson"))
#Right
cor.test(datW[datW$morph == "PI",]$march_palatine_rt, datW[datW$morph == "PI",]$may_palatine_rt, method=c("pearson"))
#Left
cor.test(datW[datW$morph == "PI",]$march_palatine_lt, datW[datW$morph == "PI",]$may_palatine_lt, method=c("pearson"))
#PL
#Mean
cor.test(datW[datW$morph == "PL",]$march_palatine_mean, datW[datW$morph == "PL",]$may_palatine_mean, method=c("pearson"))
#Right
cor.test(datW[datW$morph == "PL",]$march_palatine_rt, datW[datW$morph == "PL",]$may_palatine_rt, method=c("pearson"))
#Left
cor.test(datW[datW$morph == "PL",]$march_palatine_lt, datW[datW$morph == "PL",]$may_palatine_lt, method=c("pearson"))
#SB
#Mean
cor.test(datW[datW$morph == "SB",]$march_palatine_mean, datW[datW$morph == "SB",]$may_palatine_mean, method=c("pearson"))
#Right
cor.test(datW[datW$morph == "SB",]$march_palatine_rt, datW[datW$morph == "SB",]$may_palatine_rt, method=c("pearson"))
#Left
cor.test(datW[datW$morph == "SB",]$march_palatine_lt, datW[datW$morph == "SB",]$may_palatine_lt, method=c("pearson"))


#### vomer
cor.test(datW$march_vomer, datW$may_vomer, method=c("pearson"))
##### By morph
#LB
cor.test(datW[datW$morph == "LB",]$march_vomer, datW[datW$morph == "LB",]$may_vomer, method=c("pearson"))
#PI
cor.test(datW[datW$morph == "PI",]$march_vomer, datW[datW$morph == "PI",]$may_vomer, method=c("pearson"))
#PL
cor.test(datW[datW$morph == "PL",]$march_vomer, datW[datW$morph == "PL",]$may_vomer, method=c("pearson"))
#SB
cor.test(datW[datW$morph == "SB",]$march_vomer, datW[datW$morph == "SB",]$may_vomer, method=c("pearson"))


#### glossohyal
cor.test(datW$march_glossohyal, datW$may_glossohyal, method=c("pearson"))
##### By morph
#LB
cor.test(datW[datW$morph == "LB",]$march_glossohyal, datW[datW$morph == "LB",]$may_glossohyal, method=c("pearson"))
#PI
cor.test(datW[datW$morph == "PI",]$march_glossohyal, datW[datW$morph == "PI",]$may_glossohyal, method=c("pearson"))
#PL
cor.test(datW[datW$morph == "PL",]$march_glossohyal, datW[datW$morph == "PL",]$may_glossohyal, method=c("pearson"))
#SB
cor.test(datW[datW$morph == "SB",]$march_glossohyal, datW[datW$morph == "SB",]$may_glossohyal, method=c("pearson"))


## Linear regression
#### dentary
#Mean
tidy(lm(may_dentary_mean ~ march_dentary_mean, data= datW))
glance(lm(may_dentary_mean ~ march_dentary_mean, data= datW))
#Right
tidy(lm(may_dentary_rt ~ march_dentary_rt, data= datW))
glance(lm(may_dentary_rt ~ march_dentary_rt, data= datW))
#Left
tidy(lm(may_dentary_lt ~ march_dentary_lt, data= datW))
glance(lm(may_dentary_lt ~ march_dentary_lt, data= datW))

##### By morph
by_morph <- group_by(datW, morph)
#Mean
do(by_morph,tidy(lm(may_dentary_mean ~ march_dentary_mean, data = .)))
do(by_morph,glance(lm(may_dentary_mean ~ march_dentary_mean, data = .)))
#Right
do(by_morph,tidy(lm(may_dentary_rt ~ march_dentary_rt, data = .)))
do(by_morph,glance(lm(may_dentary_rt ~ march_dentary_rt, data = .)))
#Left
do(by_morph,tidy(lm(may_dentary_lt ~ march_dentary_lt, data = .)))
do(by_morph,glance(lm(may_dentary_lt ~ march_dentary_lt, data = .)))


#####maxilla
#Mean
tidy(lm(may_maxilla_mean ~ march_maxilla_mean, data= datW))
glance(lm(may_maxilla_mean ~ march_maxilla_mean, data= datW))
#Right
tidy(lm(may_maxilla_rt ~ march_maxilla_rt, data= datW))
glance(lm(may_maxilla_rt ~ march_maxilla_rt, data= datW))
#Left
tidy(lm(may_maxilla_lt ~ march_maxilla_lt, data= datW))
glance(lm(may_maxilla_lt ~ march_maxilla_lt, data= datW))

##### By morph
by_morph <- group_by(datW, morph)
#Mean
do(by_morph,tidy(lm(may_maxilla_mean ~ march_maxilla_mean, data = .)))
do(by_morph,glance(lm(may_maxilla_mean ~ march_maxilla_mean, data = .)))
#Right
do(by_morph,tidy(lm(may_maxilla_rt ~ march_maxilla_rt, data = .)))
do(by_morph,glance(lm(may_maxilla_rt ~ march_maxilla_rt, data = .)))
#Left
do(by_morph,tidy(lm(may_maxilla_lt ~ march_maxilla_lt, data = .)))
do(by_morph,glance(lm(may_maxilla_lt ~ march_maxilla_lt, data = .)))


####premaxilla
#Mean
tidy(lm(may_premaxilla_mean ~ march_premaxilla_mean, data= datW))
glance(lm(may_premaxilla_mean ~ march_premaxilla_mean, data= datW))
#Right
tidy(lm(may_premaxilla_rt ~ march_premaxilla_rt, data= datW))
glance(lm(may_premaxilla_rt ~ march_premaxilla_rt, data= datW))
#Left
tidy(lm(may_premaxilla_lt ~ march_premaxilla_lt, data= datW))
glance(lm(may_premaxilla_lt ~ march_premaxilla_lt, data= datW))

##### By morph
by_morph <- group_by(datW, morph)
#Mean
do(by_morph,tidy(lm(may_premaxilla_mean ~ march_premaxilla_mean, data = .)))
do(by_morph,glance(lm(may_premaxilla_mean ~ march_premaxilla_mean, data = .)))
#Right
do(by_morph,tidy(lm(may_premaxilla_rt ~ march_premaxilla_rt, data = .)))
do(by_morph,glance(lm(may_premaxilla_rt ~ march_premaxilla_rt, data = .)))
#Left
do(by_morph,tidy(lm(may_premaxilla_lt ~ march_premaxilla_lt, data = .)))
do(by_morph,glance(lm(may_premaxilla_lt ~ march_premaxilla_lt, data = .)))


###### palatine
#Mean
tidy(lm(may_palatine_mean ~ march_palatine_mean, data= datW))
glance(lm(may_palatine_mean ~ march_palatine_mean, data= datW))
#Right
tidy(lm(may_palatine_rt ~ march_palatine_rt, data= datW))
glance(lm(may_palatine_rt ~ march_palatine_rt, data= datW))
#Left
tidy(lm(may_palatine_lt ~ march_palatine_lt, data= datW))
glance(lm(may_palatine_lt ~ march_palatine_lt, data= datW))

##### By morph
by_morph <- group_by(datW, morph)
#Mean
do(by_morph,tidy(lm(may_palatine_mean ~ march_palatine_mean, data = .)))
do(by_morph,glance(lm(may_palatine_mean ~ march_palatine_mean, data = .)))
#Right
do(by_morph,tidy(lm(may_palatine_rt ~ march_palatine_rt, data = .)))
do(by_morph,glance(lm(may_palatine_rt ~ march_palatine_rt, data = .)))
#Left
do(by_morph,tidy(lm(may_palatine_lt ~ march_palatine_lt, data = .)))
do(by_morph,glance(lm(may_palatine_lt ~ march_palatine_lt, data = .)))

#### vomer
tidy(lm(may_vomer ~ march_vomer, data= datW))
glance(lm(may_vomer ~ march_vomer, data= datW))
##### By morph
by_morph <- group_by(datW, morph)

do(by_morph,tidy(lm(may_vomer ~ march_vomer, data = .)))
do(by_morph,glance(lm(may_vomer ~ march_vomer, data = .)))

#####glossohyal
tidy(lm(may_glossohyal ~ march_glossohyal, data= datW))
glance(lm(may_glossohyal ~ march_glossohyal, data= datW))
##### By morph
by_morph <- group_by(datW, morph)

do(by_morph,tidy(lm(may_glossohyal ~ march_glossohyal, data = .)))
do(by_morph,glance(lm(may_glossohyal ~ march_glossohyal, data = .)))

#Image preparation
max2<-ggplot(data= datW, aes(x=march_maxilla_mean, y=may_maxilla_mean, color=morph)) + geom_point() + geom_smooth(method=lm, se=T, aes(fill=morph)) +labs(x="maxilla mean (March)", y="maxilla mean (May)", col ="Morph", fill="Morph", title="maxilla" ) + guides(fill="none", col="none") + scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))+ scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")) + theme_bw()
den2<-ggplot(data= datW, aes(x=march_dentary_mean, y=may_dentary_mean, color=morph)) + geom_point() + geom_smooth(method=lm, se=T, aes(fill=morph)) +labs(x="dentary mean (March)", y="dentary mean (May)",col="Morph", fill="Morph", title="dentary" ) + guides(fill="none", col= "none") + scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))+ scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")) + theme_bw()
pre2<-ggplot(data= datW, aes(x=march_premaxilla_mean, y=may_premaxilla_mean, color=morph)) + geom_point() + geom_smooth(method=lm, se=T, aes(fill=morph)) +labs(x="premaxilla mean (March)", y="premaxilla mean (May)", col="Morph",fill="Morph", title="premaxilla") + guides(fill = "none", col= "none") + scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))+ scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")) + theme_bw()
pal2<-ggplot(data= datW, aes(x=march_palatine_mean, y=may_palatine_mean, color=morph)) + geom_point() + geom_smooth(method=lm, se=T, aes(fill=morph)) +labs(x="palatine mean (March)", y="palatine mean (May)", col="Morph",fill="Morph", title="palatine" ) + guides(fill="none", col= "none") + scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))+ scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")) + theme_bw()
vom2<-ggplot(data= datW, aes(x=march_vomer, y=may_vomer, color=morph)) + geom_point() + geom_smooth(method=lm, se=T, aes(fill=morph)) +labs(x="vomer (March)", y="vomer (May)", col="Morph", fill="Morph",title="vomer" ) + guides(fill="none", col="none") + scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))+ scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")) + theme_bw()
glos2<-ggplot(data= datW, aes(x=march_glossohyal, y=may_glossohyal, color=morph)) + geom_point() + geom_smooth(method=lm, se=T, aes(fill=morph)) +labs(x="glossohyal (March)", y="glossohyal (May)", col="Morph", fill="Morph", title="glossohyal") + guides(fill="none", col="none") + scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))+ scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")) + theme_bw()

pdf("MM_newcol.pdf")
multiplot(den2, max2, pre2, pal2, vom2,glos2, layout = matrix(c(1,2,3,4,5,6) ,nrow = 3,byrow = T))
dev.off()

###############################################################
############# Deviation, left and right side ########## 

#All morphs together
tidy(lm(dentary_lt ~ dentary_rt, data= datW))
glance(lm(dentary_lt ~ dentary_rt, data= datW))

tidy(lm(maxilla_lt ~ maxilla_rt, data= datW))
glance(lm(maxilla_lt ~ maxilla_rt, data= datW))

tidy(lm(premaxilla_lt ~ premaxilla_rt, data= datW))
glance(lm(premaxilla_lt ~ premaxilla_rt, data= datW))

tidy(lm(palatine_lt ~ palatine_rt, data= datW))
glance(lm(palatine_lt ~ palatine_rt, data= datW))

#By morph
by_morph <- group_by(datW, morph)
do(by_morph,tidy(lm(dentary_lt ~ dentary_rt, data= .)))
do(by_morph,glance(lm(dentary_lt ~ dentary_rt, data= .)))

do(by_morph,tidy(lm(maxilla_lt ~ maxilla_rt, data= .)))
do(by_morph,glance(lm(maxilla_lt ~ maxilla_rt, data= .)))

do(by_morph,tidy(lm(premaxilla_lt ~ premaxilla_rt, data= .)))
do(by_morph,glance(lm(premaxilla_lt ~ premaxilla_rt, data= .)))

do(by_morph,tidy(lm(palatine_lt ~ palatine_rt, data= .)))
do(by_morph,glance(lm(palatine_lt ~ palatine_rt, data= .)))

#Image preparation

ggplot(data= datW, aes(x=dentary_lt, y=dentary_rt)) + geom_point() + geom_smooth(method=lm, se=T) +labs(x="dentary Left", y="dentary Right", title="dentary") + theme_bw()
denLR<-ggplot(data= datW, aes(x=dentary_lt, y=dentary_rt, 
                              color=morph)) + geom_point() + geom_smooth(method=lm, se=T, aes(fill=morph))+labs(x="dentary Left", y="dentary Right",col="Morph", fill="Morph", title="dentary" ) + guides(fill="none", col="none") + scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))+ scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")) + theme_bw()

ggplot(data= datW, aes(x=maxilla_lt, y=maxilla_rt)) + geom_point() + geom_smooth(method=lm, se=T) +labs(x="maxilla Left", y="maxilla Right", title="maxilla") + theme_bw()
maxLR<-ggplot(data= datW, aes(x=maxilla_lt, y=maxilla_rt, 
                              color=morph)) + geom_point() + geom_smooth(method=lm, se=T, aes(fill=morph))+labs(x="maxilla Left", y="maxilla Right",col="Morph", fill="Morph", title="maxilla" ) + guides(fill="none", col="none") + scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))+ scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")) + theme_bw()

ggplot(data= datW, aes(x=premaxilla_lt, y=premaxilla_rt)) + geom_point() + geom_smooth(method=lm, se=T) +labs(x="premaxilla Left", y="premaxilla Right", title="premaxilla") + theme_bw()
preLR<-ggplot(data= datW, aes(x=premaxilla_lt, y=premaxilla_rt, 
                              color=morph)) + geom_point() + geom_smooth(method=lm, se=T, aes(fill=morph))+labs(x="premaxilla Left", y="premaxilla Right",col="Morph", fill="Morph", title="premaxilla" ) + guides(fill="none", col="none") + scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))+ scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")) + theme_bw()

ggplot(data= datW, aes(x=palatine_lt, y=palatine_rt)) + geom_point() + geom_smooth(method=lm, se=T) +labs(x="palatine Left", y="palatine Right", title="palatine") + theme_bw()
palLR<-ggplot(data= datW, aes(x=palatine_lt, y=palatine_rt, 
                              color=morph)) + geom_point() + geom_smooth(method=lm, se=T, aes(fill=morph))+labs(x="palatine Left", y="palatine Right",col="Morph", fill="Morph", title="palatine" ) + guides(fill="none", col="none") + scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))+ scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")) + theme_bw()


pdf("VH_newcol.pdf")
multiplot(denLR, maxLR, preLR, palLR, layout = matrix(c(1,2,3,4) ,nrow = 2,byrow = T))
dev.off()

# make deviations for wide format
datW$MaxD<-as.numeric(datW$maxilla_lt-datW$maxilla_rt)
datW$DenD<-as.numeric(datW$dentary_lt-datW$dentary_rt)
datW$PreD<-as.numeric(datW$premaxilla_lt-datW$premaxilla_rt)
datW$PalD<-as.numeric(datW$palatine_lt-datW$palatine_rt)

datW$MaxDM1 <- as.numeric(datW$march_maxilla_lt - datW$march_maxilla_rt)
datW$DenDM1 <- as.numeric(datW$march_dentary_lt - datW$march_dentary_rt)
datW$PreDM1 <- as.numeric(datW$march_premaxilla_lt - datW$march_premaxilla_rt)
datW$PalDM1 <- as.numeric(datW$march_palatine_lt - datW$march_palatine_rt)

datW$MaxDM2 <- as.numeric(datW$may_maxilla_lt - datW$may_maxilla_rt)
datW$DenDM2 <- as.numeric(datW$may_dentary_lt - datW$may_dentary_rt)
datW$PreDM2 <- as.numeric(datW$may_premaxilla_lt - datW$may_premaxilla_rt)
datW$PalDM2 <- as.numeric(datW$may_palatine_lt - datW$may_palatine_rt)


# Mean
temp<-as.data.frame(table(datW$DenD))
sum(temp[,2])
sum(datW$DenD, na.rm = T)
35/sum(temp[,2])

temp<-as.data.frame(table(datW$MaxD))
sum(temp[,2])
sum(datW$MaxD, na.rm = T)
34/sum(temp[,2])

temp<-as.data.frame(table(datW$PreD))
sum(temp[,2])
sum(datW$PreD, na.rm = T)
45/sum(temp[,2])

temp<-as.data.frame(table(datW$PalD))
sum(temp[,2])
sum(datW$PalD, na.rm = T)
40/sum(temp[,2])

#March
temp<-as.data.frame(table(datW$DenDM1))
sum(temp[,2])
sum(datW$DenDM1, na.rm = T)
63/sum(temp[,2])

temp<-as.data.frame(table(datW$MaxDM1))
sum(temp[,2])
sum(datW$MaxDM1, na.rm = T)
54/sum(temp[,2])

temp<-as.data.frame(table(datW$PreDM1))
sum(temp[,2])
sum(datW$PreDM1, na.rm = T)
74/sum(temp[,2])

temp<-as.data.frame(table(datW$PalDM1))
sum(temp[,2])
sum(datW$PalDM1, na.rm = T)
54/sum(temp[,2])

#May
temp<-as.data.frame(table(datW$DenDM2))
sum(temp[,2])
sum(datW$DenDM2, na.rm = T)
80/sum(temp[,2])

temp<-as.data.frame(table(datW$MaxDM2))
sum(temp[,2])
sum(datW$MaxDM2, na.rm = T)
44/sum(temp[,2])

temp<-as.data.frame(table(datW$PreDM2))
sum(temp[,2])
sum(datW$PreDM2, na.rm = T)
73/sum(temp[,2])

temp<-as.data.frame(table(datW$PalDM2))
sum(temp[,2])
sum(datW$PalDM2, na.rm = T)
65/sum(temp[,2])

## Are mean deviations different from zero?
#Mean
t.test(datW$DenD) # 0.2
t.test(datW$MaxD) # 0.3
t.test(datW$PreD) # 0.027
t.test(datW$PalD) # 0.036

#March
t.test(datW$DenDM1) # 0.002
t.test(datW$MaxDM1) # 0.09
t.test(datW$PreDM1) # 0.04957
t.test(datW$PalDM1) # 0.1096

#May
t.test(datW$DenDM2) # 0.398
t.test(datW$MaxDM2) # 0.9093
t.test(datW$PreDM2) # 0.03939
t.test(datW$PalDM2) # 0.02071

## Skewness estimates.
skewness(datW$DenD, na.rm = T)
skewness(datW$MaxD, na.rm = T)
skewness(datW$PreD, na.rm = T)
skewness(datW$PalD, na.rm = T) 

kurtosis(datW$DenD, na.rm = T)
kurtosis(datW$MaxD, na.rm = T)
kurtosis(datW$PreD, na.rm = T)
kurtosis(datW$PalD, na.rm = T)

jarque.test(as.vector(na.omit(datW$DenD)))
jarque.test(as.vector(na.omit(datW$MaxD)))
jarque.test(as.vector(na.omit(datW$PreD)))
jarque.test(as.vector(na.omit(datW$PalD)))


#Image preparations
hist(datW$MaxD)
#MaxDM2
#All morph: coord_cartesian(y = c(0,50), x=c(-6,6)), mean
#All morphs: coord_cartesian(y = c(0,80), x=c(-6,6)), May
#By morph: coord_cartesian(y = c(0,20), x=c(-6,6)), mean

denVH_hist<-ggplot(datW, aes(DenDM2)) + geom_bar(fill = "#cc663aff", col = "black") + labs(x= "Difference, Left - Right", y="Count", title ="dentary")+ coord_cartesian(y = c(0,80), x=c(-6,6)) + theme_bw()
denVH_hist<-ggplot(datW, aes(DenDM2, fill = morph,color = morph)) + geom_bar(position="identity", alpha=0.3) +
  facet_grid(morph ~ .) + labs(x= "Difference, Left - Right", y="Count", col="Morph", fill="Morph", title ="dentary") + guides(col = "none", fill = "none") + coord_cartesian(y = c(0,31), x=c(-6,6)) +
  scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))+ scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")) + theme_bw()

maxVH_hist<-ggplot(datW, aes(MaxDM2)) + geom_bar(fill= "#398ebcff", col ="black") + labs(x= "Difference, Left - Right", y="Count", title ="maxilla")+ coord_cartesian(y = c(0,80), x=c(-6,6)) + theme_bw()
maxVH_hist<- ggplot(datW, aes(MaxDM2, fill = morph,color = morph)) + geom_bar(position="identity", alpha=0.3) +
  facet_grid(morph ~ .) + labs(x= "Difference, Left - Right", y="Count", col="Morph", fill="Morph", title ="maxilla") + guides(col = "none", fill = "none") + coord_cartesian(y = c(0,31), x=c(-6,6)) + 
  scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))+ scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")) + theme_bw()

preVH_hist<-ggplot(datW, aes(PreDM2)) + geom_bar(fill = "#0ee2feff", col = "black") + labs(x= "Difference, Left - Right", y="Count", title ="premaxilla")+ coord_cartesian(y = c(0,80), x=c(-6,6)) + theme_bw()
preVH_hist<-ggplot(datW, aes(PreDM2, fill = morph,color = morph)) + geom_bar(position="identity", alpha=0.3) +
  facet_grid(morph ~ .) + labs(x= "Difference, Left - Right", y="Count", col="Morph", fill="Morph", title ="premaxilla") + guides(col = "none", fill = "none") + coord_cartesian(y = c(0,31), x=c(-6,6)) +
  scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))+ scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")) + theme_bw()

palVH_hist<-ggplot(datW, aes(PalDM2)) + geom_bar(fill = "#0abbb5ff", col = "black") + labs(x= "Difference, Left - Right", y="Count", title ="palatine")+  coord_cartesian(y = c(0,80), x=c(-6,6)) + theme_bw()
palVH_hist<-ggplot(datW, aes(PalDM2, fill = morph,color = morph)) + geom_bar(position="identity", alpha=0.3) +
  facet_grid(morph ~ .) + labs(x= "Difference, Left - Right", y="Count", col="Morph", fill="Morph", title ="palatine") + guides(col = "none", fill = "none") + coord_cartesian(y = c(0,31), x=c(-6,6)) +
  scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))+ scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")) + theme_bw()


pdf("VH_hist.pdf")
multiplot(denVH_hist, maxVH_hist, preVH_hist, palVH_hist,  layout = matrix(c(1,2,3,4) ,nrow = 2,byrow = T))
dev.off()

## Do deviations differ by morph or are they effected by length

anova(glm(datW$DenD~datW$morph*log(datW$length_cm)), test = "Chisq")
summary(glm(datW$DenD~datW$morph*log(datW$length_cm)))

anova(glm(datW$MaxD~datW$morph*log(datW$length_cm)),test = "Chisq")
summary(glm(datW$MaxD~datW$morph*log(datW$length_cm)))

anova(glm(datW$PreD~datW$morph*log(datW$length_cm)),test = "Chisq")
summary(glm(datW$PreD~datW$morph*log(datW$length_cm)))

anova(glm(datW$PalD~datW$morph*log(datW$length_cm)),test = "Chisq")
summary(glm(datW$PalD~datW$morph*log(datW$length_cm)))

#Images preparations, Deviation between left and right
p1<-ggplot(data=datW, aes(x=log(length_cm),y=DenD)) + geom_point() + geom_smooth(method="lm", se=T) + labs(x="Log Length (cm)", y="Variation in teeth", title="dentary teeth deviation") + theme_bw()
p2<-ggplot(data=datW, aes(x=log(length_cm),y=DenD,color=morph)) + geom_point() + geom_smooth(method="lm", se=T, aes(fill=morph)) +labs(x="Log Length (cm)", y="Variation in teeth", color = "Morph",fill= "Morph",  title = "dentary teeth deviation by Morph")+guides(scale = "none") + theme_bw()
p1+(p2+ scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695")) + scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")))

p1<-ggplot(data=datW, aes(x=log(length_cm),y=MaxD)) + geom_point() + geom_smooth(method="lm", se=T) + labs(x="Log Length (cm)", y="Variation in teeth", title="maxilla teeth deviation") + theme_bw()
p2<-ggplot(data=datW, aes(x=log(length_cm),y=MaxD,color=morph)) + geom_point() + geom_smooth(method="lm", se=T, aes(fill=morph)) +labs(x="Log Length (cm)", y="Variation in teeth", color = "Morph", fill="Morph", title = "maxilla teeth deviation by Morph")+guides(scale = "none", fill="none", col="none") + theme_bw()
p1+(p2+ scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695")) + scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")))

p1<-ggplot(data=datW, aes(x=log(length_cm),y=PalD)) + geom_point() + geom_smooth(method="lm", se=T) + labs(x="Log Length (cm)", y="Variation in teeth", title="palatine teeth deviation") + theme_bw()
p2<-ggplot(data=datW, aes(x=log(length_cm),y=PalD,color=morph)) + geom_point() + geom_smooth(method="lm", se=T, aes(fill=morph)) +labs(x="Log Length (cm)", y="Variation in teeth", color = "Morph", fill="Morph", title = "palatine teeth deviation by Morph")+guides(scale = "none", fill="none", col="none") + theme_bw()
p1+(p2+ scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695")) + scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")))

p1<-ggplot(data=datW, aes(x=log(length_cm),y=PreD)) + geom_point() + geom_smooth(method="lm", se=T) + labs(x="Log Length (cm)", y="Variation in teeth", title="premaxilla teeth deviation") + theme_bw()
p2<-ggplot(data=datW, aes(x=log(length_cm),y=PreD,color=morph)) + geom_point() + geom_smooth(method="lm", se=T, aes(fill=morph)) +labs(x="Log Length (cm)", y="Variation in teeth", color = "Morph",fill="Morph", title = "premaxilla teeth deviation by Morph")+guides(scale = "none", fill="none", col="none") + theme_bw()
p1+(p2+ scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695")) + scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")))

## Do deviations correlate with each other
ggpairs(datW[,c(3,63:66)])

res <- rcorr(as.matrix(datW[,c(3,63:66)]))
round(res$P, 3)

##################################################
######## Data Examination #############

datW<- dplyr::group_by(datW, morph, sex)

kable(dplyr::summarise(filter(datW,!is.na(dentary), !is.na(morph)), 
                       mean = round(mean(dentary, na.rm = T),2), 
                       median= round((median(dentary, na.rm=T)),2), 
                       sd = round(sd(dentary, na.rm = T),2), count=n()))

datW<-dplyr::group_by(datW)

## Check for all bones: dentary, maxilla, premaxilla, palatine, vomer and glossohyal

###Visualizing the data

## Violin plots
#Non size corrected 
IMG<-ggplot(data=datW, aes(morph, y=dentary, fill=morph)) + geom_violin(trim=FALSE) +labs(x="Morph", y="Number of teeth", fill="Morph", col="Sex") + ggtitle("Dentary teeth") + geom_boxplot(width=0.3, fill="white",aes(col=sex))+ guides(fill="none", col="none")  +theme_bw()
IMG + scale_fill_manual(values=c("#74add1","#f46d43","#a50026","#313695")) + scale_color_manual(values = c("#af8dc3", "#7fbf7b")) 

## Run for all bones: dentary, maxilla, premaxilla, palatine, vomer and glossohyal

#Size corrected
IMG<-ggplot(data=datW, aes(morph, y=dentary/log(length_cm), fill=morph)) + geom_violin(trim=FALSE) +labs(x="Morph", y="Number of teeth, corrected for length", fill="Morph", col = "Sex") + ggtitle("Dentary teeth") + geom_boxplot(width=0.3, fill="white",aes(col=sex))+ guides(fill ="none", col = "none") + theme_bw()
IMG + scale_fill_manual(values=c("#74add1","#f46d43","#a50026","#313695")) + scale_color_manual(values = c("#af8dc3", "#7fbf7b")) 

## Run for all bones: dentary, maxilla, premaxilla, palatine, vomer and glossohyal

#Point/scatter plots
p1<-ggplot(data=datW, aes(x=log(length_cm),y=dentary)) + geom_point() + geom_smooth(method="lm", se=T) + labs(x="Log Length (cm)", y="Number of teeth", title="Dentary teeth") +theme_bw()
p2<-ggplot(data=datW, aes(x=log(length_cm),y=dentary,color=morph)) + geom_point() + geom_smooth(method="lm", se=T, aes(fill=morph)) +labs(x="Log Length (cm)", y="Number of teeth", color = "Morph", fill = "Morph", title = "Dentary teeth by Morph")  + theme_bw()
p1+(p2+ scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695")) + scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")))

## Run for all bones: dentary, maxilla, premaxilla, palatine, vomer and glossohyal

########################################################
######## Analyses of the long data #############

######################### glmmTMB #################

datL$MorSex<-paste(datL$morph,datL$sex)
datL$LogLen<-log(datL$length_cm)

datLX <- datL[!is.na(datL$sex),]
dim(datL)
dim(datLX)
datL<-datLX
head(datL)
str(datL)

# All Bones

BonesNr_mod2 <- glmmTMB(tooth_count ~ (bone + morph +sex + log2(length_cm))^2 
                      + (1 | morph:sex:id_me/replica),
                      data = datL)
summary(emtrends(BonesNr_mod2, pairwise~bone, var = "log2(length_cm)"), infert = T)
emmip(BonesNr_mod2, bone~ log2(length_cm), cov.reduce = range)

BonesNr_mod2 <- glmmTMB(tooth_count ~ (bone + morph +sex + length_cm)^2 
                      + (1 | morph:sex:id_me/replica),
                      data = datL)
summary(emtrends(BonesNr_mod2, pairwise~bone, var = "length_cm"), infert = T)
emmip(BonesNr_mod2, bone~ length_cm, cov.reduce = range)

BonesNr_mod2 <- glmmTMB(tooth_count ~ (bone + morph +sex + logCsize_dentary)^2 
                      + (1 | morph:sex:id_me/replica),
                      data = datL)
summary(emtrends(BonesNr_mod2, pairwise~bone, var = "logCsize_dentary"), infert = T)
emmip(BonesNr_mod2, bone~ logCsize_dentary, cov.reduce = range)


##Dentary analysis
#also run for other bones
DenNr_mod2 <- glmmTMB(tooth_count ~ (morph +sex + log2(length_cm))^2 + side + side:sex + side:morph
                      + (1 + side| morph:sex:id_me/replica),
                      data = datL[datL$bone=="Dentary",])
#OR

DenNr_mod2 <- glmmTMB(tooth_count ~ (morph +sex + logCsize_dentary)^2 + side + side:sex + side:morph
                      + (1 + side| morph:sex:id_me/replica),
                      data = datL[datL$bone=="Maxilla" & !is.na(datL$logCsize_dentary),])

summary(DenNr_mod2)
car::Anova(DenNr_mod2)

DenpairsM<- emmeans(DenNr_mod2, specs =~ morph)
DenpairsS<- emmeans(DenNr_mod2, specs =~ sex)
DenpairsMS<- emmeans(DenNr_mod2, specs =~ morph | sex)
DenpairsSM<- emmeans(DenNr_mod2, specs = ~sex | morph)
Denside<- emmeans(DenNr_mod2, specs =~ side)
DensideM<- emmeans(DenNr_mod2, specs =~ side |morph)
DensideS<- emmeans(DenNr_mod2, specs =~ side | sex)

pairs(DenpairsM)
confint(pairs(DenpairsM))
plot(DenpairsM,
     xlab = "model estimates, teeth count") +
  theme_bw() +
  theme(text = element_text(size = 16))

plot(pairs(DenpairsM)) +
  geom_vline(xintercept = 0, lty = 2 , alpha = 0.5) +
  xlab("Estimated difference in teeth count")

### Run for all other bones: Maxilla, Premaxilla, Palatine, Vomer and Glossohyal 

###################################################################
################# Teeth angle analysis ############################

###Reading data and general data preparation

angmethod<-read.csv2("Teeth_angle_method.csv", na.strings = "NA", stringsAsFactors = TRUE)

angdatW<-read.csv2("Teeth_angleW.csv", na.strings = "NA", stringsAsFactors = TRUE)
angdatL<-read.csv2("Teeth_angleL.csv", na.strings = "NA", sep=",", stringsAsFactors = TRUE)

angdatL$length_cm<-as.numeric(angdatL$length_cm)
angdatL$weight_g<-as.numeric(angdatL$weight_g)
angdatL$logCsize_M<- as.numeric(angdatL$logCsize_M)
angdatL$top<-as.numeric(angdatL$top)
angdatL$side<-as.numeric(angdatL$side)
angdatL$angle<-as.numeric(angdatL$angle)

##########################################
# Checking on the repeatability of the methods

angmethod$maxDvertleft<-as.numeric(angmethod$angle_vert_left_1-angmethod$angle_vert_left_2)
angmethod$maxDvertright<-as.numeric(angmethod$angle_vert_right_1-angmethod$angle_vert_right_2)
angmethod$maxDtsleft<-as.numeric(angmethod$angle_ts_left_1-angmethod$angle_ts_left_2)
angmethod$maxDtsright<-as.numeric(angmethod$angle_ts_right_1-angmethod$angle_ts_right_2)
angmethod$maxcombinevl<-as.numeric((angmethod$angle_vert_left_1+angmethod$angle_vert_left_2)/2)
angmethod$maxcombinevr<-as.numeric((angmethod$angle_vert_right_1+angmethod$angle_vert_right_2)/2)
angmethod$maxcombinetsl<-as.numeric((angmethod$angle_ts_left_1+angmethod$angle_ts_left_2)/2)
angmethod$maxcombinetsr<-as.numeric((angmethod$angle_ts_right_1+angmethod$angle_ts_right_2)/2)

### deviation between methods

#Fronto-lateral method
ggplot(data=angmethod, aes(x = angle_vert_left_1, y = angle_vert_left_2)) +
  geom_point() + geom_smooth(method = 'lm') + labs(x= "Trial 1 Left", y = "Trial 2 Left",  title = "Fronto-lateral method") + theme_bw() + 
  ggplot(data=angmethod, aes(x = angle_vert_right_1, y = angle_vert_right_2))+
  geom_point() + geom_smooth(method = 'lm') + labs(x= "Trial 1 Right", y = "Trial 2 Right") + theme_bw() 

modelVL <- lm(formula = angle_vert_left_1~angle_vert_left_2, data=angmethod)
summary(modelVL)
modelVR <- lm(formula = angle_vert_right_1~angle_vert_right_2, data=angmethod)
summary(modelVR)

MaxDiffVert<-c(angmethod$maxDvertleft,angmethod$maxDvertright)
hist(MaxDiffVert)

ggpairs(angmethod[,c(3,9,4,10)])
res <- rcorr(as.matrix(angmethod[,c(3,9,4,10)]))
round(res$P, 3)

#Medial/ventral method
ggplot(data=angmethod, aes(x = angle_ts_left_1, y = angle_ts_left_2)) +
  geom_point() + geom_smooth(method = 'lm') + labs(x= "Trial 1 Left", y = "Trial 2 Left",  title = "Medial/ventral method") + theme_bw() +
  ggplot(data=angmethod, aes(x = angle_ts_right_1, y = angle_ts_right_2))+
  geom_point() + geom_smooth(method = 'lm') + labs(x= "Trial 1 Right", y = "Trial 2 Right") + theme_bw()

modelTSL <- lm(formula = angle_ts_left_1~angle_ts_left_2, data=angmethod)
summary(modelTSL)
modelTSR <- lm(formula = angle_ts_right_1~angle_ts_right_2, data=angmethod)
summary(modelTSR)

MaxDiffTS<-c(angmethod$maxDtsleft,angmethod$maxDtsright)
hist(MaxDiffTS)

ggpairs(angmethod[,c(15,16,17,18)])
res1 <- rcorr(as.matrix(angmethod[,c(15,16,17,18)]))
round(res1$P, 3)

##it shows that the TS method is more repetible (not the best graph for that)

angmethod$maxcombinevl<-as.numeric((angmethod$angle_vert_left_1+angmethod$angle_vert_left_2)/2)
angmethod$maxcombinevr<-as.numeric((angmethod$angle_vert_right_1+angmethod$angle_vert_right_2)/2)
angmethod$maxcombinetsl<-as.numeric((angmethod$angle_ts_left_1+angmethod$angle_ts_left_2)/2)
angmethod$maxcombinetsr<-as.numeric((angmethod$angle_ts_right_1+angmethod$angle_ts_right_2)/2)

angmethod$maxDevvert<-as.numeric(angmethod$maxcombinevl-angmethod$maxcombinevr)
angmethod$maxDevts<-as.numeric(angmethod$maxcombinetsl-angmethod$maxcombinetsr)

hist(angmethod$maxDevvert,breaks=10)
hist(angmethod$maxDevts,breaks=10)

#It appears that the Medial/ventral method is more repeatable then the Fronto-lateral method. All following analysis will only be done on angel measurements using the Medial/ventral method.

##############
#### Angle analysis

angdatW$mean_side_5<-as.numeric((angdatW$side_left_5+angdatW$side_right_5)/2)
angdatW$mean_top_5<-as.numeric((angdatW$top_left_5+angdatW$top_right_5)/2)
angdatW$mean_side_1<-as.numeric((angdatW$side_left_1+angdatW$side_right_1)/2)
angdatW$mean_top_1<-as.numeric((angdatW$top_left_1+angdatW$top_right_1)/2)

angdatW$meanangle5<-as.numeric((angdatW$angle_left_5+angdatW$angle_right_5)/2)
angdatW$meanangle1<-as.numeric((angdatW$angle_left_1+angdatW$angle_right_1)/2)
angdatW$meanangle<-as.numeric((angdatW$angle_left_5+angdatW$angle_right_5+angdatW$angle_left_1+angdatW$angle_right_1)/4)

angdatW<- dplyr::group_by(angdatW, morph, sex)

kable(dplyr::summarise(filter(angdatW,!is.na(meanangle5), !is.na(morph)), 
                       mean = round(mean(meanangle5, na.rm = T),2), 
                       median= round((median(meanangle5, na.rm=T)),2), 
                       sd = round(sd(meanangle5, na.rm = T),2), count=n()))

angdatW<-dplyr::group_by(angdatW)
# Also run for meanangle1 and meanangle

#Visualizations of angle data
ggplot(angdatW, aes(meanangle5, fill = morph,color = morph)) + geom_histogram(position="identity", alpha=0.3, binwidth = 2.5) +
  facet_grid(morph ~ .) + labs(x= "Angle of 5th teeth", y = "Count", col="Morph", fill = "Morph") + theme_bw() + 
  scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695")) + scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695"))
ggplot(angdatW, aes(meanangle1, fill = morph,color = morph)) + geom_histogram(position="identity", alpha=0.3) +
  facet_grid(morph ~ .) + labs(x= "Angle of 1st teeth", y = "Count", col="Morph", fill = "Morph") + theme_bw() + 
  scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695")) + scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695"))
ggplot(angdatW, aes(meanangle, fill = morph,color = morph)) + geom_histogram(position="identity", alpha=0.3, binwidth = 2.5) +
  facet_grid(morph ~ .) + labs(x= "Angle of teeth", y = "Count", col="Morph", fill = "Morph") + theme_bw() + 
  scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695")) + scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695"))


scatterplot(angdatW$length_cm,angdatW$meanangle,groups = angdatW$morph, grid = TRUE, smooth = FALSE, regLine = FALSE, col= c("#74add1","#f46d43","#a50026","#313695"), pch=c(20,20,20,20), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(8,55), ylim=c(10,65), xlab = "Length (cm)", ylab = "Angle of teeth")
scatterplot(log(angdatW$length_cm),angdatW$meanangle,groups = angdatW$morph, grid = TRUE, smooth = FALSE, regLine = FALSE, col= c("#74add1","#f46d43","#a50026","#313695"), pch=c(20,20,20,20), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(2.2,4.3), ylim=c(10,65), xlab = "Log Length (cm)", ylab = "Angle of teeth")

scatterplot(angdatW$length_cm,angdatW$meanangle,groups = angdatW$sex, grid = FALSE, smooth = FALSE, regLine = FALSE, col= c("#af8dc3", "#7fbf7b"), pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),ylim=c(10,70),xlim=c(-5,60), xlab = "Length (cm)", ylab = "Angle of teeth")
scatterplot(angdatW$length_cm,angdatW$meanangle,groups = angdatW$morph_sex, grid = TRUE, smooth = FALSE, regLine = FALSE, col= c("#74add1","#74add1","#f46d43","#f46d43","#a50026","#a50026","#313695","#313695"), pch=c(20,17,20,17,20,17,20,17), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(8,55), ylim=c(10,65), xlab = "Length (cm)", ylab = "Angle of teeth")

#########################################################
### Can we transfer data to the long data format
## Testing variation between left and right and between measurements 
#for tooth 1 and 5

ggpairs(angdatW[,c(3,9,12,15,18)])
#R1*L1 = 0.698

ttest <- function(reg, coefnum, val){
  co <- coef(summary(reg))
  tstat <- (co[coefnum,1]-val)/co[coefnum,2]
  2 * pt(abs(tstat), reg$df.residual, lower.tail = FALSE)
}

R1_L1<-glm(angdatW$angle_right_1~angdatW$angle_left_1)
R1_L5<-glm(angdatW$angle_right_1~angdatW$angle_left_5)
R1_R5<-glm(angdatW$angle_right_1~angdatW$angle_right_5)
R5_L1<-glm(angdatW$angle_right_5~angdatW$angle_left_1)
L1_L5<-glm(angdatW$angle_left_1~angdatW$angle_left_5)
R5_L5<-glm(angdatW$angle_right_5~angdatW$angle_left_5)

R1_L1$coefficients # run also for other models

#Here we will use R1_L1 (0.6642361) as a reference and test whether other models differ from it.
ttest(R1_L5, 2, 0.66424) # 0.02685032 > 0.00625
ttest(R1_R5, 2, 0.66424) # 0.866025 > 0.00625
ttest(R5_L1, 2, 0.66424) # 0.05159788 > 0.00625
ttest(L1_L5, 2, 0.66424) # 0.8605534 > 0.00625
ttest(R5_L5, 2, 0.66424) # 0.8117855 > 0.00625

#Now we will use R5_L5 (0.67524) as the reference
R5_L5<-glm(angdatW$angle_right_5~angdatW$angle_left_5)
R5_L5$coefficients

ttest(R1_L5, 2, 0.67524) # 0.01559954 > 0.00625
ttest(R1_R5, 2, 0.67524) # 0.7015725 > 0.00625
ttest(R5_L1, 2, 0.67524) # 0.031945 > 0.00625
ttest(L1_L5, 2, 0.67524) # 0.9640765 > 0.00625
ttest(R1_L1, 2, 0.67524) # 0.8160331 > 0.00625

########################################################
######## Analyses of the long data, angle data  #############

######################### glmmTMB #################

with(angdatW,table(sex,morph))
with(angdatL,table(sex,morph))

angdatL$MorSex<-paste(angdatL$morph,angdatL$sex)
angdatL$LogLen<-log(angdatL$length_cm)

Maxang_mod2 <- glmmTMB(angle ~ (morph +sex + log2(length_cm))^2 + sideLR + sideLR:sex + sideLR:morph + tooth + tooth:sex + tooth:morph + tooth:sideLR
                      + (1 + sideLR| morph:sex:id_me) + (1 + tooth| morph:sex:id_me),
                      data = angdatL)
#OR

Maxang_mod2 <- glmmTMB(angle ~ (morph +sex + logCsize_M)^2 + sideLR + sideLR:sex + sideLR:morph + tooth + tooth:sex + tooth:morph + tooth:sideLR
                       + (1 + sideLR| morph:sex:id_me) + (1 + tooth| morph:sex:id_me),
                       data = angdatL)

summary(Maxang_mod2)
car::Anova(Maxang_mod2)

MaxpairsM <- emmeans(Maxang_mod2, specs =~ morph)
MaxpairsS <- emmeans(Maxang_mod2, specs =~ sex)
MaxpairsMS <- emmeans(Maxang_mod2, specs =~ morph|sex)
MaxpairsSM <- emmeans(Maxang_mod2, specs = ~ sex | morph)
Maxside<- emmeans(Maxang_mod2, specs =~ sideLR)
MaxsideM<- emmeans(Maxang_mod2, specs =~ sideLR |morph)
MaxsideS<- emmeans(Maxang_mod2, specs =~ sideLR | sex)

pairs(MaxpairsM)
confint(pairs(MaxpairsM))
plot(MaxpairsM,
     xlab = "model estimates, teeth Angle") +
  theme_bw() +
  theme(text = element_text(size = 16))
plot(pairs(MaxpairsM)) +
  geom_vline(xintercept = 0, lty = 2 , alpha = 0.5) +
  xlab("Estimated difference in teeth Angle")


############################################
######## 1997 survey analysis ############

###Reading data and general data preparation
dat97W<-read.csv("ttal_W.csv", na.strings = "NA",sep = ",", stringsAsFactors = TRUE)
dat97L<-read.csv("ttal_L.csv", na.strings = "NA", sep = ",", stringsAsFactors = TRUE)

### Teeth count

dat97W<- dplyr::group_by(dat97W, morph, sex)

kable(dplyr::summarise(filter(dat97W,!is.na(den_total), !is.na(morph)), mean = round(mean(den_total, na.rm = T),2), median= round((median(den_total, na.rm=T)),2), 
                       sd = round(sd(den_total, na.rm = T),2), count=n()))

dat97W<-dplyr::group_by(dat97W)
#Run for all other bones, Maxilla, premaxilla, palatine, vomer and glossohyal.

dat97WCo<-cor(dat97W[,c(3,13,25,19,16,22,10,26)],use="pairwise.complete.obs")
p.mat<-cor_pmat(dat97W[,c(3,13,25,19,16,22,10,26)],use="pairwise.complete.obs")
ggcorrplot(dat97WCo, p.mat = p.mat) + ggtitle("All morphs")

ggpairs(dat97W[,c(3,13,25,19,16,22,10,26)])
res <- rcorr(as.matrix(dat97W[,c(3,13,25,19,16,22,10,26)]))
round(res$P, 3)

##Visualization of data
ggplot(data=dat97W, aes(morph, y=den_total, fill=morph)) + geom_violin(trim=FALSE) +labs(x="Morph", y="Number of teeth", fill="Morph", col="Sex") + ggtitle("Dentary teeth") + geom_boxplot(width=0.3, fill="white",aes(col=sex))+ guides(fill="none", col="none")  + theme_bw() + scale_fill_manual(values=c("#74add1","#f46d43","#a50026","#313695")) + scale_color_manual(values = c("#af8dc3", "#7fbf7b"))
# Run for other bones.

p1<-ggplot(data=dat97W, aes(x=log(length_cm),y=den_total)) + geom_point() + geom_smooth(method="lm", se=T) + labs(x="Log Length (cm)", y="Number of teeth", title="Dentary teeth") + theme_bw()
p2<-ggplot(data=dat97W, aes(x=log(length_cm),y=den_total,color=morph)) + geom_point() + geom_smooth(method="lm", se=T, aes(fill=morph)) +labs(x="Log Length (cm)", y="Number of teeth", color = "Morph",fill = "Morph", title = "Dentary teeth by Morph")+guides(scale = "none") + theme_bw()
p1+(p2+ scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695")) + scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")))
# Run for all other bones.

######## glmmTMB analysis ############

# To insure that the model behave properly different log transformation need to be used for different bones:
#Dentary = log(length_cm)
#Maxilla = log10(length_cm)
#Premaxilla = log(length_cm)
#Palatine = log(length_cm)
#Vomer = log(length_cm)
#Glossohyal = log10(length_cm)

DenNr_mod97 <- glmmTMB(teeth_count ~ (morph+sex+log(length_cm))^2 
                      + (1| morph:sex:id), data = dat97L[dat97L$bone=="Dentary",])


summary(DenNr_mod97)
car::Anova(DenNr_mod97)

Denpairs97M<- emmeans(DenNr_mod97, specs =~ morph)
Denpairs97S<- emmeans(DenNr_mod97, specs =~ sex)
Denpairs97MS<- emmeans(DenNr_mod97, specs =~ morph | sex)
Denpairs97SM<- emmeans(DenNr_mod97, specs = ~sex | morph)

pairs(Denpairs97M)
confint(pairs(Denpairs97M))
plot(Denpairs97M,
     xlab = "model estimates, teeth count") +
  theme_bw() +
  theme(text = element_text(size = 16))

plot(pairs(Denpairs97M)) +
  geom_vline(xintercept = 0, lty = 2 , alpha = 0.5) +
  xlab("Estimated difference in teeth count")
# Also run for all other bones, remember:
#Dentary = log(length_cm)
#Maxilla = log10(length_cm)
#Premaxilla = log(length_cm)
#Palatine = log(length_cm)
#Vomer = log(length_cm)
#Glossohyal = log10(length_cm)

########################
#### Angle analysis 97 survey data ##

dat97W<- dplyr::group_by(dat97W, morph, sex)

kable(dplyr::summarise(filter(dat97W,!is.na(angle), !is.na(morph)), mean = round(mean(angle, na.rm = T),2), median = round((median(angle, na.rm=T)),2), 
                       sd = round(sd(angle, na.rm = T),2), count=n()))

dat97W<-dplyr::group_by(dat97W)

### Data Visualization 

scatterplot(dat97W$length_cm,dat97W$angle,groups = dat97W$morph, grid = TRUE, smooth = FALSE, regLine = FALSE, col= c("#74add1","#f46d43","#a50026","#313695"), pch=c(20,20,20,20), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(8,55), ylim=c(-15,80), xlab = "Length (cm)", ylab = "Angle of teeth")
scatterplot(dat97W$length_cm,dat97W$angle,groups = dat97W$morph_sex, grid = TRUE, smooth = FALSE, regLine = FALSE, col= c("#74add1","#74add1","#f46d43","#f46d43","#a50026","#a50026","#313695","#313695"), pch=c(20,17,20,17,20,17,20,17), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(8,55), ylim=c(-15,80), xlab = "Length (cm)", ylab = "Angle of teeth")

## glmmTMB

Maxang_mod97 <- glmmTMB(angle ~ (morph+sex+log10(length_cm))^2 
                        + (1| morph:sex:id), data = dat97W)
summary(Maxang_mod97)
car::Anova(Maxang_mod97)

Maxpairs97M<- emmeans(Maxang_mod97, specs =~ morph)
Maxpairs97S<- emmeans(Maxang_mod97, specs =~ sex)
Maxpairs97MS<- emmeans(Maxang_mod97, specs =~ morph | sex)
Maxpairs97SM<- emmeans(Maxang_mod97, specs = ~sex | morph)

pairs(Maxpairs97M)
confint(pairs(Maxpairs97M))
plot(Maxpairs97M,
     xlab = "model estimates, teeth count") +
  theme_bw() +
  theme(text = element_text(size = 16))

plot(pairs(Maxpairs97M)) +
  geom_vline(xintercept = 0, lty = 2 , alpha = 0.5) +
  xlab("Estimated difference in teeth count")



