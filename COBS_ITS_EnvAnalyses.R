#Elizabeth Bach
#mixed models of environmental variables with COBS ITS dataset
#22 July 2014

library(lme4)
library(lmerTest)
library(bbmle)


metadata<-read.csv("COBS_ITS_metadata.csv")
str(metadata)
metadata$ExtC<-as.numeric(levels(metadata$ExtC))[metadata$ExtC]
metadata$MBC<-as.numeric(levels(metadata$MBC))[metadata$MBC]
metadata$MBC.MBN<-as.numeric(levels(metadata$MBC.MBN))[metadata$MBC.MBN]
metadata$TC<-as.numeric(levels(metadata$TC))[metadata$TC]
metadata$TN<-as.numeric(levels(metadata$TN))[metadata$TN]
metadata$CN<-as.numeric(levels(metadata$CN))[metadata$CN]
str(metadata)
head(metadata)

#Focusing on pH, water_content, and roots first, highest correlations (see AMFcol_analyses.R for AMF)
hist(metadata$ph)
pH.null<-lmer(ph~1+(1|Block),data=metadata,REML=FALSE)
pH.model.full<-lmer(ph~Date*Crop+(1|Block), data=metadata, REML=FALSE)
pH.model.main<-lmer(ph~Date+Crop+(1|Block), data=metadata, REML=FALSE)
AICtab(pH.null,pH.model.full,pH.model.main)
anova(pH.null,pH.model.full,pH.model.main)
#Main model lowet AIC, but only by 0.3
anova(pH.model.full)
#No interactions, so proceed with main effects
anova(pH.model.main)
difflsmeans(pH.model.main, ddf="Satterthwaite",type=3,method.grad="simple")
#Very strong effect of crop, P=PF>>CC

hist(metadata$water_content_soil)
GWC.null<-lmer(water_content_soil~1+(1|Block),data=metadata,REML=FALSE)
GWC.model.full<-lmer(water_content_soil~Date*Crop+(1|Block), data=metadata, REML=FALSE)
GWC.model.main<-lmer(water_content_soil~Date+Crop+(1|Block), data=metadata, REML=FALSE)
AICtab(GWC.null,GWC.model.full,GWC.model.main)
anova(GWC.null,GWC.model.full,GWC.model.main)
#main model lowest AIC, by 4
anova(GWC.model.full)
#No interactions
anova(GWC.model.main)
difflsmeans(GWC.model.main, ddf="Satterthwaite",type=3,method.grad="simple")
ddply(metadata, .(Crop), summarise, mean=mean(water_content_soil)) 
#Significant effect of date and crop, Oct>July, P=PF>CC

hist(metadata$RootBiomass)
#Three distict "peaks" representing the cropping systems
Root.null<-lmer(RootBiomass~1+(1|Block),data=metadata,REML=FALSE)
Root.model.full<-lmer(RootBiomass~Date*Crop+(1|Block), data=metadata, REML=FALSE)
Root.model.main<-lmer(RootBiomass~Date+Crop+(1|Block), data=metadata, REML=FALSE)
AICtab(Root.null,Root.model.full,Root.model.main)
anova(Root.null,Root.model.full,Root.model.main)
#Main model lowest AIC, by 4
anova(Root.model.main)
#Super strong crop response (of course)

TN.data<-cbind(metadata[,1:5],metadata[24])
TN.data2<-subset(TN.data, TN.data$TN>0)
TN.data2$TN
hist(TN.data2$TN)
TN.null<-lmer(TN~1+(1|Block),data=TN.data2,REML=FALSE)
TN.model.full<-lmer(TN~Date*Crop+(1|Block), data=TN.data2, REML=FALSE)
TN.model.main<-lmer(TN~Date+Crop+(1|Block), data=TN.data2, REML=FALSE)
AICtab(TN.null,TN.model.full,TN.model.main)
anova(TN.null,TN.model.full,TN.model.main)
#Main model best fit by 4
anova(TN.model.full)
#No interactions
anova(TN.model.main)
difflsmeans(TN.model.main, ddf="Satterthwaite",type=3,method.grad="simple")
ddply(TN.data2, .(Crop), summarise, mean=mean(TN)) 
#Very strong crop effect, PF>P>CC

TC.data<-cbind(metadata[,1:5],metadata[23])
TC.data2<-subset(TC.data, TC.data$TC>0)
TC.data2$TC
hist(TC.data2$TC)
TC.null<-lmer(TC~1+(1|Block),data=TC.data2,REML=FALSE)
TC.model.full<-lmer(TC~Date*Crop+(1|Block), data=TC.data2, REML=FALSE)
TC.model.main<-lmer(TC~Date+Crop+(1|Block), data=TC.data2, REML=FALSE)
AICtab(TC.null,TC.model.full,TC.model.main)
anova(TC.null,TC.model.full,TC.model.main)
#Main model best fit by ~4
anova(TC.model.full)
#No interactions
anova(TC.model.main)
difflsmeans(TC.model.main, ddf="Satterthwaite",type=3,method.grad="simple")
ddply(TC.data2, .(Crop), summarise, mean=mean(TC)) 
#Very strong crop effect, PF>P>CC

ExtN.data<-cbind(metadata[,1:5],metadata[11])
ExtN.data2<-subset(ExtN.data, ExtN.data$ExtN>0)
ExtN.data2$ExtN
hist(ExtN.data2$ExtN)
#skewed left
ExtN.null<-lmer(ExtN~1+(1|Block),data=ExtN.data2,REML=FALSE)
ExtN.model.full<-lmer(ExtN~Date*Crop+(1|Block), data=ExtN.data2, REML=FALSE)
ExtN.model.main<-lmer(ExtN~Date+Crop+(1|Block), data=ExtN.data2, REML=FALSE)
AICtab(ExtN.null,ExtN.model.full,ExtN.model.main)
anova(ExtN.null,ExtN.model.full,ExtN.model.main)
#Main model best fit by 3
anova(ExtN.model.full)
#No interactions
anova(ExtN.model.main)
difflsmeans(ExtN.model.main, ddf="Satterthwaite",type=3,method.grad="simple")
ddply(ExtN.data2, .(Crop), summarise, mean=mean(ExtN)) 
#Very strong crop effect, CC>PF=P

#Now look at variable least likely to have crop effect ?CN
CN.data<-cbind(metadata[,1:5],metadata[25])
CN.data2<-subset(CN.data, CN.data$CN>0)
CN.data2$CN
hist(CN.data2$CN)
CN.null<-lmer(CN~1+(1|Block),data=CN.data2,REML=FALSE)
CN.model.full<-lmer(CN~Date*Crop+(1|Block), data=CN.data2, REML=FALSE)
CN.model.main<-lmer(CN~Date+Crop+(1|Block), data=CN.data2, REML=FALSE)
AICtab(CN.null,CN.model.full,CN.model.main)
anova(CN.null,CN.model.full,CN.model.main)
#Full model lowest AIC, by 0.5
anova(CN.model.full)
#Date*Crop interaction, P=0.1
anova(CN.model.main)
difflsmeans(CN.model.main, ddf="Satterthwaite",type=3,method.grad="simple")
ddply(CN.data2, .(Crop), summarise, mean=mean(CN)) 
#Strong crop effect, PF>P=CC

#Enzymes: BX & NAG
BX.data<-cbind(metadata[,1:6],metadata[18])
BX.data2<-subset(BX.data, BX.data$BX>0)
#removing 1 outlier
BX.data3<-subset(BX.data2, BX.data2$BX<400)
BX.data3$BX
hist(BX.data3$BX)
head(BX.data3)
#Skewed somewhat left, log transform is good idea with all the enzymes
BX.null<-lmer(log(BX+1)~1+(1|Block),data=BX.data3,REML=FALSE)
BX.model.full<-lmer(log(BX+1)~Date*Crop*SoilFrac+(1|Block), data=BX.data3, REML=FALSE)
BX.model.main<-lmer(log(BX+1)~Date+Crop+SoilFrac+(1|Block), data=BX.data3, REML=FALSE)
AICtab(BX.null,BX.model.full,BX.model.main)
anova(BX.null,BX.model.full,BX.model.main)
#main model lowest AIC, by 0.9
anova(BX.model.full)
#Date*Crop interaction, building model to include that interaction only
BX.model<-lmer(log(BX+1)~Date+Crop+SoilFrac+Date*Crop+(1|Block), data=BX.data3, REML=FALSE)
AICtab(BX.null,BX.model.full,BX.model.main,BX.model)
anova(BX.null,BX.model.full,BX.model.main,BX.model)
#BX.model lowest AIC by 23
anova(BX.model)
#Crop*Date interaction, P<0.0001, but no significant effect of soil aggregate

NAG.data<-cbind(metadata[,1:6],metadata[20])
NAG.data2<-subset(NAG.data, NAG.data$NAG>0)
hist(NAG.data2$NAG)
#removing 1 outlier
NAG.data3<-subset(NAG.data2, NAG.data2$NAG<400)
NAG.data3$NAG
hist(NAG.data3$NAG)
head(NAG.data3)
#Skewed somewhat left, log transform is good idea with all the enzymes
NAG.null<-lmer(log(NAG+1)~1+(1|Block),data=NAG.data3,REML=FALSE)
NAG.model.full<-lmer(log(NAG+1)~Date*Crop*SoilFrac+(1|Block), data=NAG.data3, REML=FALSE)
NAG.model.main<-lmer(log(NAG+1)~Date+Crop+SoilFrac+(1|Block), data=NAG.data3, REML=FALSE)
AICtab(NAG.null,NAG.model.full,NAG.model.main)
anova(NAG.null,NAG.model.full,NAG.model.main)
#full model best fit by 30
anova(NAG.model.full)
#Date*Crop*Soil 3-way interaction P=0.05, Date*Crop interaction very significant
#Model 1 with both Date*Crop and 3-way
NAG.model1<-lmer(log(NAG+1)~Date+Crop+SoilFrac+Date*Crop+Date*Crop*SoilFrac+(1|Block), data=NAG.data3, REML=FALSE)
#Model 2 with only Date*Crop
NAG.model2<-lmer(log(NAG+1)~Date+Crop+SoilFrac+Date*Crop+(1|Block), data=NAG.data3, REML=FALSE)
AICtab(NAG.null,NAG.model.full,NAG.model.main,NAG.model1,NAG.model2)
anova(NAG.null,NAG.model.full,NAG.model.main,NAG.model1,NAG.model2)
#NAG.model2 has lowest AIC, by 16, better fit
anova(NAG.model2) 
#Date*Crop interaction significant, P<0.0001; SoilFrac main effect also significant P=0.002
