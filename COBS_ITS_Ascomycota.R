#Elizabeth Bach
#COBS ITS data 2012
#Analysis of Ascomycota only
# 15 July 2014

#Use Asco.data generated in COBS_ITS_PhylaSubset.R
#Graphs can befound in COBS_ITS_PhylaSubset.R as well
#will need these libraries for mixed models, data manipulation
library(reshape)
library(lme4)
library(lmerTest)
library(bbmle)

#Some orders show distinct SoilFrac and Agg, not much date effect
head(Asco.data)

Asco.null<-lmer(value~1+(1|Block), data=Asco.data, REML=FALSE)
Asco.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Asco.data, REML=FALSE)
Asco.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Asco.data, REML=FALSE)
AICtab(Asco.null,Asco.model.full,Asco.model.main)
anova(Asco.null,Asco.model.full,Asco.model.main)
#Main model best fit by 5.6
anova(Asco.model.main)
anova(Asco.model.full)
#No interactions in full model, proceed with main
difflsmeans(Asco.model.main, ddf="Satterthwaite",type=3,method.grad="simple")
#Significant effect of Date (Oct>July, P=0.001) and Crop (P=PF>CC, P=0.007)

#unidentified=Orbiliales
Asco.unid<-subset(Asco.data, Asco.data$Order=="o__unidentified")
head(Asco.unid)
str(Asco.unid)
Asco.orbi<-subset(Asco.unid, Asco.unid$Species=="s__Orbiliomycetessp")
str(Asco.orbi)
#All "unidentified order are identified as species Orbiliomycetes sp., without identification of higher taxanomic status
#Orbiliomycetes is class, Orbiliales is order, Orbiliaceae is family
Asco.null<-lmer(value~1+(1|Block), data=Asco.orbi, REML=FALSE)
Asco.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Asco.orbi, REML=FALSE)
Asco.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Asco.orbi, REML=FALSE)
AICtab(Asco.null,Asco.model.full,Asco.model.main)
anova(Asco.null,Asco.model.full,Asco.model.main)
#Main model best fit by 4.5
anova(Asco.model.main)
anova(Asco.model.full)
#No significant interactions, proceed with main effects model
difflsmeans(Asco.model.main, ddf="Satterthwaite",type=3,method.grad="simple")
#Significant effect of Date: Oct>July (P=0.02)
#Significant effect of Crop: PF>P=CC (P=0.0005)
Asco.orbiC<-ddply(Asco.orbi, .(Crop, Genus), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.orbiC)
ggplot(Asco.orbiC)+geom_pointrange(aes(x=Species,y=mean,ymax=high95,ymin=low95, color=Crop),position=position_dodge(width=1))+coord_flip()+scale_y_log10())

Asco.orbiD<-ddply(Asco.orbi, .(Date, Species), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.orbiD)
ggplot(Asco.orbiD)+geom_pointrange(aes(x=Species,y=mean,ymax=high95,ymin=low95, color=Date),position=position_dodge(width=1))+coord_flip()+scale_y_log10()


#Sordariales
Asco.sord<-subset(Asco.data, Asco.data$Order=="o__Sordariales")
head(Asco.sord)
str(Asco.sord)
#Only detected in 2 samples, P46.MM Oct. 2012 and PF41.micro Oct. 2012, very low abundance
#single species, Humicola nigescens
#Note, this order is Andy Miller (new INHS boss) area of expertise

#Pleosporales
Asco.pleo<-subset(Asco.data, Asco.data$Order=="o__Pleosporales")
head(Asco.pleo)
str(Asco.pleo)
Asco.null<-lmer(value~1+(1|Block), data=Asco.pleo, REML=FALSE)
Asco.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Asco.pleo, REML=FALSE)
Asco.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Asco.pleo, REML=FALSE)
AICtab(Asco.null,Asco.model.full,Asco.model.main)
anova(Asco.null,Asco.model.full,Asco.model.main)
#Null model best fit, by 5, main effects next best fit
anova(Asco.model.main)
anova(Asco.model.full)
#No significant effects
Asco.pleoC<-ddply(Asco.pleo, .(Crop, Genus), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.pleoC)
ggplot(Asco.pleoC)+geom_pointrange(aes(x=Genus,y=mean,ymax=high95,ymin=low95, color=Crop),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

Asco.pleoD<-ddply(Asco.pleo, .(Date, Genus), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.pleoD)
ggplot(Asco.pleoD)+geom_pointrange(aes(x=Genus,y=mean,ymax=high95,ymin=low95, color=Date),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

Asco.pleoSF<-ddply(Asco.pleo, .(SoilFrac, Genus), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.pleoSF)
ggplot(Asco.pleoSF)+geom_pointrange(aes(x=Genus,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()
#Genus Leptosphaerulina is most abundant, but doesn't differ by crop, date, or SoilFrac
#Genus Drechslera only found in 3 CC samples, abundance ranges 16-90, all Drechslera sp. BAFC3419
#Genus Alternaria very low abundance in only found in 2 samples, all Alternaria sp. 3MU_2012
Asco.drec<-subset(Asco.data, Asco.data$Genus=="g__Drechslera")
head(Asco.drec)

Asco.alt<-subset(Asco.data, Asco.data$Genus=="g__Alternaria")
head(Asco.alt)


#Pezizales
Asco.pez<-subset(Asco.data, Asco.data$Order=="o__Pezizales")
head(Asco.pez)
str(Asco.pez)
Asco.null<-lmer(value~1+(1|Block), data=Asco.pez, REML=FALSE)
Asco.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Asco.pez, REML=FALSE)
Asco.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Asco.pez, REML=FALSE)
AICtab(Asco.null,Asco.model.full,Asco.model.main)
anova(Asco.null,Asco.model.full,Asco.model.main)
#Null model lowest AIC by 1, main effects next best fit
anova(Asco.model.main)
anova(Asco.model.full)
#No significant interactions, so using main effects model
#Significant effect of date, Oct>July P=0.047
#SoilFrac P=0.08
#All samples either Peziza varia or unkn, both show date effect, unkn species in genus Peziza
Asco.pezD<-ddply(Asco.pez, .(Date, Species), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.pezD)
ggplot(Asco.pezD)+geom_pointrange(aes(x=Species,y=mean,ymax=high95,ymin=low95, color=Date),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

Asco.pezSF<-ddply(Asco.pez, .(SoilFrac, Species), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.pezSF)
ggplot(Asco.pezSF)+geom_pointrange(aes(x=Species,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()
#Peziza varia is common wood-decomposition cup fungus
#Only found in CC and PF plots, again no P indicating fertilizer effect or fluke?  No diff in CC and PF
Asco.pezC<-ddply(Asco.pez, .(Crop, Species), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.pezC)
ggplot(Asco.pezC)+geom_pointrange(aes(x=Species,y=mean,ymax=high95,ymin=low95, color=Crop),position=position_dodge(width=1))+coord_flip()+scale_y_log10()


#Incertaesedis
Asco.inc<-subset(Asco.data, Asco.data$Order=="o__Incertaesedis")
head(Asco.inc)
str(Asco.inc)
#15 observations
Asco.null<-lmer(value~1+(1|Block), data=Asco.inc, REML=FALSE)
Asco.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Asco.inc, REML=FALSE)
Asco.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Asco.inc, REML=FALSE)
AICtab(Asco.null,Asco.model.full,Asco.model.main)
anova(Asco.null,Asco.model.full,Asco.model.main)
#Data suffeciently unbalanced to make contrast models not possible
Asco.incD<-ddply(Asco.inc, .(Date, Species), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.incD)
ggplot(Asco.incD)+geom_pointrange(aes(x=Species,y=mean,ymax=high95,ymin=low95, color=Date),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

Asco.incSF<-ddply(Asco.inc, .(SoilFrac, Species), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.incSF)
ggplot(Asco.incSF)+geom_pointrange(aes(x=Species,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

Asco.incC<-ddply(Asco.inc, .(Crop, Species), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.incC)
ggplot(Asco.incC)+geom_pointrange(aes(x=Species,y=mean,ymax=high95,ymin=low95, color=Crop),position=position_dodge(width=1))+coord_flip()+scale_y_log10()
#Three species:  Vermispora spAS60291 (appears to be plant pathogen?), Scolecobasidium humicola (endophyte, promotes growth in organic tomato production), Polyscytalum pustulans (a potato plant pathogen)
#Also a couple of unkw species
#Very sporadic distribution among samples, no clear pattern, except all species consistently appear in PF, but not necessarily every plot


#Hypocreales
Asco.hypo<-subset(Asco.data, Asco.data$Order=="o__Hypocreales")
head(Asco.hypo)
str(Asco.hypo)
#150 observations!
Asco.null<-lmer(value~1+(1|Block), data=Asco.hypo, REML=FALSE)
Asco.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Asco.hypo, REML=FALSE)
Asco.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Asco.hypo, REML=FALSE)
AICtab(Asco.null,Asco.model.full,Asco.model.main)
anova(Asco.null,Asco.model.full,Asco.model.main)
#Null model and main model essentially same AIC
anova(Asco.model.main)
anova(Asco.model.full)
#Date*Crop interaction is significant in full model, P=0.03, Crop is significant in main effects model (P=0.01)
#Truncated model to include main factors + Date*Crop interaction
Hypo.model<-lmer(value~Date+Crop+SoilFrac+Date*Crop+(1|Block), data=Asco.hypo, REML=FALSE)
AICtab(Asco.null,Asco.model.full,Asco.model.main, Hypo.model)
#Best AIC fit
anova(Hypo.model)
#Date*Crop interaction significant P=0.04
Asco.hypoCD<-ddply(Asco.hypo, .(Date, Crop, Order), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.hypoCD)
ggplot(Asco.hypoCD)+geom_pointrange(aes(x=Order,y=mean,ymax=high95,ymin=low95, group=Crop, color=Date),position=position_dodge(width=1))+coord_flip()+scale_y_log10()
#Interaction is from similar abundance in P and PF at both sampling dates, but greater presence in CC in Oct
#Includes three families, including unk, and there are some different patterns within families, so will look more closely at that
Asco.hypoCD<-ddply(Asco.hypo, .(Date, Crop, Family), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.hypoCD)
ggplot(Asco.hypoCD)+geom_pointrange(aes(x=Family,y=mean,ymax=high95,ymin=low95, group=Crop, color=Date),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

#Family Hypocreaceae
Asco.hypoc<-subset(Asco.data, Asco.data$Family=="f__Hypocreaceae")
head(Asco.hypoc)
str(Asco.hypoc)
#64 observations
Asco.null<-lmer(value~1+(1|Block), data=Asco.hypoc, REML=FALSE)
Asco.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Asco.hypoc, REML=FALSE)
Asco.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Asco.hypoc, REML=FALSE)
AICtab(Asco.null,Asco.model.full,Asco.model.main)
anova(Asco.null,Asco.model.full,Asco.model.main)
#Null model lowest AIC by 4
anova(Asco.model.main)
anova(Asco.model.full)
#No interactions, Crop significant (P=0.04), PF=CC>P, all observations are Trichoderma citrinoviride, except 1 unk
difflsmeans(Asco.model.main, ddf="Satterthwaite",type=3,method.grad="simple")
Asco.hypoC<-ddply(Asco.hypoc, .(Crop, Species), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.hypoC)
ggplot(Asco.hypoC)+geom_pointrange(aes(x=Species,y=mean,ymax=high95,ymin=low95, color=Crop),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

#Family Bionectriaceae
Asco.bion<-subset(Asco.data, Asco.data$Family=="f__Bionectriaceae")
head(Asco.bion)
str(Asco.bion)
#33 observations
Asco.null<-lmer(value~1+(1|Block), data=Asco.bion, REML=FALSE)
Asco.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Asco.bion, REML=FALSE)
Asco.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Asco.bion, REML=FALSE)
AICtab(Asco.null,Asco.model.full,Asco.model.main)
anova(Asco.null,Asco.model.full,Asco.model.main)
#Main effects model has lowest AIC
anova(Asco.model.main)
anova(Asco.model.full)
#No interactions, SoilFrac highly significant: P=0.006, highly abundant in LM and micro, small presence in MM, only 1 in SM and WS
difflsmeans(Asco.model.main, ddf="Satterthwaite",type=3,method.grad="simple")
Asco.bionSF<-ddply(Asco.bion, .(SoilFrac, Family), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.bionSF)
ggplot(Asco.bionSF)+geom_pointrange(aes(x=Family,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()
#Species:  Clonostachysroseaf.catenulata + unkw, family-level relationship is same for both species
Asco.bionSF<-ddply(Asco.bion, .(SoilFrac, Species), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.bionSF)
ggplot(Asco.bionSF)+geom_pointrange(aes(x=Species,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

#unk Family
Asco.unk<-subset(Asco.hypo, Asco.hypo$Family=="")
head(Asco.unk)
str(Asco.unk)
#53 observations, abundances ~10-50
Asco.null<-lmer(value~1+(1|Block), data=Asco.unk, REML=FALSE)
Asco.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Asco.unk, REML=FALSE)
Asco.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Asco.unk, REML=FALSE)
AICtab(Asco.null,Asco.model.full,Asco.model.main)
anova(Asco.null,Asco.model.full,Asco.model.main)
#Null model is lowest AIC by 3, main model is next best fit
anova(Asco.model.full)
anova(Asco.model.main)
#Full model shows significant interaction of Date*Crop (P=0.004), main effects of Date and Crop both marginal at P=0.06
#Probably not worth interpreting this, just leave at order level

Asco.unkC<-ddply(Asco.unk, .(Crop, Family), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.unkC)
ggplot(Asco.unkC)+geom_pointrange(aes(x=Family,y=mean,ymax=high95,ymin=low95, color=Crop),position=position_dodge(width=1))+coord_flip()+scale_y_log10()


#Helotiales
Asco.helo<-subset(Asco.data, Asco.data$Order=="o__Helotiales")
head(Asco.helo)
str(Asco.helo)
#13 observations, low abundances, probably not contributing much
Asco.heloD<-ddply(Asco.helo, .(Date, Species), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.heloD)
ggplot(Asco.heloD)+geom_pointrange(aes(x=Species,y=mean,ymax=high95,ymin=low95, color=Date),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

Asco.heloSF<-ddply(Asco.helo, .(SoilFrac, Species), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.heloSF)
ggplot(Asco.heloSF)+geom_pointrange(aes(x=Species,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

Asco.heloC<-ddply(Asco.helo, .(Crop, Species), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.heloC)
ggplot(Asco.heloC)+geom_pointrange(aes(x=Species,y=mean,ymax=high95,ymin=low95, color=Crop),position=position_dodge(width=1))+coord_flip()+scale_y_log10()
#Generally low abundance, consistently present in Oct., but not July.  Not consistent across crop or aggregate fraction


#Capnodiales
Asco.cap<-subset(Asco.data, Asco.data$Order=="o__Capnodiales")
head(Asco.cap)
str(Asco.cap)
#Only 6 observations, all Cladosporium sp4MU_2012, value=1, no consistency in samples, super common air-borne fungus
Asco.cap

#unk
Asco.unk<-subset(Asco.data, Asco.data$Order=="")
head(Asco.unk)
str(Asco.unk)
#113 observations, can't do much interpretation, but looking a patterns just incase something pops out
Asco.null<-lmer(value~1+(1|Block), data=Asco.unk, REML=FALSE)
Asco.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Asco.unk, REML=FALSE)
Asco.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Asco.unk, REML=FALSE)
AICtab(Asco.null,Asco.model.full,Asco.model.main)
anova(Asco.null,Asco.model.full,Asco.model.main)
#Full model best fit by 5
anova(Asco.model.full)
#Significant Crop*SoilFrac interaction, truncated model
Asco.model.unk<-lmer(value~Date+Crop+SoilFrac+Crop*SoilFrac+(1|Block), data=Asco.unk, REML=FALSE)
AICtab(Asco.null,Asco.model.full,Asco.model.main,Asco.model.unk)
#Much better fit
anova(Asco.model.unk)
difflsmeans(Asco.model.unk, ddf="Satterthwaite",type=3,method.grad="simple")
#Highly significnat Crop*SoilFrac interaction, P<0.0001, 
Asco.unkCS<-ddply(Asco.unk, .(Crop, SoilFrac, Family), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.unkCS)
ggplot(Asco.unkCS)+geom_pointrange(aes(x=Family,y=mean,ymax=high95,ymin=low95, group=Crop, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()
#Each cropping system has unique rank of aggregate fractions, lots of overlap
#Again, no power to pull this apart further, but interesting to note amount of unidentified Ascomycota

