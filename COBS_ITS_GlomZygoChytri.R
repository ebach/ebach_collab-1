#Elizabeth Bach
#COBS ITS data 2012
#Analysis of "minor" phyla: Glomeromycota, Zygomycota, Chytridiomycota
# 17 July 2014

#Glomeromycota
#Use Glom.data generate in COBS_ITS_PhylaSubset.R + graphs
head(Glom.data)
str(Glom.data)
#Only 16 observations, abundance 1-2
Glom.data
#Most are identified only to family Glomeraceae
#Two species identified: Entrophopora sp. (P13WS, PF23LM), and Paraglomus laccatum (P46 WS)
#These very few/low reads from Glomeromycota reflect primary bias, as root data indicates lots of AMF in these systems
#Or AMF, moslty associated with root tissue, not so much soil?
#Only 2 observations in PF, both times PF23
#Only 1 micro presence

#Zygomycota
#Use Zygo.data from COBS_ITS_PhylaSubset.R
head(Zygo.data)
str(Zygo.data)
#51 observations, 3 families
#start by exploring phyla-level relationships
Zygo.null<-lmer(value~1+(1|Block), data=Zygo.data, REML=FALSE)
Zygo.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Zygo.data, REML=FALSE)
Zygo.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Zygo.data, REML=FALSE)
AICtab(Zygo.null,Zygo.model.full,Zygo.model.main)
anova(Zygo.null,Zygo.model.full,Zygo.model.main)
#null model lowest AIC by 3, full model next best
anova(Zygo.model.main)
anova(Zygo.model.full)
#Significant Crop*SoilFrac and Crop*Date interactions
#modified model
Zygo.model<-lmer(value~Date+Crop+SoilFrac+Crop*SoilFrac+Date*SoilFrac+(1|Block), data=Zygo.data, REML=FALSE)
AICtab(Zygo.null,Zygo.model.full,Zygo.model.main,Zygo.model)
#null still lower AIC, but Zygo.model only 0.6 greater
anova(Zygo.model)
#Crop*SoilFrac interaction, P=0.02; Date*SoilFrac interaction, P=0.002
#Graph of differences
Zygo.phylaCS<-ddply(Zygo.data, .(Crop, SoilFrac), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Zygo.phylaCS)
ggplot(Zygo.phylaCS)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95, group=Crop, color=SoilFrac),position=position_jitter())+coord_flip()+scale_y_log10()
#high variability, PF, CC slightly higher than P, CC seems to have higher Zygo presence in SM and micro, PF mor in LM, WS, P not much difference
Zygo.phylaDS<-ddply(Zygo.data, .(Date, SoilFrac), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Zygo.phylaDS)
ggplot(Zygo.phylaDS)+geom_pointrange(aes(x=Date,y=mean,ymax=high95,ymin=low95, group=Date, color=SoilFrac),position=position_jitter())+coord_flip()+scale_y_log10()
#July sampling has greater Zygp presence in LM. MM, Oct. greater presence in micro, WS
difflsmeans(Zygo.model, ddf="Satterthwaite",type=3,method.grad="simple")
#In CC: micro>LM, MM>LM, WS>LM, micro<MM, MM>SM
#In P:  LM>micro, WS>micro
#In July: LM>micro, MM>micro, SM>micro, WS>micro
#In Oct:  micro>LM, MM>LM, WS>LM, micro>SM
#Intersting shift in Zygomycota presence from large to small aggregate fractions...
#Looking at families:
#Family Mucoraceae
Zygo.mucor<-subset(Zygo.data, Zygo.data$Family=="f__Mucoraceae")
head(Zygo.mucor)
str(Zygo.mucor)
#Only 7 observations
Zygo.mucor
#All are Mucor hiemalis F. hiemalis, only occur in PF, but all plots represented, all fractions
#M. hiemalis is a common plant pathogen

#Family Mortierellaceae
Zygo.mort<-subset(Zygo.data, Zygo.data$Family=="f__Mortierellaceae")
head(Zygo.mort)
str(Zygo.mort)
#9 observations
Zygo.mort
#All are from genus Mortierella, either unidentified species or Mortierella horticola
#Mostly in CC and PF samples, only one P sample
#micros, LM, and SM fractions, both sampling dates
#Common soil fungus, especially on root tissue, easily culturable

#Family Kickxellaceae
Zygo.kick<-subset(Zygo.data, Zygo.data$Family=="f__Kickxellaceae")
head(Zygo.kick)
str(Zygo.kick)
#35 Observations, the most abundant family
Zygo.null<-lmer(value~1+(1|Block), data=Zygo.kick, REML=FALSE)
Zygo.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Zygo.kick, REML=FALSE)
Zygo.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Zygo.kick, REML=FALSE)
AICtab(Zygo.null,Zygo.model.full,Zygo.model.main)
anova(Zygo.null,Zygo.model.full,Zygo.model.main)
#Full model lowest AIC by 22
anova(Zygo.model.full)
#no 3-way interaction, but all 2-ways are signifcant
Kick.model<-lmer(value~Date+Crop+SoilFrac+Date*Crop+Date*SoilFrac+Crop*SoilFrac+(1|Block), data=Zygo.kick, REML=FALSE)
AICtab(Zygo.null,Zygo.model.full,Zygo.model.main,Kick.model)
#Full model still better fit
Zygo.kickCS<-ddply(Zygo.kick, .(Crop, SoilFrac), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Zygo.kickCS)
ggplot(Zygo.kickCS)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95, group=Crop, color=SoilFrac),position=position_jitter())+coord_flip()+scale_y_log10()
Zygo.kickDS<-ddply(Zygo.kick, .(Date, SoilFrac), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Zygo.kickDS)
ggplot(Zygo.kickDS)+geom_pointrange(aes(x=Date,y=mean,ymax=high95,ymin=low95, group=Date, color=SoilFrac),position=position_jitter())+coord_flip()+scale_y_log10()
Zygo.kickCD<-ddply(Zygo.kick, .(Date, Crop), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Zygo.kickCD)
ggplot(Zygo.kickCD)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95, group=Crop, color=Date),position=position_jitter())+coord_flip()+scale_y_log10()
#With all these interactions, not sure there really is much of a story here, Kickxellaceae is just spotty, but present
#Common soil fungus, easily culturable


#Chytridiomycota, use Chytri.data from COBS_ITS_PhylaSubset.R
head(Chytri.data)
str(Chytri.data)
#42 observations
#order
Chytri.order<-ddply(Chytri.data, .(Order), summarise,.progress="text",
mean=mean(value),
high95=boot.high(value),
low95=boot.low(value)
)
head(Chytri.order)

ggplot(Chytri.order)+geom_pointrange(aes(x=Order,y=mean,ymax=high95,ymin=low95),position=position_dodge(width=1))+coord_flip()+scale_y_log10()
#5 orders, fairly low abundance in each order, look at phyla-level responses first
Chytri.null<-lmer(value~1+(1|Block), data=Chytri.data, REML=FALSE)
Chytri.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Chytri.data, REML=FALSE)
Chytri.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Chytri.data, REML=FALSE)
AICtab(Chytri.null,Chytri.model.full,Chytri.model.main)
anova(Chytri.null,Chytri.model.full,Chytri.model.main)
#Full model best fit, but 7.6
anova(Chytri.model.full)
#Significant 3-way interaction, P=0.002
#Taking a look at main responses first
Chytri.phylaC<-ddply(Chytri.data, .(Crop), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Chytri.phylaC)
ggplot(Chytri.phylaC)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95, color=Crop),position=position_jitter())+coord_flip()+scale_y_log10()
#not really any crop differences
Chytri.phylaSF<-ddply(Chytri.data, .(SoilFrac), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Chytri.phylaSF)
ggplot(Chytri.phylaSF)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_jitter())+coord_flip()+scale_y_log10()
#Only 1 occurance in SM, other fractions are similar
Chytri.phylaD<-ddply(Chytri.data, .(Date), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Chytri.phylaD)
ggplot(Chytri.phylaD)+geom_pointrange(aes(x=Date,y=mean,ymax=high95,ymin=low95, color=Date),position=position_jitter())+coord_flip()+scale_y_log10()
#No date differences
#Worth interpreting a 3-way interaction?  Is that driven by sporadic sample presence? Not significant in any other spotty groups, so maybe there is something in it...
#Chytridiomycota are "Water molds" so I wouldn't expect that they're big players, also primers biased against

#Order unidentified
Chytri.unid<-subset(Chytri.data, Chytri.data$Order=="o__unidentified")
head(Chytri.unid)
Chytri.unid
#Only 2 sampels, abundance 1-4
#Order Spizellomycetales
Chytri.spiz<-subset(Chytri.data, Chytri.data$Order=="o__Spizellomycetales")
head(Chytri.spiz)
Chytri.spiz
#Only 2 samples, genus Spizellomyces
#Order Rhizophydiales
Chytri.rhiz<-subset(Chytri.data, Chytri.data$Order=="o__Rhizophydiales")
head(Chytri.rhiz)
str(Chytri.rhiz)
#19 observations
Chytri.rhizG<-ddply(Chytri.rhiz, .(Genus), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Chytri.rhizG)
#All from genus Operculomyces
Chytri.rhiz
#Only found in P and PF!  Looks like balance of all fractions
#Abundance only 1-4

#Order Rhizophlyctidales, new order described in 2008, common, cellulose-degrading isolated from soils see Letcher et al. 2008, Mycological Research
Chytri.rhyc<-subset(Chytri.data, Chytri.data$Order=="o__Rhizophlyctidales")
head(Chytri.rhyc)
str(Chytri.rhyc)
#10 observations, all in genus Rhizophylctis, distributed among all cropping systems , LM, SM, micro fractions
Chytri.rhyc

#Order Chytridiales
Chytri.chyt<-subset(Chytri.data, Chytri.data$Order=="o__Chytridiales")
head(Chytri.chyt)
str(Chytri.chyt)
#Only 3 observations, low abundance