#Elizabeth Bach
#COBS ITS data 2012
#Analysis for Basidiomycota only
# 14 July 2014

#Use Basidio.data generated in COBS_ITS_PhylaSubset.R
#Graphs can befound in COBS_ITS_PhylaSubset.R as well
#will need these libraries for mixed models, data manipulation
library(reshape)
library(lme4)
library(lmerTest)
library(bbmle)

#Analysis of Basiodiomycota sub-set data (COBS_ITS_PhylAnalyses.R)
#Main effects model best, Crop is only significant factor P=0.01 for fully Phylum
#Note: Basidio.order from COBS_ITS_PhlyaSubset.R is mean values, for graphing, use totals generated here for analysis

#Basidio.orderT<-ddply(Basidio.data, .(SampleName, Date, Block, Crop, SoilFrac, Order), summarise,.progress="text",value=sum(value))
#head(Basidio.orderT)

Basidio.null<-lmer(value~1+(1|Block), data=Basidio.orderT, REML=FALSE)
Basidio.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Basidio.orderT, REML=FALSE)
Basidio.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Basidio.orderT, REML=FALSE)
AICtab(Basidio.null,Basidio.model.full,Basidio.model.main)
anova(Basidio.null,Basidio.model.full,Basidio.model.main)
anova(Basidio.model.main)
anova(Basidio.model.full)
#Null model best fit,  Crop is marginally siginifcant in main model (next best fit), nothing significant in full model
#Order as a variable is higly significant
#Zoom-in, focus on Corticiales, Agaricales (most abundant), and Hymenochaetales, Geastrales, Cantharellaes (most different among fractions)

#Order corticiales
Basidio.corti<-subset(Basidio.data, Basidio.data$Order=="o__Corticiales")
head(Basidio.corti)
str(Basidio.corti)
#Only 18 observations
Basidio.corti
#Only identified one species, Limonomyces sp.
#Only present in P, but found in all plots, mix of aggregate fractions, sampling dates
#Looking at SoilFrac and Date responses
Basidio.null<-lmer(value~1+(1|Block), data=Basidio.agari, REML=FALSE)
Basidio.model.full<-lmer(value~Date*SoilFrac+(1|Block), data=Basidio.agari, REML=FALSE)
Basidio.model.main<-lmer(value~Date+SoilFrac+(1|Block), data=Basidio.agari, REML=FALSE)
AICtab(Basidio.null,Basidio.model.full,Basidio.model.main)
#Null model best fit, main is next by 6.2
anova(Basidio.model.main)
anova(Basidio.model.full)
#No factors or interactions significant

#Order Agaricales
Basidio.agari<-subset(Basidio.data, Basidio.data$Order=="o__Agaricales")
head(Basidio.agari)
Basidio.null<-lmer(value~1+(1|Block), data=Basidio.agari, REML=FALSE)
Basidio.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Basidio.agari, REML=FALSE)
Basidio.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Basidio.agari, REML=FALSE)
AICtab(Basidio.null,Basidio.model.full,Basidio.model.main)
anova(Basidio.null,Basidio.model.full,Basidio.model.main)
#null model best fit by 6, main is next best
anova(Basidio.model.main)
anova(Basidio.model.full)
#No significant effects of Date, Crop, or SoilFrac within the order
Basidio.agariSF<-ddply(Basidio.agari, .(SoilFrac, Family), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Basidio.agariSF)
ggplot(Basidio.agariSF)+geom_pointrange(aes(x=Family,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()
#Families Bolbitiaceae, Entolomataceae, Typhulaceae, and Psathyrellaceae merit follow-up analysis
#Inocybaceae only present in micros!
Basidio.agariC<-ddply(Basidio.agari, .(Crop, Family), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Basidio.agariC)
ggplot(Basidio.agariC)+geom_pointrange(aes(x=Family,y=mean,ymax=high95,ymin=low95, color=Crop),position=position_dodge(width=1))+coord_flip()+scale_y_log10()
#Families Strophariaceae, Psathyrellaceae, Marasmiaceae, Inocybaceae merit further investigation from crop perspective
Basidio.agariD<-ddply(Basidio.agari, .(Date, Family), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Basidio.agariD)
ggplot(Basidio.agariD)+geom_pointrange(aes(x=Family,y=mean,ymax=high95,ymin=low95, color=Date),position=position_dodge(width=1))+coord_flip()+scale_y_log10()
#Inocybaceae only present in July, Marasmiaceae much higher in July, Niaceae only present in Oct

#Family Bolbitiaceae
Basidio.bolbi<-subset(Basidio.data, Basidio.data$Family=="f__Bolbitiaceae")
head(Basidio.bolbi)
Basidio.null<-lmer(value~1+(1|Block), data=Basidio.bolbi, REML=FALSE)
Basidio.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Basidio.bolbi, REML=FALSE)
Basidio.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Basidio.bolbi, REML=FALSE)
AICtab(Basidio.null,Basidio.model.full,Basidio.model.main)
anova(Basidio.null,Basidio.model.full,Basidio.model.main)
#null best fit by 7, main effects next best fit
anova(Basidio.model.main)
anova(Basidio.model.full)
#no significance

#Family Entolomataceae
Basidio.entol<-subset(Basidio.data, Basidio.data$Family=="f__Entolomataceae")
head(Basidio.entol)
Basidio.null<-lmer(value~1+(1|Block), data=Basidio.entol, REML=FALSE)
Basidio.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Basidio.entol, REML=FALSE)
Basidio.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Basidio.entol, REML=FALSE)
AICtab(Basidio.null,Basidio.model.full,Basidio.model.main)
anova(Basidio.null,Basidio.model.full,Basidio.model.main)
#null model lowest AIC, by 3, main effects next best fit
anova(Basidio.model.main)
anova(Basidio.model.full)
#NS effects, SoilFrac P=0.1

#Family Typhulaceae
Basidio.typhul<-subset(Basidio.data, Basidio.data$Family=="f__Typhulaceae")
head(Basidio.typhul)
Basidio.null<-lmer(value~1+(1|Block), data=Basidio.typhul, REML=FALSE)
Basidio.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Basidio.typhul, REML=FALSE)
Basidio.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Basidio.typhul, REML=FALSE)
AICtab(Basidio.null,Basidio.model.full,Basidio.model.main)
anova(Basidio.null,Basidio.model.full,Basidio.model.main)
#null best fit by 8.5, full model next best fit
anova(Basidio.model.full)
#Crop*SoilFrac interaction, P=0.05
Basidio.typhulFig<-ddply(Basidio.typhul, .(SoilFrac, Crop, Species), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Basidio.typhulFig)
ggplot(Basidio.typhulFig)+geom_pointrange(aes(x=Species,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()
#1 genus, 2 speceis present, not present in any LM, and only 1 MM, abundance in SM and Micro varies by crop

#Family Psathyrellaceae -> This one has interesting crop and SoilFrac response
Basidio.psath<-subset(Basidio.data, Basidio.data$Family=="f__Psathyrellaceae")
head(Basidio.psath)
Basidio.null<-lmer(value~1+(1|Block), data=Basidio.psath, REML=FALSE)
Basidio.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Basidio.psath, REML=FALSE)
Basidio.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Basidio.psath, REML=FALSE)
AICtab(Basidio.null,Basidio.model.full,Basidio.model.main)
anova(Basidio.null,Basidio.model.full,Basidio.model.main)
#null model and main essentially same AIC
anova(Basidio.model.main)
anova(Basidio.model.full)
#SoilFrac P=0.009, Crop P=0.07
difflsmeans(Basidio.model.main, ddf="Satterthwaite",type=3,method.grad="simple")
#P<PF=CC, LM<WS, micro<WS, MM<WS
#lower presence in Prairie ->responding to fertilizer?
#Enriched in whole soil compaired to aggregate fractions, except SM

#Family Inocybaceae
Basidio.inocy<-subset(Basidio.data, Basidio.data$Family=="f__Inocybaceae")
head(Basidio.inocy)
#This family is only present in P46micro_July2012 and PF15micro_July2012, abundance 3 and 6, probably not really contributing, certianly not a pattern
#Inocybaceae is one of the most widely distributed and specious families within Agaricales, some species form extomycorrihzal associations
#For further information see: http://inocybaceae.org/

#Family Strophariaceae
Basidio.stroph<-subset(Basidio.data, Basidio.data$Family=="f__Strophariaceae")
head(Basidio.stroph)
Basidio.null<-lmer(value~1+(1|Block), data=Basidio.stroph, REML=FALSE)
Basidio.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Basidio.stroph, REML=FALSE)
Basidio.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Basidio.stroph, REML=FALSE)
AICtab(Basidio.null,Basidio.model.full,Basidio.model.main)
anova(Basidio.null,Basidio.model.full,Basidio.model.main)
#null model best fit by 2, main effects next best fit
anova(Basidio.model.main)
anova(Basidio.model.full)
#Full model detects significant Date*Crop*SoilFrac interaction P=0.002, 
#main effects, crop highly significant, P=0.003
Basidio.strophFig<-ddply(Basidio.stroph, .(Crop, Species), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Basidio.strophFig)
ggplot(Basidio.strophFig)+geom_pointrange(aes(x=Species,y=mean,ymax=high95,ymin=low95, color=Crop),position=position_dodge(width=1))+coord_flip()+scale_y_log10()
#Only 2 species identified, unidentified species most abundant, more abundant in P, PF than CC, very high abundance in P
#Might be worth mentioning due to high abundance of unidentified species

#Family Marasmiaceae
Basidio.mara<-subset(Basidio.data, Basidio.data$Family=="f__Marasmiaceae")
head(Basidio.mara)
Basidio.null<-lmer(value~1+(1|Block), data=Basidio.mara, REML=FALSE)
Basidio.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Basidio.mara, REML=FALSE)
Basidio.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Basidio.mara, REML=FALSE)
AICtab(Basidio.null,Basidio.model.full,Basidio.model.main)
#null model best fit by 8.3, main effects next best fit
anova(Basidio.model.main)
anova(Basidio.model.full)
#nothing significant

#Family Niaceae
Basidio.niac<-subset(Basidio.data, Basidio.data$Family=="f__Niaceae")
head(Basidio.niac)
Basidio.corti$Block<-as.factor(Basidio.corti$Block)
Basidio.null<-lmer(value~1+(1|Block), data=Basidio.niac, REML=FALSE)
Basidio.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Basidio.niac, REML=FALSE)
Basidio.model.main<-lmer(value~Dat+Crop+SoilFrac+(1|Block), data=Basidio.niac, REML=FALSE)
AICtab(Basidio.null,Basidio.model.full,Basidio.model.main)
#Date and Crop not recognized as factors here!!
str(Basidio.niac)
#This indicates R thinks both Date and Crop contains factors with 2 and 3 levels, respectively


#Order Hymenochaetales has only one representing family:  Schizoporaceae, no identifiedn species, found only in:
#P24.SM.Oct2012, PF15.MM.Oct2012, PF41.MM.Oct2012, P31.WS.Oct2012, and PF41.WS.July2012
#abundace values in 10's and 20's
Basidio.hyme<-subset(Basidio.data, Basidio.data$Order=="o__Hymenochaetales")
Basidio.hyme

#Order Geastrales
Basidio.geast<-subset(Basidio.data, Basidio.data$Order=="o__Geastrales")
Basidio.geast
#Only 1 species in this order found, Sphaerobolus sp.
#Only found int PF15.micro and PF15.SM in July 2012, but abundance is in the 40's, meaning very abundant in these 2 fractions
#That's a bit interesting it only shows up there, but not really enought to merit discussion

#Order Cantharellales, single family Ceratobasidiaceae
Basidio.canth<-subset(Basidio.data, Basidio.data$Order=="o__Cantharellales")
head(Basidio.canth)
str(Basidio.canth)
Basidio.null<-lmer(value~1+(1|Block), data=Basidio.canth, REML=FALSE)
Basidio.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Basidio.canth, REML=FALSE)
Basidio.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Basidio.canth, REML=FALSE)
AICtab(Basidio.null,Basidio.model.full,Basidio.model.main)
anova(Basidio.null,Basidio.model.full,Basidio.model.main)
#Null model best by 10, main effects next best model
anova(Basidio.model.main)
anova(Basidio.model.full)
#No significant effects or interactions, looking at graph
Basidio.canthSF<-ddply(Basidio.canth, .(SoilFrac, Family), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Basidio.canthSF)
ggplot(Basidio.canthSF)+geom_pointrange(aes(x=Family,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()
#Seems like aggregate micro<LM would be significant, block effect?
Basidio.canthC<-ddply(Basidio.canth, .(Crop, Family), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Basidio.canthC)
ggplot(Basidio.canthC)+geom_pointrange(aes(x=Family,y=mean,ymax=high95,ymin=low95, color=Crop),position=position_dodge(width=1))+coord_flip()+scale_y_log10()
#Seems like PF>CC, P intermediary
Basidio.canthD<-ddply(Basidio.canth, .(Date, Family), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Basidio.canthD)
ggplot(Basidio.canthD)+geom_pointrange(aes(x=Family,y=mean,ymax=high95,ymin=low95, color=Date),position=position_dodge(width=1))+coord_flip()+scale_y_log10()
#No date affect
#Only 2 genera, Thanatephorus and unknown, seems like there would be differences in aggregate fraction and crop,
Basidio.canthSF<-ddply(Basidio.canth, .(SoilFrac, Genus), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Basidio.canthSF)
ggplot(Basidio.canthSF)+geom_pointrange(aes(x=Genus,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

Basidio.canthC<-ddply(Basidio.canth, .(Crop, Genus), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Basidio.canthC)
ggplot(Basidio.canthC)+geom_pointrange(aes(x=Genus,y=mean,ymax=high95,ymin=low95, color=Crop),position=position_dodge(width=1))+coord_flip()+scale_y_log10()
#Thanatephorus not present in CC

Basidio.canthD<-ddply(Basidio.canth, .(Date, Genus), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Basidio.canthD)
ggplot(Basidio.canthD)+geom_pointrange(aes(x=Genus,y=mean,ymax=high95,ymin=low95, color=Date),position=position_dodge(width=1))+coord_flip()+scale_y_log10()
#slightly higher in Oct, but no major date difference

#There are differences among Thantephorus, but not unknown species
Basidio.than<-subset(Basidio.data, Basidio.data$Genus=="g__Thanatephorus")
head(Basidio.than)
str(Basidio.than)
Basidio.null<-lmer(value~1+(1|Block), data=Basidio.than, REML=FALSE)
Basidio.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Basidio.than, REML=FALSE)
Basidio.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Basidio.than, REML=FALSE)
AICtab(Basidio.null,Basidio.model.full,Basidio.model.main)
#main effects model best fit by 2.0
anova(Basidio.model.main)
anova(Basidio.model.full)
#Insuficient presence in samples to generate interaction terms (e.g. lots of 0's)
#Main effect of Date (P<0.0001), Crop (P<0.0001), and SoilFrac(P=0.02)
difflsmeans(Basidio.model.main, ddf="Satterthwaite",type=3,method.grad="simple")

#Order Atheliales
Basidio.athe<-subset(Basidio.data, Basidio.data$Order=="o__Atheliales")
head(Basidio.athe)
str(Basidio.athe)
Basidio.athe
#All identified as single species, Athelia bombacina, which is a sapbrob on Norway Spruce litter?
#Only present in P13, P31, and PF23, mostly SM, but also a couple of LM, micro, and 1 WS
#Abundance 6-52