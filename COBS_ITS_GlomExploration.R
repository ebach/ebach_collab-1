#Elizabeth Bach
#COBS ITS:  Further Taxa exploration:  SoilFrac
#26 July 2014

#Looking at non-rarified Glomeromycota
head(data.nosing.reads[,1:10])
data_melt<-melt(data.nosing.reads, id=c("SampleName","Date","Block","Crop","SoilFrac"))
head(data_melt)
taxonomy<-read.csv("Ebach_ITS_COBS_taxonomy.csv")
head(taxonomy)
head(data_melt)
data_taxa<-merge(data_melt,taxonomy,by.x="variable",by.y="X.OTU.ID")
head(data_taxa)
Glom<-subset(data_taxa, data_taxa$Phylum=="p__Glomeromycota")
Glom$Phylum<-factor(Glom$Phylum)
Glom$Order<-factor(Glom$Order)
Glom$Class<-factor(Glom$Class)
Glom$Family<-factor(Glom$Family)
Glom$Genus<-factor(Glom$Genus)
Glom$Species<-factor(Glom$Species)
str(Glom)
#Remove 0's
Glom2<-Glom[which(Glom$value>0),]
str(Glom2)
Glom.SoilFrac<-ddply(Glom2, .(SoilFrac), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Glom.SoilFrac)
ggplot(Glom.SoilFrac)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))+coord_flip()+scale_y_log10()
#mixed models
Glom.null<-lmer(value~1+(1|Block), data=Glom2, REML=FALSE)
Glom.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Glom2, REML=FALSE)
Glom.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Glom2, REML=FALSE)
AICtab(Glom.null,Glom.model.full,Glom.model.main)
anova(Glom.null,Glom.model.full,Glom.model.main)
#null model lowest AIC by 8, main model next best
anova(Glom.model.full)
#Date*Crop and Crop*SoilFrac interactions
Glom.model<-lmer(value~Date+Crop+SoilFrac+Crop*SoilFrac+Date*Crop+(1|Block), data=Glom2, REML=FALSE)
AICtab(Glom.null,Glom.model.full,Glom.model.main,Glom.model)
#Glom.model is high AIC, only Date*Crop interaction remains significant, data imbalance likely driving this
anova(Glom.model.main)
anova(Glom.model)
list(unique(Glom2$Crop))
list(unique(Glom2$SoilFrac))
list(unique(Glom2$Date))
Glom.micro<-subset(Glom2, Glom2$SoilFrac=="Micro")
Glom.micro

Glom.order<-ddply(Glom, .(SoilFrac, Order), summarise,.progress="text",mean=mean(value),high95=boot.high(value),
low95=boot.low(value)
)
head(Glom.order)
ggplot(Glom.order)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95, group=Order, color=Order),position=position_jitter())+coord_flip()+scale_y_log10()

Glomerales<-subset(Glom2, Glom2$Order=="o__Glomerales")
Diversisporales<-subset(Glom2, Glom2$Order=="o__Diversisporales")
Paraglomerales<-subset(Glom2, Glom2$Order=="o__Paraglomerales")

Glomerales.null<-lmer(value~1+(1|Block), data=Glomerales, REML=FALSE)
Glomerales.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Glomerales, REML=FALSE)
Glomerales.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Glomerales, REML=FALSE)
AICtab(Glomerales.null,Glomerales.model.full,Glomerales.model.main)
anova(Glomerales.null,Glomerales.model.full,Glomerales.model.main)
#null model lowest AIC, but only by 0.9
anova(Glomerales.model.main)
anova(Glomerales.model.full)
#Date*Crop interaction
Glomerales.model<-lmer(value~Date+Crop+SoilFrac+Date*Crop+(1|Block), data=Glomerales, REML=FALSE)
AICtab(Glomerales.null,Glomerales.model.full,Glomerales.model.main,Glomerales.model)
#main model still a better fit, crop and date significant only