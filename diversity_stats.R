# this needs data.nosing.rar from data_wrangling.R

head(data.nosing.rar[,1:10])

#calculating richness, shannons, and evenness

richness<-fisher.alpha(data.nosing.rar[,-c(1:5)],1)
shannons<-diversity(data.nosing.rar[,-c(1:5)])
evenness<-shannons/log(richness)
hist(richness)
div_stats<-data.frame(data.nosing.rar[,1:5],richness,shannons,evenness)


head(div_stats)

#looking at shannons first
ggplot(div_stats)+geom_histogram(aes(shannons^2.7))

summary(test<-aov(richness~Date+Crop+SoilFrac, data=div_stats))
TukeyHSD(test)
?diversity

ggplot(div_stats)+geom_boxplot(aes(x=SoilFrac, y=richness))

head(data.nosing.rar[,1:10])
data_melt<-melt(data.nosing.rar.2, id=c("SampleName","Date","CropBlock","Crop","SoilFrac","trt"))
taxonomy<-read.csv(file.choose())
head(taxonomy)
head(data_melt)
data_taxa<-merge(data_melt,taxonomy,by.x="variable",by.y="X.OTU.ID")
head(data_taxa)

library(plyr)
stats<-ddply(data_taxa, .(SoilFrac,order), summarise,.progress="text",
mean=mean(value),
high95=boot.high(value),
low95=boot.low(value)
)
?geom_pointrange
ggplot(stats)+geom_pointrange(aes(x=phylum,y=mean,ymax=high95,ymin=low95,colour=Crop,shape=Date),position=position_dodge(width=1))+coord_flip()+scale_y_log10()
stats$SoilFrac<-factor(stats$SoilFrac, levels=c("Micro","SM","MM","LM","WS"))
ggplot(stats)+geom_pointrange(aes(x=order,y=mean,ymax=high95,ymin=low95,colour=SoilFrac),position=position_dodge(width=.5))+coord_flip()+scale_y_log10()