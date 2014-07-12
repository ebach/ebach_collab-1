# this needs data.nosing.rar from data_wrangling.R

head(data.nosing.rar[,1:10])

#calculating richness, shannons, and evenness

richness<-fisher.alpha(data.nosing.rar[,-c(1:5)],1)
head(richness)
shannons<-diversity(data.nosing.rar[,-c(1:5)])
evenness<-shannons/log(richness)
hist(richness)
div_stats<-data.frame(data.nosing.rar[,1:5],richness,shannons,evenness)


head(div_stats)

#looking at data distribution
ggplot(div_stats)+geom_histogram(aes(shannons^2.7))
ggplot(div_stats)+geom_histogram(aes(richness^2.7))
ggplot(div_stats)+geom_histogram(aes(evenness^2.7))

#richness is skewed heavily to right
ggplot(div_stats)+geom_histogram(aes(richness))
#keeping log transformation helps normalize data

#Testing main effects of date, crop, soil frac on diversity measures
summary(test<-aov(richness~Date+Crop+SoilFrac, data=div_stats))
TukeyHSD(test)
#SoilFrac P=0.001, micro > LM, WS

summary(test<-aov(evenness~Date+Crop+SoilFrac, data=div_stats))
TukeyHSD(test)
#no effects of any main factor

summary(test<-aov(shannons~Date+Crop+SoilFrac, data=div_stats))
TukeyHSD(test)
#SoilFrac P=0.02, micro>LM

?diversity

ggplot(div_stats)+geom_boxplot(aes(x=SoilFrac, y=richness))

head(data.nosing.rar[,1:10])
data_melt<-melt(data.nosing.rar, id=c("SampleName","Date","CropBlock","Crop","SoilFrac"))

taxonomy<-read.csv(file.choose())
head(taxonomy)
head(data_melt)
data_taxa<-merge(data_melt,taxonomy,by.x="variable",by.y="X.OTU.ID")
head(data_taxa)
data_taxa2<-data_taxa[ which(data_taxa$value>0),]
head(data_taxa2)

#Bootstrap functions from R. Williams
boot.high<-function(XX){
boot.mean<-numeric(1000)
for (i in 1:1000){
 boot.mean[i]<-mean(sample(XX,replace=T))
}
return(quantile(boot.mean,(0.975)))
}

boot.low<-function(XX){
boot.mean<-numeric(1000)
for (i in 1:1000){
 boot.mean[i]<-mean(sample(XX,replace=T))
}
return(quantile(boot.mean,(0.025)))
}

library(plyr)
Phyla.data<-ddply(data_taxa2, .(SoilFrac, Phylum), summarise,.progress="text",
mean=mean(value),
high95=boot.high(value),
low95=boot.low(value)
)
head(Phyla.data)

?geom_pointrange
ggplot(Phyla.data)+geom_pointrange(aes(x=Phylum,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()
stats$SoilFrac<-factor(stats$SoilFrac, levels=c("Micro","SM","MM","LM","WS"))
ggplot(stats)+geom_pointrange(aes(x=Phylum,y=mean,ymax=high95,ymin=low95,colour=SoilFrac),position=position_dodge(width=.5))+coord_flip()+scale_y_log10()

