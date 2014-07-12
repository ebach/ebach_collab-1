#Elizabeth Bach
#COBS Aggregates Fungal ITS analysis
#taxanomic-level analysis, looking among aggregates, date and crop
#11 July 2014

#Use data_taxa2 from diversity_stats.R and Basidio.data, Asco.data, Glom.data, Zygo.data, and Chytri.data from COBS_ITS_PhylaSubset
#reshape library contains colsplit function used to generate crop, SoilFrac, and Date in summerized files

library(reshape)
#Phyla-level diversity, summarizes total counts for each phylum
#Phyla.data2<-ddply(data_taxa2, .(SampleName, Date, Crop, SoilFrac, Phylum), summarise, .drop=FALSE, .progress="text",total=sum(value))
#head(Phyla.data2)
Phyla.totals<-cast(data_taxa2, SampleName + Date + Crop + SoilFrac ~ Phylum, sum)
head(Phyla.totals)

richness<-fisher.alpha(Phyla.totals[,-c(1:5)])
head(richness)
#Calculates Fisher's "alpha" value for Fisher's logrithmic series, not true total count
hist(richness)
shannons<-diversity(Phyla.totals[,-c(1:5)])
head(shannons)
hist(shannons)
evenness<-shannons/log(richness)
hist(evenness)
Phyla_stats<-data.frame(Phyla.totals[,1:5],richness,shannons,evenness)
head(Phyla_stats)

#Main effects
summary(test<-aov(shannons~Date+Crop+SoilFrac, data=Phyla.totals))
TukeyHSD(test)
#Crop effect significant
summary(test<-aov(richness~Date+Crop+SoilFrac, data=Phyla.totals))
TukeyHSD(test)
#SoilFrac effect significant
summary(test<-aov(evenness~Date+Crop+SoilFrac, data=Phyla.totals))
TukeyHSD(test)
#NS

#Full model
#Figure out how to get block column back to run mixed models!
summary(test<-aov(shannons~Date*Crop*SoilFrac, data=Phyla.totals))
TukeyHSD(test)

#need count data file
Pdiv_stats<-data.frame(Phyla.data[,1:2],shannons)
hist(log(shannons+1))
#data skewed right, log helps
summary(test<-aov(shannons~SoilFrac, data=Pdiv_stats))

evenness<-shannons/log(richness)
hist(richness)
div_stats<-data.frame(data.nosing.rar[,1:5],richness,shannons,evenness)

