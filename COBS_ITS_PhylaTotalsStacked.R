#Elizabeth Bach
#COBS ITS Stacked Bar graphs of phyla
#July 26, 2014

#Phyla-level diversity, summarizes total counts for each phylum
Phyla.data2<-ddply(data_taxa2, .(SampleName, Date, Crop, Block, SoilFrac, Phylum), summarise, .drop=FALSE, .progress="text",total=sum(value))
head(Phyla.data2)
#Phyla.totals<-cast(data_taxa2, SampleName + Date + Crop + Block + SoilFrac ~ Phylum, sum)
#head(Phyla.totals)

#Mean "values"
Phyla.Crop<-ddply(Phyla.data2, .(Crop,Phylum), summarise,.progress="text",mean=mean(total),high95=boot.high(total),
low95=boot.low(total)
)
head(Phyla.Crop)

ggplot(Phyla.Crop, aes(Crop, mean, fill=Phylum))+geom_bar(stat="identity")

#Total"values"
Phyla.CropT<-ddply(Phyla.data2, .(Crop,Phylum), summarise,.progress="text",total=sum(total))
head(Phyla.CropT)
ggplot(Phyla.CropT, aes(Crop, total, fill=Phylum))+geom_bar(stat="identity")

#SoilFrac
#Total"values"
Phyla.SoilFracT<-ddply(Phyla.data2, .(SoilFrac,Phylum), summarise,.progress="text",total=sum(total))
head(Phyla.SoilFracT)
ggplot(Phyla.SoilFracT, aes(SoilFrac, total, fill=Phylum))+geom_bar(stat="identity")
