#Elizabeth Bach, Ryan Williams
#COBS Aggregate ITS data
#Subsetting by Phyla and looking within
#10 July 2014

library(plyr)
library(ggplot2)

#Uses "data_taxa" file generated in "diversity_stats.R" file

head(data_taxa)
data_taxa2<-data_taxa[ which(data_taxa$value>0),]
head(data_taxa2)

Zygo.data<-subset(data_taxa2, data_taxa2$Phylum=="p__Zygomycota")
head(Zygo.data)

Glom.data<-subset(data_taxa2, data_taxa2$Phylum=="p__Glomeromycota")
head(Glom.data)

Asco.data<-subset(data_taxa2, data_taxa2$Phylum=="p__Ascomycota")
head(Asco.data)

Basidio.data<-subset(data_taxa2, data_taxa2$Phylum=="p__Basidiomycota")
head(Basidio.data)

Chytri.data<-subset(data_taxa2, data_taxa2$Phylum=="p__Chytridiomycota")
head(Chytri.data)

Unk.data<-subset(data_taxa2, data_taxa2$Phylum=="")
head(Unk.data)

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


#Basidiomycota
#order
#SoilFrac
Basidio.order<-ddply(Basidio.data, .(SoilFrac, Order), summarise,.progress="text",
mean=mean(value),
high95=boot.high(value),
low95=boot.low(value)
)
head(Basidio.order)

ggplot(Basidio.order)+geom_pointrange(aes(x=Order,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

#Crop
Basidio.order<-ddply(Basidio.data, .(Crop, Order), summarise,.progress="text",
mean=mean(value),
high95=boot.high(value),
low95=boot.low(value)
)
head(Basidio.order)

ggplot(Basidio.order)+geom_pointrange(aes(x=Order,y=mean,ymax=high95,ymin=low95, color=Crop),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

#Date
#Crop
Basidio.order<-ddply(Basidio.data, .(Date, Order), summarise,.progress="text",
mean=mean(value),
high95=boot.high(value),
low95=boot.low(value)
)
head(Basidio.order)

ggplot(Basidio.order)+geom_pointrange(aes(x=Order,y=mean,ymax=high95,ymin=low95, color=Date),position=position_dodge(width=1))+coord_flip()+scale_y_log10()


#Family
Basidio.family<-ddply(Basidio.data, .(SoilFrac, Family), summarise,.progress="text",
mean=mean(value),
high95=boot.high(value),
low95=boot.low(value)
)
head(Basidio.family)

ggplot(Basidio.family)+geom_pointrange(aes(x=Family,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

#Species
Basidio.species<-ddply(Basidio.data, .(SoilFrac, Species), summarise,.progress="text",
mean=mean(value),
high95=boot.high(value),
low95=boot.low(value)
)
head(Basidio.species)

ggplot(Basidio.species)+geom_pointrange(aes(x=Species,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

#Ascomycota
#order
Asco.order<-ddply(Asco.data, .(SoilFrac, Order), summarise,.progress="text",
mean=mean(value),
high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.order)

ggplot(Asco.order)+geom_pointrange(aes(x=Order,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

#Family
Asco.family<-ddply(Asco.data, .(SoilFrac, Family), summarise,.progress="text",
mean=mean(value),
high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.family)

ggplot(Asco.family)+geom_pointrange(aes(x=Family,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

#Species
Asco.species<-ddply(Asco.data, .(SoilFrac, Species), summarise,.progress="text",
mean=mean(value),
high95=boot.high(value),
low95=boot.low(value)
)
head(Asco.species)

ggplot(Asco.species)+geom_pointrange(aes(x=Species,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

#Glomeromycota
#Order
Glom.order<-ddply(Glom.data, .(SoilFrac, Order), summarise,.progress="text",
mean=mean(value),
high95=boot.high(value),
low95=boot.low(value)
)
head(Glom.order)

ggplot(Glom.order)+geom_pointrange(aes(x=Order,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

#Family
Glom.family<-ddply(Glom.data, .(SoilFrac, Family), summarise,.progress="text",
mean=mean(value),
high95=boot.high(value),
low95=boot.low(value)
)
head(Glom.family)

ggplot(Glom.family)+geom_pointrange(aes(x=Family,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

#Species
Glom.species<-ddply(Glom.data, .(SoilFrac, Species), summarise,.progress="text",
mean=mean(value),
high95=boot.high(value),
low95=boot.low(value)
)
head(Glom.species)

ggplot(Glom.species)+geom_pointrange(aes(x=Species,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

#Zygomycota
#order
Zygo.order<-ddply(Zygo.data, .(SoilFrac, Order), summarise,.progress="text",
mean=mean(value),
high95=boot.high(value),
low95=boot.low(value)
)
head(Zygo.order)

ggplot(Zygo.order)+geom_pointrange(aes(x=Order,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

#Family
Zygo.family<-ddply(Zygo.data, .(SoilFrac, Family), summarise,.progress="text",
mean=mean(value),
high95=boot.high(value),
low95=boot.low(value)
)
head(Zygo.family)

ggplot(Zygo.family)+geom_pointrange(aes(x=Family,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

#Species
Zygo.species<-ddply(Zygo.data, .(SoilFrac, Species), summarise,.progress="text",
mean=mean(value),
high95=boot.high(value),
low95=boot.low(value)
)
head(Zygo.species)

ggplot(Zygo.species)+geom_pointrange(aes(x=Species,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

#Chytridiomycota
#order
Chytri.order<-ddply(Chytri.data, .(SoilFrac, Order), summarise,.progress="text",
mean=mean(value),
high95=boot.high(value),
low95=boot.low(value)
)
head(Chytri.order)

ggplot(Chytri.order)+geom_pointrange(aes(x=Order,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

#Family
Chytri.family<-ddply(Chytri.data, .(SoilFrac, Family), summarise,.progress="text",
mean=mean(value),
high95=boot.high(value),
low95=boot.low(value)
)
head(Chytri.family)

ggplot(Chytri.family)+geom_pointrange(aes(x=Family,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()

#Species
Chytri.species<-ddply(Chytri.data, .(SoilFrac, Species), summarise,.progress="text",
mean=mean(value),
high95=boot.high(value),
low95=boot.low(value)
)
head(Chytri.species)

ggplot(Chytri.species)+geom_pointrange(aes(x=Species,y=mean,ymax=high95,ymin=low95, color=SoilFrac),position=position_dodge(width=1))+coord_flip()+scale_y_log10()



