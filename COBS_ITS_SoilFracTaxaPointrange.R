#Elizabeth Bach
#COBS ITS:  Aggregate Target Taxa Figure
#28 July 2014

#Use taxa.interest from COBS_ITS_TargetTaxaNMDS.R
#family Bionectriaceae (Ascomycota), genus Peziza (Ascomycota), family Psathyrellaceae (Basidiomycota), genus Thanatephorus (Basidiomycota)
#show significant effect of soil aggregate fraction

taxa.long<-subset(melt(taxa.interest, id=c("SampleName","Crop","SoilFrac","Date","Block")), value>0)
head(taxa.long)
taxa.SoilFrac<-taxa.long[taxa.long$variable %in% c("f__Bionectriaceae","g__Peziza","f__Psathyrellaceae","g__Thanatephorus"),]
head(taxa.SoilFrac)
taxa.summary<-ddply(taxa.SoilFrac, .(SoilFrac, variable), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(taxa.summary) 

ggplot(taxa.summary)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95))+facet_wrap(~variable, scales="free",ncol=2)+theme_bw()