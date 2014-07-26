#Elizabeth Bach
#COBS ITS Presence/Absence community NMDS + target taxa
#25 July 2014

library(reshape)
library(grid)
library(ggplot2)

#use data.nosing rar generated in "data_wrangling.R"
#taxonomy code from diversity_stats.R
head(data.nosing.rar[,1:10])
data_melt<-melt(data.nosing.rar, id=c("SampleName","Date","Block","Crop","SoilFrac"))
head(data_melt)

#taxonomy<-read.csv(file.choose())
taxonomy<-read.csv("Ebach_ITS_COBS_taxonomy.csv")
head(taxonomy)
head(data_melt)
data_taxa<-merge(data_melt,taxonomy,by.x="variable",by.y="X.OTU.ID")
head(data_taxa)
data_taxa2<-data_taxa[ which(data_taxa$value>0),]
head(data_taxa2)
data_phyla<-data.frame(cast(data_taxa2, SampleName~Phylum, value="value", fun.aggregate=sum, add.missing=TRUE))
head(data_phyla)
head(data.nosing.rar[,1:10])
merged_taxa<-merge(data_phyla, data.nosing.rar, by="SampleName")
dim(merged_taxa)
dim(data.nosing.rar)
names(merged_taxa[,1:30])

Basidio.data<-subset(data_taxa2, data_taxa2$Phylum=="p__Basidiomycota")


ggplot.NMDS<-function(XX,ZZ,COLORS){
	library(ggplot2)
MDS1<-data.frame(scores(XX))$NMDS1
MDS2<-data.frame(scores(XX))$NMDS2
Treatment<-ZZ

NMDS<-data.frame(MDS1,MDS2,Treatment)

NMDS.mean=aggregate(NMDS[,1:2],list(group=Treatment),mean)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }

  df_ell <- data.frame()
  for(g in levels(NMDS$Treatment)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Treatment==g,],
                    veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                    ,group=g))
  }

X1<-ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = Treatment),size=3,alpha=0.75) +
    geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=2, linetype=5)+theme_bw()+theme(aspect.ratio=1)+scale_color_manual(values=COLORS)+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+theme(legend.title=element_text(size=15),legend.text=element_text(size=15))
X1    
}

mds.pa<-metaMDS(decostand(merged_taxa[,-c(1:11)],"pa" ),k=6,autotransform=FALSE, na.rm=TRUE)

#Figure summarizing target taxa in presence/absence NMDS
Limonomyces<-subset(data_taxa2, data_taxa2$Genus=="g__Limonomyces")
head(Limonomyces)
Limonomyces.value<-data.frame(cast(Limonomyces, SampleName~Genus, value="value", fun.aggregate=sum, add.mising=TRUE))
head(Limonomyces.value)
Atheliales<-subset(data_taxa2, data_taxa2$Order=="o__Atheliales")
head(Atheliales)
Atheliales.value<-data.frame(cast(Atheliales, SampleName~Order, value="value", fun.aggregate=sum, add.mising=TRUE))
head(Atheliales.value)
UnkBasidio<-subset(Basidio.data, Basidio.data$Order=="") #Note see Basidio code above for "Basidio.data"
UnkBasidio.value<-data.frame(cast(UnkBasidio, SampleName~Order, value="value", fun.aggregate=sum, add.mising=TRUE))
head(UnkBasidio.value)
Thanatephorus<-subset(data_taxa2, data_taxa2$Genus=="g__Thanatephorus")
Thanatephorus.value<-data.frame(cast(Thanatephorus, SampleName~Genus, value="value", fun.aggregate=sum, add.mising=TRUE))
head(Thanatephorus.value)
Psathyrellaceae<-subset(data_taxa2, data_taxa2$Family=="f__Psathyrellaceae")
Psathyrellaceae.value<-data.frame(cast(Psathyrellaceae, SampleName~Family, value="value", fun.aggregate=sum, add.mising=TRUE))
head(Psathyrellaceae.value)
Strophariaceae<-subset(data_taxa2, data_taxa2$Family=="f__Strophariaceae")
Strophariaceae.value<-data.frame(cast(Strophariaceae, SampleName~Family, value="value", fun.aggregate=sum, add.mising=TRUE))
head(Strophariaceae.value)
Peziza<-subset(data_taxa2, data_taxa2$Genus=="g__Peziza")
head(Peziza)
Peziza.value<-data.frame(cast(Peziza, SampleName~Genus, value="value", fun.aggregate=sum, add.mising=TRUE))
head(Peziza.value)
Bionectriaceae<-subset(data_taxa2, data_taxa2$Family=="f__Bionectriaceae")
Bionectriaceae.value<-data.frame(cast(Bionectriaceae, SampleName~Family, value="value", fun.aggregate=sum, add.mising=TRUE))
head(Bionectriaceae.value)
Glomerales<-subset(data_taxa2, data_taxa2$Order=="o__Glomerales")
head(Glomerales)
Glomerales.value<-data.frame(cast(Glomerales, SampleName~Order, value="value", fun.aggregate=sum, add.mising=TRUE))
head(Glomerales.value)
Operculomyces<-subset(data_taxa2, data_taxa2$Genus=="g__Operculomyces")
Operculomyces.value<-data.frame(cast(Operculomyces, SampleName~Genus, value="value", fun.aggregate=sum, add.mising=TRUE))
head(Operculomyces.value)
SampleInfo<-data.nosing.rar[1:5]
taxa.interest1<-merge(SampleInfo, Limonomyces.value, by="SampleName", all=TRUE)
taxa.interest2<-merge(taxa.interest1, Atheliales.value, by="SampleName", all=TRUE)
taxa.interest3<-merge(taxa.interest2, UnkBasidio.value, by="SampleName", all=TRUE)
taxa.interest4<-merge(taxa.interest3, Thanatephorus.value, by="SampleName", all=TRUE)
taxa.interest5<-merge(taxa.interest4, Psathyrellaceae.value, by="SampleName", all=TRUE)
taxa.interest6<-merge(taxa.interest5, Strophariaceae.value, by="SampleName", all=TRUE)
taxa.interest7<-merge(taxa.interest6, Peziza.value, by="SampleName", all=TRUE)
taxa.interest8<-merge(taxa.interest7, Bionectriaceae.value, by="SampleName", all=TRUE)
taxa.interest9<-merge(taxa.interest8, Glomerales.value, by="SampleName", all=TRUE)
taxa.interest<-merge(taxa.interest9, Operculomyces.value, by="SampleName", all=TRUE)
taxa.interest[is.na(taxa.interest)]<-0
head(taxa.interest)
str(taxa.interest)
dim(taxa.interest)
dim(data.nosing.rar)

IntVectors1<-envfit(mds.pa, taxa.interest[,6:15], na.rm=TRUE)
IntVectors1
vectors<-data.frame(IntVectors1$vectors[1:4])
vectors
names<-c("Limonomyces","Atheliales","UnkBasidio","Thanatephorus","Psathyrellaceae","Strophariaceae","Peziza","Bionectriaceae","Glomerales","Operculomyces")
IntVectors2<-data.frame(names, vectors)

ggplot.NMDS(mds.pa, (taxa.interest$SoilFrac), rainbow(5))+geom_point(data=IntVectors2, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="grey",inherit_aes=FALSE)+
geom_text(data=IntVectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=4)

ggplot.NMDS(mds.pa, (taxa.interest$Crop), rainbow(3))+geom_point(data=IntVectors2, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="grey",inherit_aes=FALSE)+
geom_text(data=IntVectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=4)


