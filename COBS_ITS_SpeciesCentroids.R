#Elizabeth Bach
#COBS ITS 2012:  Community + taxonomy to generate species centroids
#23 July 2014

library(reshape)
library(lme4)
library(lmerTest)
library(bbmle)
library(grid)

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


#Presence/Absence
mds.pa<-metaMDS(decostand(merged_taxa[,-c(1:11)],"pa" ),k=6,autotransform=FALSE, na.rm=TRUE)

envectors1<-envfit(mds.pa, merged_taxa[,2:7], na.rm=TRUE)
envectors1
vectors<-data.frame(envectors1$vectors[1:4])
names<-c("V1","Ascomycota","Basidiomycota","Chytridiomycota","Glomeromycota","Zygomycota")
vectors2<-subset(data.frame(names,vectors), pvals<0.05)
vectors2

ggplot.NMDS(mds.pa, (merged_taxa$SoilFrac), rainbow(5))+geom_point(data=vectors2, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

ggplot.NMDS(mds.pa, (merged_taxa$Crop), rainbow(3))+geom_point(data=vectors2, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

#Abundance
mds.ab<-metaMDS(decostand(merged_taxa[,-c(1:11)],"total" ),k=6,autotransform=FALSE, na.rm=TRUE)

envectors2<-envfit(mds.ab, merged_taxa[,2:7], na.rm=TRUE)
envectors2
vectors<-data.frame(envectors2$vectors[1:4])
vectors
names<-c("Basidiomycota")
vectors3<-data.frame(names, (subset(vectors, pvals<0.06))) #Zygomycota, P=0.095
vectors3

ggplot.NMDS(mds.ab, (merged_taxa$Crop), rainbow(3))+geom_point(data=vectors3, aes(x=arrows.NMDS1,y=arrows.NMDS2),arrow=arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors3,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

ggplot.NMDS(mds.ab, (merged_taxa$SoilFrac), rainbow(5))+geom_point(data=vectors3, aes(x=arrows.NMDS1,y=arrows.NMDS2),arrow=arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors3,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

#Diving into Phyla to see relationships with lower taxanomic levels
#Basidiomycota
Basidio.data<-subset(data_taxa2, data_taxa2$Phylum=="p__Basidiomycota")
head(Basidio.data)
Basidio.data$Order<-factor(Basidio.data$Order)
Basidio.data$Family<-factor(Basidio.data$Family)
Basidio.data$Class<-factor(Basidio.data$Class)
Basidio.data$Genus<-factor(Basidio.data$Genus)
Basidio.data$Species<-factor(Basidio.data$Species)
subset(Basidio.data, Basidio.data$Order=="o__unidentified")
subset(Basidio.data, Basidio.data$Order=="")
#By class:  Presence/absence
class_Basidio<-data.frame(cast(Basidio.data, SampleName~Class, value="value", fun.aggregate=sum, add.missing=TRUE))
head(class_Basidio)
Basidio_class<-merge(class_Basidio, data.nosing.rar, by="SampleName")
dim(Basidio_class)
dim(data.nosing.rar)
names(Basidio_class[,1:10])
#usde full dataset MDS, only fit new vectors!

Classvectors1<-envfit(mds.pa, Basidio_class[,2:5], na.rm=TRUE)
Classvectors1
vectors<-data.frame(Classvectors1$vectors[1:4])
names<-c("Agaricomycetes","unidentified")
vectors2<-data.frame(names, subset(vectors, pvals<0.05))
vectors2

ggplot.NMDS(mds.pa, (Basidio_class$SoilFrac), rainbow(5))+geom_point(data=vectors2, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

ggplot.NMDS(mds.pa, (Basidio_class$Crop), rainbow(3))+geom_point(data=vectors2, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

#Agaricomycetes
Agarico_data<-subset(Basidio.data, Basidio.data$Class=="c__Agaricomycetes")
head(Agarico_data)
Agarico_data$Order<-factor(Agarico_data$Order)
list(unique(Agarico_data$Order))
class_Agarico<-data.frame(cast(Agarico_data, SampleName~Order, value="value", fun.aggregate=sum, add.missing=TRUE))
head(class_Agarico)
Agarico_class<-merge(class_Agarico, data.nosing.rar, by="SampleName")
dim(Agarico_class)
dim(data.nosing.rar)
names(Agarico_class[,1:15])

Classvectors1<-envfit(mds.pa, Agarico_class[,2:11], na.rm=TRUE)
Classvectors1
vectors<-data.frame(Classvectors1$vectors[1:4])
names<-c("unidentified","Corticiales")
vectors2<-data.frame(names, subset(vectors, pvals<0.06))
vectors2

ggplot.NMDS(mds.pa, (Agarico_class$SoilFrac), rainbow(5))+geom_point(data=vectors2, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

ggplot.NMDS(mds.pa, (Agarico_class$Crop), rainbow(3))+geom_point(data=vectors2, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

#Order Corticiales, has only 1 genus, Limnomyces, and unidentifyable species within Limnomyces
Corti.data<-subset(Agarico_data, Agarico_data$Order=="o__Corticiales")
head(Corti.data)
Corti.data$Genus<-factor(Corti.data$Genus)
Corti.data$Species<-factor(Corti.data$Species)
list(unique(Corti.data$Genus))
str(Corti.data)
genera_Corti<-data.frame(cast(Corti.data, SampleName~Genus, value="value", fun.aggregate=sum, add.missing=TRUE))
head(genera_Corti)
Corti_genera<-merge(genera_Corti, data.nosing.rar, by="SampleName")
dim(Corti_genera)
dim(data.nosing.rar)
names(Corti_genera[,1:10])


Genvectors1<-envfit(mds.pa, Corti_genera[,2], na.rm=TRUE)
Genvectors1
vectors<-data.frame(Genvectors1$vectors[1:4])
vectors
names<-c("Limonomyces")
vectors2<-data.frame(names, subset(vectors, pvals<0.06))
vectors2

ggplot.NMDS(mds.pa, (Corti_genera$SoilFrac), rainbow(5))+geom_point(data=vectors2, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

ggplot.NMDS(mds.pa, (Corti_genera$Crop), rainbow(3))+geom_point(data=vectors2, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

#Unidentified Agaricomycetes..any way to subset or find particular association?
Ukn.agari<-subset(Agarico_data, Agarico_data$Order=="")
head(Ukn.agari)
list(unique(Ukn.agari$Genus))
No resolution below Agaricomyctes

#Abundance:  Basidiomycota only phyla significantly corrrelated with a community shift, directionally between CC and PF
list(unique(Basidio.data$Class))
subset(Basidio.data, Basidio.data$Class=="")
subset(Basidio.data, Basidio.data$Class=="c__unidentified") #only 4 sequences in this category, 3 from P46, 1 from PF32

Classvectors1ab<-envfit(mds.ab, Basidio_class[,2:5], na.rm=TRUE)
Classvectors1ab
vectorsab<-data.frame(Classvectors1ab$vectors[1:4])
names<-c("Unknown","unidentified")
vectors2ab<-data.frame(names, subset(vectorsab, pvals<0.05))
vectors2ab

ggplot.NMDS(mds.ab, (Basidio_class$SoilFrac), rainbow(5))+geom_point(data=vectors2ab, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors2ab,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

ggplot.NMDS(mds.ab, (Basidio_class$Crop), rainbow(3))+geom_point(data=vectors2ab, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors2ab,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

#Unknown Phyla: Presence/Absence
Unk.data<-subset(data_taxa2, data_taxa2$Phylum=="")
Unk.data$Phylum<-factor(Unk.data$Phylum)
head(Unk.data)
list(unique(Unk.data$Kingdom))
length(Unk.data)
#all positively id-ed as Fungi, no further classification, highly present in Prairie

class_unk<-data.frame(cast(Unk.data, SampleName~Phylum, value="value", fun.aggregate=sum, add.missing=TRUE))
head(class_unk)
str(class_unk)

#Chytridiomycota
Chytri.data<-subset(data_taxa2, data_taxa2$Phylum=="p__Chytridiomycota")
Chytri.data$Phylum<-factor(Chytri.data$Phylum)
Chytri.data$Class<-factor(Chytri.data$Class)
Chytri.data$Order<-factor(Chytri.data$Order)
Chytri.data$Family<-factor(Chytri.data$Family)
Chytri.data$Genus<-factor(Chytri.data$Genus)
Chytri.data$Species<-factor(Chytri.data$Species)
head(Chytri.data)
list(unique(Chytri.data$Class))   #Only 1 class
list(unique(Chytri.data$Order))   #6 orders, including "" and "unidentified"
list(unique(Chytri.data$Family))  #5 faimlies
list(unique(Chytri.data$Genus))   #6 genera
list(unique(Chytri.data$Species)) #5 species

#look at orders first and work down, mostly likely only 1 identified family/genus per order
order_Chytri<-data.frame(cast(Chytri.data, SampleName~Order, value="value", fun.aggregate=sum, add.missing=TRUE))
head(order_Chytri)
Chytri_order<-merge(order_Chytri, data.nosing.rar, by="SampleName")
dim(Chytri_order)
dim(data.nosing.rar)
names(Chytri_order[,1:12])
#use full dataset MDS, only fit new vectors!

ChytriOrders1<-envfit(mds.pa, Chytri_order[,2:7], na.rm=TRUE)
ChytriOrders1
ChytriVectors<-data.frame(ChytriOrders1$vectors[1:4])
names<-c("unidentified")
ChytriVectors2<-data.frame(names, subset(ChytriVectors, pvals<0.05))
ChytriVectors2

ggplot.NMDS(mds.pa, (Chytri_order$SoilFrac), rainbow(5))+geom_point(data=ChytriVectors2, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="grey",inherit_aes=FALSE)+
geom_text(data=ChytriVectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

ggplot.NMDS(mds.pa, (Chytri_order$Crop), rainbow(3))+geom_point(data=ChytriVectors2, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="grey",inherit_aes=FALSE)+
geom_text(data=ChytriVectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

#unidentified Chytri is only siginifcant correlation, but only 2 samples, most sequences from Rhizophydiales
#Rhyzophydiales sequences only found in P & PF samples, low abundance, all from genus Operculomyces

#Glomeromycota
Glom.data<-subset(data_taxa2, data_taxa2$Phylum=="p__Glomeromycota")
Glom.data$Phylum<-factor(Glom.data$Phylum)
Glom.data$Class<-factor(Glom.data$Class)
Glom.data$Order<-factor(Glom.data$Order)
Glom.data$Family<-factor(Glom.data$Family)
Glom.data$Genus<-factor(Glom.data$Genus)
Glom.data$Species<-factor(Glom.data$Species)
head(Glom.data)
list(unique(Glom.data$Class))   #Only 1 class
list(unique(Glom.data$Order))   #3 orders
list(unique(Glom.data$Family))  #3 faimlies
list(unique(Glom.data$Genus))   #4 genera, including ""
list(unique(Glom.data$Species)) #3 species, including ""

#look at orders first and work down, mostly likely only 1 identified family/genus per order
order_Glom<-data.frame(cast(Glom.data, SampleName~Order, value="value", fun.aggregate=sum, add.missing=TRUE))
head(order_Glom)
Glom_order<-merge(order_Glom, data.nosing.rar, by="SampleName")
dim(Glom_order)
dim(data.nosing.rar)
names(Glom_order[,1:10])
#use full dataset MDS, only fit new vectors!

GlomOrders1<-envfit(mds.pa, Glom_order[,2:4], na.rm=TRUE)
GlomOrders1
GlomVectors<-data.frame(GlomOrders1$vectors[1:4])
names<-c("Glomerales","Paraglomerales")
GlomVectors2<-data.frame(names, subset(GlomVectors, pvals<0.05))
GlomVectors2

ggplot.NMDS(mds.pa, (Glom_order$SoilFrac), rainbow(5))+geom_point(data=GlomVectors2, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="grey",inherit_aes=FALSE)+
geom_text(data=GlomVectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

ggplot.NMDS(mds.pa, (Glom_order$Crop), rainbow(3))+geom_point(data=GlomVectors2, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="grey",inherit_aes=FALSE)+
geom_text(data=GlomVectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

#How many observatoins in Paraglomeromycota?  Only 2, both in P46
subset(Glom.data, Glom.data$Order=="o__Paraglomerales")

