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
#??For Presence/Absence should be taking out all the 0's??  Doesn't seem to change story, but keep thinking about this.
data_phyla<-data.frame(cast(data_taxa2, SampleName~Phylum, value="value", fun.aggregate=sum, add.missing=TRUE))
head(data_phyla)
head(data.nosing.rar[,1:10])
merged_taxa<-merge(data_phyla, data.nosing.rar, by="SampleName")
dim(merged_taxa)
dim(data.nosing.rar)
names(merged_taxa[,1:30])

Basidio.data<-subset(data_taxa2, data_taxa2$Phylum=="p__Basidiomycota")

#NMDS plotting function, from R. Williams
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

#Adding environmental vectors to NMDS
#use data.metadata2 from COBS_ITS_ComMetadataAnalysis.R
#Presence/Absence 
envectors1<-envfit(mds.pa, data.metadata2[,7:15], na.rm=TRUE)
head(envectors1)
envectors1
vectors<-data.frame(envectors1$vectors[1:4])
names<-c("water_content","AP","BG","BX","CB","NAG","TC","TN","CN")
vectors2<-subset(data.frame(names,vectors), pvals<0.05)
vectors2

ggplot.NMDS(mds.pa, (taxa.interest$Crop), rainbow(3))+geom_point(data=IntVectors2, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="grey",inherit_aes=FALSE)+
geom_text(data=IntVectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=4)+geom_segment(data=vectors2, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

#BX, TC, TN, water content correlated in same direction as
#Thanatephorus, Psathyrellaceae

ggplot.NMDS(mds.pa, (taxa.interest$SoilFrac), rainbow(5))+geom_point(data=IntVectors2, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="grey",inherit_aes=FALSE)+
geom_text(data=IntVectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=4)+geom_segment(data=vectors2, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

#Correlations between BX, Thanatephorus and Psathyrellaceae
taxa_enzy<-merge(taxa.interest, data.metadata2[1:15], by="SampleName")
head(taxa_enzy)
BX.model<-lmer(log(BX+1)~g__Thanatephorus+(1|Block),data=taxa_enzy,REML=FALSE)
anova(BX.model)
BX.model2<-lmer(log(BX+1)~g__Thanatephorus+Crop+Date+(1|Block),data=taxa_enzy,REML=FALSE)
anova(BX.model2)
#Thanatephorus has high P value (>0.8) in any model

BX.model<-lmer(log(BX+1)~f__Psathyrellaceae+(1|Block),data=taxa_enzy,REML=FALSE)
anova(BX.model)
BX.model2<-lmer(log(BX+1)~f__Psathyrellaceae+Crop+Date+(1|Block),data=taxa_enzy,REML=FALSE)
anova(BX.model2)
#Psathyrellaceae is NS (P~0.24) in any model

BX.model<-lmer(log(BX+1)~TC+(1|Block),data=taxa_enzy,REML=FALSE)
anova(BX.model)
BX.model2<-lmer(log(BX+1)~TC+Crop+Date+(1|Block),data=taxa_enzy,REML=FALSE)
anova(BX.model2)
#Highly significantly correlated with TC, any model

BX.model<-lmer(log(BX+1)~water_content_soil+(1|Block),data=taxa_enzy,REML=FALSE)
anova(BX.model)
BX.model2<-lmer(log(BX+1)~water_content_soil+Crop+Date+(1|Block),data=taxa_enzy,REML=FALSE)
anova(BX.model2)
#Highly correlated with water content in any model as well

#Subset to get rid of 0s in Thant and Psath
Than.data<-subset(taxa_enzy, taxa_enzy$g__Thanatephorus>0)
head(Than.data)
Psath.data<-subset(taxa_enzy, taxa_enzy$f__Psathyrellaceae>0)
head(Psath.data)

BX.model<-lmer(log(BX+1)~g__Thanatephorus+(1|Block),data=Than.data,REML=FALSE)
anova(BX.model)
BX.model2<-lmer(log(BX+1)~g__Thanatephorus+Crop+Date+(1|Block),data=Than.data,REML=FALSE)
anova(BX.model2)
BX.model3<-lmer(log(BX+1)~g__Thanatephorus*Crop*Date+(1|Block),data=Than.data,REML=FALSE)
anova(BX.model3)
BX.model4<-lmer(log(BX+1)~g__Thanatephorus+Crop+Date+Crop*Date+(1|Block),data=Than.data,REML=FALSE)
anova(BX.model4)
AICtab(BX.model,BX.model2,BX.model3,BX.model4)
#model 3 is lowest AIC, although not enough Than presence to do full factorial comparisons
#Than. marginally significant, P=0.07

BX.model<-lmer(log(BX+1)~f__Psathyrellaceae+(1|Block),data=Psath.data,REML=FALSE)
anova(BX.model)
BX.model2<-lmer(log(BX+1)~f__Psathyrellaceae+Crop+Date+(1|Block),data=Psath.data,REML=FALSE)
anova(BX.model2)
BX.model3<-lmer(log(BX+1)~f__Psathyrellaceae*Crop*Date+(1|Block),data=Psath.data,REML=FALSE)
anova(BX.model3)
BX.model4<-lmer(log(BX+1)~f__Psathyrellaceae+Crop+Date+Crop*Date+(1|Block),data=Psath.data,REML=FALSE)
anova(BX.model4)
AICtab(BX.model,BX.model2,BX.model3,BX.model4)
#BX model is lowest AIC, no factor of Psath significant, except Crop*Date interaction in model3, may be limited by unbalanced contrasts

#Abundance
mds.ab<-metaMDS(decostand(data.metadata2[,-c(1:19)],"total" ),k=6,autotransform=FALSE, na.rm=TRUE)

envectors2<-envfit(mds.ab, data.metadata2[,7:15], na.rm=TRUE)
envectors2
vectors<-data.frame(envectors2$vectors[1:4])
vectors
names<-c("water_content","AP","BG","BX","CB","NAG","TC","TN","CN")
vectors.ab2<-subset(data.frame(names,vectors), pvals<0.05)
vectors.ab2

IntVectors1ab<-envfit(mds.ab, taxa.interest[,6:15], na.rm=TRUE)
IntVectors1ab
vectors<-data.frame(IntVectors1ab$vectors[1:4])
vectors
names<-c("Limonomyces","Atheliales","UnkBasidio","Thanatephorus","Psathyrellaceae","Strophariaceae","Peziza","Bionectriaceae","Glomerales","Operculomyces")
IntVectors2ab<-data.frame(names, vectors)

ggplot.NMDS(mds.ab, (data.metadata2$Crop.x), rainbow(3))+geom_segment(data=vectors.ab2, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors.ab2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)+geom_point(data=IntVectors2ab, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="grey",inherit_aes=FALSE)+
geom_text(data=IntVectors2ab,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=4)

#NAG correlated in same direction as Limonomyces, both correlated in Prairie direction
#TC, TN, CN, water_content correlated in same direction as Strophariaceae and Bionectriaceae (Bion. a NS correlation) in PF direction
#Correlation models

NAG.model<-lmer(log(NAG+1)~g__Limonomyces+(1|Block),data=taxa_enzy,REML=FALSE)
anova(NAG.model)
NAG.model2<-lmer(log(NAG+1)~g__Limonomyces+Crop+Date+(1|Block),data=taxa_enzy,REML=FALSE)
anova(NAG.model2)
NAG.model3<-lmer(log(NAG+1)~g__Limonomyces*Crop*Date+(1|Block),data=taxa_enzy,REML=FALSE)
anova(NAG.model3)
NAG.model4<-lmer(log(NAG+1)~g__Limonomyces+Crop+Date+Crop*Date+(1|Block),data=taxa_enzy,REML=FALSE)
anova(NAG.model4)
AICtab(NAG.model,NAG.model2,NAG.model3,NAG.model4)
#model 4 best fit, only slightly better than model 3, in both models Limonomyces NS

#Subset to remove 0s in Limonomyces
Limono.data<-subset(taxa_enzy, taxa_enzy$g__Limonomyces>0)
head(Limono.data)

NAG.model<-lmer(log(NAG+1)~g__Limonomyces+(1|Block),data=Limono.data,REML=FALSE)
anova(NAG.model)
NAG.model2<-lmer(log(NAG+1)~g__Limonomyces+Date+(1|Block),data=Limono.data,REML=FALSE)
anova(NAG.model2)
NAG.model3<-lmer(log(NAG+1)~g__Limonomyces*Date+(1|Block),data=Limono.data,REML=FALSE)
anova(NAG.model3)
AICtab(NAG.model,NAG.model2,NAG.model3)
#Model 3 lowest AIC, Limonomyces P=0.10, can't compare Crop because only appears in P

