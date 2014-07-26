#Elizabeth Bach
#COBS ITS 2012:  Community + ecological co-variates
#21 July 2014


library(reshape)
library(lme4)
library(lmerTest)
library(bbmle)
library(grid)

#use data.nosing rar generated in "data_wrangling.R"

metadata<-read.csv("COBS_ITS_metadata.csv")
str(metadata)
metadata$ExtC<-as.numeric(levels(metadata$ExtC))[metadata$ExtC]
metadata$MBC<-as.numeric(levels(metadata$MBC))[metadata$MBC]
metadata$MBC.MBN<-as.numeric(levels(metadata$MBC.MBN))[metadata$MBC.MBN]
metadata$TC<-as.numeric(levels(metadata$TC))[metadata$TC]
metadata$TN<-as.numeric(levels(metadata$TN))[metadata$TN]
metadata$CN<-as.numeric(levels(metadata$CN))[metadata$CN]
str(metadata)

#Use envfit from vegan package to fit environmental vectors in ordination
#generate single data frame with original community data + "metadata"
data.metadata<-merge(metadata, data.nosing.rar, by="SampleName")
head(data.metadata[,1:30])
str(data.metadata[,6:26])

#subset to remove environmental variables only measured at WS level
metadata.agg<-cbind(data.metadata[1:6],data.metadata[8],data.metadata[16:20],data.metadata[23:25])
head(metadata.agg)
str(metadata.agg)

data.metadata2<-merge(metadata.agg, data.nosing.rar, by="SampleName")
head(data.metadata2[,1:20])

#Now run mds function on the community data only
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
mds.pa<-metaMDS(decostand(data.metadata2[,-c(1:19)],"pa" ),k=6,autotransform=FALSE, na.rm=TRUE)

envectors1<-envfit(mds.pa, data.metadata2[,7:15], na.rm=TRUE)
head(envectors1)
envectors1
vectors<-data.frame(envectors1$vectors[1:4])
names<-c("water_content","AP","BG","BX","CB","NAG","TC","TN","CN")
vectors2<-subset(data.frame(names,vectors), pvals<0.06) #BG P=0.057, so I decided to include it
vectors2

ggplot.NMDS(mds.pa, (data.metadata2$Crop.x), rainbow(3))+geom_segment(data=vectors2, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

ggplot.NMDS(mds.pa, (data.metadata2$Date.x), rainbow(2))

ggplot.NMDS(mds.pa, (data.metadata2$SoilFrac.x), rainbow(5))+geom_segment(data=vectors2, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)
#SoilFrac does support different communities through p/a data

#Total abundance
mds.ab<-metaMDS(decostand(data.metadata2[,-c(1:19)],"total" ),k=6,autotransform=FALSE, na.rm=TRUE)

envectors2<-envfit(mds.ab, data.metadata2[,7:15], na.rm=TRUE)
envectors2
vectors<-data.frame(envectors2$vectors[1:4])
vectors
names<-c("water_content","AP","BG","BX","CB","NAG","TC","TN","CN")
vectors.ab2<-subset(data.frame(names,vectors), pvals<0.05)
vectors.ab2

ggplot.NMDS(mds.ab, (data.metadata2$Crop.x), rainbow(3))+geom_segment(data=vectors.ab2, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors.ab2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

ggplot.NMDS(mds.ab, (data.metadata2$SoilFrac.x), rainbow(5))+geom_segment(data=vectors.ab2, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors.ab2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)
#note, SoilFrac does not support different communities by abundance measure
ggplot.NMDS(mds.ab, (data.metadata2$Date.x), rainbow(2))

#To Do:
#Look at mixed models to see which env. variables are changing with crop, date, soil frac
#all have crop effect
#emphasize ones without treatment responses, as these may be good drivers of communities (e.g. driving communities beyond treatment effects)
#Perhaps cropping system differences will be better with WS, and use SoilFrac effects on aggregate-specific measures
#Bullet point results

#Aggregate Fraction (includes whole soil measures)
#Presence/Absence
mds.pa<-metaMDS(decostand(data.metadata[,-c(1:29)],"pa" ),k=6,autotransform=FALSE)

envectors1<-envfit(mds.pa, data.metadata[,7:25], na.rm=TRUE)
envectors1
vectors<-data.frame(envectors1$vectors[1:4])
names<-c("ph","water_content","MBN","ExtC","ExtN","BD","TP","AMFcol","MBC","AP","BG","BX","CB","NAG","RootBiomass","MBC.MBN","TC","TN","CN")
vectors2<-subset(data.frame(names,vectors), pvals<0.05)
vectors2


ggplot.NMDS(mds.pa, (data.metadata$SoilFrac.x), rainbow(5))+geom_segment(data=vectors2, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

#Abundance
mds.ab<-metaMDS(decostand(data.metadata[,-c(1:29)],"total" ),k=6,autotransform=FALSE, na.rm=TRUE)

envectors2<-envfit(mds.ab, data.metadata[,7:25], na.rm=TRUE)
envectors2
vectors.ab<-data.frame(envectors2$vectors[1:4])
vectors.ab2<-subset(data.frame(names,vectors.ab), pvals<0.05)
vectors.ab2

ggplot.NMDS(mds.ab, (data.metadata$SoilFrac.x), rainbow(5))+geom_segment(data=vectors.ab2, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors.ab2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)
