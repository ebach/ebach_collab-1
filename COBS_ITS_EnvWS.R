#Elizabeth Bach
#COBS ITS 2012 data:  Environmental varibales, subset for whole soil only-measures
#22 July 2014

library(vegan)
library(ggplot2)
library(grid)

#Subset "data.metadata" from COBS_ITS_ComMetadataAnalysis.R for whole soil only

metadata.WS<-subset(data.metadata, SoilFrac.x=="WS")
#remove NAs
metadata.WS2<-subset(metadata.WS, !is.na(metadata.WS[,7:25]))
head(metadata.WS2[,1:26])
#metadata.WS2$TC
str(metadata.WS2)

#Now run community NMDS analysis on WS only
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
mds.pa<-metaMDS(decostand(metadata.WS[,-c(1:29)],"pa" ),k=6,autotransform=FALSE, na.rm=TRUE)

envectors1<-envfit(mds.pa, metadata.WS[,7:25], na.rm=TRUE)
head(envectors1)
envectors1
vectors<-data.frame(envectors1$vectors[1:4])
names<-c("ph","water_content","MBN","ExtC","ExtN","BD","TP","AMFcol","MBC","AP","BG","BX","CB","NAG","RootBiomass","MBC.MBN","TC","TN","CN")
vectors2<-subset(data.frame(names,vectors), pvals<0.05)
vectors2

ggplot.NMDS(mds.pa, (metadata.WS$Crop.x), rainbow(3))+geom_segment(data=vectors2, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

#Total abundance
mds.ab<-metaMDS(decostand(metadata.WS[,-c(1:29)],"total" ),k=6,autotransform=FALSE, na.rm=TRUE)

envectors2<-envfit(mds.ab, metadata.WS[,7:25], na.rm=TRUE)
envectors2
vectors.ab<-data.frame(envectors2$vectors[1:4])
vectors.ab2<-subset(data.frame(names,vectors.ab), pvals<0.05)
vectors.ab2

ggplot.NMDS(mds.ab, (metadata.WS$Crop.x), rainbow(3))+geom_segment(data=vectors.ab2, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors.ab2,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=5)

