#Elizabeth Bach
#COBS ITS 2012:  Community + ecological co-variates
#21 July 2014


library(reshape)
library(lme4)
library(lmerTest)
library(bbmle)
library(grid)

#use "mds" generated in multivariate_tests.R, and data.nosing rar

metadata<-read.csv("COBS_ITS_metadata.csv")
str(metadata)
metadata$ExtC<-as.numeric(levels(metadata$ExtC))[metadata$ExtC]
metadata$MBC<-as.numeric(levels(metadata$MBC))[metadata$MBC]
metadata$MBC.MBN<-as.numeric(levels(metadata$MBC.MBN))[metadata$MBC.MBN]
metadata$TC<-as.numeric(levels(metadata$TC))[metadata$TC]
metadata$TN<-as.numeric(levels(metadata$TN))[metadata$TN]
metadata$CN<-as.numeric(levels(metadata$CN))[metadata$CN]
str(metadata)

metadata.mds<-merge(mds.points12, metadata, by="SampleName")
head(metadata.mds)
str(metadata.mds)

#Use envfit from vegan package to fit environmental vectors in ordination
#generate single data frame with original community data + "metadata"
data.metadata<-merge(metadata, data.nosing.rar, by="SampleName")
head(data.metadata[,1:30])
str(data.metadata[,6:26])

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
mds.pa<-metaMDS(decostand(data.metadata[,-c(1:29)],"pa" ),k=6,autotransform=FALSE)

envectors1<-envfit(mds.pa, data.metadata[,7:25], na.rm=TRUE)
head(envectors1)
envectors1
vectors<-scores(envectors1, "vectors")
names<-c("ph","water_content","MBN","ExtC","ExtN","BD","TP","AMFcol","MBC","AP","BG","BX","CB","NAG","RootBiomass","MBC.MBN","TC","TN","CN")
vectors2<-data.frame(names,vectors)
str(vectors2)
vectors2
ggplot.NMDS(mds.pa, (data.metadata$Crop.x), rainbow(3))+geom_segment(data=vectors2, aes(x=0,xend=NMDS1,y=0,yend=NMDS2),arrow=arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE)+
geom_text(data=vectors2,aes(x=NMDS1,y=NMDS2,label=names),size=5)
ggplot.NMDS(mds.pa, (data.metadata$Date.x), rainbow(2))

#Total abundance
mds.ab<-metaMDS(decostand(data.metadata[,-c(1:29)],"total" ),k=6,autotransform=FALSE)
ggplot.NMDS(mds.ab, (data.metadata$Crop.x), rainbow(3))

envectors2<-envfit(mds.ab, data.metadata[,7:25], na.rm=TRUE)
head(envectors2)
envectors2
vectors<-scores(envectors2, "vectors")


#To Do:
#pull out NS variables
#Look at mixed models to see which env. variables are changing with crop, date, soil frac
#emphasize ones without treatment responses, as these may be good drivers of communities (e.g. driving communities beyond treatment effects)
#do a WS subset, look at MBC, pH, AMF, etc. with WS only, as I do not have aggregate-specific measures of those
#Perhaps cropping system differences will be better with WS, and use SoilFrac effects on aggregate-specific measures
#Bullet point results
