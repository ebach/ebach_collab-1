# This starts from data.nosing.rar as a result from data_wrangling.R

#first transforming to scaled data within each sample (between 0 and 1) for multivariate tests.  This uses decostand() from the vegan pacakge
head(data.nosing.rar[,1:5])
data.trans.rar<-cbind(data.nosing.rar[,1:5],decostand(data.nosing.rar[,-c(1:5)],"total"))

#using adonis for the test; I pulled out non-significant interactions
adonis(data.trans.rar[,-c(1:5)]~data.trans.rar$Date*data.trans.rar$Crop+data.trans.rar$SoilFrac, permutations=9999)
                                        # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# data.trans.rar$Date                      1     0.767 0.76746  2.2292 0.01970 0.0004 ***
# data.trans.rar$Crop                      2     4.442 2.22111  6.4515 0.11405 0.0001 ***
# data.trans.rar$SoilFrac                  4     1.523 0.38065  1.1057 0.03909 0.1535    
# data.trans.rar$Date:data.trans.rar$Crop  2     1.231 0.61566  1.7883 0.03161 0.0004 ***
# Residuals                               90    30.985 0.34428         0.79554           
# Total                                   99    38.949                 1.00000           
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#now transforming to presence/absence (0 and 1) for multivariate tests.  This uses decostand() from the vegan pacakge
head(data.nosing.rar[,1:5])
data.trans.rar<-cbind(data.nosing.rar[,1:5],decostand(data.nosing.rar[,-c(1:5)],"pa"))

#using adonis for the test; I pulled out non-significant interactions and jaccard for presence/absence
adonis(data.trans.rar[,-c(1:5)]~data.trans.rar$Date*data.trans.rar$Crop+data.trans.rar$SoilFrac, permutations=9999)

# when you rarefy by 1000 you get
                                        # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# data.trans.rar$Date                      1    0.4082 0.40824  1.5080 0.01459 0.0009 ***
# data.trans.rar$Crop                      2    1.3746 0.68732  2.5388 0.04913 0.0001 ***
# data.trans.rar$SoilFrac                  4    1.2217 0.30541  1.1281 0.04366 0.0224 *  
# data.trans.rar$Date:data.trans.rar$Crop  2    0.6126 0.30628  1.1314 0.02189 0.0658 .  
# Residuals                               90   24.3649 0.27072         0.87074           
# Total                                   99   27.9820                 1.00000           
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# when you rarefy by 2000 you get

                                        # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# data.trans.rar$Date                      1    0.2648 0.26482  1.3149 0.01739 0.0118 *  
# data.trans.rar$Crop                      2    0.8664 0.43319  2.1509 0.05691 0.0001 ***
# data.trans.rar$SoilFrac                  4    0.9494 0.23735  1.1785 0.06236 0.0046 ** 
# data.trans.rar$Date:data.trans.rar$Crop  2    0.4560 0.22800  1.1321 0.02995 0.0726 .  
# Residuals                               63   12.6880 0.20140         0.83339           
# Total                                   72   15.2246                 1.00000           
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# run this code
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

# doing the mds for crop*date based on composition (or structure depending on who is defining it)
mds<-metaMDS(decostand(subset(data.nosing.reads, reads>100)[,-c(1:6)],"pa" ),k=6,autotransform=FALSE)
ggplot.NMDS(mds, (subset(data.nosing.reads, reads>100)$Crop), rainbow(3))
fig_16S<-ggplot.NMDS(mds, (data.nosing.rar$Crop), rainbow(3))
# doing the mds for soilfrac based on p/a 
mds<-metaMDS(decostand(data.nosing.rar[,-c(1:5)],"pa"),k=2,autotransform=FALSE)
ggplot.NMDS(mds,data.nosing.rar$SoilFrac, rainbow(5))


