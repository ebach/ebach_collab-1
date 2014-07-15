#Elizabeth Bach
#COBS Aggregates Fungal ITS analysis
#taxanomic-level analysis, looking among aggregates, date and crop
#11 July 2014

#Use data_taxa2 from diversity_stats.R and Basidio.data, Asco.data, Glom.data, Zygo.data, and Chytri.data from COBS_ITS_PhylaSubset
#reshape library contains colsplit function used to generate crop, SoilFrac, and Date in summerized files
library(reshape)
library(lme4)
library(lmerTest)
library(bbmle)

#Phyla-level diversity, summarizes total counts for each phylum
Phyla.data2<-ddply(data_taxa2, .(SampleName, Date, Crop, Block, SoilFrac, Phylum), summarise, .drop=FALSE, .progress="text",total=sum(value))
head(Phyla.data2)
Phyla.totals<-cast(data_taxa2, SampleName + Date + Crop + Block + SoilFrac ~ Phylum, sum)
head(Phyla.totals)

#mixed model for differences across phyla
Phyla.null<-lmer(total~1+(1|Block), data=Phyla.data2, REML=FALSE)
Phyla.model.full<-lmer(total~Phylum*Date*Crop*SoilFrac+(1|Block), data=Phyla.data2, REML=FALSE)
Phyla.model.main<-lmer(total~Phylum+Date+Crop+SoilFrac+(1|Block), data=Phyla.data2, REML=FALSE)
AICtab(Phyla.null,Phyla.model.full,Phyla.model.main)
anova(Phyla.null,Phyla.model.full,Phyla.model.main)
anova(Phyla.model.main)
anova(Phyla.model.full)
anova(Phyla.null)
#Phylum has very strong effect, so best to analyize crop, date, frac, within each phyla

#Ascomycota analysis
Asco.null<-lmer(value~1+(1|Block), data=Asco.data, REML=FALSE)
Asco.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Asco.data, REML=FALSE)
Asco.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Asco.data, REML=FALSE)
AICtab(Asco.null,Asco.model.full,Asco.model.main)
anova(Asco.null,Asco.model.full,Asco.model.main)
anova(Asco.model.main)
anova(Asco.model.full)
#Main effects model is best fit, no significant interactions in full model anyway
#Date and Crop both significant, SoilFrac not

#Basidiomycota analysis
Basidio.null<-lmer(value~1+(1|Block), data=Basidio.data, REML=FALSE)
Basidio.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Basidio.data, REML=FALSE)
Basidio.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Basidio.data, REML=FALSE)
AICtab(Basidio.null,Basidio.model.full,Basidio.model.main)
anova(Basidio.null,Basidio.model.full,Basidio.model.main)
anova(Basidio.model.main)
anova(Basidio.model.full)
#Null model is lowest AIC, but only by 2.  Main effects is next best, significant crop effect
#No interactions
#Null may be significant due to 0s in some date/crop/soilFrac combos

#Glomeromycota analysis
Glom.null<-lmer(value~1+(1|Block), data=Glom.data, REML=FALSE)
Glom.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Glom.data, REML=FALSE)
Glom.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Glom.data, REML=FALSE)
AICtab(Glom.null,Glom.model.full,Glom.model.main)
anova(Glom.null,Glom.model.full,Glom.model.main)
anova(Glom.model.main)
anova(Glom.model.full)
#Null model lowest AIC, but likely driven by zeros in some clases (e.g. "unbalanced") Full model is next best fit
#no interactions, but Crop and SoilFrac are significant

#Zygomycota analysis
Zygo.null<-lmer(value~1+(1|Block), data=Zygo.data, REML=FALSE)
Zygo.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Zygo.data, REML=FALSE)
Zygo.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Zygo.data, REML=FALSE)
AICtab(Zygo.null,Zygo.model.full,Zygo.model.main)
anova(Zygo.null,Zygo.model.full,Zygo.model.main)
anova(Zygo.model.main)
anova(Zygo.model.full)
#Full model best fit, by three
#Crop*SoilFrac and Date*SoilFrac interactions significant
#Interesting that Zygo is showing differences across fractions!

#Chytridiomycota analysis
Chytri.null<-lmer(value~1+(1|Block), data=Chytri.data, REML=FALSE)
Chytri.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Chytri.data, REML=FALSE)
Chytri.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Chytri.data, REML=FALSE)
AICtab(Chytri.null,Chytri.model.full,Chytri.model.main)
anova(Chytri.null,Chytri.model.full,Chytri.model.main)
anova(Chytri.model.main)
anova(Chytri.model.full)
#Null model is lowest AIC, in general Chytridiomycota has low abundance, may be contributing
#Full model is next best choice
#Date*Crop, Crop*SoilFrac, Date*SoilFrac all significant interactions

#unk analysis
Unk.null<-lmer(value~1+(1|Block), data=Unk.data, REML=FALSE)
Unk.model.full<-lmer(value~Date*Crop*SoilFrac+(1|Block), data=Unk.data, REML=FALSE)
Unk.model.main<-lmer(value~Date+Crop+SoilFrac+(1|Block), data=Unk.data, REML=FALSE)
AICtab(Unk.null,Unk.model.full,Unk.model.main)
anova(Unk.null,Unk.model.full,Unk.model.main)
anova(Unk.model.main)
anova(Unk.model.full)
#main model is best fit, but full model shows significant Date*SoilFrac interaction


