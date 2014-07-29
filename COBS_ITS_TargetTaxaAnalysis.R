#Elizabeth Bach
#COBS ITS Target Taxa:  Point Range Graphs + Crop/Agg stats
#28 July, 2014

#Use taxa.interest from COBS_ITS_TargetTaxaNMDS.R
#subset so can evaluate July and Oct independently
taxa.July<-subset(taxa.interest, taxa.interest$Date=="Jul-12")
taxa.Oct<-subset(taxa.interest, taxa.interest$Date=="Oct-12")

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

#Limonomyces only found in P samples, represented in 18 samples, removing Crop from model to look at SoilFrac and Date
Limono.null<-lmer(g__Limonomyces~1+(1|Block), data=taxa.interest, REML=FALSE)
Limono.model.full<-lmer(g__Limonomyces~Date*SoilFrac+(1|Block), data=taxa.interest, REML=FALSE)
Limono.model.main<-lmer(g__Limonomyces~Date+SoilFrac+(1|Block), data=taxa.interest, REML=FALSE)
AICtab(Limono.null,Limono.model.full,Limono.model.main)
#null model lowest AIC by 2
anova(Limono.model.full)
anova(Limono.model.main)
difflsmeans(Limono.model.main, ddf="Satterthwaite",type=3,method.grad="simple")
#Date is significant
#Looking at pointrange:
Limono.Crop<-ddply(taxa.interest, .(Crop), summarise, .progress="text",mean=mean(g__Limonomyces),
high95=boot.high(g__Limonomyces), low95=boot.low(g__Limonomyces))
ggplot(Limono.Crop)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))
#Crop is most inresting diffrence

Limono.Date<-ddply(taxa.interest, .(Date), summarise, .progress="text",mean=mean(g__Limonomyces),
high95=boot.high(g__Limonomyces), low95=boot.low(g__Limonomyces))
ggplot(Limono.Date)+geom_pointrange(aes(x=Date,y=mean,ymax=high95,ymin=low95))

Limono.SoilFrac<-ddply(taxa.interest, .(SoilFrac), summarise, .progress="text",mean=mean(g__Limonomyces),
high95=boot.high(g__Limonomyces), low95=boot.low(g__Limonomyces))
ggplot(Limono.SoilFrac)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))
#No clear statistical differences, slightly lower in MM

#Great abundance in July, looking within July for aggregate differences:
Limono.JulySoilFrac<-ddply(taxa.July, .(SoilFrac), summarise, .progress="text",mean=mean(g__Limonomyces),
high95=boot.high(g__Limonomyces), low95=boot.low(g__Limonomyces))
ggplot(Limono.JulySoilFrac)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))
#Doesn't really change relationships

#Atheliales, all are Athelia bombacina, not well represented across samples, but abundant when found
Atheli.null<-lmer(o__Atheliales~1+(1|Block), data=taxa.interest, REML=FALSE)
Atheli.model.full<-lmer(o__Atheliales~Date*SoilFrac*Crop+(1|Block), data=taxa.interest, REML=FALSE)
Atheli.model.main<-lmer(o__Atheliales~Date+SoilFrac+Crop+(1|Block), data=taxa.interest, REML=FALSE)
AICtab(Atheli.null,Atheli.model.full,Atheli.model.main)
#null model best fit by 3, main model next lowest AIC
anova(Atheli.model.main)
anova(Atheli.model.full)
#Crop is significant main effect, full model shows Date*Crop interaction
Atheli.Crop<-ddply(taxa.interest, .(Crop), summarise, .progress="text",mean=mean(o__Atheliales),
high95=boot.high(o__Atheliales), low95=boot.low(o__Atheliales))
ggplot(Atheli.Crop)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))

Atheli.Date<-ddply(taxa.interest, .(Date), summarise, .progress="text",mean=mean(o__Atheliales),
high95=boot.high(o__Atheliales), low95=boot.low(o__Atheliales))
ggplot(Atheli.Date)+geom_pointrange(aes(x=Date,y=mean,ymax=high95,ymin=low95))

Atheli.SoilFrac<-ddply(taxa.interest, .(SoilFrac), summarise, .progress="text",mean=mean(o__Atheliales),
high95=boot.high(o__Atheliales), low95=boot.low(o__Atheliales))
ggplot(Atheli.SoilFrac)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))


#UnkBasidio, highly abundant, espeically in prairies
UnkBasidio.null<-lmer(V1~1+(1|Block), data=taxa.interest, REML=FALSE)
UnkBasidio.model.full<-lmer(V1~Date*SoilFrac*Crop+(1|Block), data=taxa.interest, REML=FALSE)
UnkBasidio.model.main<-lmer(V1~Date+SoilFrac+Crop+(1|Block), data=taxa.interest, REML=FALSE)
AICtab(UnkBasidio.null,UnkBasidio.model.full,UnkBasidio.model.main)
#main model has lowest AIC by 19
anova(UnkBasidio.model.main)
anova(UnkBasidio.model.full)
#Crop & Date significant, no interactions, so main model is best
difflsmeans(UnkBasidio.model.main, ddf="Satterthwaite",type=3,method.grad="simple")
UnkBasidio.Crop<-ddply(taxa.interest, .(Crop), summarise, .progress="text",mean=mean(V1),
high95=boot.high(V1), low95=boot.low(V1))
ggplot(UnkBasidio.Crop)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))

UnkBasidio.Date<-ddply(taxa.interest, .(Date), summarise, .progress="text",mean=mean(V1),
high95=boot.high(V1), low95=boot.low(V1))
ggplot(UnkBasidio.Date)+geom_pointrange(aes(x=Date,y=mean,ymax=high95,ymin=low95))

UnkBasidio.SoilFrac<-ddply(taxa.interest, .(SoilFrac), summarise, .progress="text",mean=mean(V1),
high95=boot.high(V1), low95=boot.low(V1))
ggplot(UnkBasidio.SoilFrac)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))
#Looking within dates
#July
UnkBasidio.July.null<-lmer(V1~1+(1|Block), data=taxa.July, REML=FALSE)
UnkBasidio.July.full<-lmer(V1~SoilFrac*Crop+(1|Block), data=taxa.July, REML=FALSE)
UnkBasidio.July.main<-lmer(V1~SoilFrac+Crop+(1|Block), data=taxa.July, REML=FALSE)
AICtab(UnkBasidio.July.null,UnkBasidio.July.full,UnkBasidio.July.main)
#main modle lowest AIC by 7,
anova(UnkBasidio.July.main) #Crop significant, no SoilFrac effect
anova(UnkBasidio.July.full) #No interactions
difflsmeans(UnkBasidio.July.main, ddf="Satterthwaite",type=3,method.grad="simple")

UnkBasidio.Crop<-ddply(taxa.July, .(Crop), summarise, .progress="text",mean=mean(V1),
high95=boot.high(V1), low95=boot.low(V1))
ggplot(UnkBasidio.Crop)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))
#Same as full dataset

UnkBasidio.SoilFrac<-ddply(taxa.July, .(SoilFrac), summarise, .progress="text",mean=mean(V1),
high95=boot.high(V1), low95=boot.low(V1))
ggplot(UnkBasidio.SoilFrac)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))
#LM lower than rest, but still not statistically significant

#October
UnkBasidio.Oct.null<-lmer(V1~1+(1|Block), data=taxa.Oct, REML=FALSE)
UnkBasidio.Oct.full<-lmer(V1~SoilFrac*Crop+(1|Block), data=taxa.Oct, REML=FALSE)
UnkBasidio.Oct.main<-lmer(V1~SoilFrac+Crop+(1|Block), data=taxa.Oct, REML=FALSE)
AICtab(UnkBasidio.Oct.null,UnkBasidio.Oct.full,UnkBasidio.Oct.main)
#main effects lowest AIC by 6
anova(UnkBasidio.Oct.main) #crop only significant factor

UnkBasidio.Crop<-ddply(taxa.Oct, .(Crop), summarise, .progress="text",mean=mean(V1),
high95=boot.high(V1), low95=boot.low(V1))
ggplot(UnkBasidio.Crop)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))
#Same as full dataset, just tighter errors on P, PF (PF nearly 0)

UnkBasidio.SoilFrac<-ddply(taxa.Oct, .(SoilFrac), summarise, .progress="text",mean=mean(V1),
high95=boot.high(V1), low95=boot.low(V1))
ggplot(UnkBasidio.SoilFrac)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))
#lots of spread, somewhat higher in WS, but nothing statistically robust

#Remove 0s
UnkBasidio.data<-subset(taxa.interest, taxa.interest$V1>0)
head(UnkBasidio.data)

UnkBasidio.SoilFrac2<-ddply(UnkBasidio.data, .(SoilFrac), summarise, .progress="text",mean=mean(V1),
high95=boot.high(V1), low95=boot.low(V1))
ggplot(UnkBasidio.SoilFrac2)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))


#Thanatephorus, Main effect of Date (P<0.0001), Crop (P<0.0001), and SoilFrac(P=0.02) (see Basidiomycota code)
#Including all 0s (presence/absence?)
Than.Crop<-ddply(taxa.interest, .(Crop), summarise, .progress="text",mean=mean(g__Thanatephorus),
high95=boot.high(g__Thanatephorus), low95=boot.low(g__Thanatephorus))
ggplot(Than.Crop)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))

Than.Date<-ddply(taxa.interest, .(Date), summarise, .progress="text",mean=mean(g__Thanatephorus),
high95=boot.high(g__Thanatephorus), low95=boot.low(g__Thanatephorus))
ggplot(Than.Date)+geom_pointrange(aes(x=Date,y=mean,ymax=high95,ymin=low95))

Than.SoilFrac<-ddply(taxa.interest, .(SoilFrac), summarise, .progress="text",mean=mean(g__Thanatephorus),
high95=boot.high(g__Thanatephorus), low95=boot.low(g__Thanatephorus))
ggplot(Than.SoilFrac)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))
#Dropping 0's (abundance?)
Than.data<-subset(taxa.interest, taxa.interest$g__Thanatephorus>0)
head(Than.data)
Than.Crop2<-ddply(Than.data, .(Crop), summarise, .progress="text",mean=mean(g__Thanatephorus),
high95=boot.high(g__Thanatephorus), low95=boot.low(g__Thanatephorus))
ggplot(Than.Crop2)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))

Than.SoilFrac2<-ddply(Than.data, .(SoilFrac), summarise, .progress="text",mean=mean(g__Thanatephorus),
high95=boot.high(g__Thanatephorus), low95=boot.low(g__Thanatephorus))
ggplot(Than.SoilFrac2)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))

#Psathyrellaceae (inky caps)
#SoilFrac P=0.009, Crop P=0.07, P<PF=CC, LM<WS, micro<WS, MM<WS (see Basidiomycota code)
#Including all 0s (presence/absence)
Psath.Crop<-ddply(taxa.interest, .(Crop), summarise, .progress="text",mean=mean(f__Psathyrellaceae),
high95=boot.high(f__Psathyrellaceae), low95=boot.low(f__Psathyrellaceae))
ggplot(Psath.Crop)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))

Psath.Date<-ddply(taxa.interest, .(Date), summarise, .progress="text",mean=mean(f__Psathyrellaceae),
high95=boot.high(f__Psathyrellaceae), low95=boot.low(f__Psathyrellaceae))
ggplot(Psath.Date)+geom_pointrange(aes(x=Date,y=mean,ymax=high95,ymin=low95))

Psath.SoilFrac<-ddply(taxa.interest, .(SoilFrac), summarise, .progress="text",mean=mean(f__Psathyrellaceae),
high95=boot.high(f__Psathyrellaceae), low95=boot.low(f__Psathyrellaceae))
ggplot(Psath.SoilFrac)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))

#Removing 0s to see stat. differences
Psath.data<-subset(taxa.interest, taxa.interest$f__Psathyrellaceae>0)
head(Psath.data)

Psath.Crop2<-ddply(Psath.data, .(Crop), summarise, .progress="text",mean=mean(f__Psathyrellaceae),
high95=boot.high(f__Psathyrellaceae), low95=boot.low(f__Psathyrellaceae))
ggplot(Psath.Crop2)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))

Psath.SoilFrac2<-ddply(Psath.data, .(SoilFrac), summarise, .progress="text",mean=mean(f__Psathyrellaceae),
high95=boot.high(f__Psathyrellaceae), low95=boot.low(f__Psathyrellaceae))
ggplot(Psath.SoilFrac2)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))

#Strophariaceae, strong crop effect from models (see Basidiomycota.R), also significant Date*Crop*SoilFrac intercation, worth pulling apart
#looking at full set (includes 0's)
Stroph.Crop<-ddply(taxa.interest, .(Crop), summarise, .progress="text",mean=mean(f__Strophariaceae),
high95=boot.high(f__Strophariaceae), low95=boot.low(f__Strophariaceae))
ggplot(Stroph.Crop)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))

Stroph.Date<-ddply(taxa.interest, .(Date), summarise, .progress="text",mean=mean(f__Strophariaceae),
high95=boot.high(f__Strophariaceae), low95=boot.low(f__Strophariaceae))
ggplot(Stroph.Date)+geom_pointrange(aes(x=Date,y=mean,ymax=high95,ymin=low95))

Stroph.SoilFrac<-ddply(taxa.interest, .(SoilFrac), summarise, .progress="text",mean=mean(f__Strophariaceae),
high95=boot.high(f__Strophariaceae), low95=boot.low(f__Strophariaceae))
ggplot(Stroph.SoilFrac)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))

#Subsetting to remove 0's
Stroph.data<-subset(taxa.interest, taxa.interest$f__Strophariaceae>0)
head(Stroph.data)
Stroph.Crop2<-ddply(Stroph.data, .(Crop), summarise, .progress="text",mean=mean(f__Strophariaceae),
high95=boot.high(f__Strophariaceae), low95=boot.low(f__Strophariaceae))
ggplot(Stroph.Crop2)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))

Stroph.SoilFrac2<-ddply(Stroph.data, .(SoilFrac), summarise, .progress="text",mean=mean(f__Strophariaceae),
high95=boot.high(f__Strophariaceae), low95=boot.low(f__Strophariaceae))
ggplot(Stroph.SoilFrac2)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))

#Looking at each sampling date
StrophJuly.Crop<-ddply(taxa.July, .(Crop), summarise, .progress="text",mean=mean(f__Strophariaceae),
high95=boot.high(f__Strophariaceae), low95=boot.low(f__Strophariaceae))
ggplot(StrophJuly.Crop)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))
#Cropping system difference less strong

StrophJuly.SoilFrac<-ddply(taxa.July, .(SoilFrac), summarise, .progress="text",mean=mean(f__Strophariaceae),
high95=boot.high(f__Strophariaceae), low95=boot.low(f__Strophariaceae))
ggplot(StrophJuly.SoilFrac)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))
#Run model to test SoilFrac in July, micro and SM seem higher than others
StrophJuly.null<-lmer(f__Strophariaceae~1+(1|Block), data=taxa.July, REML=FALSE)
StrophJuly.model.full<-lmer(f__Strophariaceae~Crop*SoilFrac+(1|Block), data=taxa.July, REML=FALSE)
StrophJuly.model.main<-lmer(f__Strophariaceae~Crop+SoilFrac+(1|Block), data=taxa.July, REML=FALSE)
AICtab(StrophJuly.null,StrophJuly.model.full,StrophJuly.model.main)
#null model has lowest AIC
anova(Stroph.model.main) #Only Cropping system significant

StrophOct.Crop<-ddply(taxa.Oct, .(Crop), summarise, .progress="text",mean=mean(f__Strophariaceae),
high95=boot.high(f__Strophariaceae), low95=boot.low(f__Strophariaceae))
ggplot(StrophOct.Crop)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))
#Croppping systme difference strong, similar to full dataset

StrophOct.SoilFrac<-ddply(taxa.Oct, .(SoilFrac), summarise, .progress="text",mean=mean(f__Strophariaceae),
high95=boot.high(f__Strophariaceae), low95=boot.low(f__Strophariaceae))
ggplot(StrophOct.SoilFrac)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))
#very little diff between fractions (slightly lower in WS)

#Without 0s in July:  13 observations, 2 CCs, rest P and PF, very high values for SM, micros
StrophJuly.data<-subset(taxa.July, taxa.July$f__Strophariaceae>0)

StrophJuly.SoilFrac<-ddply(StrophJuly.data, .(SoilFrac), summarise, .progress="text",mean=mean(f__Strophariaceae),
high95=boot.high(f__Strophariaceae), low95=boot.low(f__Strophariaceae))
ggplot(StrophJuly.SoilFrac)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))
#run model without 0s
StrophJuly.null<-lmer(f__Strophariaceae~1+(1|Block), data=StrophJuly.data, REML=FALSE)
StrophJuly.model.full<-lmer(f__Strophariaceae~Crop*SoilFrac+(1|Block), data=StrophJuly.data, REML=FALSE)
StrophJuly.model.main<-lmer(f__Strophariaceae~Crop+SoilFrac+(1|Block), data=StrophJuly.data, REML=FALSE)
AICtab(StrophJuly.null,StrophJuly.model.full,StrophJuly.model.main)
#full model lowest AIC
anova(StrophJuly.model.full)
#Significant crop*SoilFrac interacation, not present in enough samples to pull apart by cropping system

#Peziza, Date is significant (P=0.047), SoilFrac marginal (P=0.08), see Ascomycota.R
#looking at full set (includes 0's)
Pez.Crop<-ddply(taxa.interest, .(Crop), summarise, .progress="text",mean=mean(g__Peziza),
high95=boot.high(g__Peziza), low95=boot.low(g__Peziza))
ggplot(Pez.Crop)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))

Pez.Date<-ddply(taxa.interest, .(Date), summarise, .progress="text",mean=mean(g__Peziza),
high95=boot.high(g__Peziza), low95=boot.low(g__Peziza))
ggplot(Pez.Date)+geom_pointrange(aes(x=Date,y=mean,ymax=high95,ymin=low95))

Pez.SoilFrac<-ddply(taxa.interest, .(SoilFrac), summarise, .progress="text",mean=mean(g__Peziza),
high95=boot.high(g__Peziza), low95=boot.low(g__Peziza))
ggplot(Pez.SoilFrac)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))

#Dropping 0s, 16 observations, mostly from Oct.
Pez.data<-subset(taxa.interest, taxa.interest$g__Peziza>0)
head(Pez.data)
Pez.Crop2<-ddply(Pez.data, .(Crop), summarise, .progress="text",mean=mean(g__Peziza),
high95=boot.high(g__Peziza), low95=boot.low(g__Peziza))
ggplot(Pez.Crop2)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))

Pez.SoilFrac2<-ddply(Pez.data, .(SoilFrac), summarise, .progress="text",mean=mean(g__Peziza),
high95=boot.high(g__Peziza), low95=boot.low(g__Peziza))
ggplot(Pez.SoilFrac2)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))
#Running models to look at SoilFrac effect without 0s
Pez.null<-lmer(g__Peziza~1+(1|Block), data=Pez.data, REML=FALSE)
Pez.model.full<-lmer(g__Peziza~Date*Crop*SoilFrac+(1|Block), data=Pez.data, REML=FALSE)
Pez.model.main<-lmer(g__Peziza~Date+Crop+SoilFrac+(1|Block), data=Pez.data, REML=FALSE)
AICtab(Pez.null,Pez.model.full,Pez.model.main)
#full model best fit
anova(Pez.model.full)
#Date interactions, take a look at Oct seperately
Pez.Oct<-subset(taxa.Oct, taxa.Oct$g__Peziza>0)
Pez.Oct
Pez.SoilFracOct<-ddply(Pez.Oct, .(SoilFrac), summarise, .progress="text",mean=mean(g__Peziza),
high95=boot.high(g__Peziza), low95=boot.low(g__Peziza))
ggplot(Pez.SoilFracOct)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))
Pez.null<-lmer(g__Peziza~1+(1|Block), data=Pez.Oct, REML=FALSE)
Pez.model.full<-lmer(g__Peziza~Crop*SoilFrac+(1|Block), data=Pez.Oct, REML=FALSE)
Pez.model.main<-lmer(g__Peziza~Crop+SoilFrac+(1|Block), data=Pez.Oct, REML=FALSE)
AICtab(Pez.null,Pez.model.full,Pez.model.main)
#main model and full model essentially same AIC
anova(Pez.model.full) #no interaction
anova(Pez.model.main) #Crop and SoilFrac both significant
difflsmeans(Pez.model.main, ddf="Satterthwaite",type=3,method.grad="simple")

Pez.CropOct<-ddply(Pez.Oct, .(Crop), summarise, .progress="text",mean=mean(g__Peziza),
high95=boot.high(g__Peziza), low95=boot.low(g__Peziza))
ggplot(Pez.CropOct)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))

#Bionectriaceae, strong SoilFrac effect, P=0.006 (see Ascomycota.R)
#looking at full set (includes 0's)
Bion.Crop<-ddply(taxa.interest, .(Crop), summarise, .progress="text",mean=mean(f__Bionectriaceae),
high95=boot.high(f__Bionectriaceae), low95=boot.low(f__Bionectriaceae))
ggplot(Bion.Crop)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))

Bion.Date<-ddply(taxa.interest, .(Date), summarise, .progress="text",mean=mean(f__Bionectriaceae),
high95=boot.high(f__Bionectriaceae), low95=boot.low(f__Bionectriaceae))
ggplot(Bion.Date)+geom_pointrange(aes(x=Date,y=mean,ymax=high95,ymin=low95))

Bion.SoilFrac<-ddply(taxa.interest, .(SoilFrac), summarise, .progress="text",mean=mean(f__Bionectriaceae),
high95=boot.high(f__Bionectriaceae), low95=boot.low(f__Bionectriaceae))
ggplot(Bion.SoilFrac)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))

#Dropping 0s, clarifies, stregthens diffs, LM=micro for highest abundance, oddly, note low abundance across the board
Bion.data<-subset(taxa.interest, taxa.interest$f__Bionectriaceae>0)
head(Bion.data)

Bion.SoilFrac2<-ddply(Bion.data, .(SoilFrac), summarise, .progress="text",mean=mean(f__Bionectriaceae),
high95=boot.high(f__Bionectriaceae), low95=boot.low(f__Bionectriaceae))
ggplot(Bion.SoilFrac2)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))

Bion.Crop2<-ddply(Bion.data, .(Crop), summarise, .progress="text",mean=mean(f__Bionectriaceae),
high95=boot.high(f__Bionectriaceae), low95=boot.low(f__Bionectriaceae))
ggplot(Bion.Crop2)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))

#Glomerales
#looking at full set (includes 0's)
Glom.Crop<-ddply(taxa.interest, .(Crop), summarise, .progress="text",mean=mean(o__Glomerales),
high95=boot.high(o__Glomerales), low95=boot.low(o__Glomerales))
ggplot(Glom.Crop)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))

Glom.Date<-ddply(taxa.interest, .(Date), summarise, .progress="text",mean=mean(o__Glomerales),
high95=boot.high(o__Glomerales), low95=boot.low(o__Glomerales))
ggplot(Glom.Date)+geom_pointrange(aes(x=Date,y=mean,ymax=high95,ymin=low95))

Glom.SoilFrac<-ddply(taxa.interest, .(SoilFrac), summarise, .progress="text",mean=mean(o__Glomerales),
high95=boot.high(o__Glomerales), low95=boot.low(o__Glomerales))
ggplot(Glom.SoilFrac)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))

#Dropping 0s, 
Glom.data<-subset(taxa.interest, taxa.interest$o__Glomerales>0)
head(Glom.data)

Glom.SoilFrac2<-ddply(Glom.data, .(SoilFrac), summarise, .progress="text",mean=mean(o__Glomerales),
high95=boot.high(o__Glomerales), low95=boot.low(o__Glomerales))
ggplot(Glom.SoilFrac2)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))

Glom.Crop2<-ddply(Glom.data, .(Crop), summarise, .progress="text",mean=mean(o__Glomerales),
high95=boot.high(o__Glomerales), low95=boot.low(o__Glomerales))
ggplot(Glom.Crop2)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))

Glom.July<-subset(taxa.July, taxa.July$o__Glomerales>0)
Glom.July
Glom.SoilFracJuly<-ddply(Glom.July, .(SoilFrac), summarise, .progress="text",mean=mean(o__Glomerales),
high95=boot.high(o__Glomerales), low95=boot.low(o__Glomerales))
ggplot(Glom.SoilFracJuly)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))

#Operculomyces
#looking at full set (includes 0's)
Operc.Crop<-ddply(taxa.interest, .(Crop), summarise, .progress="text",mean=mean(g__Operculomyces),
high95=boot.high(g__Operculomyces), low95=boot.low(g__Operculomyces))
ggplot(Operc.Crop)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))

Operc.Date<-ddply(taxa.interest, .(Date), summarise, .progress="text",mean=mean(g__Operculomyces),
high95=boot.high(g__Operculomyces), low95=boot.low(g__Operculomyces))
ggplot(Operc.Date)+geom_pointrange(aes(x=Date,y=mean,ymax=high95,ymin=low95))

Operc.SoilFrac<-ddply(taxa.interest, .(SoilFrac), summarise, .progress="text",mean=mean(g__Operculomyces),
high95=boot.high(g__Operculomyces), low95=boot.low(g__Operculomyces))
ggplot(Operc.SoilFrac)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))
#running models
Operc.null<-lmer(g__Operculomyces~1+(1|Block), data=taxa.interest, REML=FALSE)
Operc.model.full<-lmer(g__Operculomyces~Date*Crop*SoilFrac+(1|Block), data=taxa.interest, REML=FALSE)
Operc.model.main<-lmer(g__Operculomyces~Date+Crop+SoilFrac+(1|Block), data=taxa.interest, REML=FALSE)
AICtab(Operc.null,Operc.model.full,Operc.model.main)
#null model best fit, main second lowest AIC
anova(Operc.model.main) #Crop is only significat effect, P=0.02
anova(Operc.model.full) #Significant Crop*SoilFrac interaction

#Dropping 0s, 
Operc.data<-subset(taxa.interest, taxa.interest$g__Operculomyces>0)
head(Operc.data)
Operc.data

Operc.SoilFrac2<-ddply(Operc.data, .(SoilFrac), summarise, .progress="text",mean=mean(g__Operculomyces),
high95=boot.high(g__Operculomyces), low95=boot.low(g__Operculomyces))
ggplot(Operc.SoilFrac2)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))

Operc.null2<-lmer(g__Operculomyces~1+(1|Block), data=Operc.data, REML=FALSE)
Operc.model.full2<-lmer(g__Operculomyces~Date*Crop*SoilFrac+(1|Block), data=Operc.data, REML=FALSE)
Operc.model.main2<-lmer(g__Operculomyces~Date+Crop+SoilFrac+(1|Block), data=Operc.data, REML=FALSE)
AICtab(Operc.null2,Operc.model.full2,Operc.model.main2)
#Full model best fit, very similar AIC to main
anova(Operc.model.full2) #unbalanced, gives Date*SoilFrac interaction, but dropping values

#most from July, subset to just look at July
Operc.SoilFracJuly<-ddply(taxa.July, .(SoilFrac), summarise, .progress="text",mean=mean(g__Operculomyces),
high95=boot.high(g__Operculomyces), low95=boot.low(g__Operculomyces))
ggplot(Operc.SoilFracJuly)+geom_pointrange(aes(x=SoilFrac,y=mean,ymax=high95,ymin=low95))

Operc.CropJuly<-ddply(taxa.July, .(Crop), summarise, .progress="text",mean=mean(g__Operculomyces),
high95=boot.high(g__Operculomyces), low95=boot.low(g__Operculomyces))
ggplot(Operc.CropJuly)+geom_pointrange(aes(x=Crop,y=mean,ymax=high95,ymin=low95))
