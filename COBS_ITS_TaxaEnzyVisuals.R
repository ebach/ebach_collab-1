#Elizabeth Bach
#COBS ITS:  Scatterplots of key taxa with enzy activity
#Chosen from NMDS work with env. vectors + species centroids
#30 July 2014

#Statistically, none of these regressions are significant,
#I want to look at the visually to see if uneven sample distribution is playing a role

#BX significantly correltated with prairie communities
#Regressing against some key taxa associated with prairies
#use taxa_enzy from COBS_ITS_TargetTaxaNMDS.R

Than.data<-subset(taxa_enzy, taxa_enzy$g__Thanatephorus>0)
Than.data2<-subset(Than.data, Than.data$g__Thanatephorus<100)
head(Than.data)
Psath.data<-subset(taxa_enzy, taxa_enzy$f__Psathyrellaceae>0)
head(Psath.data)
Stroph.data<-subset(taxa_enzy, taxa_enzy$f__Strophariaceae>0)
head(Stroph.data)

ggplot(Than.data2)+geom_point(aes(x=BX, y=g__Thanatephorus))
ggplot(Psath.data)+geom_point(aes(x=BX, y=f__Psathyrellaceae))
ggplot(Stroph.data)+geom_point(aes(x=BX, y=f__Strophariaceae))
#Strophariaceae only one with trend toward a correlation
#Is there a better way to look at presence/absence?

#NAG
Limono.data<-subset(taxa_enzy, taxa_enzy$g__Limonomyces>0)
head(Limono.data)
Limono.data2<-subset(Limono.data, Limono.data$g__Limonomyces<300)

ggplot(Limono.data)+geom_point(aes(x=NAG, y=g__Limonomyces))
ggplot(Limono.data2)+geom_point(aes(x=NAG, y=g__Limonomyces))
#With outlier removed, there may be postive correlation with NAG and Limonomyces