# This starts from data.nosing.rar as a result from data_wrangling.R

#first transforming to scaled data within each sample (between 0 and 1) for multivariate tests.  This uses decostand() from the vegan pacakge
head(data.nosing.rar[,1:5])
data.trans.rar<-cbind(data.nosing.rar[,1:5],decostand(data.nosing.rar[,-c(1:5)],"total"))

#using adonis for the test; I pulled out non-significant interactions
adonis(data.trans.rar[,-c(1:5)]~data.trans.rar$Date*data.trans.rar$Crop+data.trans.rar$SoilFrac, permutations=9999)
