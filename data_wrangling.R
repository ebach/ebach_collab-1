rm(list=ls())

#all useful packages, install.packages(...package name...) if you don't have them
library(labdsv)
library(vegan)
library(plyr)
library(reshape)
library(ggplot2)



# setwd("/Users/ryanjwtx/Desktop/bach_collab/")
dataset<-read.csv("EBach_COBS_ITS_data.csv")


# Read in the data
dataset<-read.csv(file.choose())

# First we remove the singletons using the dropspc() function form the labdsv package.  In the line below I bind metadata (columns 1 through 5)
# and the singleton-removed dataset.  Note that for the dropspc function I only use columns 6 through 9538 as these are all the OTUs (no metadata),
# and I remove species that occur 1 or less times


data.nosing<-cbind(dataset[,1:5],dropspc(dataset[,6:9538],1))

# Now its time to figure out how many reads to rarefy by...I added a column to our dataset of the total number of reads per sample (row)

reads<-rowSums(data.nosing[,-c(1:5)])
data.nosing.reads<-cbind(reads,data.nosing)
head(data.nosing.reads[,1:10])

hist(reads)
# the distribution of reads isn't that skewed...this makes it difficult to decide how many samples to rarefy by...

#we can see how many samples we have by subsetting by number of reads

dim(subset(data.nosing.reads, reads > 999))

# we maintain a high level of samples with a rarefication at 1000

# lets create rarefaction curves for each sample starting around 1000
rared<-rarefy(subset(data.nosing.reads, reads > 1499)[,-c(1:6)],sample=c(1,10,25,50,75,100,250,500,750,1000,1250,1500),se=FALSE )
rared_melt<-melt(rared)
names(rared_melt)<-c("sample","sample_size","OTUs")
rared_melt$sample_size<-c(rep(1,97),rep(10,97),rep(25,97),rep(50,97),rep(75,97),rep(100,97),rep(250,97),rep(500,97),rep(750,97),rep(1000,97),rep(1250,97),rep(1500,97))
head(rared_melt)

ggplot(rared_melt)+geom_line(aes(x=sample_size,y=OTUs,colour=sample,group=sample))+theme(aspect.ratio=1)+theme_bw()

# I will decide with 1000 (100 samples) and do a sensitivity analysis to see if moving this number changes the results
head(data.nosing[,1:5])
data.nosing.rar<-cbind(subset(data.nosing, reads > 999)[,1:5],rrarefy(subset(data.nosing,reads > 999)[,-c(1:5)],1000))
head(data.nosing[,1:10])