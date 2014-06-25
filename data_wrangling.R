rm(list=ls())

#all useful packages, install.packages(...package name...) if you don't have them
library(labdsv)
library(vegan)
library(plyr)
library(reshape)
library(ggplot2)



# setwd("/Users/ryanjwtx/Desktop/bach_collab/")
# dataset<-read.csv("EBach_COBS_ITS_data.csv")


# Read in the data
dataset<-read.csv(file.choose())

# First we remove the singletons using the dropspc() function form the labdsv package.  In the line below I bind metadata (columns 1 through 5)
# and the singleton-removed dataset.  Note that for the dropspc function I only use columns 6 through 9538 as these are all the OTUs (no metadata),
# and I remove species that occur 1 or less times


data.nosing<-cbind(dataset[,1:5],dropspc(dataset[,6:9538],1))

#figure out how many samples to rarefy by
summary(rowSums(data.nosing[,-c(1:5)]))
reads<-rowSums(data.nosing[,-c(1:5)])
data.nosing.reads<-cbind(reads,data.nosing)
hist(reads)
rared<-rarefy(subset(data.nosing.reads, reads > 1999)[,-c(1:6)],sample=c(1,10,100,250,500,750,1000,2000),se=FALSE )
rared_melt<-melt(rared)
summary(rared_melt)
dim(subset(data.nosing.reads, reads > 999))