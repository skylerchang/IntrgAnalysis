library(tidyverse)
library(here)
library(RColorBrewer)
library(reshape2)
library(gsubfn)
library(gridExtra)
library(ggplot2)


#*************************** adjust settings ************************
#Is SampleNo coding? values: T/F; (e.g. MRD - yes: jake-01.1 (pbmcs), jake-01.2 (plasma); Histiocytoma - no)
sampleNoCoding<-1

#********************************************************************
setwd(here())
getwd()

targetDir<-'../Results/'

files<-list.files(targetDir,pattern = "^locus.*xlsx")
#pick which file to use
# 1 = runs 21, 24, 26
# 2 = run 28

j<-1
j<-2

  L<-readxl::read_excel(paste0(targetDir,files[j]))
  
# adding the means 
L$mean.IGH<-(L$IGH.rep1+L$IGH.rep2)/2
L$mean.TRB<-(L$TRB.rep1+L$TRB.rep2)/2

#calculate the ratio 
L$locus.Ratio<-L$mean.IGH/L$mean.TRB

#plot ratios agaisnt each other 
#no log
pdf("../Results/locus ratio.pdf")
ggplot(L,aes(locus.Ratio,submission))+geom_point()
dev.off()
#log axis 
pdf("../Results/locus ratio(log-axis).pdf")
ggplot(L,aes(locus.Ratio,submission))+geom_point()+scale_x_log10()+ geom_vline(xintercept = 1)
dev.off()
#log transformation
pdf("../Results/locus ratio(log).pdf")
ggplot(L,aes(log(locus.Ratio),submission))+geom_point()+ geom_vline(xintercept = 0)
dev.off()
#summary
summary(L$mean.IGH)
summary(L$mean.TRB)
summary(L$locus.Ratio)

