#script calculates diversity for split 'Clntab_RDS' files: 'Clntab_RDS_noTemplate' & 'Clntab_RDS_templateOnly'
#v14: process data from multiple sequencing runs

library(plyr)
library(tidyverse)
library(here)
#library(reshape2)
library(RColorBrewer)
#library(ggpubr)
library(openxlsx)
#library(ggExtra)
library(MASS)

setwd(here())
getwd()

source("2_postInterrogate_R/functions.R")

#=========== adjust =============

sampleNoCodesForFraction<-T        #T or F

rdsPath<-'../Data/Clntab_RDS/clntab_vAndJ.rds'
rdsPaths<-c(
  '../../Run13_k9MRD-Reilly-AllTimePoints/Data/Clntab_RDS/clntab_vAndJ_filtered.rds',
  '../../Run17_Bella_Marishka/Data/Clntab_RDS/clntab_vAndJ_filtered.rds',
  '../../Run22_MRD_Daisy_Sona_UntilRelapse/Data/Clntab_RDS/clntab_vAndJ_filtered.rds',
  '../../Run27_MRD_Eddie_Jake_UntilRelapse/Data/Clntab_RDS/clntab_vAndJ_filtered.rds',
  '../../Run29_MRD_Bentley_Cheech_Gypsy_UntilRelapse/Data/Clntab_RDS/clntab_vAndJ_filtered.rds'
  )
runs<-c(13,17,22,27,29)

outPathTile<-'../Results/TrackClones/Tileplots/'
outPathBar<-'../Results/TrackClones/Barplots/'
outPathDensity<-'../Results/TrackClones/DensityPlots/'
#===================================
dir.create(outPathTile,recursive=T)
dir.create(outPathBar,recursive=T)
dir.create(outPathDensity,recursive=T)

indexSummary<-tibble()
templateSummary<-tibble()
diversity<-tibble()

for (k in 1:length(rdsPaths)){
  print(rdsPaths[k])
  
  t<-read_rds(rdsPaths[k])
  
  #datalist contains one tibble per sample
  datalist<-t[[1]]
  print(paste0("Number of files: ",length(datalist)))
  #sample names are stored separately in a vector
  files_short<-t[[2]]
  files_short
  
  for (i in 1:length(datalist)){
    a<-datalist[[i]][,c("vGene","jGene","aaSeq","aaLength","readCount","completeNtSeq","locus")]
    if(is.null(a)){next}
    if(nrow(a)==0){next}
    
    #============ check for index sequences =======
    a$index<-'no'
    a$index[grepl("CARADYYDSFWAAFGYW",a$aaSeq)]<-'bella'
    a$index[grepl("CVTSVFSPW",a$aaSeq)]<-'bentley'
    a$index[grepl("CAPEIGWAAEYW",a$aaSeq)]<-'cheech'
    a$index[grepl("CAKDYYGTRNIYSLDYW",a$aaSeq)]<-'daisy'
    a$index[grepl("CVRSRRGYPETV#GMDYW",a$aaSeq)]<-'eddie1'
    a$index[grepl("CAKDGLVVPTGSIEYW",a$aaSeq)]<-'eddie2'
    a$index[grepl("CAKEWSSSNWYGDVGMDYW",a$aaSeq)]<-'gypsy'
    a$index[grepl("CARDYHGISWLEYW",a$aaSeq)]<-'jake'
    a$index[grepl("CVPFSPYGSWFADDQW",a$aaSeq)]<-'marishka'
    a$index[grepl("CVRGANSWFPDFDYW",a$aaSeq)]<-'reilly'
    a$index[grepl("CAKELYYDSYSVDYW",a$aaSeq)]<-'sona'
    
    #========== summarize index sequences ==============
    temp<-ddply(a[,c("readCount","index")],"index",na.rm=T, numcolwise(sum))
    if (nrow(temp)==0){next}
    temp$totalReads<-sum(temp$readCount)
    temp$percentage<-temp$readCount/temp$totalReads*100
    temp$run<-runs[k]
    temp$filename<-files_short[i]
    indexSummary<-bind_rows(indexSummary,temp)
  }
}
c<-indexSummary

#remove controls for now -> revisit
remove<-c("ntc","k9","pc")
c<-c[!grepl(paste0(remove,collapse='|'),c$filename),]

#remove non-index clones
c<-c[c$index!='no',]

#split filename
c<-splitFilename(c,sampleNoCodesForFraction)

#=============================================================
#========= index clone analysis - tile plots =================
#=============================================================

#++++++++++++ functions +++++++++++++++++
#plot tiles - facet: patient~fraction
#x: data, y: patient index; z: readCount or percentage (implement)
plotTile<-function(x,y){
  p<-ggplot(x,aes(submission,replicate,fill=readCount))+
    geom_tile(colour="darkgrey",size=0.25)+
    facet_grid(patient~fraction)+
    scale_fill_gradient(high='red',low='white')+
    scale_x_discrete(limits=(c(paste0(0,seq(1,9)),10,11,12)))+
    coord_equal()+
    guides(fill = guide_colourbar(title=NULL,barheight=2))
  if (y==1){
    p+theme(axis.title.x=element_blank(),strip.text.x = element_blank(),strip.background = element_rect(fill="red"))
  }else{
    p+theme(axis.title.x=element_blank(),strip.text.x = element_blank(),strip.background = element_rect(fill="lightgreen"))
  }
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#========= generate empty tibble ===================
#create empty matrix with the following variables to create the default plotting 'area'
#readCount,percentage,replicate,fraction,run,patient,submission
patient<-c("bella","bentle","cheech","daisy","eddie","gypsy","jake","marish","sona","reilly")
timepoint<-c(7,3,6,12,8,11,13,9,6,10)   #the number of submissions per patient
run<-c(17,29,29,22,27,29,27,17,22,13)   #the run number in which it was sequenced
empty<-tibble()
for (i in 1:length(timepoint)){
  temp<-tibble(
    readCount=rep(0,(4*timepoint[i])),
    percentage=rep(0,(4*timepoint[i])),
    replicate=rep(c("rep1","rep2"),(2*timepoint[i])),
    fraction=rep(rep(c("PBMCs","Plasma"),each=2),timepoint[i]),
    run=rep(run[i],4*timepoint[i]),
    patient=rep(patient[i],4*timepoint[i]),
    submission=formatC(rep(1:timepoint[i],each=4),width=2,flag="0")
  )
  empty<-bind_rows(empty,temp)
}

#plot to check layout
plotTile(empty,1)
empty[rowSums(is.na(empty))>0,]
#========= generate index tibble ===================
clones<-c("bella","bentle","cheech","daisy","eddie1","eddie2","gypsy","jake","marishka","sona","reilly")

for (k in 1:length(clones)){
  #subset by index clone
  cc<-c[c$index==clones[k],]
  cc<-subset(cc,select=c(readCount,percentage,replicate,fraction,run,ownerPatient,submission))

  #combine empty and x and sum up
  data<-as_tibble(bind_rows(empty,cc))
  data<-as_tibble(ddply(data,c("replicate","fraction","run","patient","submission"),numcolwise(sum)))
  
  #split by patient
  data<-split(data,data$patient)
  names(data)
  
  pl<-list()
  #values 'data': "bella"  "bentle" "cheech" "daisy"  "eddie"  "gypsy"  "jake"   "marish" "reilly" "sona" 
  for (l in 1:length(data)){
    pl[[l]]<-plotTile(data[[l]],l)
  }
  temp<-data[[l]]
  temp[temp$patient=='bella',]
  
  library(ggpubr)
  pdf(paste0(outPathTile,"tilePlot_clone-",clones[k],".pdf"))
  print(ggarrange(plotlist=pl,ncol=1))
  dev.off()
}


#================================================
#========= generate bar plots ==================
#================================================

#+++++++++++++++++++++ function +++++++++++++++++++++
plotCount<-function(x,y,z){
  ggplot(x,aes(submission,readCount,fill=replicate))+
    geom_col(position=position_dodge(preserve = "single"))+
    facet_grid(index~fraction)+
    scale_fill_brewer(palette="Paired")+
    scale_x_discrete(limits=y,labels=z)+
    xlab("Time till relapse [weeks]")+ylab("Read count neoplastic clone")
}
plotCount_freeScale<-function(x,y,z){
  ggplot(x,aes(submission,readCount,fill=replicate))+
    geom_col(position=position_dodge(preserve = "single"))+
    facet_grid(index~fraction,scale='free')+
    scale_fill_brewer(palette="Paired")+
    scale_x_discrete(limits=y,labels=z)+
    xlab("Time till relapse [weeks]")+ylab("Read count neoplastic clone")
}
plotPercentage<-function(x,y,z){
  ggplot(x,aes(submission,percentage,fill=replicate))+
    geom_col(position=position_dodge(preserve = "single"))+
    facet_grid(index~fraction)+
    scale_fill_brewer(palette="Paired")+
    scale_x_discrete(limits=y,labels=z)+
    xlab("Time till relapse [weeks]")+ylab("Percentage of neoplastic clone [%]")
}
plotPercentage_freeScale<-function(x,y,z){
  ggplot(x,aes(submission,percentage,fill=replicate))+
    geom_col(position=position_dodge(preserve = "single"))+
    facet_grid(index~fraction,scale='free')+
    scale_fill_brewer(palette="Paired")+
    scale_x_discrete(limits=y,labels=z)+
    xlab("Time till relapse [weeks]")+ylab("Percentage of neoplastic clone [%]")
}
#===================
submissions<-c("bella","bentle","cheech","daisy","eddie","eddie","gypsy","jake","marish","sona")
limits<-list(
  paste0("bella-0",seq(0,6,1)),
  paste0("bentle-0",seq(0,2,1)),
  paste0("cheech-0",seq(0,5,1)),
  c(paste0("daisy-0",seq(0,9,1)),paste0("daisy-",seq(10,11,1))),
  paste0("eddie-0",seq(0,7,1)),
  paste0("eddie-0",seq(0,7,1)),
  c(paste0("gypsy-0",seq(0,9,1)),"gypsy-10"),
  c(paste0("jake-0",seq(0,9,1)),paste0("jake-",seq(10,12,1))),
  paste0("marish-0",seq(0,8,1)),
  paste0("sona-0",seq(0,5,1))
)
labels<-list(
  c(17,16,13,10,7,4,0),    #bella           
  c(4,2,0),                #bentley
  c(14,12,9,7,4,0),        #cheech
  c(38,35,34,32,30,27,23,18,14,12,4,0), #daisy
  c(18,17,15,12,9,6,3,0),               #eddie1
  c(18,17,15,12,9,6,3,0),               #eddie2
  c(31,29,26,24,21,17,13,9,7,2,0),      #gypsy
  c(49,47,43,41,38,34,30,26,24,12,6,1,0),  #jake
  c(28,26,23,21,18,14,8,6,0),           #marishka
  c(8,7,6,5,2,0)                      #sona
  
)

i<-indexSummary
i<-splitFilename(i,sampleNoCodesForFraction)

for (j in 1:length(submissions)){
  ii<-i[grepl(submissions[j],i$submission) & i$index!='no',]
  
  pdf(paste0(outPathBar,submissions[j],"_count.pdf"))
  print(plotCount(ii,limits[[j]],labels[[j]]))
  dev.off()
  
  pdf(paste0(outPathBar,submissions[j],"_countFreeScale.pdf"))
  print(plotCount_freeScale(ii,limits[[j]],labels[[j]]))
  dev.off()
  
  pdf(paste0(outPathBar,submissions[j],"_percent.pdf"))
  print(plotPercentage(ii,limits[[j]],labels[[j]]))
  dev.off()
  
  pdf(paste0(outPathBar,submissions[j],"_percentFreeScale.pdf"))
  print(plotPercentage_freeScale(ii,limits[[j]],labels[[j]]))
  dev.off()
}

#print to xlsx
submissions<-unique(submissions)
wb<-createWorkbook()
for (j in 1:length(submissions)){
  temp<-i[i$patient==submissions[j],]
  temp<-temp[temp$index!='no',]
  addWorksheet(wb,submissions[j])
  writeData(wb,submissions[j],temp)
}
saveWorkbook(wb,'../Results/cloneAbundanceByPatient.xlsx',overwrite=T)

#still to do:
#jake
#13 timpoints in timeline sheet, only 12 sequenced -> no relapse sample; last sample is 1w before relapse
#not included in loop version:
#  +annotate(geom="text",x=13,y=15,label="No sample",angle=90)

#reilly - including post relapse
#13 samples sequenced (00-12); samples to & including replapse: 10 => 3 additional time points post relapse
labels<-c(32,30,27,25,22,16,12,8,6,0,-2,-5,-7)
limits<-c(paste0("reilly-0",seq(0,9,1)),paste0("reilly-",seq(10,12,1)))
pdf("../Results/Clones/reilly_inclPostRelapse.pdf")
ggplot(i.reilly,aes(submission,percentage,fill=replicate))+
  geom_col(position="dodge")+
  facet_wrap(~fraction,ncol=1)+
  scale_fill_brewer(palette="Paired")+
  scale_x_discrete(limits=limits,labels=labels)+
  xlab("Time till relapse [weeks]")+ylab("Reads of neoplastic clone [%]")
dev.off()

#reilly
#13 samples sequenced (00-12); samples to & including replapse: 10 => 3 additional time points post relapse
labels<-c(32,30,27,25,22,16,12,8,6,0)
limits<-paste0("reilly-0",seq(0,9,1))
pdf("../Results/Clones/reilly.pdf")
ggplot(i.reilly,aes(submission,percentage,fill=replicate))+
  geom_col(position="dodge")+
  facet_wrap(~fraction,ncol=1)+
  scale_fill_brewer(palette="Paired")+
  scale_x_discrete(limits=limits,labels=labels)+
  xlab("Time till relapse [weeks]")+ylab("Reads of neoplastic clone [%]")
dev.off()

#=========================================================================================
#================ sample cross-contamination - read count vs. percentage =================
#=========================================================================================
i<-indexSummary
i<-splitFilename(i,sampleNoCodesForFraction)
#remove lines that are not index sequences
i<-i[i$index!='no',]

i$indexId<-substring(i$index,1,4)
i$dogId<-substring(i$submission,1,4)

i$indexType<-ifelse(i$indexId==i$dogId,'mrd','contamination')

#======== determine which seqs occur in both replicates =====
#create additional col with fused index, submission and fraction
i<-i %>% unite(index_submission_fraction,index,submission,fraction,remove=F)
i$index_submission_fraction<-as.factor(i$index_submission_fraction)
#determine frequency
n_occur<-data.frame(table(i$index_submission_fraction))
#make vector with all that occur >1 time
inBothReplicates<-n_occur$Var1[n_occur$Freq>1]
#create new column in i
i$inBothReplicates<-F
i$inBothReplicates[i$index_submission_fraction %in% inBothReplicates]<-T

#======== plot =======
point<-function(x){
  ggplot(x,aes(readCount,percentage,fill=inBothReplicates))+geom_point(pch=21)
}
pointLog<-function(x){
  ggplot(x,aes(log10(readCount),log10(percentage),fill=inBothReplicates))+geom_point(pch=21)
}
plotDensity<-function(x,y){
  density<-kde2d(log10(x$readCount),log10(x$percentage),n=50)
  filled.contour(density,color.palette=colorRampPalette(y))
}
#define palettes
palette.mrd<-c('darkgrey','lightgrey','white','blue','darkblue')
palette.contamination<-c('darkgrey','lightgrey','white','lightyellow','orange','red','darkred')

#plot all
point(i)
pointLog(i)
plotDensity(i.mrd,palette.mrd)


#subset by mrd
i.mrd<-i[i$indexType=='mrd',]
ggplot(i.mrd,aes(inBothReplicates,percentage))+geom_boxplot()
point(i.mrd)
pointLog(i.mrd)
plotDensity(i.mrd,palette.mrd)
#subset by mrd + in both
i.mrd.inBoth<-i.mrd[i.mrd$inBothReplicates==T,]
point(i.mrd.inBoth)
pointLog(i.mrd.inBoth)
plotDensity(i.mrd.inBoth,palette.mrd)
#subset by mrd + in one
i.mrd.inOne<-i.mrd[i.mrd$inBothReplicates==F,]
point(i.mrd.inOne)
pointLog(i.mrd.inOne)
plotDensity(i.mrd.inOne,palette.mrd)

#subset by contamination
i.contamination<-i[i$indexType=='contamination',]
point(i.contamination)
pointLog(i.contamination)
plotDensity(i.contamination,palette.contamination)
#subset by contamination + in both
i.contamination.inBoth<-i.contamination[i.contamination$inBothReplicates==T,]
point(i.contamination.inBoth)
pointLog(i.contamination.inBoth)
plotDensity(i.contamination.inBoth,palette.contamination)
#subset by contamination + in one
i.contamination.inOne<-i.contamination[i.contamination$inBothReplicates==F,]
point(i.contamination.inOne)
pointLog(i.contamination.inOne)
plotDensity(i.contamination.inOne,palette.contamination)


#subset by replicate
i.both<-i[i$inBothReplicates==T,] %>% 
  group_by(index_submission_fraction) %>%
  mutate(percentage.mean=mean(percentage),percentage.sd=sd(percentage),readCount.mean=mean(readCount),n=n())

ggplot(i.both,aes(percentage.mean,percentage.sd))+geom_point()

#readCount vs. percentage - no facet
ggplot(i.both,aes(log10(readCount),log10(percentage),fill=indexType))+geom_point(pch=21)+geom_density_2d(aes(color=indexType))

#readCount vs. percentage
pdf(paste0(outPathDensity,"density_facetConVsMrd.pdf"))
ggplot(i.both,aes(log10(readCount),log10(percentage),fill=indexType))+facet_wrap(~indexType)+geom_density_2d(aes(color=indexType),bins=40)+geom_point(pch=21)
dev.off()

#readCount vs. percentage - SD as size
ggplot(i.both,aes(log10(readCount),log10(percentage),fill=indexType,size=percentage.sd))+geom_point(pch=21)+facet_wrap(~indexType)+geom_density_2d(aes(color=indexType))

#readCount vs. percentage - with line connecting replicates
ggplot(i.both,aes(log10(readCount),log10(percentage),fill=indexType))+geom_point(pch=21)+facet_wrap(~indexType)+geom_density_2d(aes(color=indexType))+geom_line(group=i.both$index_submission_fraction)

#focus on outliers in contaminated
problemCases<-i.both[log10(i.both$readCount)>log10(100),] #& i.both$indexType=='contamination'
pcVector<-problemCases$index_submission_fraction
problemCasesBoth<-i.both[i.both$index_submission_fraction %in% pcVector,]

ggplot(problemCasesBoth,aes(log10(readCount),log10(percentage),fill=indexType))+geom_point(pch=21)+facet_wrap(~indexType)+geom_density_2d(aes(color=indexType))+geom_line(group=problemCasesBoth$index_submission_fraction)+geom_text(aes(label=problemCasesBoth$index_submission_fraction),hjust=0, vjust=0)

problemCasesBothCont<-ungroup(problemCasesBothCont)
ggplot(problemCasesBothCont,aes(log10(readCount),log10(percentage),fill=indexType))+
  geom_point(pch=21)+
  facet_wrap(~indexType)+
  #geom_line(group=problemCasesBothCont$index_submission_fraction)+
  geom_text(aes(label=problemCasesBothCont$index_submission_fraction),hjust=0, vjust=0,size=5)

problemCasesBothCont$index_submission_fraction<-as.factor(problemCasesBothCont$index_submission_fraction)
str(problemCasesBothCont)
problemCasesBothCont<-problemCasesBoth[problemCasesBoth$indexType=='contamination',]


#================== hexagonal binning ====================
hexbin<-function(x){
  ggplot(x,aes(readCount,percentage))+geom_hex()+
  labs(title="Index clone cross-contamination",x="Number of reads [n]",y="Percentage of reads [%]")
}

pdf("../Results/SampleCrossContamination/contamination_hexbin.pdf")
hexbin(ii)
dev.off()

pdf("../Results/SampleCrossContamination/contamination_hexbin_readCountUnder1000.pdf")
hexbin(ii[ii$readCount<1000,])
dev.off()

#plot as points with labels
ii$label<-paste0(ii$index," in ",ii$submission," ",ii$fraction," ",ii$replicate)

pdf("../Results/SampleCrossContamination/contamination_pointwithLabel.pdf")
ggplot(ii,aes(readCount,percentage))+geom_point()+
  labs(title="Index clone cross-contamination",x="Read coverage",y="Percentage of reads [%]")+geom_text(aes(label=ii$label),hjust=0,vjust=0,readCount=2)+xlim(0,9000)
dev.off()

pdf("../Results/SampleCrossContamination/contamination_pointwithLabel_readCountUnder1000.pdf")
ggplot(ii[ii$readCount<1000,],aes(readCount,percentage))+geom_point()+
  labs(title="Index clone cross-contamination",x="Read coverage",y="Percentage of reads [%]")+geom_text(aes(label=ii[ii$readCount<1000,]$label),hjust=0,vjust=0,readCount=2)+xlim(0,600)
dev.off()

#marginal plots don't work with hexagonal binning because they only consider the number of bins and not the number of dots within a bin



#================== point plot ====================
ii$run<-as.factor(ii$run)
point<-function(x){
  ggplot(x,aes(readCount,percentage))+geom_point(pch=21,aes(fill=run))+
    labs(title="Index clone cross-contamination",x="Number of reads [n]",y="Percentage of reads [%]")
}

p2<-point(ii)

pdf("../Results/SampleCrossContamination/contamination_point_wMarginalViolin.pdf")
ggMarginal(p2+theme(legend.position="left"),type="violin",readCount=10)
dev.off()
pdf("../Results/SampleCrossContamination/contamination_point_wMarginalBoxplot.pdf")
ggMarginal(p2+theme(legend.position="left"),type="boxplot",readCount=10)
dev.off()


#summary statistics
summary(ii$readCount)
summary(ii$percentage)

temp<-subset(ii,select=c(readCount,percentage))
temp<-gather(temp)
pdf("../Results/SampleCrossContamination/contamination_boxplots_readAndPercentage.pdf")
ggplot(temp,aes(key,value))+geom_boxplot()+facet_wrap(~key,scales="free")
dev.off()

#======= plot by clone =====
i<-indexSummary
ii<-i[substring(i$index,1,4)!=substring(i$submission,1,4),]
ii<-ii[ii$index!='no',]

colorCount<-length(unique(ii$patient))
getPalette<-colorRampPalette(brewer.pal(9,"Set1"))
pdf("../Results/SampleCrossContamination/contamination_point_coverageAsreadCount.pdf")
ggplot(ii,aes(index,percentage,fill=patient))+geom_jitter(pch=21,aes(readCount=ii$readCount))+scale_readCount(range = c(0,10))+scale_fill_manual(values=getPalette(colorCount))
dev.off()

pdf("../Results/SampleCrossContamination/contamination_point_coverageAsreadCount_facetIndex.pdf")
ggplot(ii,aes(patient,percentage,fill=patient))+geom_jitter(pch=21,width=0.2,aes(readCount=ii$readCount))+scale_readCount(range = c(0,10))+facet_wrap(~index)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="top")+scale_fill_manual(values=getPalette(colorCount))
dev.off()

