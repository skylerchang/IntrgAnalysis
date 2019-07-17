#script calculates diversity for split 'Clntab_RDS' files: 'Clntab_RDS_noTemplate' & 'Clntab_RDS_templateOnly'
#v14: process data from multiple sequencing runs

library(tidyverse)
library(vegan)
library(plyr)
library(here)
library(reshape2)
library(RColorBrewer)
library(ggpubr)
library(openxlsx)
library(ggExtra)

setwd(here())
getwd()

#=========== adjust =============

sampleNoCodesForFraction<-T        #T or F

rdsPath<-'../Data/Clntab_RDS/clntab_vAndJ.rds'
rdsPaths<-c(
  '../../Run13_k9MRD-Reilly-AllTimePoints/Data/Clntab_RDS/clntab_vAndJ.rds',
  '../../Run17_Bella_Marishka/Data/Clntab_RDS/clntab_vAndJ.rds',
  '../../Run22_MRD_Daisy_Sona_UntilRelapse/Data/Clntab_RDS/clntab_vAndJ.rds',
  '../../Run27_MRD_Eddie_Jake_UntilRelapse/Data/Clntab_RDS/clntab_vAndJ.rds',
  '../../Run29_MRD_Bentley_Cheech_Gypsy_UntilRelapse/Data/Clntab_RDS/clntab_vAndJ.rds'
  )
runs<-c(13,17,22,27,29)

loci<-c("IGH")
#===================================
indexSummary<-tibble()
templateSummary<-tibble()
diversity<-tibble()

for (k in 1:length(rdsPaths)){
  print(rdsPaths[k])
  
  ifelse(k==1,thisIsReilly<-T,thisIsReilly<-F)
  
  t<-read_rds(rdsPaths[k])
  
  #datalist contains one tibble per sample
  datalist<-t[[1]]
  print(paste0("Number of files: ",length(datalist)))
  #sample names are stored separately in a vector
  files_short<-t[[2]]
  files_short
  
  for (i in 1:length(datalist)){
    a<-datalist[[i]][,c("vGene","jGene","aaSeq","aaLength","size","completeNtSeq","vAndJchainSimplified")]
    if(is.null(a)){next}
    if(nrow(a)==0){next}
    #============ standardize sample name =======
    #assume the following default format:
    # [1] "16-088089-3D1P3C1P1C1_Mcdermott-Gibbie"      
    # [2] "16-088089-3D1P4C1P1C1_Mcdermott-Gibbie"   
    #deviations:
    #[92] "k9-pc_D16-07-2FD4P59C1P1C1_OWN33-PAT1_S92"      
    #[93] "k9-pc-st_D16-07-2FD1P32C1P1C1_OWN33-PAT1_S93"  
    #<sampleId>_<sampleName>
    file<-sub("_S[0-9]+$","",files_short[i])
    file<-sub("ntc","NTC",file)
    #if this is reilly
    if (thisIsReilly==T){
      file<-sub("[0-9]+_","",file)
      file<-paste0(file,"_McCormick-Reilly")
    }
    #tranform samples to fit default format -> only one '_' which needs to be between sampleId & sampleName
    #run13: "35_NTC-D1P122_S35","42_NTC-D1P123_S42","49_D16-D1P133_S49","56_D16-D1P134_S56"
    #run17: "D16-07-1FD1P267C1P1C1_OWN33-PAT1_S71","D16-07-1FD1P268C1P1C1_OWN33-PAT1_S72"
    #       "NTC-can-0D1P182C1P1C1_NTC-NTC-canine_S69","NTC-can-0D1P183C1P1C1_NTC-NTC-canine_S70"
    #run22: "D16-07-1FD1P337C1P1C1","D16-07-1FD1P338C1P1C1","D16-07-ST-1FD1P335C1P1C1","D16-07-ST-1FD1P336C1P1C1"
    file<-sub("NTC","ntc-ntc",file)
    file<-sub("D16","pc-pc",file)
    file<-sub("k9-pc-st_D16-07","pcST-pcST",file)
    file<-sub("k9-pc_D16-07","pc-pc",file)
    file<-sub("k9-nt-st_H2O","ntc-ntc",file)
    file<-sub("k9-nt_H2O","ntcST-ntcST",file)
    sampleId<-sub("_.*$","",file)
    sampleName<-sub("^.*_","",file)
    
    if (thisIsReilly==T){
      sampleId<-paste0(sampleId,"C1P1C1")
      file<-paste(sampleId,sampleName,sep="_")
    }
    #============ extract replicate & sample =======
    replicate<-ifelse(grepl("D[0-9]+P[0-9]*[02468]C",sampleId),"rep1","rep2")
    if(sampleNoCodesForFraction==T){
      fraction<-ifelse(grepl("-1D[0-9]+P[0-9]+C[0-9]+P[0-9]+C[0-9]+$",sampleId),"PBMCs","Plasma")
    }else{
      fraction<-"fraction1"
    }
    sampleIdShort<-sub("D[0-9]+P[0-9]+C[0-9]+P[0-9]+C[0-9]+$","",sampleId)
    submission<-sub("-[0-9a-zA-Z]+$","",sampleIdShort)
    patient<-sub("-.*","",submission)
    
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
    
    #============ check for templates =======
    a$template<-'no'
    #Akash‘s templates
    a$template[grepl("TGTGCATCACGACACAGTGGTCTGG",a$completeNtSeq)]<-'t1'
    a$template[grepl("TGTGCATCACGACCAGATCCACAGATCCATTGGTTACTGG",a$completeNtSeq)]<-'t2'
    #Tamara‘s templates
    a$template[grepl("TGTTCGCCTTATCGCCTTATGG",a$completeNtSeq)]<-'IGHV1-30_IGHJ4'
    a$template[grepl("TGTCTAGTACGCCTCTCTGCCTCTCTGCTAGTACGTGG",a$completeNtSeq)]<-'IGHV1-30_IGHJ6'
    a$template[grepl("TGTGCTTCTGCCTTTCTGCCTGCTCAGGATTCTGCCTGCTCAGGAGCTCAGGATTCTCTGG",a$completeNtSeq)]<-'IGHV3-1_IGHJ4'
    a$template[grepl("TGTGAGGAGTCCGTAGAGAGAGGAGTCCAGCGTAGCCATGCCTAAGGAGTCCCAGCCTCGGTAGAGAGAGCGCTGG",a$completeNtSeq)]<-'IGHV3-1_IGHJ6'
    a$template[grepl("TGTAAGGAGTAACTGCATAACTGCATACTAAGCCTAAGGAGTATGG",a$completeNtSeq)]<-'IGHV4-1_IGHJ4'
    a$template[grepl("TGTGTGCCTCTTTCCTCTACTAGATCGCCTCTCTATTATCCTCTAGAGTAGAGTAAGGAGTAGATCGCTATCCTCTGTAAGGAGTCCTCTACCTGG",a$completeNtSeq)]<-'IGHV4-1_IGHJ6'
    
    #split in separate tibbles
    a.withoutTemplate<-a[a$template=='no',]
    a.templateOnly<-a[a$template!='no',]
    
    #========== summarize index sequences ==============
    temp<-ddply(a.withoutTemplate[,c("size","index")],"index",na.rm=T, numcolwise(sum))
    if (nrow(temp)==0){next}
    temp$totalReads<-sum(temp$size)
    temp$percentage<-temp$size/temp$totalReads*100
    temp$submission<-submission
    temp$replicate<-replicate
    temp$fraction<-fraction
    temp$run<-runs[k]
    temp$patient<-patient
    indexSummary<-bind_rows(indexSummary,temp)
    
    #============ summarize templates for sample and store in tibble =======
    temp<-ddply(a.templateOnly,"template",numcolwise(sum))
    if(nrow(temp)==0){next}
    temp<-as_tibble(temp[,c("template","size")])
    temp$sample<-files_short[i]
    temp$totalReads<-sum(temp$size)
    temp$percentage<-temp$size/temp$totalReads*100
    temp$submission<-submission
    temp$replicate<-replicate
    temp$fraction<-fraction
    temp$run<-runs[k]
    temp$patient<-patient
    templateSummary<-bind_rows(templateSummary,temp)
  }
}
c<-indexSummary

#remove controls for now -> revisit
remove<-c("ntc","k9","pc")
c<-c[!grepl(paste0(remove,collapse='|'),c$submission),]
#remove non-index clones
c<-c[c$index!='no',]

#create col 'submissionNo'
c$submissionNo<-sub(".*-","",c$submission)

#========= templates - quick peak ==================
ggplot(templateSummary,aes(template,size,fill=patient))+geom_col()+coord_flip()+theme(legend.position="none")
ggplot(templateSummary,aes(fraction,size,fill=template))+geom_col()+scale_fill_brewer(palette="Dark2")+coord_flip()+theme(legend.position="top")+facet_wrap(~patient)


#=============================================================
#========= index clone analysis - tile plots =================
#=============================================================

#++++++++++++ functions +++++++++++++++++
#plot tiles - facet: patient~fraction
#x: data, y: patient index; z: size or percentage (implement)
plotTile<-function(x,y){
  p<-ggplot(x,aes(submissionNo,replicate,fill=size))+
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
#size,percentage,replicate,fraction,run,patient,submissionNo
patient<-c("bella","bentle","cheech","daisy","eddie","gypsy","jake","marish","sona","reilly")
timepoint<-c(7,3,6,12,8,11,13,9,6,10)   #the number of submissions per patient
run<-c(17,29,29,22,27,29,27,17,22,13)   #the run number in which it was sequenced
empty<-tibble()
for (i in 1:length(timepoint)){
  temp<-tibble(
    size=rep(0,(4*timepoint[i])),
    percentage=rep(0,(4*timepoint[i])),
    replicate=rep(c("rep1","rep2"),(2*timepoint[i])),
    fraction=rep(rep(c("PBMCs","Plasma"),each=2),timepoint[i]),
    run=rep(run[i],4*timepoint[i]),
    patient=rep(patient[i],4*timepoint[i]),
    submissionNo=formatC(rep(1:timepoint[i],each=4),width=2,flag="0")
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
  cc<-subset(cc,select=c(size,percentage,replicate,fraction,run,patient,submissionNo))

  #combine empty and x and sum up
  data<-as_tibble(bind_rows(empty,cc))
  data<-as_tibble(ddply(data,c("replicate","fraction","run","patient","submissionNo"),numcolwise(sum)))
  
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
  pdf(paste0("../Results/Clones/TilePlots/tilePlot_clone-",clones[k],".pdf"))
  print(ggarrange(plotlist=pl,ncol=1))
  dev.off()
}


#================================================
#========= generate bar plots ==================
#================================================

#+++++++++++++++++++++ function +++++++++++++++++++++
plotCount<-function(x,y,z){
  ggplot(x,aes(submission,size,fill=replicate))+
    geom_col(position=position_dodge(preserve = "single"))+
    facet_grid(index~fraction)+
    scale_fill_brewer(palette="Paired")+
    scale_x_discrete(limits=y,labels=z)+
    xlab("Time till relapse [weeks]")+ylab("Read count neoplastic clone")
}
plotCount_freeScale<-function(x,y,z){
  ggplot(x,aes(submission,size,fill=replicate))+
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

for (j in 1:length(submissions)){
  ii<-i[grepl(submissions[j],i$submission) & i$index!='no',]
  
  pdf(paste0("../Results/Clones/",submissions[j],"_count.pdf"))
  print(plotCount(ii,limits[[j]],labels[[j]]))
  dev.off()
  
  pdf(paste0("../Results/Clones/",submissions[j],"_countFreeScale.pdf"))
  print(plotCount_freeScale(ii,limits[[j]],labels[[j]]))
  dev.off()
  
  pdf(paste0("../Results/Clones/",submissions[j],"_percent.pdf"))
  print(plotPercentage(ii,limits[[j]],labels[[j]]))
  dev.off()
  
  pdf(paste0("../Results/Clones/",submissions[j],"_percentFreeScale.pdf"))
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


#================ sample cross-contamination - read count vs. percentage =================
i<-indexSummary
#filter for lines where $index != $submission => cross-contamination
ii<-i[substring(i$index,1,4)!=substring(i$submission,1,4),]
ii<-ii[ii$index!='no',]

#================== hexagonal binning ====================
hexbin<-function(x){
  ggplot(x,aes(size,percentage))+geom_hex()+
  labs(title="Index clone cross-contamination",x="Number of reads [n]",y="Percentage of reads [%]")
}

pdf("../Results/SampleCrossContamination/contamination_hexbin.pdf")
hexbin(ii)
dev.off()

pdf("../Results/SampleCrossContamination/contamination_hexbin_sizeUnder1000.pdf")
hexbin(ii[ii$size<1000,])
dev.off()

#plot as points with labels
ii$label<-paste0(ii$index," in ",ii$submission," ",ii$fraction," ",ii$replicate)

pdf("../Results/SampleCrossContamination/contamination_pointwithLabel.pdf")
ggplot(ii,aes(size,percentage))+geom_point()+
  labs(title="Index clone cross-contamination",x="Read coverage",y="Percentage of reads [%]")+geom_text(aes(label=ii$label),hjust=0,vjust=0,size=2)+xlim(0,9000)
dev.off()

pdf("../Results/SampleCrossContamination/contamination_pointwithLabel_sizeUnder1000.pdf")
ggplot(ii[ii$size<1000,],aes(size,percentage))+geom_point()+
  labs(title="Index clone cross-contamination",x="Read coverage",y="Percentage of reads [%]")+geom_text(aes(label=ii[ii$size<1000,]$label),hjust=0,vjust=0,size=2)+xlim(0,600)
dev.off()

#marginal plots don't work with hexagonal binning because they only consider the number of bins and not the number of dots within a bin



#================== point plot ====================
ii$run<-as.factor(ii$run)
point<-function(x){
  ggplot(x,aes(size,percentage))+geom_point(pch=21,aes(fill=run))+
    labs(title="Index clone cross-contamination",x="Number of reads [n]",y="Percentage of reads [%]")
}

p2<-point(ii)

pdf("../Results/SampleCrossContamination/contamination_point_wMarginalViolin.pdf")
ggMarginal(p2+theme(legend.position="left"),type="violin",size=10)
dev.off()
pdf("../Results/SampleCrossContamination/contamination_point_wMarginalBoxplot.pdf")
ggMarginal(p2+theme(legend.position="left"),type="boxplot",size=10)
dev.off()


#summary statistics
summary(ii$size)
summary(ii$percentage)

temp<-subset(ii,select=c(size,percentage))
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
pdf("../Results/SampleCrossContamination/contamination_point_coverageAsSize.pdf")
ggplot(ii,aes(index,percentage,fill=patient))+geom_jitter(pch=21,aes(size=ii$size))+scale_size(range = c(0,10))+scale_fill_manual(values=getPalette(colorCount))
dev.off()

pdf("../Results/SampleCrossContamination/contamination_point_coverageAsSize_facetIndex.pdf")
ggplot(ii,aes(patient,percentage,fill=patient))+geom_jitter(pch=21,width=0.2,aes(size=ii$size))+scale_size(range = c(0,10))+facet_wrap(~index)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="top")+scale_fill_manual(values=getPalette(colorCount))
dev.off()

