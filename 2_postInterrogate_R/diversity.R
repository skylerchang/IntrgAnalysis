#script calculates diversity for split 'Clntab_RDS' files: 'Clntab_RDS_noTemplate' & 'Clntab_RDS_templateOnly'

library(tidyverse)
library(vegan)
library(plyr)
library(here)
library(reshape2)
library(openxlsx)
library(RColorBrewer)
library(ggforce)

source("2_postInterrogate_R/functions.R")

setwd(here())
getwd()

#=========== adjust =============
run<-30
sampleNoCodesForFraction<-T        #T or F

rdsPath<-'../Data/Clntab_RDS/clntab_vAndJ.rds'
resultsPathDiversity<-'../Results/Diversity/'
resultsPathLocus<-'../Results/LocusStats/'
resultsPathReads<-'../Results/ReadCounts/'
resultsPathTemplate<-'../Results/Template/'

#===============================
dir.create(resultsPathDiversity,recursive=T)
dir.create(resultsPathLocus,recursive=T)
dir.create(resultsPathReads,recursive=T)
dir.create(resultsPathTemplate,recursive=T)

t<-read_rds(rdsPath)

#datalist contains one tibble per sample
datalist<-t[[1]]
length(datalist)
#sample names are stored separately in a vector
files_short<-t[[2]]
length(files_short)

#filter out WGA samples
includedLogical<-!grepl("WGA",files_short)
files_short<-files_short[includedLogical]
datalist<-datalist[includedLogical]

#initialize tibbles
templateSummary<-tibble()
readCountSummary<-tibble()
locusSummary<-tibble()
diversitySummary<-tibble()

for (i in 1:length(files_short)){  
  print(files_short[i])
  a<-datalist[[i]][,c("vGene","jGene","aaSeq","aaLength","size","completeNtSeq","vAndJchainSimplified")]
  colnames(a)[7]<-'locus'
  count.all<-nrow(a)
  if(nrow(a)==0){next}
  a$sample<-files_short[i]
  
  #============ split dataset into +/- templates =======
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
  
  #create template summary
  temp<-ddply(a[,c("template","size")],"template",numcolwise(sum))
  temp$filename<-files_short[i]
  templateSummary<-bind_rows(templateSummary,temp)
  
  #count template reads 
  count.templates<-nrow(a[a$template!='no',])
  
  #continue with dataset w/o templates
  a<-a[a$template=='no',]
  
  #split by locus
  igh<-a[a$locus=='IGH',]
  trb<-a[a$locus=='TRB',]
  count.igh.before<-nrow(igh)
  count.trb.before<-nrow(trb)
  
  #============ filter seqs by seqs & ^C.{3,30}[FW]$ =======
  igh<-igh[grepl("^C.{3,30}[W]$",igh$aaSeq),]
  trb<-trb[grepl("^C.{3,30}[F]$",trb$aaSeq),]
  count.igh.after<-nrow(igh)
  count.trb.after<-nrow(trb)
  count.igh.anchorFiltered<-count.igh.before-count.igh.after
  count.trb.anchorFiltered<-count.trb.before-count.trb.after
  
  #============= create read count summary =============
  #checksum
  count.all-(count.templates+
               count.igh.anchorFiltered+count.igh.after+
               count.trb.anchorFiltered+count.trb.after
             )
  #create tibble
  temp<-tibble(
    filename=files_short[i],
    total=count.all,
    templates=count.templates,
    igh=count.igh.after,
    falseAnchor.igh=count.igh.anchorFiltered,
    trb=count.trb.after,
    falseAnchor.trb=count.trb.anchorFiltered
    )
  readCountSummary<-bind_rows(readCountSummary,temp)
  
  #============ create locus summary =======
  totalReads<-count.igh.after+count.trb.after
  temp<-tibble(
    filename=files_short[i],
    totalReads=totalReads,
    igh.count=count.igh.after,
    igh.percent=count.igh.after/totalReads*100,
    trb.count=count.trb.after,
    trb.percent=count.trb.after/totalReads*100,
    locusRatio=igh.count/trb.count
    )
  locusSummary<-bind_rows(locusSummary,temp)
  
  #============ create diversity summary =======
  calculateDiversity<-function(x,y){
    #aggregate all lines with identical aaSeq
    x<-as_tibble(ddply(x[,c("aaSeq","size")],"aaSeq",numcolwise(sum)))
    #create tibble with diversity indexes
    tibble(
      filename=y,
      shannon=diversity(x$size,index = "shannon",base = exp(1)),
      effectiveSpecies=exp(shannon),
      simpson=diversity(x$size,index = "simpson"),
      readCount=sum(x$size),
      clonotypeCount=nrow(x)
      )
  }
  #calculate diversity indexes and store in 'diversity'
  temp<-calculateDiversity(igh,files_short[i])
  temp$locus<-'IGH'
  diversitySummary<-bind_rows(diversitySummary,temp)
  
  temp<-calculateDiversity(trb,files_short[i])
  temp$locus<-'TRB'
  diversitySummary<-bind_rows(diversitySummary,temp)
}

#===============================================================
#===================== read count summary ======================
#===============================================================
#print xlsx
wb<-createWorkbook()
addWorksheet(wb,"summary")
writeData(wb,"summary",readCountSummary)
saveWorkbook(wb,"../Results/ReadCounts/readCountSummary.xlsx",overwrite=T)

readCount.long<-gather(subset(readCountSummary,select=-(total)),variable,value,-filename)
pdf('../Results/ReadCounts/readCountSummary.pdf')
ggplot(readCount.long,aes(filename,value,fill=variable))+geom_col()+coord_flip()+scale_fill_brewer(palette='Dark2')
dev.off()

#===============================================================
#========================== template summary ===================
#===============================================================
templateSummary
#plot stratified by file
pdf(paste0(resultsPathTemplate,'templatesBytemplate.pdf'))
ggplot(templateSummary,aes(template,size,fill=filename))+geom_col()+coord_flip()+theme(legend.position="none")
dev.off()

#exclude 'no' template
templateSummary2<-templateSummary[templateSummary$template!='no' & templateSummary$template!='t2',]

#plot stratified by template - bar plot
pdf(paste0(resultsPathTemplate,'templatesBySample_barplot.pdf'))
ggplot(templateSummary,aes(filename,size,fill=template))+geom_col()+coord_flip()
dev.off()
#plot stratified by template - line plot
pdf(paste0(resultsPathTemplate,'templatesBySample_lineplot.pdf'))
ggplot(templateSummary2,aes(filename,size,group=template))+geom_point()+geom_line(aes(color=template))+coord_flip()
dev.off()

#===============================================================
#========================== locus summary ======================
#===============================================================
l<-locusSummary

#function splits file name into components creating a new col for each; starting format:
#<id>_<ownerPatient>_<idNumber>
#strsplit by '_' -> 3 parts: 1) id, 2) ownerPatient, 3) idNumber
#1) id is then broken down further into: pcr, sample, submission, replicate, fraction
splitFilename<-function(x){
  #transform old control names to required format
  x$filename<-x$filename %>%
    sub("k9-pc-st_D16-07","D16-07-st",.) %>% 
    sub("k9-pc_D16-07","D16-07",.) %>% 
    sub("k9-nt-st_H2O","ntc",.) %>% 
    sub("k9-nt_H2O","ntc-st",.)
  splitFilename<-strsplit(x$filename,"_")
  x$id<-laply(splitFilename, '[[', 1)
  x$ownerPatient<-laply(splitFilename, '[[', 2)
  x$idNumber<-laply(splitFilename, '[[', 3)
  x$pcr<-sub("C[0-9]+P[0-9]+C[0-9]+$","",x$id)
  x$sample<-sub("-D[0-9]+P[0-9]+$","",x$pcr)
  x$submission<-sub("-[0-9a-zA-Z]+$","",x$sample)
  x$replicate<-ifelse(grepl("D[0-9]+P[0-9]*[13579]$",x$pcr),"rep1","rep2")
  if(sampleNoCodesForFraction==T){
    x$fraction<-ifelse(grepl("-1$",sample),"fraction1","fraction2")
  }else{
    x$fraction<-rep("fraction1",length(sample))
  }
  x
}
l<-splitFilename(l)

#=============== save as xlsx & rds =================
#save as rds
saveRDS(l,paste0(resultsPathLocus,"locusSummary.rds"))

#save as xlsx
wb<-createWorkbook()
addWorksheet(wb,"summary")
writeData(wb,"summary",l)
saveWorkbook(wb,paste0(resultsPathLocus,"locusSummary.xlsx"),overwrite = T)

#=============== plot locus count data =================
#subset count data, reshape to long format, split file name
l.count.long<-subset(l,select=c(filename,igh.count,trb.count)) %>%
  gather(locus,value,-filename) %>%
  splitFilename()
#create 4 replicates by incorporating the locus information -> plot loci in different colors
l.count.long$replicate2<-paste(l.count.long$locus,l.count.long$replicate,sep='-')

pdf(paste0(resultsPathLocus,"ighVsTrb_withReplicates_barplot.pdf"))
ggplot(l.count.long,aes(locus,value,fill=replicate2))+geom_col(position=position_dodge())+facet_wrap(~sample)+scale_fill_brewer(palette="Paired")
dev.off()

#================= plot locus ratio data =================
l.splitFilename<-splitFilename(l)

#total reads vs locus ratio
pdf(paste0(resultsPathLocus,"readsVsLocusRatio_all.pdf"))
ggplot(l.splitFilename,aes(totalReads,locusRatio,fill=sample))+geom_point(pch=21)+theme(legend.position="top")
dev.off()

#total reads vs locus ratio - with label
pdf(paste0(resultsPathLocus,"readsVsLocusRatio_all_withLabel.pdf"))
ggplot(l.splitFilename,aes(totalReads,log(locusRatio),fill=sample))+geom_point(pch=21)+geom_label(aes(label=sample))+theme(legend.position="none")
dev.off()

pdf(paste0(resultsPathLocus,"readsVsSample_all_boxplot.pdf"))
ggplot(l.splitFilename,aes(totalReads,log(locusRatio),group=sample))+geom_boxplot(aes(fill=sample))
dev.off()

pdf(paste0(resultsPathLocus,"readsVslocusRatio_all_smoothedCondMeans.pdf"))
ggplot(l.splitFilename,aes(totalReads,log(locusRatio)))+geom_point()+geom_smooth()
dev.off()

#remove outliers 
#exclude<-c("k9-nt-H2O-0","k9-nt-st-H2O-1","09-020952-16")
#l.splitFilename<-l.splitFilename[!l.splitFilename$sample %in% exclude,]

#===============================================================
#======================= diversity summary =====================
#===============================================================
d<-diversitySummary

d.split1<-splitFilename(d)
d.split2<-split(d.split1,d.split1$sample)
str(d.split2)

#print to xlsx
wb<-createWorkbook()
addWorksheet(wb,"summary")
writeData(wb,"summary",d.split1)
for (i in 1:length(d.split2)){
  addWorksheet(wb,d.split2[[i]]$sample[1])
  writeData(wb,d.split2[[i]]$sample[1],d.split2[i],colNames = TRUE)
}
saveWorkbook(wb,paste0(resultsPathDiversity,"diversitySummary.xlsx"),overwrite = T)

#====== transpose igh/trb lines to columns -> one pcrId per line =============
d2<-d %>%
  gather(key, value, -filename, -locus) %>%
  unite(indexByLocus, key, locus) %>%
  spread(indexByLocus,value)

d2<-splitFilename(d2)

#save as rds file
saveRDS(d2,paste0(resultsPathDiversity,"diversitySummary.rds"))

#transform to long format
d.long<-gather(d,index,value,-filename,-locus,-readCount,-clonotypeCount)
d.long<-splitFilename(d.long)

d.long$ownerPatient<-as.factor(d.long$ownerPatient)
colorCount<-nlevels(d.long$ownerPatient)
getPalette<-colorRampPalette(brewer.pal(9,'Set1'))

#jitter all indexes
pdf(paste0(resultsPathDiversity,"diversityPlot_allIndexes.pdf"))
ggplot(d.long,aes(locus,value))+
  geom_jitter(aes(shape=locus,color=ownerPatient),size=3)+
  facet_wrap(~index, scales="free")+
  scale_color_manual(values=getPalette(colorCount))
dev.off()

#point - by size - all indexes 
pdf(paste0(resultsPathDiversity,"readCountVsDiversity_facetLocusIndex.pdf"))
ggplot(d.long,aes(readCount,value))+
  geom_point(mapping=aes(shape=locus,color=ownerPatient),size=3)+
  facet_grid(index~locus,scales="free")+
  scale_color_manual(values=getPalette(colorCount))
dev.off()


#================= merge replicates ================
d2<-splitFilename(d)
d2<-subset(d2,select=c(submission,sample,locus,replicate,shannon,effectiveSpecies,simpson))

d2<-d2 %>% gather(index,value,-submission,-sample,-locus,-replicate) %>%
  unite(index2,index,replicate) %>%
  spread(index2,value)

d2$shannon.mean<-(d2$shannon_rep1+d2$shannon_rep2)/2
d2$simpson.mean<-(d2$simpson_rep1+d2$simpson_rep2)/2
d2$effectiveSpecies.mean<-(d2$effectiveSpecies_rep1+d2$effectiveSpecies_rep2)/2

d2<-subset(d2,select=c(submission,sample,locus,shannon.mean,simpson.mean,effectiveSpecies.mean))

d2<-d2 %>% gather(index,value,-submission,-sample,-locus) %>%
  unite(index3,index,locus) %>%
  spread(index3,value)

#save as xlsx
wb<-createWorkbook()
addWorksheet(wb,'summary')
writeData(wb,'summary',d2)
saveWorkbook(wb,paste0(resultsPathDiversity,'diversity_withMeans.xlsx'),overwrite = T)

#save as rds
write_rds(d2,paste0(resultsPathDiversity,'diversity_withMeans.rds'))
