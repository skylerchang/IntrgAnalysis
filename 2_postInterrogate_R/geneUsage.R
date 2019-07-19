
library(plyr)
library(tidyverse)
library(gridExtra)
library(here)
library(openxlsx)
library(RColorBrewer)

setwd(here())

#========== adjust the following variables ================
t<-read_rds('../Data/clntab_RDS/clntab_vAndJ_filtered.rds')
outpath<-'../Results/GeneUsage/'
loci<-c("IGH","TRB")

#==========================================================
dir.create(outpath)

#datalist contains one tibble per sample
datalist<-t[[1]]
#sample names are stored separately in a vector
files_short<-t[[2]]

#plotlist_histogramFacet<-list()
#plotlist_density<-list()

v<-tibble()
j<-tibble()

for (i in 1:length(datalist)){
  t<-datalist[[i]][,c("vGene","jGene","aaSeq","aaLength","readCount","locus")]
  for (k in 1:length(loci)){
    tt<-t[t$locus==loci[k],]
    tt<-tt[rowSums(is.na(tt))==0,]
    
    #================ vGene subgroup usage =================================
    tt.v<-tt
    #determine subgroup
    tt.v$vSubgroup<-NA
    temp<-gsub("-[0-9]+","",tt.v$vGene)
    temp<-strsplit(temp,split="=")
    temp<-lapply(temp,unique)
    tt.v$vSubgroup<-sapply(temp,paste,collapse="=")
    tt.v<-tt.v[!grepl("=",tt.v$vSubgroup),]

    #aggregate & sort rows with identical V subgroup
    ttt.v<-as_tibble(aggregate(tt.v$readCount,by=list(vSubgroup=tt.v$vSubgroup),FUN=sum))
    ttt.v<-ttt.v[order(-ttt.v$x),]
    ttt.v<-subset(ttt.v,!is.na(vSubgroup))
    ttt.v$vSubgroup<-factor(ttt.v$vSubgroup,levels=ttt.v$vSubgroup)
    ttt.v$x.perc=ttt.v$x/sum(ttt.v$x,na.rm=T)*100
    ttt.v$rank<-1:nrow(ttt.v)
    
    #Add to main table
    ttt.v$locus<-loci[k]
    ttt.v$animal<-files_short[i]
    v<-bind_rows(v,ttt.v)
    
    #================ jGene usage =================================
    #aggregate & sort rows with identical j gene
    tt.j<-as_tibble(aggregate(tt$readCount,by=list(jGene=tt$jGene),FUN=sum))
    tt.j<-tt.j[order(-tt.j$x),]
    tt.j<-subset(tt.j,!is.na(jGene))
    tt.j$jGene<-factor(tt.j$jGene,levels=tt.j$jGene)
    tt.j$x.perc=tt.j$x/sum(tt.j$x,na.rm=T)*100
    tt.j$rank<-1:nrow(tt.j)
    
    #Add to main table
    tt.j$locus<-loci[k]
    tt.j$animal<-files_short[i]
    j<-bind_rows(j,tt.j)
    
  }
}

#================ vGene subgroup usage =================================
v<-v[v$vSubgroup!="NA",]
v[rowSums(is.na(v))>0,]

#calculate median percentage for subgroups; used to order levels for plotting based on frequency
medianlist<-list()
for (i in 1:length(loci)){
  #subset by locus
  vv<-v[v$locus==loci[i],]
  
  #======== by sample ==========
  #convert 'long' -> 'wide' to get summary table by sample
#  vv.wide<-dcast(vv,vSubgroup~animal,value.var = "x.perc")
#  wb<-createWorkbook()
#  worksheet<-paste0("V-",loci[i],"-BySample")
#  addWorksheet(wb,worksheet)
#  writeData(wb,worksheet,vv.wide)
  
  #======== summary ==========
  #calculate summary stats
  vv.long<-as_tibble(ddply(vv,~vSubgroup,summarize,median=median(x.perc),mean=mean(x.perc),sd=sd(x.perc),min=min(x.perc),max=max(x.perc)))
  #write to xlsx file
  wb<-createWorkbook()
  worksheet<-paste0("V-",loci[i],"-Summary")
  addWorksheet(wb,worksheet)
  writeData(wb,worksheet,vv.long)
  #sort descendingly by median
  vv.long<-vv.long[order(-vv.long$median),]
  #store in list
  medianlist[[i]]<-vv.long$vSubgroup[1:10]
}
saveWorkbook(wb,"geneUsageSummary_vGene.xlsx",overwrite = T)
medianlistVector<-unlist(medianlist)

v$locus<-as.factor(v$locus)

#subset v to include the 10 most commonly used genes only (by median)
s<-v[v$vSubgroup %in% medianlistVector,]
s<-s[!s$vSubgroup=="NA",]

#split dataset to allow for individual legends per facet
ss<-split(s,f=s$locus)
str(ss)
#adjust level order in decreasing order of median
for (i in 1:length(loci)){
  ss[[i]]$vSubgroup<-factor(ss[[i]]$vSubgroup,levels=medianlist[[i]])
}
str(ss)

#boxplot
p1<-ggplot(ss$IGH,aes(vSubgroup,x.perc,fill=vSubgroup))+geom_boxplot(outlier.shape=NA)+geom_point(pch = 21,position = position_jitterdodge(jitter.width=2))+facet_wrap(~locus)+scale_fill_brewer(palette="Set3")+scale_x_discrete(name="Subgroup",limits = rev(levels(ss$TRA$vSubgroup)))+scale_y_continuous(labels = function(x) paste0(x, "%"))+labs(color='Subgroup')+coord_flip()+theme(legend.position = "none")+theme(axis.title.x = element_blank())

p2<-ggplot(ss$TRB,aes(vSubgroup,x.perc,fill=vSubgroup))+geom_boxplot(outlier.shape=NA)+geom_point(pch = 21,position = position_jitterdodge(jitter.width=2))+facet_wrap(~locus)+scale_fill_brewer(palette="Set3")+scale_x_discrete(name="Subgroup",limits = rev(levels(ss$TRB$vSubgroup)))+scale_y_continuous(labels = function(x) paste0(x, "%"))+labs(color='Subgroup')+coord_flip()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
  

#p2<-p1 %+% ss$TRB
#p3<-p1 %+% ss$TRD
#p4<-p1 %+% ss$TRG
pdf(paste0(outpath,"usage_vGene_subgroup.pdf"))
grid.arrange(p1,p2,ncol=2)
dev.off()

#================ jGene usage =================================
str(j)

wb<-createWorkbook()
#calculate median percentage for subgroups; used to order levels for plotting based on frequency
medianlist<-list()
for (i in 1:length(loci)){
  #subset by locus
  jj<-j[j$locus==loci[i],]
  
  #======== by sample ==========
  #convert 'long' -> 'wide' to get summary table by sample
#  jj.wide<-dcast(jj,jGene~animal,value.var = "x.perc")
#  worksheet<-paste0("J-",loci[i],"-BySample")
#  addWorksheet(wb,worksheet)
#  writeData(wb,worksheet,jj.wide)
  
  #======== summary ==========
  #calculate summary stats
  jj.long<-as_tibble(ddply(jj,~jGene,summarize,median=median(x.perc),mean=mean(x.perc),sd=sd(x.perc),min=min(x.perc),max=max(x.perc)))
  #write to xlsx file
  worksheet<-paste0("J-",loci[i],"-Summary")
  addWorksheet(wb,worksheet)
  writeData(wb,worksheet,jj.long)
  #sort descendingly by median
  jj.long<-jj.long[order(-jj.long$median),]
  #store in list
  medianlist[[i]]<-jj.long$jGene[1:10]
}
saveWorkbook(wb,"geneUsageSummary_jGene.xlsx",overwrite = T)

medianlistVector<-unlist(medianlist)

j$locus<-as.factor(j$locus)

#subset v to include the 10 most commonly used genes only (by median)
s<-j[j$jGene %in% medianlistVector,]
s<-s[!s$jGene=="NA",]


#split dataset to allow for individual legends per facet
ss<-split(s,f=s$locus)
#adjust level order in decreasing order of median
for (i in 1:length(loci)){
  ss[[i]]$jGene<-factor(ss[[i]]$jGene,levels=medianlist[[i]])
}
str(ss)

#boxplot
p1<-ggplot(ss$IGH,aes(jGene,x.perc,fill=jGene))+geom_boxplot(outlier.shape=NA)+geom_point(pch = 21,position = position_jitterdodge(jitter.width=2))+facet_wrap(~locus)+scale_fill_brewer(palette="Set3")+scale_x_discrete(name="Gene",limits = rev(levels(ss$TRA$jGene)))+scale_y_continuous(labels = function(x) paste0(x, "%"))+coord_flip()+theme(legend.position = "none")+theme(axis.title.x = element_blank())

p2<-ggplot(ss$TRB,aes(jGene,x.perc,fill=jGene))+geom_boxplot(outlier.shape=NA)+geom_point(pch = 21,position = position_jitterdodge(jitter.width=2))+facet_wrap(~locus)+scale_fill_brewer(palette="Set3")+scale_x_discrete(name="Gene",limits = rev(levels(ss$TRB$jGene)))+scale_y_continuous(labels = function(x) paste0(x, "%"))+coord_flip()+theme(legend.position = "none")+theme(axis.title.x = element_blank())

pdf(paste0(outpath,"usage_jGene_subgroup.pdf"))
grid.arrange(p1,p2,ncol=2)
dev.off()



