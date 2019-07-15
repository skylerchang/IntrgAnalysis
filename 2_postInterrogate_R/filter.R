#script calculates diversity for split 'Clntab_RDS' files: 'Clntab_RDS_noTemplate' & 'Clntab_RDS_templateOnly'

library(tidyverse)
library(plyr)
library(here)
library(openxlsx)
library(RColorBrewer)


setwd(here())
getwd()

#=========== adjust =============
run<-30
sampleNoCodesForFraction<-F        #T or F

rdsPath<-'../Data/Clntab_RDS/clntab_vAndJ.rds'
resultsPathReads<-'../Results/ReadCounts/'
resultsPathTemplate<-'../Results/Template/'

#===============================
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

#initialize things
templateSummary<-tibble()
readCountSummary<-tibble()
filteredData<-list()
filteredDataForDB<-tibble()

for (i in 1:length(files_short)){  
  print(files_short[i])
  clUp2Id<-files_short[i] %>% strsplit(.,"_") %>% unlist() %>% .[1]
  a<-datalist[[i]][,c("vGene","jGene","aaSeq","aaLength","size","completeNtSeq","vAndJchainSimplified")]
  colnames(a)[5]<-'readCount'
  colnames(a)[7]<-'locus'

  count.all<-nrow(a)
  if(nrow(a)==0){next}
  a$filename<-files_short[i]
  
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
  temp<-ddply(a[,c("filename","template","readCount")],c("template","filename"),numcolwise(sum))
  templateSummary<-bind_rows(templateSummary,temp)
  
  #count template reads 
  count.templates<-nrow(a[a$template!='no',])
  
  #continue with dataset w/o templates
  a<-a[a$template=='no',]
  a<-subset(a,select=-(template))
  
  #split by locus
  igh<-a[a$locus=='IGH',]
  trb<-a[a$locus=='TRB',]
  count.igh.before<-nrow(igh)
  count.trb.before<-nrow(trb)
  
  #============ filter seqs by seqs & ^C.{3,30}[FW]$ =======
  igh<-igh[grepl("^C.{3,30}[W]$",igh$aaSeq),]
  trb<-trb[grepl("^C.{3,30}[F]$",trb$aaSeq),]
  
  temp<-bind_rows(igh,trb)
  filteredData[[i]]<-temp
  #replace col 'sample' by col 'clUp2Id' -> for import in DB
  temp<-subset(temp,select=-(filename))
  if (nrow(temp)!=0){
    temp$clUp2Id<-clUp2Id
    filteredDataForDB<-bind_rows(filteredDataForDB,temp)
  }

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
}


#=================== save data as rds ==========================
saveRDS(list(filteredData,files_short),'../Data/Clntab_RDS/clntab_vAndJ_filtered.rds')
dir.create('../Data/Clntab_delim/')
write_delim(filteredDataForDB,'../Data/Clntab_delim/clntab_vAndJ_filtered.txt')

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
ggplot(templateSummary,aes(template,readCount,fill=filename))+geom_col()+coord_flip()+theme(legend.position="none")
dev.off()

#exclude 'no' template
templateSummary2<-templateSummary[templateSummary$template!='no' & templateSummary$template!='t2',]

#plot stratified by template - bar plot
pdf(paste0(resultsPathTemplate,'templatesBySample_barplot.pdf'))
ggplot(templateSummary,aes(filename,readCount,fill=template))+geom_col()+coord_flip()
dev.off()
#plot stratified by template - line plot
pdf(paste0(resultsPathTemplate,'templatesBySample_lineplot.pdf'))
ggplot(templateSummary2,aes(filename,readCount,group=template))+geom_point()+geom_line(aes(color=template))+coord_flip()
dev.off()

