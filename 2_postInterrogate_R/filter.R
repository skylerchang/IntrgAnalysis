#script calculates diversity for split 'Clntab_RDS' files: 'Clntab_RDS_noTemplate' & 'Clntab_RDS_templateOnly'

library(tidyverse)
library(plyr)
library(here)
library(openxlsx)
library(RColorBrewer)
library(tcltk)

setwd(here())
getwd()

source("2_postInterrogate_R/functions.R")

#=========== adjust =============
msgBox <- tkmessageBox(title = "Alert",
                       message = "Don't forget to adjust 'run' and 'sampleNoCodesForFraction'", icon = "info", type = "ok")

run<-25
sampleNoCodesForFraction<-F       #T or F

rdsPath<-'../Data/Clntab_RDS/clntab_vAndJ.rds'
resultsPathReads<-'../Results/ReadCounts/'
resultsPathTemplate<-'../Results/Templates/'

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

#for run '26-24-21' 
#problem: multiple ntc & D 16 controls result in >2 replicates
#solution: convert ntc & D16 controls from different runs into different submissions
exclude<-c(
  #Run 21
  "NTC-can-0D1P146C1P1C1_NTC-NTC-canine_S29",
  "NTC-can-0D1P147C1P1C1_NTC-NTC-canine_S30",
  #Run 24
  "k9-nt_H2O-0D1P454C1L1P1_NTC_S25",
  "k9-nt_H2O-0D1P455C1L1P1_NTC_S26",
  "k9-nt_H2O-0D1P456C1L1P1_NTC-Template_S27",
  "k9-nt_H2O-0D1P457C1L1P1_NTC-Template_S28",
  #Run 26
  "k9-nt_H2O-0D1P460C1L1P1_NTC_S57",
  "k9-nt_H2O-0D1P461C1L1P1_NTC_S58",
  "k9-nt_H2O-0D1P393C1L1P1_NTC-Template_S59",
  "k9-nt_H2O-0D1P394C1L1P1_NTC-Template_S60",
  #Run21:
  "D16-07-1FD1P250C1P1C1_OWN33-PAT1_S23",
  "D16-07-1FD1P251C1P1C1_OWN33-PAT1_S24",
  #Run24:
  "k9-pc_D16-07-6D1P89C1L1P1_D16_S29",
  "k9-pc_D16-07-6D1P90C1L1P1_D16_S30",           
  "k9-pc_D16-07-6D1P91C1L1P1_D16-Template_S31",
  "k9-pc_D16-07-6D1P92C1L1P1_D16-Template_S32"    
)
if (run=='26-24-21'){
  #filter out WGA samples (152=>104)
  includedLogical<-!grepl("WGA",files_short)
  files_short<-files_short[includedLogical]
  datalist<-datalist[includedLogical]
  length(files_short)
  
  #exclude redundant controls (104-16=88)
  includedLogical<-!grepl(paste(exclude,collapse="|"),files_short)
  files_short<-files_short[includedLogical]
  datalist<-datalist[includedLogical]
  
  #rename controls with aberrant format
  files_short<-files_short %>%
    #Run 21
    sub("NTC-can-0D1P171C1P1C1_NTC-NTC-canine_S21","21-H2O-0D1P171C1P1C1_21ntc-21ntc_S21",.) %>%
    sub("NTC-can-0D1P172C1P1C1_NTC-NTC-canine_S22","21-H2O-0D1P172C1P1C1_21ntc-21ntc_S22",.) %>%
    #Run 24
    sub("k9-nt_H2O-0D1P391C1L1P1_NTC_S57","24-H2O-0D1P391C1L1P1_24ntc-24ntc_S57",.) %>%
    sub("k9-nt_H2O-0D1P392C1L1P1_NTC_S58","24-H2O-0D1P392C1L1P1_24ntc-24ntc_S58",.) %>%
    sub("k9-nt-st_H2O-0D1P393C1L1P1_NTC-Template_S59","24-H2O-st-0D1P393C1L1P1_24ntcST-24ntcST_S59",.) %>%
    sub("k9-nt-st_H2O-0D1P394C1L1P1_NTC-Template_S60","24-H2O-st-0D1P394C1L1P1_24ntcST-24ntcST_S60",.) %>%
    #Run 26
    sub("k9-nt_H2O-0D1P459C1L1P1_NTC_S25","26-H2O-0D1P459C1L1P1_26ntc-26ntc_S25",.) %>%
    sub("k9-nt_H2O-0D1P458C1L1P1_NTC_S26","26-H2O-0D1P458C1L1P1_26ntc-26ntc_S26",.) %>%
    sub("k9-nt-st_H2O-0D1P456C1L1P1_NTC-Template_S27","26-H2O-st-0D1P456C1L1P1_26ntcST-26ntcST_S27",.) %>%
    sub("k9-nt-st_H2O-0D1P457C1L1P1_NTC-Template_S28","26-H2O-st-0D1P457C1L1P1_26ntcST-26ntcST_S28",.) %>%
    #Run21:
    sub("D16-07-1FD1P158C1P1C1_OWN33-PAT1_S31","21-D16-07-1FD1P158C1P1C1_21pcc-21pcc_S31",.) %>%
    sub("D16-07-1FD1P159C1P1C1_OWN33-PAT1_S32","21-D16-07-1FD1P159C1P1C1_21pcc-21pcc_S32",.) %>%
    #Run24:
    sub("k9-pc_D16-07-6D1P29C1L1P1_D16_S61","24-D16-07-6D1P29C1L1P1_24pcc-24pcc_S61",.) %>%            
    sub("k9-pc_D16-07-6D1P30C1L1P1_D16_S62","24-D16-07-6D1P30C1L1P1_24pcc-24pcc_S62",.) %>%      
    sub("k9-pc_D16-07-6D1P31C1L1P1_D16-Template_S63","24-D16-07-st-6D1P31C1L1P1_24pccST-24pccST_S63",.) %>%  
    sub("k9-pc_D16-07-6D1P32C1L1P1_D16-Template_S64","24-D16-07-st-6D1P32C1L1P1_24pccST-24pccST_S64",.)
  length(files_short)
}
str(files_short)
#initialize things
templateSummary<-tibble()
readCountSummary<-tibble()
filteredData<-list()
filteredDataForDB<-tibble()

for (i in 1:length(files_short)){  
  print(files_short[i])
  clUp2Id<-files_short[i] %>% strsplit(.,"_") %>% unlist() %>% .[1]
  a<-datalist[[i]][,c("vGene","jGene","aaSeq","aaLength","readCount","completeNtSeq","locus")]
  if(is.null(a)){
    print("no data - NULL")
    next
  }
  if(nrow(a)==0){
    print("no data")
    next
  }
  count.all<-sum(a$readCount)
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
  templates<-a[a$template!='no',]
  count.templates<-sum(templates$readCount)
  
  #continue with dataset w/o templates
  a<-a[a$template=='no',]
  count.noTemplates<-sum(a$readCount)
  a<-subset(a,select=-(template))
  
  #split by locus
  igh<-a[a$locus=='IGH',]
  trb<-a[a$locus=='TRB',]
  trg<-a[a$locus=='TRG',]
  count.igh.before<-sum(igh$readCount)
  count.trb.before<-sum(trb$readCount)
  count.trg.before<-sum(trg$readCount)
  count.all.before.sum<-count.igh.before+count.trb.before+count.trg.before
  
  #============ filter seqs by seqs & ^C.{3,30}[FW]$ =======
  igh<-igh[grepl("^C.{3,30}[W]$",igh$aaSeq),]
  trb<-trb[grepl("^C.{3,30}[F]$",trb$aaSeq),]
  trg<-trg[grepl("^C.{3,30}[FAC]$",trg$aaSeq),]
  
  temp<-bind_rows(igh,trb,trg)
  filteredData[[i]]<-temp
  #replace col 'sample' by col 'clUp2Id' -> for import in DB
  temp<-subset(temp,select=-(filename))
  if (nrow(temp)!=0){
    temp$clUp2Id<-clUp2Id
    filteredDataForDB<-bind_rows(filteredDataForDB,temp)
  }

  count.igh.after<-sum(igh$readCount)
  count.trb.after<-sum(trb$readCount)
  count.trg.after<-sum(trg$readCount)
  count.all.after.sum<-count.igh.after+count.trb.after+count.trg.after
  count.igh.anchorFiltered<-count.igh.before-count.igh.after
  count.trb.anchorFiltered<-count.trb.before-count.trb.after
  count.trg.anchorFiltered<-count.trg.before-count.trg.after
  
  #============= create read count summary =============
  #checksum
  count.all-(count.templates+
               count.igh.anchorFiltered+count.igh.after+
               count.trb.anchorFiltered+count.trb.after+
               count.trg.anchorFiltered+count.trg.after
             )
  #create tibble
  temp<-tibble(
    filename=files_short[i],
    readCount.total=count.all,
    
    readCount.template=count.templates,
    readCount.noTemplate=count.noTemplates,
    readCount.all.before.sum=count.all.before.sum,
    
    readCount.igh=count.igh.after,
    falseAnchor.igh=count.igh.anchorFiltered,
    readCount.trb=count.trb.after,
    falseAnchor.trb=count.trb.anchorFiltered,
    readCount.trg=count.trg.after,
    falseAnchor.trg=count.trg.anchorFiltered
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
readCountSummary2<-splitFilename(readCountSummary,sampleNoCodesForFraction,run)
readCountSummary2$fraction

#save as xlsx
wb<-createWorkbook()
addWorksheet(wb,"summary")
writeData(wb,"summary",readCountSummary2)
saveWorkbook(wb,paste0(resultsPathReads,"readCountSummary_repsAsRows.xlsx"),overwrite=T)

#save as rds
saveRDS(readCountSummary2,paste0(resultsPathReads,"readCountSummary_repsAsRows.rds"))

#transform to long format and split filename
readCount.long<-readCountSummary %>%
  select(-c(readCount.total,readCount.noTemplate,readCount.all.before.sum)) %>%
  gather(variable,value,-filename) %>%
  splitFilename(.,sampleNoCodesForFraction,run)

#plot by submission for 'sampleNoCodesForFraction<-T'
if (sampleNoCodesForFraction==T){
  pdf('../Results/ReadCounts/readCountSummary_bySubmission.pdf')
  ggplot(readCount.long,aes(replicate,value,fill=variable))+geom_col()+coord_flip()+scale_fill_brewer(palette='Dark2')+facet_grid(submission~fraction)+theme(strip.text.y=element_text(angle=360))
  dev.off()
  
  pdf('../Results/ReadCounts/readCountSummary_byOwnerPatient.pdf')
  ggplot(readCount.long,aes(replicate,value,fill=variable))+geom_col()+coord_flip()+scale_fill_brewer(palette='Dark2')+facet_grid(ownerPatient~fraction)+theme(strip.text.y=element_text(angle=360))
  dev.off()
}else{
  #plot by sample for 'sampleNoCodesForFraction<-T'
  pdf('../Results/ReadCounts/readCountSummary_bySample.pdf')
  ggplot(readCount.long,aes(replicate,value,fill=variable))+geom_col()+coord_flip()+scale_fill_brewer(palette='Dark2')+facet_grid(sample~fraction)+theme(strip.text.y=element_text(angle=360))
  dev.off()
}

#======= transpose reps to cols ===========
readCountSummary3<-readCountSummary2 %>%
  select(-c(filename,idNumber,id,pcr)) %>%
  gather(key,value,readCount.total:falseAnchor.trg) %>%
  unite(key_rep,key,replicate) %>%
  spread(key_rep,value)

wb<-createWorkbook()
addWorksheet(wb,'summary')
writeData(wb,'summary',readCountSummary3)
saveWorkbook(wb,paste0(resultsPathReads,"readCountSummary_repsAsCols.xlsx"),overwrite = T)

#===============================================================
#========================== template summary ===================
#===============================================================
templateSummary
#exclude 'no' template
templateSummary2<-templateSummary[templateSummary$template!='no',]

#transformt to wide format
templateSummary2.wide<-spread(templateSummary2,template,readCount)

#save as xlsx
wb<-createWorkbook()
addWorksheet(wb,'summary')
writeData(wb,'summary',templateSummary2.wide)
saveWorkbook(wb,paste0(resultsPathTemplate,'templateSummary.xlsx'),overwrite=T)

#save as rds
saveRDS(templateSummary2.wide,paste0(resultsPathTemplate,'templateSummary.rds'))

#plot stratified by file
pdf(paste0(resultsPathTemplate,'templatesBytemplate.pdf'))
ggplot(templateSummary,aes(template,readCount,fill=filename))+geom_col()+coord_flip()+theme(legend.position="none")
dev.off()

#plot stratified by template - bar plot
pdf(paste0(resultsPathTemplate,'templatesBySample_barplot.pdf'))
ggplot(templateSummary,aes(filename,readCount,fill=template))+geom_col()+coord_flip()
dev.off()
#plot stratified by template - line plot
pdf(paste0(resultsPathTemplate,'templatesBySample_lineplot.pdf'))
ggplot(templateSummary2,aes(filename,readCount,group=template))+geom_point()+geom_line(aes(color=template))+coord_flip()
dev.off()

