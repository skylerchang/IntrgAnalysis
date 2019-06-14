#script calculates diversity for split 'Clntab_RDS' files: 'Clntab_RDS_noTemplate' & 'Clntab_RDS_templateOnly'

library(tidyverse)
library(vegan)
library(plyr)
library(here)
library(reshape2)
library(openxlsx)

setwd(here())
getwd()

#=========== adjust =============

sampleNoCodesForFraction<-F        #T or F

rdsPath<-'../Data/Clntab_RDS/clntab_vAndJ.rds'
resultsPath<-'../Results/Diversity/'

dir.create(resultsPath)

loci<-c("TRB","IGH")
#===================================
t<-read_rds(rdsPath)

#datalist contains one tibble per sample
datalist<-t[[1]]
#sample names are stored separately in a vector
files_short<-t[[2]]

templateSummary<-tibble()
diversity<-tibble()

for (i in 1:length(datalist)){
  a<-datalist[[i]][,c("vGene","jGene","aaSeq","aaLength","size","completeNtSeq","vAndJchainSimplified")]
  if(nrow(a)==0){next}
  #============ standardize sample name =======
  #assume the following default format:
  # [1] "16-088089-3D1P3C1P1C1_Mcdermott-Gibbie"      
  # [2] "16-088089-3D1P4C1P1C1_Mcdermott-Gibbie"   
  #deviations:
  #[92] "k9-pc_D16-07-2FD4P59C1P1C1_OWN33-PAT1_S92"      
  #[93] "k9-pc-st_D16-07-2FD1P32C1P1C1_OWN33-PAT1_S93"  
  #<sampleId>_<sampleName>
  #tranform samples to fit default format -> only one '_' which needs to be between sampleId & sampleName
  file<-sub("_S[0-9]+$","",files_short[i])
  file<-sub("k9-pc-st_D16-07","D16-07-st",file)
  file<-sub("k9-pc_D16-07","D16-07",file)
  file<-sub("k9-nt-st_H2O","ntc",file)
  file<-sub("k9-nt_H2O","ntc-st",file)
  sampleId<-sub("_.*$","",file)
  sampleName<-sub("^.*_","",file)
  
  #============ extract replicate & sample =======
  replicate<-ifelse(grepl("D[0-9]+P[0-9]*[02468]C",sampleId),"rep1","rep2")
  if(sampleNoCodesForFraction==T){
    fraction<-ifelse(grepl("-1D[0-9]+P[0-9]+C[0-9]+P[0-9]+C[0-9]+$",sampleId),"fraction1","fraction2")
  }else{
    fraction<-"fraction1"
  }
  sampleIdShort<-sub("D[0-9]+P[0-9]+C[0-9]+P[0-9]+C[0-9]+$","",sampleId)
  submission<-sub("-[0-9a-zA-Z]+$","",sampleIdShort)
  
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
  a.woTemp<-a[a$template=='no',]
  a.tempOnly<-a[a$template!='no',]
  
  #============ summarize templates for sample and store in tibble =======
  temp<-as_tibble(ddply(a.tempOnly[,c("template","size")],"template",numcolwise(sum)))
  if(nrow(temp)==0){next}
  temp$sample<-files_short[i]
  templateSummary<-bind_rows(templateSummary,temp)
  
  #============ diversity =======
  for (j in 1:length(loci)){
    aa<-as_tibble(ddply(a.woTemp[,c("aaSeq","size")],"aaSeq",numcolwise(sum)))
    aa.shannon<-diversity(aa$size, index = "shannon",base = exp(1))
    aa.simpson<-diversity(aa$size, index = "simpson")
    aa.invsimpson<-diversity(aa$size, index = "invsimpson")
    
    temp<-tibble(shannon=aa.shannon,simpson=aa.simpson,invsimpson=aa.invsimpson,sample=files_short[i],locus=loci[j],replicate=replicate,submission=submission,sampleIdShort=sampleIdShort)
    diversity<-bind_rows(diversity,temp)
  }
}

#============= diversity summary ===========
d<-diversity#[rowSums(is.na(diversity)),]
d.split<-split(d,d$submission)
str(d.split)

#print to xlsx
wb<-createWorkbook()
addWorksheet(wb,"summary")
writeData(wb,"summary",d)
for (i in 1:length(d.split)){
  addWorksheet(wb,d.split[[i]]$submission[1])
  writeData(wb,d.split[[i]]$submission[1],d.split[i],colNames = TRUE)
}
saveWorkbook(wb,paste0(resultsPath,"diversitySummary.xlsx"),overwrite = T)

#transform to long format
d.long<-as_tibble(melt(d,id.vars=c("sample","locus","replicate","submission","sampleIdShort")))
d.long<-na.omit(d.long)
d.long.subset<-d.long[d.long$value!=Inf,]

pdf(paste0(resultsPath,"diversityPlot_compareIndexes.pdf"))
ggplot(d.long.subset,aes(variable,value))+
  geom_boxplot()+
  geom_jitter(mapping=aes(shape=d.long.subset$locus,color=d.long.subset$sampleIdShort),size=5)+
dev.off()

#======== template summary ===========
templateSummary
ggplot(templateSummary,aes(template,size,fill=sample))+geom_col()+coord_flip()
ggplot(templateSummary,aes(sample,size,fill=template))+geom_col()+coord_flip()
