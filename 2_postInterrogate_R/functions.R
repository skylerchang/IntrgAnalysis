#functions that are used by the following scripts: diversity.R, filter.R, locus.R

#x: data, y: sampleNoCodesForFraction (T/F), z: run
splitFilename<-function(x,y,z){
  #define fractions dependent on sequencing run
  fraction.blood<-c(13,17,22,27,29)
  fraction.csf1<-c(21,24,26)
  fraction.csf2<-c(28)
  fraction.ln<-c(30)
  if(z %in% fraction.blood){fraction1<-'Blood';fraction2<-'Plasma'}
  else if(z %in% fraction.csf1){fraction1<-'CSF pellet';fraction2<-'CSF supernatant'}
  else if(z %in% fraction.csf2){fraction<-'CSF'}
  else if(z %in% fraction.ln){fraction<-'Lymph node'}
  else{fraction<-'noFractionDefined'}
  
  #split filename
  splitFilename<-strsplit(x$filename,"_")
  x$id<-laply(splitFilename, '[[', 1)
  x$ownerPatient<-laply(splitFilename, '[[', 2)
  x$idNumber<-laply(splitFilename, '[[', 3)
  x$pcr<-sub("C[0-9]+P[0-9]+C[0-9]+","",x$id)
  x$pcr<-sub("C[0-9]+L[0-9]+P[0-9]+","",x$pcr)
  x$sample<-sub("-*D[0-9]+P[0-9]+$","",x$pcr)
  x$submission<-sub("-[0-9a-zA-Z]+$","",x$sample)
  x$replicate<-ifelse(grepl("D[0-9]+P[0-9]*[13579]$",x$pcr),"rep1","rep2")
  if(y==T){
    x$fraction<-ifelse(grepl("-1$",x$sample),fraction1,fraction2)
  }else{
    x$fraction<-rep(fraction,length(sample))
  }
  x
}