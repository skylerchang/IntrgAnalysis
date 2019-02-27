#!/bin/bash

#structure before:
#../../Data/Clntab/
#	sample1.fastq.processed.junctioned.profiled.clntab
#	sample2.fastq.processed.junctioned.profiled.clntab

#structure after:
#../../Data/Clntab/
#	sample1/
#		sample1.fastq.processed.junctioned.profiled.clntab
#	sample2/
#		sample2.fastq.processed.junctioned.profiled.clntab


dir1="../../Data/Clntab/"

subs=`ls $dir1`

for i in $subs; do
  echo $i
  foldername=${i/.clntab/}
  echo $foldername
  mkdir $dir1$foldername
  mv $dir1$i $dir1$foldername"/"$i
#  mv $dir1/$i/* $dir1/ 
done