#!/bin/bash

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