#!/bin/bash

dir1="/Users/SKeller/Downloads/all/"
dir1="FASTQ_Generation_2018-08-10_18_02_24Z-114998115"
dir1="../../Data/Original/"
destinationFolder='CombinedFastqFiles'

subs=`ls $dir1`

for i in $subs; do
  echo $i
#  mv $dir1/$i/* $dir1/ 
	cp $dir1/$i/* $destinationFolder/ 
done