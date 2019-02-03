#!/bin/bash

dir1="/Users/SKeller/Downloads/all/"

subs=`ls $dir1`

for i in $subs; do
  echo $i
  mv $dir1/$i/* $dir1/ 
done