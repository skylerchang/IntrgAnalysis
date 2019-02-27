#!/bin/bash

echo "Target folder: '$1'"

loci="TRA TRB TRD TRG noju"		#TRD TRG noju

for dir in $1/Data/*
do
	echo $dir
	for locus in $loci 
	do 
		fname=`basename $dir`
		awk -F'\t' '{print "\""$15"\"""\t"$21"\t"$23"\t"$27}' $dir/${fname}_${locus}.clntab > $dir/${fname}_${locus}_col15-21-23-27.clntab
	done
done


