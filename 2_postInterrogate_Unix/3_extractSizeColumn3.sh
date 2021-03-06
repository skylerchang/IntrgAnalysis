#!/bin/bash
#v3: extracts col from <locus>.clntab AND <locus>_uniqueCDR3s.clntab
echo "Target folder: '$1'"

loci="TRA TRB TRD TRG noju"		#TRD TRG noju

for dir in $1/Data/*
do
	echo $dir
	for locus in $loci 
	do 
		fname=`basename $dir`
		awk -F'\t' '{print $27}' $dir/${fname}_${locus}.clntab > $dir/${fname}_${locus}_col27.clntab
		awk -F'\t' '{print $27}' $dir/${fname}_${locus}_uniqueCDR3s.clntab > $dir/${fname}_${locus}_uniqueCDR3s_col27.clntab
	done
done


