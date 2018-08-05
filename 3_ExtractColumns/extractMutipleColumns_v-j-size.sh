#!/bin/bash
#4	sequence.3-GENE
#8	sequence.5-GENE
#15	sequence.JUNCTION.aa seq
#21	sequence.JUNCTION.nt seq.len
#23	sequence.functionality
#27 sequence.size

echo "Target folder: '$1'"

loci="TRA TRB TRD TRG noju"		#TRD TRG noju

for dir in $1/Data/*
do
	echo $dir
	for locus in $loci 
	do 
		fname=`basename $dir`
		echo -e "v\tj\tsize" > $dir/${fname}_${locus}_v-j-size.clntab
		awk -F'\t' '{print $8"\t"$4"\t"$27}' $dir/${fname}_${locus}_uniqueCDR3s.clntab >> $dir/${fname}_${locus}_v-j-size.clntab
	done
done


