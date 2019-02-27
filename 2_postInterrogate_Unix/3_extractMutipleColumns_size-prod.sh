#!/bin/bash
#4	sequence.3-GENE
#8	sequence.5-GENE
#15	sequence.JUNCTION.aa seq
#21	sequence.JUNCTION.nt seq.len
#23	sequence.functionality
#27 sequence.size

#target folder is clntab dir (e.g. 'Clntab_2018-04-17')
echo "Target folder: '$1'"

#adjust loci
loci="IGH TRA TRB TRD TRG noju"		#TRD TRG noju

for dir in $1/Data/*
do
	echo $dir
	for locus in $loci 
	do 
		fname=`basename $dir`
		echo -e "size\tprod" > $dir/${fname}_${locus}_size-prod.clntab
		awk -F'\t' '{print $27"\t"$23}' $dir/${fname}_${locus}_uniqueCDR3s.clntab >> $dir/${fname}_${locus}_size-prod.clntab
	done
done


