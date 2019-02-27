#!/bin/bash


echo "Target folder: '$1'"

#mkdir $1/SizeColumnOnlyFiles


sizeCol=23


#************* TRA/TRB/TRD/TRG/noju ***************
#split clntab file based on J gene usage
loci="TRB TRD TRG"		#TRD TRG noju

for locus in $loci 
do 
	echo "Locus: "$locus
	mkdir $1/Clntab_${locus}_col23only
	for f in $1/Clntab_${locus}/*; do fname=`basename $f`; echo $fname; awk -F'\t' '{print $2}' $f > $1/Clntab_${locus}_col23only/$fname; done
done

