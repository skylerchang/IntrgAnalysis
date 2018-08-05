#!/bin/Bash

echo "Target folder: '$1'"

#mkdir $1/Summaries


#get columns for functionality and productivity from header of first file and store for later use
echo "Assessing header..."
#get first file in directory
d1="$(ls $1/Data/* | head -1)"; d2=${d1/:/}
#get first file in directory
f1="$(ls $d2/*1000* | head -1)"
echo "f1: "$f1
echo "d2: "$d2

#determine the column that contains the productivity information, i.e. header = 'sequence.functionality' (for use below)
#colSize=2
#colProductivity=1
colNtLength="$(awk -v RS="\t" '/^sequence.JUNCTION.nt seq.len$/{print NR;}' $f1)"
echo "Column containing junctional nt seq length: "$colNtLength
colSize="$(awk -v RS="\t" '/^sequence.size$/{print NR;}' $f1)"
echo "Column containing size: "$colSize


#split clntab file based on J gene usage
loci="TRB"		#TRD TRG noju
echo "..."
for dir in $1/Data/*
do
	echo "dir: "$dir
	for locus in $loci 
	do 
		fname="$(basename $dir)"; target=$dir"/"$fname"_"$locus".clntab"
		echo $target
		for i in {40..50}
		do
			echo $i
			awk -v var="$i" -F"\t" '$21 == var { print $21"\t"$27 }' $target > $1/Data/$fname/${fname}_${locus}_${i}.clntab;
		done
	done
done

