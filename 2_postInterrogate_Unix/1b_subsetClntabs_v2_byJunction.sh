#!/bin/bash
#based on 'subsetClntabs_TRs_<version>'
#v4: sums size column to obtain total read number per file
#v5: Clntab file folder structure: sorted by sample and not loci
#v6: grep strategy changed to eliminate nojus in TR groups
#

#Usage: ./subsetClntabs_v1.sh <ClntabFolder>

#Expected data file structure for clntab analysis:
# !!! Folder & clntab name must be identical !!!
#Clntab_yyyy-mm-dd/
#	Data/
#		Sample1/
#			sample1.clntab
#		Sample2/
#			sample2.clntab
#		...
#	Results/

#============== v2 ==============
#Expected data file structure for clntab analysis:
# !!! Folder & clntab name must be identical !!!
#Run02_feTR_normal
#	Data/
#		Clntab/
#			Sample1/
#				sample1.clntab
#			Sample2/
#				sample2.clntab
#			...
#	Results/

#Usage: ./subsetClntabs_v2.sh <RunFolder>

echo "Target folder: '$1'"

mkdir ../../Results/Clntab

echo "Determining column headers ..."
#determine the column that contains the productivity information, i.e. header = 'sequence.functionality' (for use below)
#colSize=2
#colProductivity=1
colAaSeq=16
#colProductivity="$(awk -v RS="\t" '/^sequence.functionality$/{print NR;}' $f1)"
#echo "Column containing productivity: "$colProductivity
#colSize="$(awk -v RS="\t" '/^sequence.size$/{print NR;}' $f1)"
#3echo "Column containing size: "$colSize
#colFiveGene="$(awk -v RS="\t" '/^sequence.5-GENE$/{print NR;}' $f1)"
#echo "Column containing 5-gene: "$colFiveGene
#colThreeGene="$(awk -v RS="\t" '/^sequence.3-GENE$/{print NR;}' $f1)"
#echo "Column containing 5-gene: "$colThreeGene
#colAaSeq="$(awk -v RS="\t" '/^sequence.JUNCTION.aa seq$/{print NR;}' $f1)"
#echo "Column containing aa seq: "$colAaSeq




#========= split clntab into IGH/TRA/TRB/TRD/TRG/noju =================
#split clntab file based on J gene usage

echo "splitting into junction vs. noju based on 'sequence.JUNCTION.aa seq' ..."
for dir in $1/Data/Clntab/*
do
	echo
	echo "dir: "$dir
	#get first file in directory
	f1="$(ls $dir/* | head -1)"; echo "f1: "$f1
	fname=`basename $f1`; base=${fname/.clntab/}
	echo "base: "$base
	
#	grep 'noju' $f1 | awk -F"\t" '{print $16}' > $1/Data/Clntab/$base/${base}_noju2.clntab;
	grep -v 'noju' $f1 | awk -F"\t" '{print $8"\t"$4"\t"$16"\t"$28}' > $1/Data/Clntab/$base/${base}_junc.clntab;

done
