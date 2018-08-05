#!/bin/bash
#based on 'subsetClntabs_TRs_<version>'
#v4: sums size column to obtain total read number per file
#v5: Clntab file folder structure: sorted by sample and not loci
#v6: grep strategy changed to eliminate nojus in TR groups
#

#Usage: ./subsetClntabs_TRs_v5.sh <ClntabFolder>

#Expected data file structure for clntab analysis:
#Clntab_yyyy-mm-dd/
#	Data/
#		Sample1/
#			sample1.clntab
#		Sample2/
#			sample2.clntab
#		...
#	Results/


echo "Target folder: '$1'"

mkdir $1/Results/Summaries


#get columns for functionality and productivity from header of first file and store for later use
echo "Assessing header..."
#get 1st subdirectory in directory
d1="$(ls $1/Data/* | head -1)"; d2=${d1/:/}
echo "d1: "$d1
#get first file in directory
f1="$(ls $d2/* | head -1)"
echo "f1: "$f1
echo "d2: "$d2

#determine the column that contains the productivity information, i.e. header = 'sequence.functionality' (for use below)
#colSize=2
#colProductivity=1
colProductivity="$(awk -v RS="\t" '/^sequence.functionality$/{print NR;}' $f1)"
echo "Column containing productivity: "$colProductivity
colSize="$(awk -v RS="\t" '/^sequence.size$/{print NR;}' $f1)"
echo "Column containing size: "$colSize

#write header line to separate file in 'Summaries' directory
echo "generating header file ..."
sed '1q;d' $f1 > $1/Results/Summaries/clntab_header.txt
sed '1q;d' $f1 | tr '\t' $'\n' > $1/Results/Summaries/clntab_header_transposed.txt


#************* IGH/TRA/TRB/TRD/TRG/noju ***************
#split clntab file based on J gene usage
loci="IGH TRA TRB TRD TRG noju"		#TRD TRG noju
echo "splitting into IGH/TRA/TRB/TRD/TRG/noju based on J gene ..."
for dir in $1/Data/*
do
	echo
	echo "dir: "$dir
	#get first file in directory
	f1="$(ls $dir/* | head -1)"; echo "f1: "$f1
	fname=`basename $f1`; base=${fname/.clntab/}
	for locus in $loci 
	do 
		echo $locus
		if [ "$locus" = "noju" ]; then 
			grep 'noju' $f1 > $1/Data/$base/${base}_${locus}.clntab;
		else
			target=$locus"J"
			grep -v 'noju' $f1 | grep $target > $1/Data/$base/${base}_${locus}.clntab;
		fi
	done
done




#numbers check - count lines in each file
echo "producing stats ..."

#get counts for original clntab file
touch $1/Results/clntab-unireads.txt
touch $1/Results/clntab-reads.txt
echo "Counting lines in original clntab files ..."
for dir in $1/Data/*
do
	echo "dir: "$dir
	fname=`basename $dir`; echo $fname; 
	wc -l $dir/${fname}.clntab >> $1/Results/clntab-unireads.txt; 
	awk -F'\t' -v c1=$colSize '{total = total + $c1}END{print total, FILENAME}' $dir/${fname}.clntab >> $1/Results/clntab-reads.txt; 
done


for dir in $1/Data/*
do
	echo "dir: "$dir							#format: 'Clntab_test//Data/sample1'
	fname=`basename $dir`; echo $fname; 		#format: 'sample1'
	for locus in $loci 
	do 
		#print header
		echo $locus
		#output all lines of a file
		wc -l $dir/${fname}_${locus}.clntab | sed 's/^ *//g' >> $1/Results/clntab_$locus-unireads.txt
		#output the sum of all lines of col 27 (sequence.size) => to calculate total read number
		awk -F'\t' -v c1=$colSize '{total = total + $c1}END{print total, FILENAME}' $dir/${fname}_${locus}.clntab >> $1/Results/clntab_$locus-reads.txt; 
		#if sum equals zero, i.e. does not exist, insert '0'
		sed 's/^ /0 /g' $1/Results/clntab_$locus-reads.txt > tmpfile; mv tmpfile $1/Results/clntab_$locus-reads.txt
	done
done

#summarize uniread counts
echo "Unireads cross-check: cols: IGH TRA TRB TRD TRG noju sum all sum-all; col9 should equal 0"
echo "IGH TRA TRB TRD TRG noju sum all sumMinusAll sample" > $1/Results/clntab_ALL-unireads.txt;
paste $1/Results/clntab_IGH-unireads.txt $1/Results/clntab_TRA-unireads.txt $1/Results/clntab_TRB-unireads.txt $1/Results/clntab_TRD-unireads.txt $1/Results/clntab_TRG-unireads.txt $1/Results/clntab_noju-unireads.txt $1/Results/clntab-unireads.txt | awk '{print $1,$3,$5,$7,$9,$11, ($1+$3+$5+$7+$9+$11), ($13-1), ($1+$3+$5+$7+$9+$11-$13+1), $14}' 
paste $1/Results/clntab_IGH-unireads.txt $1/Results/clntab_TRA-unireads.txt $1/Results/clntab_TRB-unireads.txt $1/Results/clntab_TRD-unireads.txt $1/Results/clntab_TRG-unireads.txt $1/Results/clntab_noju-unireads.txt $1/Results/clntab-unireads.txt | awk '{print $1,$3,$5,$7,$9,$11, ($1+$3+$5+$7+$9+$11), ($13-1), ($1+$3+$5+$7+$9-$11-$13+1), $14}' >> $1/Results/clntab_ALL-unireads.txt

echo "Reads cross-check: cols: IGH TRA TRB TRD TRG noju sum all sum-all; col9 should equal 0"
echo "IGH TRA TRB TRD TRG noju sum all sumMinusAll sample" > $1/Results/clntab_ALL-reads.txt;
paste $1/Results/clntab_IGH-reads.txt $1/Results/clntab_TRA-reads.txt $1/Results/clntab_TRB-reads.txt $1/Results/clntab_TRD-reads.txt $1/Results/clntab_TRG-reads.txt $1/Results/clntab_noju-reads.txt $1/Results/clntab-reads.txt | awk '{print $1,$3,$5,$7,$9,$11, ($1+$3+$5+$7+$9+$11), ($11), ($1+$3+$5+$7+$9+$11-$13), $14}' 
paste $1/Results/clntab_IGH-reads.txt $1/Results/clntab_TRA-reads.txt $1/Results/clntab_TRB-reads.txt $1/Results/clntab_TRD-reads.txt $1/Results/clntab_TRG-reads.txt $1/Results/clntab_noju-reads.txt $1/Results/clntab-reads.txt | awk '{print $1,$3,$5,$7,$9,$11, ($1+$3+$5+$7+$9+$11), ($11), ($1+$3+$5+$7+$9+$11-$13), $14}' >> $1/Results/clntab_ALL-reads.txt

