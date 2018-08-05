#!/bin/bash

#dir1="/Users/SKeller/Downloads/all/"
#dir1="/Users/SKeller/Documents/Projects/feTR_HTS/Lymphomas/all/Clntab_2018-04-25/Data/"
dir1="/Users/SKeller/Documents/Sequencing/Runs/k9MultiLoci3-73454401/Clntab_2018-04-26/Data/"

files=`ls $dir1`

for f in $files; do
	echo "f "$f
#	base=${f%_L006_R1_001.fastq.processed.junctioned.profiled.clntab}
	base=${f%_L001_R1_001.fastq.processed.junctioned.profiled.clntab}
#  	base=$(echo $f | cut -d'_L' -f 1)
	echo $base
	mkdir $dir1/$base
	mv $dir1/$f $dir1/$base/${base}.clntab
done
