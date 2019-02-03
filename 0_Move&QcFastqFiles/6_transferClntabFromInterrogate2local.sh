#!/bin/bash

#cd to the local directory you want the files be transferred to

while read d
do 
	echo $d
	file="/intrg_storage/users/smkeller/analyses/Run21_CSF_6samples_2018_11_20_Q30L50/last_run_res/Run21_CSF_6samples_2018_11_20_Q30L50--results/"$d"/"$d".fastq.processed.junctioned.profiled.clntab"
	echo $file
	scp -p keeper@s0.arrest.tools:$file ./
done < files2transfer2.txt