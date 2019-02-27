touch lineCount.txt
for d in FASTQ_Generation_2018-08-10_18_02_24Z-114998115/*
do
	for f in $d/*
	do
		echo $f
		wc $f >> lineCount.txt
	done
done




