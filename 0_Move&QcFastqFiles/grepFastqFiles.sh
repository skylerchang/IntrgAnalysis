touch lineCount.txt
for d in FASTQ_Generation_2018-08-10_18_02_24Z-114998115/*
do
	for f in $d/*
	do
		echo $f
		wc $f >> lineCount.txt
	done
done


touch cloneBella.txt
for d in FASTQ_Generation_2018-08-10_18_02_24Z-114998115/*
do
	for f in $d/*
	do
		echo $f
		grep 'TGTGCAAGGGCAGACTACTACGATAGTTTCTGGGCTGCCTTTGGTTACTGG' $f | wc >> cloneBella.txt
	done
done


touch cloneMarishka.txt
for d in FASTQ_Generation_2018-08-10_18_02_24Z-114998115/*
do
	for f in $d/*
	do
		echo $f
		grep 'TGTGTGCCATTTAGTCCCTACGGTAGTTGGTTCGCGGATGACCAGTGG' $f | wc >> cloneMarishka.txt
	done
done


#v44-j2
touch v44-j2.txt
for d in FASTQ_Generation_2018-08-10_18_02_24Z-114998115/*
do
	for f in $d/*
	do
		echo $f
		grep 'CAGATCCACAGATCCA' $f | wc >> cloneMarishka.txt
	done
done


#v80-j3
touch v80-j3.txt
for d in FASTQ_Generation_2018-08-10_18_02_24Z-114998115/*
do
	for f in $d/*
	do
		echo $f
		grep 'ATCACGACACAGTGGT' $f | wc >> v80-j3.txt
	done
done