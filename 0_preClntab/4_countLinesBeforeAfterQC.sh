targetFolder="../../Data"
destination=$targetFolder"/Trimmed30"

mkdir -p $targetFolder"/Counts/"

touch $targetFolder"/Counts/beforeQC-names.txt"
touch $targetFolder"/Counts/beforeQC-count.txt"
touch $targetFolder"/Counts/afterQC-names.txt"
touch $targetFolder"/Counts/afterQC-count.txt"

for d in $targetFolder/"Original/"* 
do 
	echo "d: "$d
	for f in $d/*
	do
		echo "f: "$f
		echo $f >> $targetFolder"/Counts/beforeQC-names.txt"
		zgrep -Ec "$" $f >> $targetFolder"/Counts/beforeQC-count.txt"
	done 
done

for f in $targetFolder/"Trimmed30/Paired/"* 
do 
	echo "f: "$f
	echo $f >> $targetFolder"/Counts/afterQC-names.txt"
	zgrep -Ec "$" $f >> $targetFolder"/Counts/afterQC-count.txt"

done

