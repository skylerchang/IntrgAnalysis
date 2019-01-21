targetFolder="../../Data"
destination=$targetFolder"/Trimmed30"

mkdir -p $targetFolder"/Counts/"

touch $targetFolder"/Counts/before-names.txt"
touch $targetFolder"/Counts/before-count.txt"
touch $targetFolder"/Counts/after-names.txt"
touch $targetFolder"/Counts/after-count.txt"

for d in $targetFolder/"Original/"* 
do 
#	echo "d: "$d
	for f in $d/*
	do
		echo "f: "$f
		echo $f >> $targetFolder"/Counts/before-names.txt"
		zgrep -Ec "$" $f >> $targetFolder"/Counts/before-count.txt"
	done 
done

for f in $targetFolder/"Trimmed30/Paired/"* 
do 
	echo "f: "$f
	echo $f >> $targetFolder"/Counts/after-names.txt"
	zgrep -Ec "$" $f >> $targetFolder"/Counts/after-count.txt"

done

