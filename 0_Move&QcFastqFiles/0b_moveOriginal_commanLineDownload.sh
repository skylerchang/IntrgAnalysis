
targetFolder="../../Data/"
mkdir $targetFolder"Original"




for d in $targetFolder/"Basespace"/*
do 
	echo "d: "$d
		
	for f in $d/*.json
	do
		echo "f: "$f
		filename="${f##*/}"
		echo "filename1: "$filename
		part1=$(echo $filename | perl -lne 'print $1 if /^(.*_)L001/')
		echo "part1: "$part1
 	done
 	
 	for f in $d/*.gz
 	do
 		echo $f
 		filename="${f##*/}"
 		part2=$(echo $filename | perl -lne 'print $1 if /^(.*)_001.fastq.gz/')
		echo "part2: "$part2
		echo "finalFileName: "$part1$part2
		mv $f $targetFolder"/Original/"$part1$part2
 	done
done
