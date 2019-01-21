
targetFolder="../../Data/"
mkdir $targetFolder"Original"




for d in $targetFolder/"Basespace"/* 
do 
	echo "d: "$d
	
	#get first part of foldername for part1 of final file name
	dirname="${d##*/}"
	echo "dirname: "$dirname
	part1=$(echo $dirname | perl -lne 'print $1 if /(.*)_L001/')
	echo "part1: "$part1
	mkdir $targetFolder"Original/"$part1
	
 	for f in $d/*
 	do
 		echo $f
 		filename="${f##*/}"
 		part2=$(echo $filename | perl -lne 'print $1 if /(_L001.*)/')
		echo "part2: "$part2
		echo "finalFileName: "$part1$part2
		mv $f $targetFolder"/Original/"$part1"/"$part1$part2
 	done
 	rmdir $d
done

rmdir $targetFolder/"Basespace"