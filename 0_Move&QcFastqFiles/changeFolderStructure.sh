for folderpath in FASTQ_Generation_2018-12-19_19_32_44Z-144643634/* 
do 
	echo $folderpath
	#get foldername from path
	foldername="${folderpath##*/}"
	id1=$(echo $foldername | perl -lne 'print $1 if /(.*)_L001/')
	echo $id1 
	for f in $folderpath/*
	do
		echo "f: "$f
		#get foldername from path
		filename="${f##*/}"
		id2=$(echo $filename | perl -lne 'print $1 if /.*(_L001.*)/')
		concat=$id1$id2
		mkdir "FastqGzOriginal/"$id1
		echo "origin: "$f
		echo "FastqGzOriginal/"$id1"/"$concat
		cp $f "FastqGzOriginal/"$id1"/"$concat
	done
done

