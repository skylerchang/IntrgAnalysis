#qc

targetDir="../../Test"

mkdir $targetDir"/Trimmed/"
mkdir $targetDir"/TrimmingLog/"

for d in $targetDir"/Original/"*
do 
	echo "d: "$d
	
	base="${d##*/}"
	echo "base: "$base
#	subdir="/${1##*/}"
	
	
	files=`ls $d`
	set -- $files
	echo "1: "$1
	echo "2: "$2
	
	inOne=$d"/"$1
	inTwo=$d"/"$2
	outOne=$targetDir"/Trimmed/"$1
	outTwo=$targetDir"/Trimmed/"$2
	
	logHtml=$targetDir"/TrimmingLog/"$base".html"
	logJson=$targetDir"/TrimmingLog/"$base".json"
	
	echo $inOne
	echo $inTwo
	echo $outOne
	echo $outTwo
	
	fastp -i $inOne -I $inTwo -o $outOne -O $outTwo -h $logHtml -j $logJson
done

#echo "Running MultiQC ..."
#multiqc $targetDir"/OriginalQC/"
#mkdir $targetDir"/OriginalQCMulti/"
#mv multi* $targetDir"/OriginalQCMulti/"


