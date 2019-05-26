module load trimmomatic
mkdir ../../Data/Trimmed30 



targetFolder="../../Data"
destination=$targetFolder"/Trimmed30"

mkdir -p $destination"/Paired/"
mkdir -p $destination"/Unpaired/"



for f in $targetFolder/"Original/"* 
do 
	echo "f: "$f
	
	#define originalRead1&2 files
	files=`ls $f`
	set -- $files
	
	originalRead1=$f"/"$1
	originalRead2=$f"/"$2
	trimmed_P_Read1=$destination"/Paired/"$1
	trimmed_UP_Read1=$destination"/Unpaired/"$1
	trimmed_P_Read2=$destination"/Paired/"$2
	trimmed_UP_Read2=$destination"/Unpaired/"$2
	trimlogFile=$destination"/Logfiles/"$1".txt"
	
	echo "trimming ..."
	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -phred33 \
	$originalRead1 $originalRead2 \
	$trimmed_P_Read1 $trimmed_UP_Read1 \
	$trimmed_P_Read2 $trimmed_UP_Read2 \
	ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:50
done
