#On compute Canada: java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar
#On Stefan's laptop: java -jar $TRIMMOMATIC_JAR

targetFolder="../../Data"
destination=$targetFolder"/Trimmed30"

mkdir -p $destination"/Paired/"
mkdir -p $destination"/Unpaired/"
#mkdir -p $destination"/Logfiles/"

export SOFT_DIR=/usr/local/
export TRIMMOMATIC_JAR=$SOFT_DIR/Trimmomatic-0.38/trimmomatic-0.38.jar


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
