#qc

targetDir="../../Data"

mkdir $targetDir"/OriginalQC/"

for d in $targetDir"/Original/*"
do 
	echo $d
	for f in $d/*
	do
		fastqc $f --outdir=$targetDir"/OriginalQC/"
	done
done

echo "Running MultiQC ..."
multiqc $targetDir"/OriginalQC/"
mkdir $targetDir"/OriginalQCMulti/"
mv multi* $targetDir"/OriginalQCMulti/"


