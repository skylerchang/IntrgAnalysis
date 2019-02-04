#qc

targetDir="../../Data"

mkdir $targetDir"/OriginalQC/"

for f in $targetDir"/Original/*"
do 
	fastqc $f --outdir=$targetDir"/OriginalQC/"
done

echo "Running MultiQC ..."
multiqc $targetDir"/OriginalQC/"
mkdir $targetDir"/OriginalQCMulti/"
mv multi* $targetDir"/OriginalQCMulti/"


