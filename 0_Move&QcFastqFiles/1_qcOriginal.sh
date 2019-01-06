#qc

targetDir="FASTQ_Generation_2018-12-19_19_32_44Z-144643634"

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


