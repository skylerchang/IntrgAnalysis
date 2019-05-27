
module load fastqc
module load python/3.6.3
pip install --user multiqc
mkdir ../../Results



#once installed, go from here
targetDir="../../"

mkdir $targetDir"Results/OriginalQC/"

for d in $targetDir"Data/Original/*"
do 
	echo $d
	for f in $d/*
	do
		echo "fastqc-ing "$f" ..."
		fastqc $f --outdir=$targetDir"Results/OriginalQC/"
	done
done

echo "Running MultiQC ..."
multiqc $targetDir"Results/OriginalQC/"
mkdir $targetDir"Results/OriginalQCMulti/"
mv multi* $targetDir"Results/OriginalQCMulti/"
