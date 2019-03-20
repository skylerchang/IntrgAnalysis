#qc

targetDir="../../"

mkdir $targetDir"Results/Trimmed30QC/"

for f in $targetDir"Data/Trimmed30/Paired/*"
do 
	echo $f
	fastqc $f --outdir=$targetDir"Results/Trimmed30QC/"
done

echo "Running MultiQC ..."
multiqc $targetDir"Results/Trimmed30QC/"
mkdir $targetDir"Results/Trimmed30QCMulti/"
mv multi* $targetDir"Results/Trimmed30QCMulti/"
