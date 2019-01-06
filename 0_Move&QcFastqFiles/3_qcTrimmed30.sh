#qc

targetDir="FASTQ_Generation_2018-12-19_19_32_44Z-144643634"

mkdir $targetDir"/Trimmed30QC/"

for f in $targetDir"/Trimmed30/Paired/*"
do 
	echo $f
	fastqc $f --outdir=$targetDir"/Trimmed30QC/"
done

echo "Running MultiQC ..."
multiqc $targetDir"/Trimmed30QC/"
mkdir $targetDir"/Trimmed30QCMulti/"
mv multi* $targetDir"/Trimmed30QCMulti/"
