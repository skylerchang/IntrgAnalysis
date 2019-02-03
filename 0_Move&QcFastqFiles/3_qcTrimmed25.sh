#qc

targetDir="FASTQ_Generation_2018-12-19_19_32_44Z-144643634"

mkdir $targetDir"/Trimmed25QC/"

for f in $targetDir"/Trimmed25/Paired/*"
do 
	echo $f
	fastqc $f --outdir=$targetDir"/Trimmed25QC/"
done

echo "Running MultiQC ..."
multiqc $targetDir"/Trimmed25QC/"
mkdir $targetDir"/Trimmed25QCMulti/"
mv multi* $targetDir"/Trimmed25QCMulti/"
