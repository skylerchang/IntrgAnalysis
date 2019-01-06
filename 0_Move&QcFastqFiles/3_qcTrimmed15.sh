#qc

targetDir="FASTQ_Generation_2018-12-19_19_32_44Z-144643634"

mkdir $targetDir"/Trimmed15QC/"

for f in $targetDir"/Trimmed15/Paired/*"
do 
	echo $f
	fastqc $f --outdir=$targetDir"/Trimmed15QC/"
done

echo "Running MultiQC ..."
multiqc $targetDir"/Trimmed15QC/"
mkdir $targetDir"/Trimmed15QCMulti/"
mv multi* $targetDir"/Trimmed15QCMulti/"
