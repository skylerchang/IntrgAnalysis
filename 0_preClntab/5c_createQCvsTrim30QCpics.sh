#unzip & move qc results


targetFolder2="../../Results/Trimmed30QC/*.zip"
targetFolder1="../../Results/OriginalQC/*.zip"


rm "files.txt"
touch "files.txt"


for dir1 in $targetFolder1
do
d1=`echo $dir1 | sed 's/\.zip$//'`
filename1="${d1##*/}"
path1=`echo $d1"/"$filename1"/Images/per_base_quality.png"`


for dir2 in $targetFolder2
do
d2=`echo $dir2 | sed 's/\.zip$//'`
filename2="${d2##*/}"
path2=`echo $d2"/"$filename2"/Images/per_base_quality.png"`
if [ $filename1 != $filename2 ]
then
continue
fi


echo $path1 $path2 >> files.txt
done
done

file=$(cat files.txt)

montage -label "%l" $file -font FreeMono-Bold-Oblique -tile 4x4 -page letter -geometry '360x360+4+4>' -pointsize 10 ../../Data/CountPlots/montage.pdf
convert -density 300 -trim ../../Data/CountPlots/montage.pdf -quality 100 ../../Data/CountPlots/montage_small.jpg
cd ../../Data/CountPlots/
convert `ls -v *.jpg` montage_small.pdf
