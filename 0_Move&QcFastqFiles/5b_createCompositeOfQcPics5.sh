#unzip & move qc results



targetFolder="../../Data/Trimmed30QC/*.zip"
targetFolder="../../Data/OriginalQC/*.zip"

rm "files.txt"
touch "files.txt"

rm "names.txt"
touch "names.txt"

for dir in $targetFolder
do
  d=`echo $dir | sed 's/\.zip$//'`
  echo $d
  filename="${d##*/}"
  echo $filename
  echo $filename >> names.txt
  path=`echo $d"/"$filename"/Images/per_base_quality.png"`
  echo $path >> files.txt
done

file=$(cat files.txt)
echo $file

n=$(cat names.txt)
echo $n


montage $file -tile 4x8 -page letter -geometry +4+4 -pointsize 30 ../../Data/CountPlots/QC_thumbnails_original.pdf
#montage $file -tile 4x8 -page letter -geometry +4+4 -pointsize 30 ../../Data/CountPlots/QC_thumbnails_trimmed.pdf

#montage -label $f $f -tile 4x8 -page letter -geometry +4+4 -pointsize 20 test.pdf



