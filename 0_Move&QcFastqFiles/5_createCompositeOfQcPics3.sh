#unzip & move qc results


targetFolder="../../Data/OriginalQC/*.zip"

#list=()

touch "list.txt"

for dir in $targetFolder
do
  d=`echo $dir | sed 's/\.zip$//'`
  filename="${d##*/}"
  path=`echo $d"/"$filename"/Images/per_base_quality.png"`
  echo $path >> list.txt
done

v=$(cat list.txt)
echo $v

while read snack; do
	montage -tile 4x8 -page letter -geometry +4+4 -pointsize 20 -label %t $snack test.pdf

done < $v

