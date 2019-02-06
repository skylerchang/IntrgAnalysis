#unzip & move qc results


targetFolder="../../Data/OriginalQC/*.zip"


file=$(cat files.txt)
echo $file

n=$(cat names.txt)
echo $n


montage -label $n $file -tile 4x8 -page letter -geometry +4+4 -pointsize 30 test.pdf



