#unzip & move qc results


targetFolder="../../Data/OriginalQC/*.zip"
targetFolder="../../Data/Trimmed30QC/*.zip"



for zip in $targetFolder
do
  dirname=`echo $zip | sed 's/\.zip$//'`
  if mkdir "$dirname"
  then
    if cd "$dirname"
    then
      unzip ../"$zip"
      cd ..
      # rm -f $zip # Uncomment to delete the original zip file
    else
      echo "Could not unpack $zip - cd failed"
    fi
  else
    echo "Could not unpack $zip - mkdir failed"
  fi
done
