 #!/bin/bash
 SOURCES=$(cat files.txt)
 # or SOURCES=$(find . -name "*.png")
 # Iterate over known files
 for FILENAME in files.txt
 do
    # Substring up to first "_" character
    LABEL=$(echo $FILENAME | cut -d "_" -f 1)
    # Set meta-data label
    mogrify -label $LABEL $FILENAME
 done
 montage -label %l $SOURCES show: