#!/bin/bash
#
for file in *; do
   if [ -d $file ]; then
      cd $file
      mv ${file}_AciTot.txt ${file}_AciTot.csv
      mv ${file}_AciIfTot.txt ${file}_AciIfTot.csv
	rm *.txt
	sed -i -e '$a\' ${file}_AciTot.csv
	sed -i -e '$a\' ${file}_AciIfTot.csv
      cd ..
   fi
done
