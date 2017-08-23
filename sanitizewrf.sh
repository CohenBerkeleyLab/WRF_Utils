#!/bin/sh
if [ "$1" = "-h" -o "$1" = "--help" -o "$#" -eq 0 ]; then
	echo "Usage: sanitizewrf.sh [files]"
	echo " Use to remove colons from WRF output file names."
	echo " Pass a list of files to rename. Will replace any"
	echo " colons with dashes."
	echo ""
	echo " Example: sanitizewrf.sh wrfout*"
	echo ""
	return 0
fi


for file in $@
do
	newfile=$(echo $file | sed 's/\:/-/g')
#	echo $file $newfile
	mv $file $newfile
done
