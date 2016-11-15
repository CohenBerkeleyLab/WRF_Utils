#!/bin/bash
#
# This script allows easy subsetting of wrfout files after the fact by picking out
# specific variables.
#
#   Josh Laughner (joshlaugh5@gmail.com) 14 Nov 2014

# Get the location of this script, which will use to find calculated_quantities.nco if needed
currdir=$(pwd -P)
linkpath=$(readlink "$0")
tmppath=$(dirname "$linkpath")
cd "$tmppath"
scriptdir="$(pwd -P)"
cd "$currdir"
__me=$(basename $0)
# Standard variables to include
stdvars="Times XLONG XLAT no2 U V COSALPHA SINALPHA calc"

# Help text
usage="$(basename $0) -- allows specific variables to be extracted from wrfout files into wrfout_subset files

This program looks for wrfout files in the current directory and extracts a subset of variables from them.
If given with no additional arguments, this will extract the following variables by default:
    $stdvars

Otherwise, the argument specify the variable names (as seen in using ncview or ncdump) that you want to extract.
There are two exceptions:
    std - will include the default variables listed above. This way, you can add to that list of variables.
    calc - will include calculated quantities including temperature (TT), number density of air (ndens),
        NO2 number density (no2_ndens), pressure (pres), altitude (z), and box height (zlev). Note that 
        these are included in std.

If using calculated quantities, it is important that $__me reside in the same directory as 
calculated_quantities.nco (which should be in the same Git repo). Therefore, it is usually best to put a
symbolic link to $__me in the directory with the wrfout files rather than copying it.

Examples:
    ./$__me - generate wrfout_subset files containing the standard variables.
    
    ./$__me U V - generate wrfout_subset files containing only the U and V variables.

    ./$__me std no hno3 - generate wrfout_subset files containing the standard variables plus no and hno3.
"



# Parse the input
if [[ $# == 0 ]]; then
    subvars="$stdvars"
else
    subvars=""
    while [[ $# > 0 ]]
    do
    keyin="$1"
    # Ensure input is lower case
    key=$(echo $keyin | awk '{print tolower($0)}')
        case $key in
            -h|--help)
                echo "$usage"
                exit 0
                ;;
            'std')
            for v in $stdvars; do
                regex='[[:space:]]+'"$v"'[[:space:]]|^'"$v"'[[:space:]]|[[:space:]]'"$v"'$|^'"$v"'$'
                if [[ $subvars =~ $regex ]]; then
                    echo "Duplicate variable: $v, will only output once"
                else
                    subvars="$subvars $v"
                fi
            done
            ;;  
            *) # anything else should be a variable name
            regex='[[:space:]]+'"$key"'[[:space:]]|^'"$key"'[[:space:]]|[[:space:]]'"$key"'$|^'"$key"'$'
            if [[ $subvars =~ $regex ]]; then
                echo "Duplicate variable: $key, will only output once"
            else
                subvars="$subvars $key"
            fi
            ;;  
        esac
        
        shift # shift the input arguments left by one
    done
fi

# If we are using calculated quantities, then make sure we can find the calculated_quantities.nco
# file. It is in the same directory as this file, but if this is not a link, then we won't be
# pointing to the right directory
regex='[[:space:]]+calc[[:space:]]|^calc[[:space:]]|[[:space:]]calc$|^calc$'
if [[ $subvars =~ $regex ]]; then
    if [[ ! -f "$scriptdir/calculated_quantities.nco" ]]; then
        echo "Cannot find the script calculated_quantities.nco" >&2
        echo "in the directory $scriptdir" >&2
        echo "$0 should be a link to the original in the directory" >&2
        echo "containing calculated_quantities.nco as well" >&2
        exit 1
    fi
fi

# Now loop over all wrfout files in the current directory. Create a new subset file
# that contains only the specified variables.
files=wrfout*
for f in $files; do
    echo "Subsetting $f"
    savename=$(echo "$f" | sed -e 's/wrfout/wrfout_subset/')
    for var in $subvars; do
        echo "    Variable $var"
        if [[ $var == "calc" ]]; then
            ncap2 -A -v -S $scriptdir/calculated_quantities.nco $f $savename
        else
            ncks -A -v "$var" $f $savename
        fi
    done
done
