#!/bin/bash

################################################
# Small script to create comma-separated input file
# Written by Stan Kwenda - NICD-SCF
# Date: 2020-10-08

# Change log
# Modified: 2022-11-16
# 	    2023-06-17: 
#		     - Combine 1st and 2nd script (fastq renamer) into one script
#		     - Convert 1st script into a function
#		     - See below for user input
################################################

# Notes:
# pathDir ==> Absolute path to paired-end reads directory
# outFile ==> Absolute path to output file (preferably create this file in your project directory - work directory)

set -e
set -x

# Usage information
if [ $# != 2 ];
 then
        echo "Usage: `basename $0` [path to fastq files] [path for new sample index to be created by this script]"
        exit
fi


declare -x SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
echo "$SCRIPTPATH"

# create function
create_index () {
	pathDir=$1
	outFile=$2

	#  generate input file with header
	echo "sampleID,read1,read2" > $outFile

	for name in $(ls $pathDir | grep R1)
	  do
     		#nam=$(echo $name | cut -d_ -f1)
     		nam=$(echo $name | awk -F '_R{1,2}' '{print $1}')
     		paste -d "," <(echo $nam) <(find $pathDir -name "${nam}*R1*fastq*") <(find $pathDir -name "${nam}*R2*fastq*")
	done >> $outFile

}

# Execute function
create_index "$1" "$2"


# Remove underscores
Rscript $SCRIPTPATH/rename_sample_names.R "$2"


# create final sample sheet (index)
indexDir=$(dirname "$2")
newIndex="$indexDir"/sample_index_new.csv
create_index "$1" "$newIndex"


