#!/bin/bash

################################################
# Small script to create comma-separated input file
# Written by Stan Kwenda - NICD-SCF
# Date: 2020-10-08
# Modified: 2022-11-16
################################################

# Notes:
# pathDir ==> Absolute path to paired-end reads directory
# outFile ==> Absolute path to output file

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
