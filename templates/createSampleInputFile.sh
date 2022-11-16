#!/bin/bash

################################################
# Small script to create comma-separated input file
# Written by Stan Kwenda - NICD-SCF
# Date: 2020-10-08
# Modified: 2022-11-16
################################################
# Notes:
# pathDir ==> Absolute path to paired-end reads directory

pathDir=$1
outDir=$2

#  generate input file with header
echo "sampleID,read1,read2" > $outDir/sample_file.csv 

for name in $(ls $pathDir | grep R1)
  do
     nam=$(echo $name | cut -d_ -f1)
     paste -d "," <(echo $nam) <(find $pathDir -name "${nam}*R1*fastq*") <(find $pathDir -name "${nam}*R2*fastq*")
done >> $outDir/sample_file.csv
