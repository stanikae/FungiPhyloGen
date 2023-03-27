#!/bin/bash



sampleList=/scratch/sysuser/stanford/git-repos/FungiPhyloGen/sample_file.csv
refseq=~/test-data/ref-genomes/ncbi_dataset/data/GCA_000150115.1/GCA_000150115.1_ASM15011v1_genomic.fna
gbk=~/test-data/ref-genomes/ncbi_dataset/data/GCA_000150115.1/genomic.gbff
confg=/scratch/sysuser/stanford/git-repos/FungiPhyloGen/fpg_test.config
nextflow run /scratch/sysuser/stanford/git-repos/FungiPhyloGen/main_test.nf --refseq $refseq --gbk $gbk --prjName "TEST-FASTQC" -c $confg -with-conda true
