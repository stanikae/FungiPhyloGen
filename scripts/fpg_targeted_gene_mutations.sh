#!/bin/bash

snpeffVCF=$1
geneList=$2
outVCF=$3
outFile=$4


jar_path=/home/stan/anaconda3/envs/snpsift/share/snpsift-5.1-0

cat $snpeffVCF | \
	java -jar $jar_path/SnpSift.jar filter --set $geneList "(ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE in SET[0]) " > $outVCF


conda activate fpgCallVariants

SAMPLESVCF=$(bcftools query -l $outVCF | tr '\n' '|')

echo -e "CONTIG|POSITION|TYPE|QUALITY_SCORE|REF|ALT|ALLELE|ANNOTATION|PUTATIVE_IMPACT|GENENAME|GENEID|FEATURETYPE|FEATUREID|TRANSCRIPTBIOTYPE|RANK|HGVS.c|HGVS.p|cDNAposition|CDSposition|PROTEINposition|DISTANCE2FEATURE|ERRORS|$SAMPLESVCF" > filtered_header.txt

bcftools query -f '%CHROM|%POS|%TYPE|%QUAL|%REF|%ALT|%INFO/ANN{0}|[%GT|]\n' $outVCF | grep "missense_variant" >> filtered_variants.txt

#cat filtered_header.txt filtered_variants.txt > FKS1_ERG11_c_auris.txt
cat filtered_header.txt filtered_variants.txt > $outFile
