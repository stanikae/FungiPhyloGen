#!/usr/bin/env nextflow
nextflow.enable.dsl=2



process REHEADERVCF {
// tag "$sampleId"

  cpus 8
  executor 'slurm'

  conda "$params.cacheDir/fpgCallVariants"
  publishDir "$params.ann", mode: 'copy'


  input:
    path(vcf)

  output:
    path("snpeff/snpeff_ann_reheader.vcf"), emit: reh


  script:


  """
  #!/usr/bin/env bash

  if ! [[ -d snpeff ]]; then mkdir -p snpeff; fi

  # rehead (rename) samples i.e. remove paths and suffix
  paste -d "\t" <(bcftools query -l "$vcf") <(bcftools query -l "$vcf" | sed 's/\\/.*\\///' | sed 's/_sorted.*//') > ./reheader.txt
  bcftools reheader -s ./reheader.txt -o snpeff/snpeff_ann_reheader.vcf "$vcf"


  """
}




