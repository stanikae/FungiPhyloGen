#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process mqc {

  conda "$params.cacheDir/trimReads"
  publishDir "$params.mqcOut", mode: 'copy'


  input:
   // path fq
   file fq

  output:
   //path "multiqc_report.html", emit: mqc_raw
   file "*.html"


  script:
  """
  #!/bin/env bash
  nam=$params.fileN
  #multiqc -f --no-data-dir -o "." "$fq"
  multiqc -f --no-data-dir --filename "\$nam" $fq 
 
  """
}


