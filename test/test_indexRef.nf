#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.refIndex = "$params.resultsDir/index"

process INDEXREF {
  conda "$params.cacheDir/fpgAlign"
  publishDir "$params.refIndex", mode: 'copy'

  input:
  tuple file(ref)

  output:
   path(bwa) , emit: indx

  script:
  """
  #!/usr/bin/env bash
  
  if ! [[ -d bwa ]]; then mkdir bwa; fi

  bwa index -p bwa/${ref.simpleName} $ref

  """
}

workflow {
  INDEXREF(file("$params.refseq")) 
}



