#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process fqc {
  tag "${meta.id}"
  label 'process_medium'
  errorStrategy 'ignore'

  conda "$params.cacheDir/fpgtrimReads"
  
  publishDir "$params.fqcOut", mode: 'copy'

  input:
    path cont
    path adp
    // Updated to accept [meta, [reads]] structure
    tuple val(meta), path(reads)

  output:
    path "*.html", emit: html
    path "*.zip" , emit: zip

  script:
  """
  #!/usr/bin/env bash
  
  fastqc \
    --outdir . \
    --contaminants "$cont" \
    --adapters "$adp" \
    --threads ${task.cpus} \
    ${reads}
  """
}


