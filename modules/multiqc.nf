#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process mqc {
  label 'process_medium'

  conda "$params.cacheDir/fpgtrimReads"
  publishDir "$params.mqcOut", mode: 'copy'

  input:
  path multiqc_files

  output:
  path "*.html"      , emit: report
  path "*_data"      , emit: data, optional: true

  script:
  """
  #!/usr/bin/env bash
  
  multiqc . \
    --force \
    --no-data-dir \
    --filename "${params.fileN}_multiqc_report.html"
  """
}
