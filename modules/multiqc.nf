#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process mqc {
  // Use the same executor strategy as the rest
  cpus 4
  executor 'slurm'

  // --- FIX 1: Correct Environment Name ---
  conda "$params.cacheDir/fpgtrimReads"
  publishDir "$params.mqcOut", mode: 'copy'

  input:
  // We receive a list of files from the .collect() in the main workflow
  path multiqc_files

  output:
  // --- FIX 2: Named outputs for better tracking ---
  path "*.html"      , emit: report
  path "*_data"      , emit: data, optional: true

  script:
  """
  #!/usr/bin/env bash
  
  # Run MultiQC on the current directory (.) where input files are staged
  multiqc . \
    --force \
    --no-data-dir \
    --filename "${params.fileN}_multiqc_report.html"
  """
}
