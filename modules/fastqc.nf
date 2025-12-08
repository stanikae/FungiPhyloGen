#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process fqc {
  tag "${meta.id}"
  cpus 8
  executor 'slurm'
  errorStrategy 'ignore'

  // Update to the correct environment name
  conda "$params.cacheDir/fpgtrimReads"
  
  // Use the parameter passed from the main workflow
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
  
  # Ensure output directory for samples exists if needed (FastQC handles this usually)
  # Using meta.id to handle the map structure
  
  # Run FastQC
  # Note: ${reads} expands to both files in the list (R1 R2)
  
  fastqc \
    --outdir . \
    --contaminants "$cont" \
    --adapters "$adp" \
    --threads ${task.cpus} \
    ${reads}
  """
}


