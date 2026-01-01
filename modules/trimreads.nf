#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process trimReads {
  tag "${meta.id}"
  label 'process_medium'

  conda "$params.cacheDir/fpgtrimReads" 
  publishDir "$params.trm", mode: 'copy'

  input:
  tuple val(meta), path(reads) 

  output:
  tuple val(meta), path("*.fq.gz"), emit: rds

  script:
  """
  #!/usr/bin/env bash
  
  
  trim_galore -q 30 \
   --length 50 --trim-n -o "." --gzip \
   --cores ${task.cpus} --paired "${reads[0]}" "${reads[1]}" \
   --no_report_file --basename "${meta.id}"
  
  # Renaming for consistency
  
  mv ${meta.id}_val_1.fq.gz ${meta.id}_1.trimmed.fq.gz
  mv ${meta.id}_val_2.fq.gz ${meta.id}_2.trimmed.fq.gz
  """
}

workflow trim {
  take: 
    ch_input // Expects: [ [id:sample], [read1, read2] ]

  main:
    trimReads(ch_input)

  emit:
    rds = trimReads.out.rds
}




