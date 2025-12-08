#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process trimReads {
  tag "${meta.id}"
  cpus 6
  executor 'slurm'

  // Update this path if necessary
  conda "$params.cacheDir/fpgtrimReads" 
  publishDir "$params.trm", mode: 'copy'

  input:
  // We accept a meta map and a list of reads
  tuple val(meta), path(reads) 

  output:
  // We output the same structure: [ [id:sample], [trimmed_R1, trimmed_R2] ]
  tuple val(meta), path("*.fq.gz"), emit: rds

  script:
  """
  #!/usr/bin/env bash
  
  # Note: ${reads[0]} is R1, ${reads[1]} is R2
  
  trim_galore -q 30 \
   --length 50 --trim-n -o "." --gzip \
   --cores ${task.cpus} --paired "${reads[0]}" "${reads[1]}" \
   --no_report_file --basename "${meta.id}"
  
  # Renaming for consistency (Optional but recommended)
  # trim_galore adds weird suffixes like _val_1.fq.gz. 
  # This makes them standard: sampleID_1.trimmed.fq.gz
  
  mv ${meta.id}_val_1.fq.gz ${meta.id}_1.trimmed.fq.gz
  mv ${meta.id}_val_2.fq.gz ${meta.id}_2.trimmed.fq.gz
  """
}

workflow trim {
  take: 
    ch_input // Expects: [ [id:sample], [read1, read2] ]

  main:
    // Just pass the channel directly to the process
    trimReads(ch_input)

  emit:
    rds = trimReads.out.rds
}




