#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process ALIGNBWAMEM {
  tag "${meta.id}"
  label 'process_high'  
  
  conda "$params.cacheDir/fpgAlign"
  publishDir "$params.bwaMem", mode: 'copy'

  input:
    val ready 
    // FIX: Accept meta map
    tuple val(meta), path(reads)

  output:
    // FIX: Emit meta map
    tuple val(meta), path("bwa/*.bam") , emit: aln_bam

  script:
  """
  #!/usr/bin/env bash

  if ! [[ -d bwa ]]; then mkdir bwa; fi

  index=`find -L "$workDir" -name "*.bwt" | head -n1 | sed 's/.bwt//'`
 
  # FIX: Use meta.id
  bwa mem -M -t ${task.cpus} \\
  \$index \\
  $reads \\
  | samtools sort --threads ${task.cpus} --output-fmt BAM -o bwa/${meta.id}.bam -
  """
}


process MARKDUPS {
  tag "${meta.id}"
  label 'process_high'
  
  conda "$params.cacheDir/fpgAlign"
  publishDir "$params.bwaMem", mode: 'copy'

  input:
    tuple val(meta), path(bam)

  output:
    tuple val(meta), path('*marked.bam'), emit: marked
    path('*_metrics.txt'), emit: dups_metric      
      
  script:
  """ 
   #!/usr/bin/env bash
   
   # FIX: Use meta.id for output filenames
   picard MarkDuplicates \\
   --INPUT $bam \\
   --OUTPUT ${meta.id}_marked.bam \\
   --METRICS_FILE ${meta.id}_metrics.txt \\
   --ASSUME_SORTED true \\
   --VALIDATION_STRINGENCY SILENT
  """
}


process SORTMARKED {
  tag "${meta.id}"
  label 'process_high'

  conda "$params.cacheDir/fpgAlign"
  publishDir "$params.bwaMem", mode: 'copy'

  input:
    tuple val(meta), path(bam)

  output:
    tuple val(meta), path('*sorted_marked.bam'), emit: sorted

  script:
  """
   #!/usr/bin/env bash
   
   # FIX: Use meta.id
   samtools sort --threads ${task.cpus} --output-fmt BAM "$bam" -o ${meta.id}_sorted_marked.bam
  """
}


process SAMINDEX {
  tag "${meta.id}"
  label 'process_high'

  conda "$params.cacheDir/fpgAlign"
  publishDir "$params.bwaMem", mode: 'copy'

  input:
    tuple val(meta), path(bam)

  output:
    tuple val(meta), path('*.bai'), emit: bai

  script:
  """
   #!/usr/bin/env bash

   samtools index -@ ${task.cpus} $bam -o ${bam}.bai
  """
}

// --- WORKFLOWS ---
// Bypassed in main script, but updated for consistency

workflow ALN {
  take:
    vl
    cln
    
  main:
     ALIGNBWAMEM(vl, cln)
     
  emit:
    bam = ALIGNBWAMEM.out.aln_bam
}
