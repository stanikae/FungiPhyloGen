#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process ALIGNBWAMEM {
  tag "$sampleId"
  
  cpus 8
  executor 'slurm'
  
  conda "$params.cacheDir/fpgAlign"
  publishDir "$params.bwaMem", mode: 'copy'


  input:
    val ready
    tuple val(sampleId), file(reads)

  output:
    tuple val(sampleId), path("bwa/*.bam") , emit: aln_bam

  script:
  
  """
  #!/usr/bin/env bash
  
  if ! [[ -d bwa ]]; then mkdir bwa; fi

  index=`find -L "$workDir" -name "*.bwt" | head -n1 | sed 's/.bwt//'`
 
  bwa mem -M -t ${task.cpus} \\
  \$index \\
  $reads \\
  | samtools sort --threads ${task.cpus} --output-fmt BAM -o bwa/${sampleId}.bam -

  """
}



process MARKDUPS {

  cpus 12
  executor 'slurm'
  
  tag "$sampleId"

  conda "$params.cacheDir/fpgAlign"
  publishDir "$params.bwaMem", mode: 'copy'

  input:
    tuple val(sampleId), file(bam)


  output:
    //tuple val(sampleId), path("picard/*marked.bam"), emit: marked
    //tuple val(sampleId), path("picard/*_metrics.txt"), emit: dups_metric

      path('picard/*marked.bam'), emit: marked
      path('picard/*_metrics.txt'), emit: dups_metric      
     
  script:

  """ 
   #!/usr/bin/env bash
   
    if ! [[ -d picard ]]; then mkdir picard; fi
 
   picard MarkDuplicates \\
   --INPUT $bam \\
   --OUTPUT picard/${sampleId}_sorted_marked.bam \\
   --METRICS_FILE picard/${sampleId}_metrics.txt \\
   --ASSUME_SORTED true \\
   --VALIDATION_STRINGENCY SILENT


  """
}



process SAMINDEX {
  tag "$bam"
  //tag "$sampleId"
  
  cpus 10
  executor 'slurm'

  conda "$params.cacheDir/fpgAlign"
  publishDir "$params.bwaMem", mode: 'copy'

  input:
    file(bam)


  output:
      path('picard/*.bai'), emit: bai

  script:

  """
   #!/usr/bin/env bash

    if ! [[ -d picard ]]; then mkdir picard; fi

    samtools index -@ ${task.cpus} $bam -o picard/${bam}.bai

  """

}



workflow ALN {
  take:
    vl
    cln
    
  main:
     ALIGNBWAMEM(vl,
       cln
          .map { file -> def key = file.name.toString().tokenize('_').get(0).replace('[','')
             return tuple(key, file)}
          //.groupTuple( )
     )
     
  emit:
    bam = ALIGNBWAMEM.out.aln_bam
    
}


