#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process ALIGNBWAMEM {
  //tag "$sampleId"
  
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
  
  //tag "$sampleId"

  conda "$params.cacheDir/fpgAlign"
  publishDir "$params.bwaMem", mode: 'copy'

  input:
    tuple val(sampleId), file(bam)


  output:
    //tuple val(sampleId), path("picard/*marked.bam"), emit: marked
    //tuple val(sampleId), path("picard/*_metrics.txt"), emit: dups_metric

      tuple val(sampleId), path('*marked.bam'), emit: marked
      path('*_metrics.txt'), emit: dups_metric      
     
  script:

  """ 
   #!/usr/bin/env bash
   
 
   picard MarkDuplicates \\
   --INPUT $bam \\
   --OUTPUT ${sampleId}_marked.bam \\
   --METRICS_FILE ${sampleId}_metrics.txt \\
   --ASSUME_SORTED true \\
   --VALIDATION_STRINGENCY SILENT


  """
}



process SORTMARKED {
  //tag "$bam"
  //tag "$sampleId"

  cpus 10
  executor 'slurm'

  conda "$params.cacheDir/fpgAlign"
  publishDir "$params.bwaMem", mode: 'copy'

  input:
    //file(bam)
    tuple val(sampleId), path(bam)


  output:
      path('*sorted_marked.bam'), emit: sorted

  script:

  """
   #!/usr/bin/env bash
   
   samtools sort --threads ${task.cpus} --output-fmt BAM "$bam" -o ${sampleId}_sorted_marked.bam
   
  """
   
}



process SAMINDEX {
  //tag "$bam"
  //tag "$sampleId"
  
  cpus 10
  executor 'slurm'

  conda "$params.cacheDir/fpgAlign"
  publishDir "$params.bwaMem", mode: 'copy'

  input:
    file(bam)


  output:
      path('*.bai'), emit: bai

  script:

  """
   #!/usr/bin/env bash

   samtools index -@ ${task.cpus} $bam -o ${bam}.bai

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


workflow TALN {
  take:
    vl
    ch_in

  main:
     ALIGNBWAMEM(vl,ch_in )

  emit:
    bam = ALIGNBWAMEM.out.aln_bam

}
