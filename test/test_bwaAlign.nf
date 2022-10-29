#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//params.bwaMem = "$params.resultsDir/align"


process ALIGNBWAMEM {
  tag "$sampleId"
  
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
 
  bwa mem -M -t 5 \\
  \$index \\
  $reads \\
  | samtools sort --threads 5 --output-fmt BAM -o bwa/${sampleId}.bam -

  """
}



process MARKDUPS {
  
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

    samtools index $bam -o picard/${bam}.bai

  """

}



workflow ALN {
  take:
    vl
    cln
    //dups

  main:
     ALIGNBWAMEM(vl,
       cln
          .map { file -> def key = file.name.toString().tokenize('_').get(0).replace('[','')
             return tuple(key, file)}
          //.groupTuple( )
     )
     //MARKDUPS(dups
      //           .map { file -> def key = file.name.replaceAll(/.bam|.bai$/,'')
      //                          return tuple(key, file) }
      //   )


  emit:
    bam = ALIGNBWAMEM.out.aln_bam
    //markdup = MARKDUPS.out.marked
    //dupmetric = MARKDUPS.out.dups_metric

}



workflow {
  idx_ch = Channel.fromPath("$projectDir/results/index/*.bwt", checkIfExists:true)
  clean_ch = Channel.fromPath("$projectDir/results/clean_reads/*_{val_1,val_2}.fq.gz", checkIfExists:true)
  ALN(idx_ch,clean_ch)
  //MARKDUPS(ALN.out.bam | map { file -> def key = file.name.replaceAll(".bam","") 
  //                             return tuple(key, file) } )
  MARKDUPS(ALN.out.bam)
  //            .map { file -> def key = file.name.toString().tokenize('.').get(0).replace('[','')
  //           return tuple(key, file)}
  //       )

  SAMINDEX(MARKDUPS.out.marked)
}



