#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//outdir = '/home/stan/fungal-genomics/test-data/clean'
readsDir = Channel.fromFilePairs('/home/stan/fungal-genomics/test-data/raw/*_{R1_001,R1_002}*.f*q.gz')


//Channel
//    .fromFilePairs('/home/stan/fungal-genomics/test-data/raw/*_{R1_001,R1_002}*.f*q.gz')
//    .view()



process trimReads {
   conda "/home/stan/anaconda3/envs/trimReads"

   input:
    val f1
    val f2

  output:
    path trimmedReads

  
  """
  echo ${params.CONDA_BASE}"
  """
}


workflow {
 readsDir.view()
 trimReads()
 // splitLetters | flatten | convertToUpper | view { it.trim() }
}


