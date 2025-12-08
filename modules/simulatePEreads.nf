#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// params.cacheDir = "/scratch/package/anaconda/anaconda3_2023-02-27/envs"
params.cacheDir = "/spaces/stanford/anaconda3/envs"

/*
==================================================================
        CHECK INPUTS
==================================================================
*/

// Check mandatory parameters
if (params.refsDir) {
        ref_ch = Channel.fromPath("$params.refsDir/**/G*_genomic.fna", checkIfExists: true)
} else {
        exit 1, 'Please provide path to references directory ...'
}




/*
=================================================================
        Create directories - OPTIONAL
=================================================================
*/

params.inDir = file("$params.resultsDir")
params.inDir.mkdirs() ? ! params.inDir.exists() : 'Directory already exists'




/*
==================================================================
	  RUN PROGRAM
==================================================================
*/


process SIMULATEREADS {
 
  cpus 8
  executor 'slurm'
   
  conda "$params.cacheDir/fpgCallVariants"
  publishDir "$params.inDir", mode: 'copy'
  
  input:
    path(ref)

  output:
    path('*_001.fastq*')

  script:


  """
  #!/usr/bin/env bash

  name=\$(basename -s _genomic.fna "$ref") 
  wgsim -1 300 -2 300 -r 0 -R 0 -X 0 -e 0 "$ref" \${name}_R1_001.fastq \${name}_R2_001.fastq
  gzip \${name}_R1_001.fastq \${name}_R2_001.fastq

  """
}



process PEGZIP {

 publishDir "$params.inDir", mode: 'copy'

   input:
     path(reads, stageAs: "?/*")
     // path '*_val*.fq'
     // path '*_val_1.fq'

   output:
     path('*val*.fq.gz')

   script:


   """
   #!/usr/bin/env bash
  
   gzip "$reads" 

   """
}



workflow {
    SIMULATEREADS(ref_ch)
   // PEGZIP(SIMULATEREADS.collect())
}


