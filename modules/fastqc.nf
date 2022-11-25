#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//reads_ch = Channel.fromPath("$params.readsDir/**.f*q.gz", checkIfExists: true) 
//params.readsDir = Channel.fromPath("$params.rawReads/**.f*q.gz", checkIfExists: true)
//params.mqcOut = "$params.resultsDir/multiqc/raw"  //$params.mqcTyp"
//params.fqcOut = file("$params.resultsDir/fqc_raw")
//params.threads = 2


// fastqc contam and adapter files
contaminants = "$params.cacheDir/trimReads/opt/fastqc*/Configuration/contaminant_list.txt"
adapters = "$params.cacheDir/trimReads/opt/fastqc*/Configuration/adapter_list.txt"



process fqc {
  cpus 8
  executor 'slurm'

  conda "$params.cacheDir/trimReads"
  publishDir "$params.fqcOut", mode: 'copy'

  input:
    file cont
    file adp 
    path fq
    //file fq
   

  output:
    // path("*.html"), emit: fqc_out
    // path("*.zip"), emit: fqc_zip
   file "*"

  script:
  """
  #!/usr/bin/env bash
  THREADS=${task.cpus} #//$params.threads
  fastqc --contaminants "$cont" --adapters "$adp" --threads \$THREADS $fq
  
  """
}



workflow qc_ind {
    take: 
        cont
        adp
	fq
    main:
       fqc(cont, adp, fq)
    emit:
        fqc.out

}


workflow {
    qc_ind(file("$contaminants"),file("$adapters"),reads_ch)    
}




