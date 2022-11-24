#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//params.rawReads = "$params.inputDir/raw" 
//params.readsDir = Channel.fromPath("$params.rawReads/**.f*q.gz", checkIfExists: true)
//params.mqcOut = "$params.resultsDir/multiqc/raw"  //$params.mqcTyp"
//params.fqcOut = Channel.fromPath("$params.resultsDir/fqc_raw")
//params.fqcOut = file("$params.resultsDir/fqc_raw") 

params.threads = 2


// fastqc contam and adapter files
//contaminants = "$params.cacheDir/trimReads/opt/fastqc*/Configuration/contaminant_list.txt"
//adapters = "$params.cacheDir/trimReads/opt/fastqc*/Configuration/adapter_list.txt"



process mqc {
  conda "$params.cacheDir/trimReads"
  publishDir "$params.mqcOut", mode: 'copy'


  input:
   // path fq
   file fq

  output:
   //path "multiqc_report.html", emit: mqc_raw
   file "*.html"


  script:
  """
  #!/bin/env bash
  nam=$params.fileN
  #multiqc -f --no-data-dir -o "." "$fq"
  multiqc -f --no-data-dir --filename "\$nam" $fq 
 
  """
}
//multiqc -f --no-data-dir -o "." -d "$fq"



workflow {
  mqc(file("$params.resultsDir/fqc_raw"))

}


/*
workflow qc_ind {
    take: 
        cont
        adp
	fq
      //  fqq
    main:
       // fqc(cont, adp, fq)
       // mqc(fqq)
       fqc(cont, adp, fq) | collect | mqc
    emit:
        mqc.out

}


workflow {
    // dummy_ch = Channel.empty()
    //qc_ind(file("$contaminants"),file("$adapters"),Channel.fromPath("$params.rawReads/*.f*q.gz", checkIfExists: true),"$params.fqcOut")
    qc_ind(file("$contaminants"),file("$adapters"),Channel.fromPath("$params.rawReads/*.f*q.gz", checkIfExists: true))    
    //qc_ind(file("$contaminants"),file("$adapters"),Channel.fromPath("$params.rawReads/*.f*q.gz", checkIfExists: true),dummy_ch)
}
*/











