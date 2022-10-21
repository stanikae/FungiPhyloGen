#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.rawReads = "$params.inputDir/raw" 
params.readsDir = Channel.fromPath("$params.rawReads/*.f*q.gz", checkIfExists: true)
params.mqcOut = "$params.resultsDir/multiqc/raw"  //$params.mqcTyp"
params.fqcOut = file("$params.resultsDir/fqc_raw")

//params.fqcOut.mkdirs()

params.threads = 2


// fastqc contam and adapter files
contaminants = "$params.cacheDir/trimReads/opt/fastqc*/Configuration/contaminant_list.txt"
adapters = "$params.cacheDir/trimReads/opt/fastqc*/Configuration/adapter_list.txt"



process fqc {
  conda "$params.cacheDir/trimReads"
  publishDir "$params.fqcOut", mode: 'copy'

  input:
    file cont
    file adp 
    path fq

  output:
     path("*.html"), emit: fqc_out
     path("*.zip"), emit: fqc_zip
   

  script:
  """
  #!/usr/bin/env bash
  THREADS=$params.threads
  fastqc -o "." --contaminants "$cont" --adapters "$adp" --threads \$THREADS "$fq"
  # fastqc -t 2 --outdir "." "$fq"
  """
}



process mqc {
  conda "$params.cacheDir/trimReads"
  publishDir "$params.mqcOut", mode: 'copy'


  input:
    path fq


  output:
    path "multiqc_report.html", emit: mqc_raw


  script:
  """
  multiqc -f --no-data-dir -o "." -d "$fq"
  """
}


workflow qc_ind {
    take: 
        cont
        adp
	fq
        fqq
    main:
        fqc(cont, adp, fq)
        mqc(fqq)
    emit:
        mqc.out

}


workflow {
    qc_ind(file("$contaminants"),file("$adapters"),Channel.fromPath("$params.rawReads/*.f*q.gz", checkIfExists: true),"$params.fqcOut")

}












