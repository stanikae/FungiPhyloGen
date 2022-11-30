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
    tuple val(sampleID), file(read1), file(read2)
    //path fq
    //file fq
   

  output:
    // path("*.html"), emit: fqc_out
    // path("*.zip"), emit: fqc_zip
   file "*"

  script:
  """
  #!/usr/bin/env bash
  THREADS=${task.cpus} 
  nam=\$(echo "$sampleID" | cut -d_ -f1-2 | tr "_" "-")
  if ! [[ -d \$nam ]]; then mkdir -p \$nam ; fi
  fastqc -o \$nam --contaminants "$cont" --adapters "$adp" --threads \$THREADS "$read1" "$read2"
  
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
    fq_ch = index_ch.splitCsv(header:true).map { row-> tuple(row.sampleID, file(row.read1), file(row.read2)) }
    qc_ind(file("$contaminants"),file("$adapters"),fq_ch)    
}




