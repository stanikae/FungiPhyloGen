#!/usr/bin/env nextflow
nextflow.enable.dsl=2



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



workflow FASTQC {
    take: 
        cont
        adp
	fq
    main:
       fqc(cont, adp, fq)
    emit:
        fqc.out

}

