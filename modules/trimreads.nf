#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.trm = "$params.resultsDir/clean_reads"
//params.index = "$params.inputDir/sample_file.csv"
//index = "$params.inputDir/sample_file.csv"
//index = "/home/stan/git-repos/FungiPhyloGen/test/sample_file.csv"

process trimReads {
  cpus 6
  executor 'slurm'

  conda "$params.cacheDir/trimReads"
  publishDir "$params.trm", mode: 'copy'

  input:
  tuple val(sampleID), file(read1), file(read2) 

  output:
   //tuple val(sampleID), path("*val*.fq.gz"), emit: cln_reads
   path('*val*.fq.gz')

  script:
  """
  #!/usr/bin/env bash
  
  trim_galore -q 30 \
   --length 50 --trim-n -o "." --gzip \
   --cores ${task.cpus} --paired "$read1" "$read2" \
   --no_report_file --basename "$sampleID"
 
  """
}

workflow trim {
  take: infile
  main:
    infile | splitCsv(header:true) | map { row-> tuple(row.sampleID, file(row.read1), file(row.read2)) } | trimReads
  emit:
    rds = trimReads.out
}


workflow {
   
  trim(Channel.fromPath(params.index))


}









