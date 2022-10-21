#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.rawReads = "$params.inputDir/raw" //'/home/stan/fungal-genomics/test-data/raw'
params.mqcOut = "$params.resultsDir"
params.trm = "$params.resultsDir/clean_reads"
seqReads = Channel.fromFilePairs("$params.rawReads/*{R1,R2}*.f*q.gz", checkIfExists: true)
params.index = "$params.inputDir/sample_file.csv"




process readsQC {
  conda "$params.cacheDir/qcReads"
  //publishDir "$params.mqcOut", mode: 'copy'

  input:
  tuple val(sampleID), file(read1), file(read2)

  output:
    tuple val(sampleID), path("*"), emit: reads



  """
    FaQCs -1 "$read1" -2 "$read2" \
    --avg_q 30 --phiX --prefix "$sampleID" \
    --qc_only --stats sampleID.stats.txt \
    -d "." \
    -p raw.1.fq raw.2.fq


  """


}


//FaQCs -1 ../raw/811-ND_S4_L001_R1_001.fastq.gz -2 ../raw/811-ND_S4_L001_R2_001.fastq.gz --qc_only --prefix "811-ND" -t 4 --avg_q 30 -d faqcs









process baseCounts {
  conda "$params.cacheDir/trimReads"
  //conda "/home/stan/anaconda3/envs/trimReads" //"$params.cacheDir/trimReads"

  input:
    path fq

  output:
   // path "${fq.SimpleName}.mqc.tsv"
    path "${fq.SimpleName}.txt"
  
  
 """
 seqfu stats -b --csv --gc --sort-by filename $fq > ${fq.SimpleName}.txt
 
 """
}


process mqcPrep {
  conda "$params.cacheDir/trimReads"

  input:
    path fq

  output:
    path "${fq.SimpleName}.mqc.tsv"
  
  script:
  """
  seqfu stats -b --csv --gc --sort-by filename $fq --multiqc "${fq.SimpleName}.mqc.tsv"
  """

}


process mqc {
  conda "$params.cacheDir/trimReads"
  publishDir "$params.mqcOut", mode: 'copy'


  input:
    path x


  output:
    path "raw_data/multiqc_report.html", emit: mqc_raw


  script:
  """
  multiqc -f --no-data-dir -o "./raw_data" "$x"
  """
}


process trimReads {
  conda "$params.cacheDir/trimReads"
  publishDir "$params.trm", mode: 'copy'

  input:
  tuple val(sampleID), file(read1), file(read2) 

  output:
    tuple val(sampleID), path("*"), emit: reads

  script:
  """
  #!/usr/bin/env bash
  
  trim_galore -q 30 \
   --length 50 --trim-n -o "." --gzip \
   --cores 2 --paired "$read1" "$read2" \
   --no_report_file --basename "$sampleID"
 
  """
}

/*
workflow output {
 take:
    data
 main:
   baseCounts(data) | collect | mqc(baseCounts.out[1])
 emit:
   mqc.out
 //mqc.seqfu_qc.out.view()
}

*/



workflow {
 readsDir = Channel.fromPath("$params.rawReads/*.f*q.gz", checkIfExists: true)
 baseCounts(readsDir) | collectFile(name: 'seqfu.txt', keepHeader: true, sort: true)
 mqcPrep(readsDir) | collectFile(name: 'seqfu_mqc.tsv') | mqc
 Channel.fromPath(params.index) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sampleID, file(row.read1), file(row.read2)) } \
        | trimReads 
}









