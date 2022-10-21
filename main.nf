#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.rawReads = "$params.inputDir/raw"
params.index = "$params.inputDir/sample_file.csv"
params.raw = Channel.fromPath("$params.rawReads/*.f*q.gz", checkIfExists: true)
//params.cleanReads = Channel.fromPath("$params.resultsDir/clean_reads/*val*.f*q.gz", checkIfExists: true)
params.trm = "$params.resultsDir/clean_reads"


// fqc and mqc parameters
reads_ch = Channel.fromPath("$params.rawReads/*.f*q.gz", checkIfExists: true)
params.mqcOut = "$params.resultsDir/multiqc/raw"  //$params.mqcTyp"
params.fqcOut = "$params.resultsDir/fqc_raw"
params.fqcCln = file("$params.resultsDir/fqc_clean")


params.fqcCln.mkdirs() ? ! params.fqcCln.exists() : 'Directory already exists' 


// files 
contaminants = "$params.cacheDir/trimReads/opt/fastqc*/Configuration/contaminant_list.txt"
adapters = "$params.cacheDir/trimReads/opt/fastqc*/Configuration/adapter_list.txt"


//params.x1 = file('params.resultsDir/fqc_clean')
//params.x1.mkdirs()


// load modules
include { trimReads } from './modules/fpg_trimGalore.nf'
include { fqc; mqc; fqc as cln } from './modules/fastqc.nf' 
//include { fqc as raw } from './modules/fastqc.nf' //addParams(fqcOut: "$projectDir/fqc_out") 
//include { fqc as clean } from './modules/fpg_trimGalore.nf' //addParams(fqcOut: "$projectDir/fqc_clean_out")
//include { mqc } from './modules/fastqc.nf' //addParams(mqcOut: "$params.resultsDir/raw_report")
//include { mqc as mqc_clean } from './modules/fpg_trimGalore.nf' //addParams(mqcOut: "$params.resultsDir/clean_report")



// load workflows
include { qc_ind } from './modules/fastqc.nf'
include { qc_ind as tst } from './modules/fastqc.nf' //addParams(fqcOut: "$params.resultsDir/fqc_clean")
include { trim } from './modules/fpg_trimGalore.nf'




//workflow {
// (params.flag ? bar : foo) | omega
//}



workflow {
  qc_ind(file("$contaminants"),file("$adapters"),reads_ch,"$params.fqcOut")
  trim()
  
  if ( Channel.fromPath("$params.resultsDir/clean_reads/*val*.f*q.gz", checkIfExists: true) ) {
        clean_ch = Channel.fromPath("$params.resultsDir/clean_reads/*val*.f*q.gz", checkIfExists: true)
	tst(file("$contaminants"),file("$adapters"),clean_ch,"$params.fqcOut")
  }
  //qc_clean(file("$contaminants"),file("$adapters"),"$params.cleanReads","$params.fqcOut")

}


