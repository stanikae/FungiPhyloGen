#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//params.rawReads = "$params.inputDir/raw"
params.index = "$params.inputDir/sample_file.csv"
params.rawReads = Channel.fromPath("$params.readsDir/**.f*q.gz", checkIfExists: true)
//params.cleanReads = Channel.fromPath("$params.resultsDir/clean_reads/*val*.f*q.gz", checkIfExists: true)
params.trm = "$params.resultsDir/clean_reads"



reads_ch = Channel.fromPath("$params.readsDir/**.f*q.gz", checkIfExists: true)



// fqc and mqc parameters
//reads_ch = Channel.fromPath("$params.rawReads/*.f*q.gz", checkIfExists: true)

//params.mqcOut = "$params.resultsDir/multiqc/raw"  //$params.mqcTyp"
//params.fqcOut = "$params.resultsDir/fqc_raw"
params.fqcCln = file("$params.resultsDir/fqc_clean")

params.fqcCln.mkdirs() ? ! params.fqcCln.exists() : 'Directory already exists' 


params.flag = false


// files 
contaminants = "$params.cacheDir/trimReads/opt/fastqc*/Configuration/contaminant_list.txt"
adapters = "$params.cacheDir/trimReads/opt/fastqc*/Configuration/adapter_list.txt"


//params.x1 = file('params.resultsDir/fqc_clean')
//params.x1.mkdirs()


// load modules
//include { trimReads } from './modules/fpg_trimGalore.nf'
//include { fqc; mqc; fqc as cln } from './modules/fastqc.nf' 

include { trimReads } from './test_trimGalore.nf'
include { fqc as frw } from './test_fastqc.nf' addParams(fqcOut: "$params.resultsDir/fqc_raw") 
include { fqc as fcln } from './test_fastqc.nf' addParams(fqcOut: "$params.resultsDir/fqc_clean")
include { mqc as mraw} from './multiqc.nf' //addParams(mqcOut: "$params.resultsDir/multiqc/raw")
include { mqc as mclean } from './multiqc.nf' addParams(mqcOut: "$params.resultsDir/multiqc/clean")



// load workflows
//include { qc_ind } from './test_fastqc.nf'
//include { qc_ind as tst } from './modules/fastqc.nf' //addParams(fqcOut: "$params.resultsDir/fqc_clean")

include { trim } from './test_trimGalore.nf'




//workflow {
// (params.flag ? bar : foo) | omega
//}



workflow {

  frw(file("$contaminants"),file("$adapters"),reads_ch)
  trim(Channel.fromPath(params.index))

  //clean_ch = Channel.fromPath("$params.readsDir/**.f*q.gz", checkIfExists: true)
  fcln(file("$contaminants"),file("$adapters"),trim.out.rds.collect())
  
  mraw(frw.out)
  mclean(fcln.out)

//  qc_ind(file("$contaminants"),file("$adapters"),reads_ch)
//  //qc_ind(file("$contaminants"),file("$adapters"),reads_ch,"$params.fqcOut")
//  trim(Channel.fromPath(params.index))
//  //tst(file("$contaminants"),file("$adapters"),clean_ch,trim.out.evl)
//  //tst(file("$contaminants"),file("$adapters"),trim.out.evl) 
//  if ( Channel.fromPath("$params.resultsDir/clean_reads/*val*.f*q.gz", checkIfExists: true) ) {
//        clean_ch = Channel.fromPath("$params.resultsDir/clean_reads/*val*.f*q.gz", checkIfExists: true)
//   //       tst(file("$contaminants"),file("$adapters"),trim.out)
//  	tst(file("$contaminants"),file("$adapters"),clean_ch)
//  }
//  //qc_clean(file("$contaminants"),file("$adapters"),"$params.cleanReads","$params.fqcOut")

}


