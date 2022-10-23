#!/usr/bin/env nextflow

nextflow.enable.dsl=2


// Set parameters
params.index = "$params.inputDir/sample_file.csv"
params.rawReads = Channel.fromPath("$params.readsDir/**.f*q.gz", checkIfExists: true)
params.trm = "$params.resultsDir/clean_reads"
params.prjName = "Test01"

// channels
reads_ch = Channel.fromPath("$params.readsDir/**.f*q.gz", checkIfExists: true)
index_ch = Channel.fromPath(params.index)


// create directories - optional
params.fqcCln = file("$params.resultsDir/fqc_clean")
params.fqcCln.mkdirs() ? ! params.fqcCln.exists() : 'Directory already exists' 


//params.flag = false


// Define files 
contaminants = "$params.cacheDir/trimReads/opt/fastqc*/Configuration/contaminant_list.txt"
adapters = "$params.cacheDir/trimReads/opt/fastqc*/Configuration/adapter_list.txt"


// load modules
include { trimReads } from './test_trimGalore.nf'
include { fqc as frw } from './test_fastqc.nf' addParams(fqcOut: "$params.resultsDir/fqc_raw") 
include { fqc as fcln } from './test_fastqc.nf' addParams(fqcOut: "$params.resultsDir/fqc_clean")
include { mqc as mraw} from './multiqc.nf' addParams(mqcOut: "$params.resultsDir/multiqc/raw", fileN: "${params.prjName}_raw")
include { mqc as mclean } from './multiqc.nf' addParams(mqcOut: "$params.resultsDir/multiqc/clean", fileN: "${params.prjName}_clean")



// load workflows
include { trim } from './test_trimGalore.nf'


// Run analysis
workflow {

  frw(file("$contaminants"),file("$adapters"),reads_ch)   				// fastqc on raw reads
  trim(index_ch)									// trimming raw reads
  //clean_ch = Channel.fromPath("$params.readsDir/**.f*q.gz", checkIfExists: true)
  fcln(file("$contaminants"),file("$adapters"),trim.out.rds.collect())			// fastqc on clean reads
  mraw(frw.out.collect())										// multiqc on raw reads
  mclean(fcln.out.collect())									// multiqc on clean reads


}


