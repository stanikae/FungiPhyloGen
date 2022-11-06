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
params.bwaMem = file("$params.resultsDir/align")
params.bcftl = file("$params.resultsDir/variants")
params.bcftl.mkdirs() ? ! params.bcftl.exists() : 'Directory already exists'



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
include { GETREPEATS } from './test_indexRef.nf'
include { REPEATSBED } from './test_indexRef.nf'
include { MASKREF } from './test_indexRef.nf' addParams(refMasked: "$params.resultsDir/masked")
include { INDEXREF } from './test_indexRef.nf' addParams(refIndex: "$params.resultsDir/index")
include { ALIGNBWAMEM } from './test_bwaAlign.nf' addParams(bwaMem: "$params.resultsDir/align")
include { MARKDUPS } from './test_bwaAlign.nf' addParams(bwaMem: "$params.resultsDir/align")
include { SAMINDEX } from './test_bwaAlign.nf' addParams(bwaMem: "$params.resultsDir/align")
include { CALLVARIANTS } from './test_variantCalling.nf' addParams(bcftl: "$params.resultsDir/variants") 
//include { INDEXBCF } from './test_variantCalling.nf' addParams(bcftl: "$params.resultsDir/variants")
include { SOFTFILTERVCF } from './test_variantCalling.nf' addParams(bcftl: "$params.resultsDir/variants")
include { BCFNORM } from './test_variantCalling.nf' addParams(bcftl: "$params.resultsDir/variants")
include { FILTERVCF } from './test_variantCalling.nf' addParams(bcftl: "$params.resultsDir/variants")
include { VCFSNPS2FASTA } from './test_variantCalling.nf' addParams(bcftl: "$params.resultsDir/variants")
include { VCF2PHYLIP } from './test_variantCalling.nf' addParams(bcftl: "$params.resultsDir/variants")





// load workflows
include { trim } from './test_trimGalore.nf'
include { ALN } from './test_bwaAlign.nf'
include { BCFTOOLS } from './test_variantCalling.nf'

// Run analysis
workflow {

  frw(file("$contaminants"),file("$adapters"),reads_ch)   				// fastqc on raw reads
  trim(index_ch)									// trimming raw reads
  //clean_ch = Channel.fromPath("$params.readsDir/**.f*q.gz", checkIfExists: true)
  fcln(file("$contaminants"),file("$adapters"),trim.out.rds.collect())			// fastqc on clean reads
  mraw(frw.out.collect())										// multiqc on raw reads
  mclean(fcln.out.collect())									// multiqc on clean reads
  
  // prepare and index ref file
  GETREPEATS(file("$params.refseq"))
  REPEATSBED(GETREPEATS.out.delta)
  MASKREF(file("$params.refseq"),REPEATSBED.out.rpts_bed)
  INDEXREF(MASKREF.out.masked_fa)
  ALN(INDEXREF.out.prs,trim.out.rds)						// perform ref based mapping
  MARKDUPS(ALN.out.bam)
  BCFTOOLS(file("$params.refseq"),MARKDUPS.out.marked.collect()) 
  //MARKDUPS(ALN.out.bam
  //            .map { file -> def key = file.name.replaceAll(".bam","")
  //                     return tuple(key, file) } )

  //MARKDUPS(ALN.out.bam
  //            .map { file -> def key = file.name.toString().tokenize('.').get(0).replace('[','')
  //           return tuple(key, file)}
  //       )
  
  SAMINDEX(MARKDUPS.out.marked)

}


