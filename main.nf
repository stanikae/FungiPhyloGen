#!/usr/bin/env nextflow

nextflow.enable.dsl=2


// Set parameters
//params.index = "$params.inputDir/sample_file.csv"
params.rawReads = Channel.fromPath("$params.readsDir/**.f*q.gz", checkIfExists: true)
params.trm = "$params.resultsDir/clean_reads"
//params.prjName = "Test01"

// channels
reads_ch = Channel.fromPath("$params.readsDir/**.f*q.gz", checkIfExists: true)
index_ch = Channel.fromPath(params.index, checkIfExists: true)


// create directories - optional
params.fqcCln = file("$params.resultsDir/fqc_clean")
params.fqcCln.mkdirs() ? ! params.fqcCln.exists() : 'Directory already exists' 
params.bwaMem = file("$params.resultsDir/align")
params.bcftl = file("$params.resultsDir/variants")
params.iq = file("$params.resultsDir/iqtree")
params.dist = file("$params.resultsDir/snpdists")
params.ann = file("$params.resultsDir/snpeff")
params.bcftl.mkdirs() ? ! params.bcftl.exists() : 'Directory already exists'
params.iq.mkdirs() ? ! params.iq.exists() : 'Directory already exists'


//params.flag = false


// Define files 
contaminants = "$params.cacheDir/trimReads/opt/fastqc*/Configuration/contaminant_list.txt"
adapters = "$params.cacheDir/trimReads/opt/fastqc*/Configuration/adapter_list.txt"


// load modules
include { trimReads as TRIMREADS } from './modules/trimreads.nf'
include { fqc as FASTQCRAW } from './modules/fastqc.nf' addParams(fqcOut: "$params.resultsDir/fqc_raw") 
include { fqc as FASTQCCLEAN } from './modules/fastqc.nf' addParams(fqcOut: "$params.resultsDir/fqc_clean")
include { mqc as MULTIQCRAW } from './modules/multiqc.nf' addParams(mqcOut: "$params.resultsDir/multiqc/raw", fileN: "${params.prjName}_raw")
include { mqc as MULTIQCCLEAN } from './modules/multiqc.nf' addParams(mqcOut: "$params.resultsDir/multiqc/clean", fileN: "${params.prjName}_clean")
include { GETREPEATS } from './modules/indexref.nf'
include { REPEATSBED } from './modules/indexref.nf'
include { MASKREF } from './modules/indexref.nf' addParams(refMasked: "$params.resultsDir/masked")
include { INDEXREF } from './modules/indexref.nf' addParams(refIndex: "$params.resultsDir/index")
include { ALIGNBWAMEM } from './modules/bwaAlign.nf' addParams(bwaMem: "$params.resultsDir/align")
include { MARKDUPS } from './modules/bwaAlign.nf' addParams(bwaMem: "$params.resultsDir/align")
include { SAMINDEX } from './modules/bwaAlign.nf' addParams(bwaMem: "$params.resultsDir/align")
include { CALLVARIANTS } from './modules/variantCalling.nf' addParams(bcftl: "$params.resultsDir/variants") 
//include { INDEXBCF } from './modules/variantCalling.nf' addParams(bcftl: "$params.resultsDir/variants")
include { SOFTFILTERVCF } from './modules/variantCalling.nf' addParams(bcftl: "$params.resultsDir/variants")
include { BCFNORM } from './modules/variantCalling.nf' addParams(bcftl: "$params.resultsDir/variants")
include { FILTERVCF } from './modules/variantCalling.nf' addParams(bcftl: "$params.resultsDir/variants")
include { VCFSNPS2FASTA } from './modules/variantCalling.nf' addParams(bcftl: "$params.resultsDir/variants")
include { VCF2PHYLIP } from './modules/variantCalling.nf' addParams(bcftl: "$params.resultsDir/variants")
include { RUNIQTREE } from './modules/phyloTrees.nf' addParams(iq: "$params.resultsDir/iqtree")
include { RUNSNPDISTS } from './modules/phyloTrees.nf' addParams(dist: "$params.resultsDir/snpdists")
include { RUNSNPEFF } from './modules/annotate.nf' addParams(dist: "$params.resultsDir/snpeff")



// load workflows
include { trim } from './modules/trimreads.nf'
include { ALN } from './modules/bwaAlign.nf'
include { BCFTOOLS } from './modules/variantCalling.nf'
include { WFIQTREE } from './modules/phyloTrees.nf'
include { WFANNOTATESNP } from './modules/annotate.nf'



// Run analysis
workflow {

  //FASTQCRAW(file("$contaminants"),file("$adapters"),reads_ch)
  fq_ch = index_ch.splitCsv(header:true).map { row-> tuple(row.sampleID, file(row.read1), file(row.read2)) }
  FASTQCRAW(file("$contaminants"),file("$adapters"), fq_ch)   				// fastqc on raw reads
  trim(index_ch)									// trimming raw reads
  //clean_ch = Channel.fromPath("$params.readsDir/**.f*q.gz", checkIfExists: true)
  //FASTQCCLEAN(file("$contaminants"),file("$adapters"),trim.out.rds.collect())			// fastqc on clean reads
  //MULTIQCRAW(FASTQCRAW.out.collect().unique())										// multiqc on raw reads
  //MULTIQCCLEAN(FASTQCCLEAN.out.unique().collect())									// multiqc on clean reads
  
  // prepare and index ref file
  GETREPEATS(file("$params.refseq"))
  REPEATSBED(GETREPEATS.out.delta)
  MASKREF(file("$params.refseq"),REPEATSBED.out.rpts_bed)
  INDEXREF(MASKREF.out.masked_fa)
  ALN(INDEXREF.out.prs,trim.out.rds)						// perform ref based mapping
  MARKDUPS(ALN.out.bam)
  SAMINDEX(MARKDUPS.out.marked)
  BCFTOOLS(file("$params.refseq"),MARKDUPS.out.marked.collect().unique())
  WFIQTREE(BCFTOOLS.out.msa_snp) 
  WFANNOTATESNP(BCFTOOLS.out.pass_vcf,file("$params.refseq"),file("$params.gbk"))
 
  //MARKDUPS(ALN.out.bam
  //            .map { file -> def key = file.name.replaceAll(".bam","")
  //                     return tuple(key, file) } )

  //MARKDUPS(ALN.out.bam
  //            .map { file -> def key = file.name.toString().tokenize('.').get(0).replace('[','')
  //           return tuple(key, file)}
  //       )
  

}


