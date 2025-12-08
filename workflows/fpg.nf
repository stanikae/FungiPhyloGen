#!/usr/bin/env nextflow
nextflow.enable.dsl=2



// --- WORKFLOW: FUNGIPHYLOGEN ---


/*
================================================================
        Load Modules
================================================================
*/

include { fqc as FASTQCRAW } from '../modules/fastqc.nf' addParams(fqcOut: "$params.resultsDir/fqc_raw")
include { trimReads as TRIMREADS } from '../modules/trimreads.nf'
include { fqc as FASTQCCLEAN } from '../modules/fastqc.nf' addParams(fqcOut: "$params.resultsDir/fqc_clean")
include { mqc as MULTIQCRAW } from '../modules/multiqc.nf' addParams(mqcOut: "$params.resultsDir/multiqc/raw", fileN: "${params.prjName}_raw")
include { mqc as MULTIQCCLEAN } from '../modules/multiqc.nf' addParams(mqcOut: "$params.resultsDir/multiqc/clean", fileN: "${params.prjName}_clean")
include { GETREPEATS } from '../modules/indexref.nf'
include { REPEATSBED } from '../modules/indexref.nf'
include { MASKREF } from '../modules/indexref.nf' addParams(refMasked: "$params.resultsDir/masked")
include { INDEXREF } from '../modules/indexref.nf' addParams(refIndex: "$params.resultsDir/index")
include { SPADES } from '../modules/denovo.nf' addParams(deNovo: "$params.resultsDir/assemblies")
include { ALIGNBWAMEM } from '../modules/bwaAlign.nf' addParams(bwaMem: "$params.resultsDir/align")
include { MARKDUPS } from '../modules/bwaAlign.nf' addParams(bwaMem: "$params.resultsDir/align")
include { SORTMARKED } from '../modules/bwaAlign.nf' addParams(bwaMem: "$params.resultsDir/align")
include { SAMINDEX } from '../modules/bwaAlign.nf' addParams(bwaMem: "$params.resultsDir/align")
include { CALLVARIANTS } from '../modules/callVariants.nf' addParams(bcftl: "$params.resultsDir/variants")
include { FILTERSAMPLE } from '../modules/callVariants.nf' addParams(bcftl: "$params.resultsDir/variants")
//include { CALLVARIANTSgrp } from '../modules/callVariants.nf' addParams(bcftl: "$params.resultsDir/variants")
//include { INDEXBCF } from '../modules/callVariants.nf' addParams(bcftl: "$params.resultsDir/variants")
include { REHEADERVCF } from '../modules/callVariants.nf' addParams(bcftl: "$params.resultsDir/variants")
include { SOFTFILTERVCF } from '../modules/callVariants.nf' addParams(bcftl: "$params.resultsDir/variants")
include { BCFNORM } from '../modules/callVariants.nf' addParams(bcftl: "$params.resultsDir/variants")
include { FILTERVCF } from '../modules/callVariants.nf' addParams(bcftl: "$params.resultsDir/variants")
include { VCFSNPS2FASTA } from '../modules/callVariants.nf' addParams(bcftl: "$params.resultsDir/variants")
include { VCF2PHYLIP } from '../modules/callVariants.nf' addParams(bcftl: "$params.resultsDir/variants")
include { RUNIQTREE } from '../modules/phyloTrees.nf' addParams(iq: "$params.resultsDir/iqtree")
include { RUNRAPIDNJ } from '../modules/phyloTrees.nf' addParams(nj: "$params.resultsDir/rapidnj")
include { RUNSNPDISTS } from '../modules/phyloTrees.nf' addParams(dist: "$params.resultsDir/snpdists")
include { RUNSNPEFF } from '../modules/annotate.nf' addParams(ann: "$params.resultsDir/snpeff")


/*
========================================================================================
    Import Subworkflows
========================================================================================
*/


include { trim } from '../modules/trimreads.nf'
// include { ALN } from '../modules/bwaAlign.nf'
// include { TALN } from '../modules/bwaAlign.nf'
include { BCFTOOLS } from '../modules/callVariants.nf'
include { WFIQTREE } from '../modules/phyloTrees.nf'
include { WFANNOTATESNP } from '../modules/annotate.nf'
include { DENOVO } from '../modules/denovo.nf'



/*
=================================================================
        Create directories - OPTIONAL
=================================================================
*/

params.fqcCln = file("$params.resultsDir/fqc_clean")
params.fqcCln.mkdirs() ? ! params.fqcCln.exists() : 'Directory already exists'
params.bwaMem = file("$params.resultsDir/align")
params.bcftl = file("$params.resultsDir/variants")
params.iq = file("$params.resultsDir/iqtree")
params.nj = file("$params.resultsDir/rapidnj")
params.dist = file("$params.resultsDir/snpdists")
params.ann = file("$params.resultsDir/snpeff")
params.deNovo = file("$params.resultsDir/assemblies")


/*
=================================================================
        Define files
=================================================================
*/




/*
=======================================================================================
    RUN MAIN WORKFLOW
=======================================================================================
*/


workflow FUNGIPHYLOGEN {

// --- INPUT CHANNEL FROM SAMPLESHEET ---
    if (!params.samplesheet) {
        exit 1, "Input samplesheet not specified. Please provide a path using --samplesheet <file.csv>"
    }

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            // DEFENSIVE CHECK 1: Use the correct CSV headers here!
            if (!row.sampleID || !row.read1 || !row.read2) {
                log.warn "Skipping invalid row (missing columns): $row"
                return null
            }

            // DEFENSIVE CHECK 2: Check for empty strings
            if (row.read1 == "" || row.read2 == "") {
                log.warn "Skipping row with empty paths: $row"
                return null
            }

            def meta = [id: row.sampleID] // MATCHES CSV HEADER 'sampleID'
            
            def reads = [
                // MATCHES CSV HEADERS 'read1' and 'read2'
                file(row.read1, checkIfExists: true),
                file(row.read2, checkIfExists: true)
            ]
            return [meta, reads]
        }
        .filter { it != null } 
        .set { ch_input_reads }
    
    def contaminants = file("/spaces/stanford/anaconda3/envs/fpgtrimReads/opt/fastqc-0.11.9/Configuration/contaminant_list.txt",checkIfExists: true)
    def adapters     = file("/spaces/stanford/anaconda3/envs/fpgtrimReads/opt/fastqc-0.11.9/Configuration/adapter_list.txt",checkIfExists: true)

    // OPTIONAL: DEBUG PRINT
    log.info "Contaminants found: ${contaminants}"
    log.info "Adapters:     $adapters"


    // --- REFERENCE PREPARATION ---
    GETREPEATS(file("$params.refseq"))
    REPEATSBED(GETREPEATS.out.delta)
    MASKREF(file("$params.refseq"),REPEATSBED.out.rpts_bed)
    ch_ref_masked = MASKREF.out.masked_fa
    INDEXREF(ch_ref_masked)
    ch_ref_indexed = INDEXREF.out.prs


    // --- CONDITIONAL READ TRIMMING ---
    if (params.skip_trimming) {
        log.info "Skipping read trimming as requested by --skip_trimming."
        ch_reads_for_alignment = ch_input_reads
    } else {
        log.info "Trimming raw reads."

        // --- MULTIQC ON RAW READS ---
        // FIX: Removed file("$...") wrapper. Passed variables directly.
        FASTQCRAW(contaminants, adapters, ch_input_reads)
        
        // FIX: Added .zip to collect only the zip files
        MULTIQCRAW(FASTQCRAW.out.zip.collect())

        trim(ch_input_reads)
        ch_reads_for_alignment = trim.out.rds

        // --- MULTIQC ON CLEAN READS ---
        // FIX: Removed file("$...") wrapper here too.
        FASTQCCLEAN(contaminants, adapters, ch_reads_for_alignment)
        
        // FIX: Added .zip here too
        MULTIQCCLEAN(FASTQCCLEAN.out.zip.collect())
    }
    
    // --- ALIGNMENT ---
    //ALN(ch_ref_indexed,ch_reads_for_alignment) //ch_ref_indexed.prs.collect()
    ALIGNBWAMEM(ch_ref_indexed, ch_reads_for_alignment)
    //MARKDUPS(.out.bam)
    MARKDUPS(ALIGNBWAMEM.out.aln_bam)
    SORTMARKED(MARKDUPS.out.marked)
    ch_sorted_bam = SORTMARKED.out.sorted 
    SAMINDEX(ch_sorted_bam)

    // --- VARIANT CALLING, FILTERING & CONVERSION ---
    BCFTOOLS(ch_sorted_bam, ch_ref_masked)

    // --- PHYLOGENY ---
    WFIQTREE(BCFTOOLS.out.msa_snp, BCFTOOLS.out.aln_snp)
    RUNRAPIDNJ(BCFTOOLS.out.msa_snp)

    // --- SNP ANNOTATION ---
    WFANNOTATESNP(BCFTOOLS.out.pass_vcf, file(params.refseq), file(params.gbk))


   // --- DE NOVO ASSEMBLY ---
   //DENOVO(ch_reads_for_alignment)  
  
}


