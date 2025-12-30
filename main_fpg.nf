#!/usr/bin/env nextflow

// ------------------------------------------------------
// 1. VALIDATION CHECK (Fail Fast)
// ------------------------------------------------------

if (!params.active_filter) {
    log.error """
    ================================================================
    CRITICAL ERROR: Invalid Filter Profile
    ================================================================
    The profile '${params.filter_profile}' was not found in the config.
    
    Available Profiles:
    ${params.filters.keySet().join('\n    - ')}
    
    Please check your --filter_profile argument.
    ================================================================
    """
    exit 1
}

log.info """
=========================================
DEBUGGING PARAMETERS
=========================================
RefSeq         : ${params.refseq}
GBK            : ${params.gbk}
Filter Profile : ${params.filter_profile}
Filter Logic   : ${params.active_filter}
Config         : ${workflow.configFiles}
=========================================
"""

/*
========================================================================================
    nf-FPG
========================================================================================
    Github : https://github.com/stanikae/FungiPhyloGen
    Website: https://github.com/stanikae/FungiPhyloGen
    Zenodo : 
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { FUNGIPHYLOGEN as FPG } from './workflows/fpg'

//
// WORKFLOW: Run main nf-FungiPhyloGen analysis pipeline
//
workflow FUNGIPHYLOGEN {
    FPG ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// main.nf adapted/modified from: https://github.com/nf-core/bactmap/blob/master/workflows/bactmap.nf
//
workflow {
    FUNGIPHYLOGEN ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
