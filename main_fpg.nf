#!/usr/bin/env nextflow


log.info """
=========================================
DEBUGGING PARAMETERS
=========================================
RefSeq : ${params.refseq}
GBK    : ${params.gbk}
Config : ${workflow.configFiles}
=========================================
"""

/*
========================================================================================
    nf-core/bactmap
========================================================================================
    Github : https://github.com/stanikae/FungiPhyloGen
    Website: https://github.com/stanikae/FungiPhyloGen
    Slack  : 
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

//WorkflowMain.initialise(workflow, params, log)

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
