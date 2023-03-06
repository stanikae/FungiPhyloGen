#!/usr/bin/env nextflow
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

include { FUNGIPHYLOGEN as FungiPhyloGen } from './workflows/fpg'

//
// WORKFLOW: Run main nf-FungiPhyloGen analysis pipeline
//
workflow NF_FUNGIPHYLOGEN {
    FungiPhyloGen ()
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
    NF_FUNGIPHYLOGEN ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
