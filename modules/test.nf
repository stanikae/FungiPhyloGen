#!/usr/bin/env nextflow
nextflow.enable.dsl=2



Channel
    .fromFilePairs('/home/stan/fungal-genomics/test-data/raw/*_{R1_001,R1_002}*.f*q.gz')
    .view()
