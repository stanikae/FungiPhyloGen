#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process SPADES {
  label 'process_long' 
 
  //conda ("$params.condaPath/fpgDenovo" ? "$params.cacheDir/fpgDenovo" : null)
  //conda "$params.cacheDir/fpgDenovo"
  conda "$params.condaPath/fpgDenovo"
  publishDir "$params.deNovo", mode: 'copy'

  input:
  tuple val(sampleID), file(read1), file(read2) 

  output:
   //path("corrected/*.fastq*") , emit: cor
   path("scaffolds.fasta") , emit: scaf
   //path("contigs.fasta") , emit: ctg
   path("assembly_graph_with_scaffolds.gfa") , emit: scaf_gfa
   //path("assembly_graph.fastg"), emit: gfa
   //path("contigs.paths"), emit: ctg_path
   path("scaffolds.paths"), emit: scaf_path

  script:
  """
  #!/usr/bin/env bash
  
  dirNam=`echo "$sampleID"`
  if ! [[ -d \$dirNam ]]; then mkdir \$dirNam; fi

  spades.py -o \$dirNam --careful -1 "$read1" -2 "$read2" -t ${task.cpus} --cov-cutoff auto 
  
 
  """
}



workflow DENOVO {
    take:
      cln_reads

    main:

        SPADES (
                cln_reads.
                        map { it -> def key = it.name.toString().tokenize('_').get(0).replace('[','')
                                return tuple(key, file(it[0]), file(it[1])) }
                        //.groupTuple()
                        //.view()
       )
    emit:
     //cor = SPADES.out.cor
     scaf = SPADES.out.scaf
     //ctg = SPADES.out.ctg
     scaf_gfa = SPADES.out.scaf_gfa
     //ctg_path = SPADES.out.ctg_path
     scaf_path = SPADES.out.scaf_path
}


