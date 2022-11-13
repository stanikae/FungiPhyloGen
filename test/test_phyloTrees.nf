#!/usr/bin/env nextflow
nextflow.enable.dsl=2



process RUNIQTREE {
 // tag "$sampleId"

  conda "$params.cacheDir/fpgPhylogen"
  publishDir "$params.iq", mode: 'copy'


  input:
    //val ready
    path(msa)

  output:
    path("*.treefile"), emit: tre
    path("*.iqtree"), emit: iq_log
    path("*.mldist"), emit: mldist
    path("*.bionj"), emit: bionj
    path("*.log"), emit: log
    path("*.gz"),emit: zp



  script:


  """
  #!/usr/bin/env bash

  if ! [[ -d iqtree ]]; then mkdir iqtree; fi
  
  #cd "./iqtree"
  
  pfx=`basename $projectDir`

  iqtree -s $msa --prefix \$pfx -T AUTO -m GTR+F+ASC+R4 -o "reference" --threads-max 6
 
  """

}



process RUNSNPDISTS {
 // tag "$sampleId"

  conda "$params.cacheDir/fpgPhylogen"
  publishDir "$params.dist", mode: 'copy'


  input:
    //val ready
    path(msa)

  output:
    path("snpdists/*snpdist.csv"), emit: dist


  script:


  """
  #!/usr/bin/env bash

  if ! [[ -d snpdists ]]; then mkdir -p snpdists; fi

  pfx=`basename $projectDir`

  snp-dists -j 6 -c $msa > snpdists/\${pfx}.snpdist.csv
 
  if [[ -f "snpdists/\${pfx}.snpdist.csv" ]]; then
  	
     sed -i 's/_sorted_marked.bam//g' snpdists/\${pfx}.snpdist.csv
  
  fi

  """

}


workflow WFIQTREE {
  take:
    aln

  main:
     RUNIQTREE(aln)
     RUNSNPDISTS(aln)

  emit:
   iq_tre = RUNIQTREE.out.tre
   iq_log = RUNIQTREE.out.iq_log
   iq_mldist = RUNIQTREE.out.mldist
   iq_bionj = RUNIQTREE.out.bionj
   log_file = RUNIQTREE.out.log
   snp_dist = RUNSNPDISTS.out.dist
        

}



workflow {
 WFIQTREE(file("$params.resultsDir/results/variants/SNPfasta/fpg_snp_aln.fa"))

}



