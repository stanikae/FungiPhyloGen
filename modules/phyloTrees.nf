#!/usr/bin/env nextflow
nextflow.enable.dsl=2



process RUNIQTREE {
 // tag "$sampleId"

  label 'process_high'

  conda "$params.cacheDir/fpgPhylogen"
  publishDir "$params.iq", mode: 'copy'


  input:
    //val ready
    path(msa)

  output:
    path("*.contree"), emit: con
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

  iqtree -s $msa --prefix \$pfx -T AUTO -m GTR+F+ASC+R4 -o "reference" -B 1000 --threads-max ${task.cpus} 
 
  """

}





process RUNRAPIDNJ {
   label 'process_medium'
   
   conda "$params.cacheDir/fpgPhylogen"
   publishDir "$params.nj", mode: 'copy'

   input:
     path(msa)


   output:
     path("*.dist"), emit: dist
     path("*.nwk"), emit: nwk

   script:

   """
   #!/usr/bin/env bash
   
   if ! [[ -d rapidnj ]]; then mkdir rapidnj; fi

  
   rapidnj -i fa -c ${task.cpus} -b 1000 -o m -x rapid.nj.1000bootstraps.dist $msa

   rapidnj -i fa -c ${task.cpus} -b 1000 -o t -x rapid.nj.1000bootstraps.nwk $msa


   """


}


process RUNSNPDISTS {
  label 'process_high'

  conda "$params.cacheDir/fpgPhylogen"
  publishDir "$params.dist", mode: 'copy'

  input:
    path msa
    path vmsa

  output:
    path("snpdists/*snpdist.csv"), emit: dist
    path("snpdists/*.snpdistmolten.csv"), emit: molt
    path("snpdists/*vcf2phylo.snpdist.csv"), emit: vdist
    path("snpdists/*.vcf2phylo.snpdistmolten.csv"), emit: vmolt

  script:
  """
  #!/usr/bin/env bash
  set -euo pipefail

  mkdir -p snpdists

  pfx=\$(basename "$projectDir")

  snp-dists \
    -j ${task.cpus} \
    -c \
    "$msa" \
    > snpdists/\${pfx}.snpdist.csv

  snp-dists \
    -j ${task.cpus} \
    -a \
    -m \
    -c \
    "$msa" \
    > snpdists/\${pfx}.snpdistmolten.csv

  snp-dists \
    -j ${task.cpus} \
    -c \
    "$vmsa" \
    > snpdists/\${pfx}.vcf2phylo.snpdist.csv

  snp-dists \
    -j ${task.cpus} \
    -a \
    -m \
    -c \
    "$vmsa" \
    > snpdists/\${pfx}.vcf2phylo.snpdistmolten.csv

  for file in \
    snpdists/\${pfx}.snpdist.csv \
    snpdists/\${pfx}.snpdistmolten.csv \
    snpdists/\${pfx}.vcf2phylo.snpdist.csv \
    snpdists/\${pfx}.vcf2phylo.snpdistmolten.csv
  do
    if [[ -s "\$file" ]]; then
      sed -i.bak 's/_sorted_marked\\.bam//g' "\$file"
      rm -f "\${file}.bak"
    fi
  done
  """
}




workflow WFIQTREE {
  take:
    aln
    valn

  main:
     RUNIQTREE(aln)
     RUNSNPDISTS(aln,valn)

  emit:
   con_tre = RUNIQTREE.out.con
   iq_tre = RUNIQTREE.out.tre
   iq_log = RUNIQTREE.out.iq_log
   iq_mldist = RUNIQTREE.out.mldist
   iq_bionj = RUNIQTREE.out.bionj
   log_file = RUNIQTREE.out.log
   snp_dist = RUNSNPDISTS.out.dist
   snp_molt = RUNSNPDISTS.out.molt
   snp_vdist = RUNSNPDISTS.out.vdist
   snp_vmolt = RUNSNPDISTS.out.vmolt
        

}



/*
workflow {
 WFIQTREE(file("$params.resultsDir/results/variants/SNPfasta/fpg_snp_aln.fa"))

}
*/

