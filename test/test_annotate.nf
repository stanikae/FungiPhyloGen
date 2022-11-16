#!/usr/bin/env nextflow
nextflow.enable.dsl=2



process RUNSNPEFF {
 // tag "$sampleId"

  conda "$params.cacheDir/fpgCallVariants"
  publishDir "$params.ann", mode: 'copy'


  input:
    //val ready
    path(vcf)
    file(ref)
    file(gbk)

  output:
    path("snpeff/snpeff_ann.vcf"), emit: ann
    //path("*.iqtree"), emit: iq_log
    //path("*.mldist"), emit: mldist
    //path("*.bionj"), emit: bionj
    //path("*.log"), emit: log
    //path("*.gz"),emit: zp


  script:


  """
  #!/usr/bin/env bash
  
  if ! [[ -d snpeff ]]; then mkdir -p snpeff; fi
  
  
  CONDA_BASE="$params.cacheDir/fpgCallVariants"
  pfx=`basename $projectDir`
  configPath=\$(find \${CONDA_BASE}/share/* -maxdepth 1 -name "snpEff.config")
  binDir=\$(dirname `find \${CONDA_BASE}/share/* -name "snpEff.jar"`)

  # add ref and gbk to workdir
  if ! [[ -d snpeff/ref/genomes/ref ]]; then mkdir -p snpeff/ref/genomes/ref; fi

  cp "$ref" snpeff/ref/genomes/ref/ref.fa
  cp "$gbk" snpeff/ref/genomes/ref/genes.gbk
  cp \$configPath snpeff/ref/snpEff.config
  cp "$vcf" snpeff/ref/fpg.vcf
  
  # append organism description to config
  echo -e "ref.genome : FPG Reference" >> snpeff/ref/snpEff.config 
  
  # build database
  cd ./snpeff/ref
  
  if [[ -f genomes/ref/genes.gbk ]]; then
     # build database
     java -jar \$binDir/snpEff.jar build -genbank -dataDir genomes -c snpEff.config ref
     # Annotate SNPs with snpEff
     java -jar \$binDir/snpEff.jar ann ref -noLog -nodownload -onlyProtein -dataDir genomes -c snpEff.config fpg.vcf > ../snpeff_ann.vcf
  fi

  

 # if ! [[ -f ../snpeff_ann.vcf ]]; then 
 #    cp $vcf ../snpeff_ann.vcf
 # fi
    
  
  """

}


/*
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
*/

workflow WFANNOTATESNP {
  take:
    x
    y
    z

  main:
     RUNSNPEFF(x,y,z)
     //RUNSNPDISTS(aln)

  emit:
   eff = RUNSNPEFF.out.ann

  // iq_tre = RUNIQTREE.out.tre
  // iq_log = RUNIQTREE.out.iq_log
  // iq_mldist = RUNIQTREE.out.mldist
  // iq_bionj = RUNIQTREE.out.bionj
  // log_file = RUNIQTREE.out.log
  // snp_dist = RUNSNPDISTS.out.dist
        

}



workflow {
 WFIQTREE(file("$params.resultsDir/results/variants/bcftools/fpg.filt.norm.pass.vcf"),file("$params.refseq"),file("$params.gbk"))

}



