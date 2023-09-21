#!/usr/bin/env nextflow
nextflow.enable.dsl=2



process RUNSNPEFF {
 // tag "$sampleId"

  cpus 8
  executor 'slurm'

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
     java -jar \$binDir/snpEff.jar ann ref -noLog -nodownload -onlyProtein -dataDir genomes -t ${task.cpus} -c snpEff.config fpg.vcf > ../snpeff_ann.vcf
  fi

  

 # if ! [[ -f ../snpeff_ann.vcf ]]; then 
 #    cp $vcf ../snpeff_ann.vcf
 # fi
    
  
  """

}



workflow WFANNOTATESNP {
  take:
    x
    y
    z

  main:
     RUNSNPEFF(x,y,z)

  emit:
   eff = RUNSNPEFF.out.ann
        

}


