#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//params.bwaMem = "$params.resultsDir/align"


process CALLVARIANTS {
 // tag "$sampleId"
  
  conda "$params.cacheDir/fpgCallVariants"
  publishDir "$params.bcftl", mode: 'copy'


  input:
    //val ready
    //tuple val(sampleId), file(bam)
    path(ref)
    path(bam)

  output:
    //tuple val(sampleId), path("bcftools/*.bcf") , emit: vcfs
    path("bcftools/fpg.call.bcf"), emit: vcfs
    path("bcftools/fpg.call.bcf.csi"), emit: idx


  script:
  
  """
  #!/usr/bin/env bash
  
  if ! [[ -d bcftools ]]; then mkdir bcftools; fi


  bcftools mpileup -a AD -Q 30 -f $ref $bam -Ou \\
  | bcftools call --threads "$params.threads" --ploidy "$params.ploidy" -mv -Ob -o bcftools/fpg.call.bcf

  bcftools index bcftools/fpg.call.bcf

  """
}


/*
process INDEXBCF {
  
  conda "$params.cacheDir/fpgCallVariants"
  publishDir "$params.bcftl", mode: 'copy'

  input:
     path(bcf)

  output:
    path("bcftools/*.csi"), emit: bcf_idx


  script:
  
  """
  #!/usr/bin/env bash
  
  bcftools index $bcf

  """ 


}
*/


process SOFTFILTERVCF {
   
   conda "$params.cacheDir/fpgCallVariants"
   publishDir "$params.bcftl", mode: 'copy'


   input:
     path(bcf)
     path(idx)


   output:
     path("bcftools/fpg_filt.bcf"), emit: bcf_filt
     path("bcftools/fpg_filt.bcf.csi"), emit: idx_filt


   script:

    """
    #!/usr/bin/env bash

    #bcftools filter -s 'LowQual' -i  'QUAL>=30 && AD[*:1]>=25' -g8 -G10 $bcf -o bcftools/fpg_filt.bcf
    bcftools filter -s 'LowQual' -i  'QUAL>=30 && AD[*:1]>=25 && MQ>=30 && DP>=10' -g8 -G10 $bcf -o bcftools/fpg_filt.bcf

    bcftools index bcftools/fpg_filt.bcf

    """
}



process BCFNORM {

   conda "$params.cacheDir/fpgCallVariants"
   publishDir "$params.bcftl", mode: 'copy'

   input:
     path(bcf)
     path(idx)
     path(ref)


   output:
     path("bacftools/fpg.filt.norm.bcf"), emit: bcf_norm
     path("bacftools/fpg.filt.norm.bcf.csi"), emit: idx_norm

   script:

   """
   #!/usr/bin/env bash

   bcftools norm -f $ref $bcf -o bacftools/fpg.filt.norm.bcf
   
   bcftools index bacftools/fpg.filt.norm.bcf
   """


}




process FILTERVCF {

   conda "$params.cacheDir/fpgCallVariants"
   publishDir "$params.bcftl", mode: 'copy'


   input:
     path(bcf)
     path(idx)


   output:
     path("bcftools/fpg_filt.norm.pass.bcf"), emit: bcf_pass
     path("bcftools/fpg_filt.norm.pass.bcf.csi"), emit: idx_pass


   script:

    """
    #!/usr/bin/env bash
    
    bcftools view -i 'FILTER="PASS"' $bcf -o bcftools/fpg_filt.norm.pass.bcf
    
    bcftools index bcftools/fpg_filt.norm.pass.bcf

    """
}



process VCFSNPS2FASTA {

   conda "$params.cacheDir/fpgVcf2FastaEnv"
   publishDir "$params.bcftl", mode: 'copy'


   input:
     path(bcf)
     path(idx)


   output:
     path("SNPfasta/fpg_snp_aln.fa"), emit: snp_msa
     path("SNPfasta/fpg_filt.norm.pass.vcf"), emit: vcf_filt
     path("SNPfasta/vcf_infile"), emit: infile


   script:

    """
    #!/usr/bin/env bash
    
    if ! [[ -d SNPfasta ]]; then mkdir SNPfasta; fi   

    bcftools view -i 'FILTER="PASS"' -Ov -o SNPfasta/fpg_filt.norm.pass.vcf \$bcf

    readlink -f SNPfasta/fpg_filt.norm.pass.vcf > SNPfasta/vcf_infile

    python ./templates/broad-fungalgroup/scripts/SNPs/vcfSnpsToFasta.py --max_amb_samples 1 SNPfasta/vcf_infile > SNPfasta/fpg_snp_aln.fa
    
    """

}





process VCF2PHYLIP {

   conda "$params.cacheDir/fpgCallVariants"
   publishDir "$params.bcftl", mode: 'copy'


   input:
     path(bcf)
     path(idx)


   output:
     path("vcf2phylip/vcfSNPs.min2.fasta"), emit: snp_aln
     path("vcf2phylip/vcfSNPs.min2.phy"), emit: phy
     path("vcf2phylip/vcfSNPs.min2.used_sites.tsv"), emit: sites
     path("vcf2phylip/fpg_filt.norm.pass.vcf"), emit: vcf_filt


   script:

    """
    #!/usr/bin/env bash

    if ! [[ -d vcf2phylip ]]; then mkdir vcf2phylip; fi

    bcftools view  -i 'FILTER="PASS"' -Ov -o vcf2phylip/fpg_filt.norm.pass.vcf \$bcf

    python ./templates/vcf2phylip.py -i vcf2phylip/fpg_filt.norm.pass.vcf -f -w --output-folder vcf2phylip --output-prefix vcfSNPs

    """
}



workflow BCFTOOLS {
  take:
    fa
    path(bam)

  main:
     CALLVARIANTS(fa,bam)
     INDEXBCF(CALLVARIANTS.out.vcfs)
     SOFTFILTERVCF(CALLVARIANTS.out.vcfs,INDEXBCF.out.bcf_idx)
     BCFNORM(SOFTFILTERVCF.out.bcf_filt,SOFTFILTERVCF.out.idx_filt,fa) 
     FILTERVCF(BCFNORM.out.bcf_norm,BCFNORM.out.idx_norm)
     VCFSNPS2FASTA(FILTERVCF.out.bcf_pass,FILTERVCF.out.idx_pass)
     VCF2PHYLIP(FILTERVCF.out.bcf_pass,FILTERVCF.out.idx_pass)

  emit:
    bcf_raw = CALLVARIANTS.out.vcfs
    bcf_rawIdx = INDEXBCF.out.bcf_idx
    filt_bcf = SOFTFILTERVCF.out.bcf_filt
    filt_idx = SOFTFILTERVCF.out.idx_filt
    norm_bcf = BCFNORM.out.bcf_norm
    norm_idx = BCFNORM.out.idx_norm
    pass_bcf = FILTERVCF.out.bcf_pass
    pass_idx = FILTERVCF.out.idx_pass
    msa_snp = VCFSNPS2FASTA.out.snp_msa
    vcf_z = VCFSNPS2FASTA.out.vcf_filt
    vcf_list = VCFSNPS2FASTA.out.infile
    aln_snp = VCF2PHYLIP.out.snp_aln
    aln_phy = VCF2PHYLIP.out.phy
    snp_sites = VCF2PHYLIP.out.sites
    z_vcf = VCF2PHYLIP.out.vcf_filt
        

}



workflow {

 BCFTOOLS()

}



