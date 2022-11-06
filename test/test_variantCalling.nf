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
  | bcftools call --threads $params.threads --ploidy "$params.ploidy" -mv -Ob -o bcftools/fpg.call.bcf

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
     path("bcftools/fpg.filt.bcf"), emit: bcf_filt
     path("bcftools/fpg.filt.bcf.csi"), emit: idx_filt


   script:

    """
    #!/usr/bin/env bash

    if ! [[ -d bcftools ]]; then mkdir bcftools; fi
    
    bcftools filter -s 'LowQual' -i  'QUAL>=30 && AD[*:1]>=25 && MQ>=30 && DP>=10' -g8 -G10 $bcf -o bcftools/fpg.filt.bcf

    bcftools index bcftools/fpg.filt.bcf

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
     path("bcftools/fpg.filt.norm.bcf"), emit: bcf_norm
     path("bcftools/fpg.filt.norm.bcf.csi"), emit: idx_norm
     path("bcftools/norm.log"), emit: log_norm

   script:

   """
   #!/usr/bin/env bash
   
   if ! [[ -d bcftools ]]; then mkdir bcftools; fi
   
   bcftools norm -f $ref $bcf -o bcftools/fpg.filt.norm.bcf 2> bcftools/norm.log
   
   bcftools index bcftools/fpg.filt.norm.bcf

   """


}




process FILTERVCF {

   conda "$params.cacheDir/fpgCallVariants"
   publishDir "$params.bcftl", mode: 'copy'


   input:
     path(bcf)
     path(idx)


   output:
     path("bcftools/fpg.filt.norm.pass.bcf"), emit: bcf_pass
     path("bcftools/fpg.filt.norm.pass.bcf.csi"), emit: idx_pass
     path("bcftools/fpg.filt.norm.pass.vcf"), emit: vcf_pass
    // path("bcftools/fpg.filt.norm.pass.vcf.tbi"), emit: tbi_pass


   script:

    """
    #!/usr/bin/env bash
   
    if ! [[ -d bcftools ]]; then mkdir bcftools; fi

    #bcftools index $bcf
 
    bcftools view -i 'FILTER="PASS"' $bcf -o bcftools/fpg.filt.norm.pass.bcf
    
    bcftools index bcftools/fpg.filt.norm.pass.bcf

    # create vcf file n bgzip it
    bcftools view -i 'FILTER="PASS"' -Ov -o bcftools/fpg.filt.norm.pass.vcf $bcf
    #bcftools index -t bcftools/fpg.filt.norm.pass.vcf

    """
}



process VCFSNPS2FASTA {

   conda "$params.cacheDir/fpgVcf2FastaEnv"
   publishDir "$params.bcftl", mode: 'copy'


   input:
     path(vcf)
     //path(idx)


   output:
     path("SNPfasta/fpg_snp_aln.fa"), emit: snp_msa
     //path("SNPfasta/fpg.filt.norm.pass.vcf"), emit: vcf_filt
     path("SNPfasta/vcf_infile"), emit: infile


   script:

    """
    #!/usr/bin/env bash
    
    if ! [[ -d SNPfasta ]]; then mkdir SNPfasta; fi   

    #bcftools view -i 'FILTER="PASS"' -Ov -o SNPfasta/fpg.filt.norm.pass.vcf \$bcf

    readlink -f $vcf > SNPfasta/vcf_infile

    python $projectDir/templates/broad-fungalgroup/scripts/SNPs/vcfSnpsToFasta.py --max_amb_samples 1 SNPfasta/vcf_infile > SNPfasta/fpg_snp_aln.fa
    
    """

}





process VCF2PHYLIP {

   conda "$params.cacheDir/fpgCallVariants"
   publishDir "$params.bcftl", mode: 'copy'


   input:
     path(vcf)
     //path(idx)


   output:
     path("vcf2phylip/vcfSNPs.min2.fasta"), emit: snp_aln
     path("vcf2phylip/vcfSNPs.min2.phy"), emit: phy
     path("vcf2phylip/vcfSNPs.min2.used_sites.tsv"), emit: sites
    // path("vcf2phylip/fpg.filt.norm.pass.vcf"), emit: vcf_filt
     path("vcf2phylip/vcfSNPs.min2.fold.fasta"), emit: fold_aln

   script:

    """
    #!/usr/bin/env bash

    if ! [[ -d vcf2phylip ]]; then mkdir vcf2phylip; fi

    python $projectDir/templates/vcf2phylip.py -i $vcf -f -w --output-folder vcf2phylip --output-prefix vcfSNPs

    seqtk seq -Cl60 vcf2phylip/vcfSNPs.min2.fasta > vcf2phylip/vcfSNPs.min2.fold.fasta

    """
}



workflow BCFTOOLS {
  take:
    fa
    bam

  main:
     CALLVARIANTS(fa,bam)
     //INDEXBCF(CALLVARIANTS.out.vcfs)
     SOFTFILTERVCF(CALLVARIANTS.out.vcfs,CALLVARIANTS.out.idx)
     BCFNORM(SOFTFILTERVCF.out.bcf_filt,SOFTFILTERVCF.out.idx_filt,fa) 
     FILTERVCF(BCFNORM.out.bcf_norm,BCFNORM.out.idx_norm)
     VCFSNPS2FASTA(FILTERVCF.out.vcf_pass) //,FILTERVCF.out.tbi_pass)
     VCF2PHYLIP(FILTERVCF.out.vcf_pass) //,FILTERVCF.out.tbi_pass)

  emit:
    bcf_raw = CALLVARIANTS.out.vcfs
    bcf_rawIdx = CALLVARIANTS.out.idx
    filt_bcf = SOFTFILTERVCF.out.bcf_filt
    filt_idx = SOFTFILTERVCF.out.idx_filt
    norm_bcf = BCFNORM.out.bcf_norm
    //norm_idx = BCFNORM.out.idx_norm
    norm_log = BCFNORM.out.log_norm
    pass_bcf = FILTERVCF.out.bcf_pass
    pass_idx = FILTERVCF.out.idx_pass
    pass_vcf = FILTERVCF.out.vcf_pass
    //pass_tbi = FILTERVCF.out.tbi_pass
    msa_snp = VCFSNPS2FASTA.out.snp_msa
    //vcf_z = VCFSNPS2FASTA.out.vcf_filt
    vcf_list = VCFSNPS2FASTA.out.infile
    aln_snp = VCF2PHYLIP.out.snp_aln
    aln_phy = VCF2PHYLIP.out.phy
    snp_sites = VCF2PHYLIP.out.sites
    //z_vcf = VCF2PHYLIP.out.vcf_filt
    aln_fold = VCF2PHYLIP.out.fold_aln
        

}



workflow {
 bam_ch = Channel.fromPath("$params.resultsDir/align/picard/**marked.bam", checkIfExists: true).
 BCFTOOLS(file("$params.refseq"),bam_ch.collect())

}



