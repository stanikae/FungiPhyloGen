#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process CALLVARIANTS {

  cpus 8
  executor 'slurm'

  //tag "$sampleId"

  conda "$params.cacheDir/fpgCallVariants"
  publishDir "$params.bcftl", mode: 'copy'

  input:
    path(ref)
    path(bam)
    //path(idx)
    //tuple val(sampleId), file(bam)


  output:
     path('*call.bcf'), emit: bcf
     path('*call.bcf.csi'), emit: bcf_idx

  script:

  """
   #!/usr/bin/env bash
   
   dirNam=`basename -s _sorted_marked.bam "$bam"`
   #if ! [[ -d \$dirNam ]]; then mkdir \$dirNam; fi

   
   bcftools mpileup --threads ${task.cpus} -a AD -Q 30 -f $ref $bam -Ou \\
   | bcftools call --threads ${task.cpus} -a GQ --ploidy $params.ploidy -m -Ob -o \${dirNam}_call.bcf

   bcftools index \${dirNam}_call.bcf

  """
}

// 2023-10-07: remove -v from bcftools call to address issue of missing variants
// suggestion taken from biostars: https://www.biostars.org/p/9466376/


// 2023-10-08: Filter individual samples first before merging
// As suggested at https://www.biostars.org/p/411766/ this will minimize false positives unique to a single sample
// 2024-05-14: Apply filter individual VCFs before merging across the board including for C. auris


process FILTERSAMPLE {

   conda "$params.cacheDir/fpgCallVariants"
   publishDir "$params.bcftl", mode: 'copy'

   cpus 12
   executor 'slurm'

   input:
     path(bcf)
     path(idx)


   output:
     path("bcftools/filtered/*.filtered.bcf"), emit: smpl_filt
     path("bcftools/filtered/*.filtered.bcf.csi"), emit: smpl_idx


   script:
        if(params.genus=="Candida"){

          """
          #!/usr/bin/env bash

          if ! [[ -d bcftools/filtered ]]; then mkdir -p bcftools/filtered; fi

	  sampleNam=`basename -s _call.bcf "$bcf"`

      
          #bcftools view -v 'snps' --threads ${task.cpus} $bcf | \\
           #             bcftools filter -i 'QUAL/DP>2.0 && FS<=60 && MQ>=30 && DP>=10' -g8 -G10 -Ob -o bcftools/filtered/\${sampleNam}.filtered.bcf

          bcftools view --threads ${task.cpus} $bcf -Ob -o bcftools/filtered/\${sampleNam}.filtered.bcf

          bcftools index bcftools/filtered/\${sampleNam}.filtered.bcf

    	  """
   	} else {

    	  """
    	  #!/usr/bin/env bash

    	  if ! [[ -d bcftools/filtered ]]; then mkdir -p bcftools/filtered; fi

    	  sampleNam=`basename -s _call.bcf "$bcf"`


    	  bcftools view -v 'snps' --threads ${task.cpus} $bcf | \\
       		bcftools filter -i 'MQ>=40 && DP>=10 && QUAL>=30 && (MQSBZ > -2 || MQSBZ < 2) && FMT/AD > 10' -g8 -G10 -Ob -o bcftools/filtered/\${sampleNam}.filtered.bcf


    	  bcftools index bcftools/filtered/\${sampleNam}.filtered.bcf


    """
   }   
}


/*
process CALLVARIANTSgrp {
 // tag "$sampleId"

  cpus 12
  executor 'slurm'
  
  conda "$params.cacheDir/fpgCallVariants"
  publishDir "$params.bcftl", mode: 'copy'


  input:
    //val ready
    path(ref)
    path(bam)

  output:
    path("bcftools/fpg.call.bcf"), emit: vcfs
    path("bcftools/fpg.call.bcf.csi"), emit: idx


  script:
  
  """
  #!/usr/bin/env bash
  
  if ! [[ -d bcftools ]]; then mkdir bcftools; fi


  bcftools mpileup --threads ${task.cpus} -a AD -Q 30 -f $ref $bam -Ou \\
  | bcftools call --threads ${task.cpus} --ploidy $params.ploidy -mv -Ob -o bcftools/fpg.call.bcf

  bcftools index bcftools/fpg.call.bcf

  """
}
*/





process BCFMERGE {
  
  cpus 12
  executor 'slurm'

  conda "$params.cacheDir/fpgCallVariants"
  publishDir "$params.bcftl", mode: 'copy'


  input:
    path(bcf)
    path(idx)

  output:
    path("bcfmerge/fpg.call.bcf"), emit: mge
    path("bcfmerge/fpg.call.bcf.csi"), emit: mge_idx


  script:

  """
  #!/usr/bin/env bash

  if ! [[ -d bcfmerge ]]; then mkdir bcfmerge; fi


  bcftools merge -0 $bcf --threads ${task.cpus} -Ob -o bcfmerge/fpg.call.bcf
  bcftools index bcfmerge/fpg.call.bcf
 
 
 """

}




process REHEADERVCF {
// tag "$sampleId"

  cpus 8
  executor 'slurm'

  conda "$params.cacheDir/fpgCallVariants"
  publishDir "$params.bcftl", mode: 'copy'


  input:
    path(bcf)
    path(idx)

  output:
    path("bcftools/fpg.call.reheader.bcf"), emit: reh_bcf
    path("bcftools/fpg.call.reheader.bcf.csi"), emit: reh_idx

  script:


  """
  #!/usr/bin/env bash

  if ! [[ -d bcftools ]]; then mkdir bcftools; fi

  # rehead (rename) samples i.e. remove paths and suffix
  paste -d "\t" <(bcftools query -l "$bcf") <(bcftools query -l "$bcf" | sed 's/\\/.*\\///' | sed 's/_sorted.*//') > ./reheader.txt
  bcftools reheader -s ./reheader.txt -o bcftools/fpg.call.reheader.bcf "$bcf"

  bcftools index bcftools/fpg.call.reheader.bcf  
 
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



// 2023-10-08: Add bcftools -m -any option to get split multiallelic variants

process BCFNORM {
   cpus 12
   executor 'slurm'

   conda "$params.cacheDir/fpgCallVariants"
   publishDir "$params.bcftl", mode: 'copy'

   input:
     path(bcf)
     path(idx)
     path(ref)


   output:
     path("bcftools/fpg.norm.bcf"), emit: bcf_norm
     path("bcftools/fpg.norm.bcf.csi"), emit: idx_norm
     path("bcftools/norm.log"), emit: log_norm

   script:

   """
   #!/usr/bin/env bash

   if ! [[ -d bcftools ]]; then mkdir bcftools; fi

   bcftools norm -m -any --threads ${task.cpus} -f $ref $bcf -o bcftools/fpg.norm.bcf 2> bcftools/norm.log

   bcftools index bcftools/fpg.norm.bcf

   """


}




process SOFTFILTERVCF {
   
   conda "$params.cacheDir/fpgCallVariants"
   publishDir "$params.bcftl", mode: 'copy'

   cpus 12
   executor 'slurm'

   input:
     path(bcf)
     path(idx)


   output:
     path("bcftools/fpg.filt.bcf"), emit: bcf_filt
     path("bcftools/fpg.filt.bcf.csi"), emit: idx_filt


   script:
	if(params.genus=="Candida"){

  	  """
    	  #!/usr/bin/env bash

    	  if ! [[ -d bcftools ]]; then mkdir bcftools; fi

          # get count of samples in bcf file
          nsamples=\$(bcftools query --list-samples $bcf | wc -l)
          nac=`awk "BEGIN {print (90/100)*\$nsamples}" | awk '{ print int(\$1) }'`

    
   	  #bcftools view -v 'snps' --threads ${task.cpus} $bcf \\
          # | bcftools +fill-tags -- -t FORMAT/VAF \\
          # | bcftools filter --threads ${task.cpus} -s 'LowQual' -i 'FMT/VAF>0.8 && FMT/GQ>50' -g8 -G10 -Ob \\
          # | bcftools filter --threads ${task.cpus} -s 'LowQual' -i 'MQ>=40 && DP>=10 && QUAL>=30 && (MQSBZ > -2 || MQSBZ < 2) && FMT/AD > 10' -g8 -G10 -Ob | \\
          #   bcftools filter --threads ${task.cpus} -s 'LowQual' -i "QUAL>=50 & AD[*:1]>=25 & (AC[0]+AC[1])>=\$nac & DP>=10" -g8 -G10 -Ob -o bcftools/fpg.filt.bcf #| \\

	
	    bcftools view -v 'snps' --threads ${task.cpus} $bcf \\
              | bcftools +fill-tags -- -t FORMAT/VAF \\
              | bcftools filter --threads ${task.cpus} -s 'LowQual' -i '(QUAL/MAX(AD[:1]))>2.0 && FS<60 && MQ>=40 && DP>=10' -g8 -G10 -Ob \\
              | bcftools filter --threads ${task.cpus} -s 'LowQual' -i '(MQSBZ > -2 || MQSBZ < 2)' -g8 -G10 -Ob \\
	      | bcftools filter --threads ${task.cpus} -S . -i 'FMT/GQ>20' -g8 -G10 -Ob \\
              | bcftools +fill-tags -- -t F_MISSING | bcftools filter --threads ${task.cpus} -s 'LowQual' -i 'F_MISSING<0.25' -g8 -G10 -Ob \\
              | bcftools filter -s 'LowQual' -e 'AC==0' -g8 -G10 -Ob -o bcftools/fpg.filt.bcf
	

	 # main filters/paramters for the workflow -- 2024-05-14
	 #
         #   bcftools view -v 'snps' --threads ${task.cpus} $bcf \\
         #    | bcftools +fill-tags -- -t FORMAT/VAF \\
         #    | bcftools filter --threads ${task.cpus} -s 'LowQual' -i '(QUAL/MAX(AD[:1]))>2.0 && FS<60 && MQ>=40 && DP>=10' -g8 -G10 -Ob \\
         #    | bcftools filter --threads ${task.cpus} -s 'LowQual' -i '(MQSBZ > -2 || MQSBZ < 2)' -g8 -G10 -Ob \\
         #    | bcftools filter --threads ${task.cpus} -S . -i 'FMT/GQ>20 & AD[:1]>10' -g8 -G10 -Ob \\
	 #    | bcftools +fill-tags -- -t F_MISSING | bcftools filter -s 'LowQual' -e 'AC==0' -g8 -G10 -Ob -o bcftools/fpg.filt.bcf             

	

	# | bcftools +fill-tags -- -t F_MISSING |  bcftools filter --threads ${task.cpus} -s 'LowQual' -i 'F_MISSING=0' -g8 -G10 -Ob \\
        #     | bcftools filter -s 'LowQual' -e 'AC==0' -g8 -G10 -Ob -o bcftools/fpg.filt.bcf

	 # | bcftools filter --threads ${task.cpus} -S . -i 'FMT/GQ>20 & FMT/VAF>=0.8 & AD[:1]>10' -g8 -G10 -Ob \\

	  #bcftools view -v 'snps' --threads ${task.cpus} $bcf \\
          #   | bcftools +fill-tags -- -t FORMAT/VAF \\
          #   | bcftools filter --threads ${task.cpus} -s 'LowQual' -i 'QUAL>20 && FS<60 && MQ>=40 && DP>=10' -g8 -G10 -Ob \\
          #   | bcftools filter --threads ${task.cpus} -s 'LowQual' -i 'MQ>=30 && DP>=10 && QUAL>=30 && (MQSBZ > -2 || MQSBZ < 2) && FMT/AD[:1]>=10' -g8 -G10 -Ob \\
          #   | bcftools filter --threads ${task.cpus} -S . -i 'FMT/GQ>50 && FMT/VAF>=0.8 && AD[:1]>10' -g8 -G10 -Ob \\
          #   | bcftools +fill-tags -- -t F_MISSING |  bcftools filter --threads ${task.cpus} -s 'LowQual' -i 'F_MISSING<0.25' -g8 -G10 -Ob \\
          #   | bcftools filter -s 'LowQual' -e 'AC==0' -g8 -G10 -Ob -o bcftools/fpg.filt.bcf
	
          #sAVG  
          # | bcftools filter --threads ${task.cpus} -s 'LowQual' -i '(QUAL/MAX(AD[:1]))>2.0 && FS<60 && MQ>=40 && DP>=10' -g8 -G10 -Ob \\	  
 	  # | bcftools filter --threads ${task.cpus} -s 'LowQual' -i '(QUAL/AD[:1])>2.0 && FS<60 && MQ>=40 && DP>=10' -g8 -G10 -Ob \\
	  # | bcftools filter --threads ${task.cpus} -S . -i 'FMT/GQ>50 && FMT/VAF>=0.8 && AD[:1]>10' -g8 -G10 -Ob \\
	  # F_MISSING<0.25

           #  | bcftools filter -S . -e 'FMT/VAF<=0.8 | FMT/GQ<50 | AD[:1]<10' \\
           #  | bcftools filter --threads ${task.cpus} -s 'LowQual' -i 'FMT/VAF>0.8 && FMT/GQ>50' -g8 -G10 -Ob \\
           #  | bcftools +fill-tags -- -t F_MISSING | bcftools filter -s 'LowQual' -e 'F_MISSING>=0.25' -Ob \\
           #  | bcftools filter -S . -e 'DP<10 | FMT/GQ<50' \\
           #  | bcftools filter -s 'LowQual' -e 'AC==0 || AC==AN' -g8 -G10 -Ob -o bcftools/fpg.filt.bcf


             #| bcftools +setGT -- -t q -n . -i 'GT="1" & FMT/VAF>0.8 & FMT/GQ>=50 & AD[:1]>=10' \\
           
             # | bcftools filter --threads ${task.cpus} -s 'LowQual' -i "(AC[0]/AN)*100>=10" -g8 -G10 | \\

	     #| bcftools filter --threads ${task.cpus} -s 'LowQual' -i  "QUAL>=30 & (AC[0]/AN)*100>=90 & DP>=10" -g8 -G10 -Ob -o bcftools/fpg.filt.bcf
             
	    #| bcftools +setGT -- -t q -n . -i 'GT="1" & FMT/VAF>0.8 & FMT/GQ>50 & AD[:1]<10' \\
	     #| bcftools filter --threads ${task.cpus} -s 'LowQual' -i 'FMT/VAF>0.8 && FMT/GQ>50' -g8 -G10 -Ob \\
             #| bcftools +fill-tags -- -t F_MISSING | bcftools filter -s 'LowQual' -e 'F_MISSING>=0.25' -Ob -o bcftools/fpg.filt.bcf



          #   bcftools filter --threads ${task.cpus} -s 'LowQual' -i  "(QUAL/DP)>=2 & MQ>=40 & FS<=60 & AD[*:1]>=25 & ((AC[0]+AC[1])/AN)*100>=90 & DP>=10" -g8 -G10 -Ob -o bcftools/fpg.filt.bcf
               
               
               #| bcftools filter --threads ${task.cpus} -s 'LowQual' -i '((AC[0]+AC[1])/AN)*100>=75' -g8 -G10 -Ob \\
               #bcftools filter --threads ${task.cpus} -s 'LowQual' -i  "QUAL>=50 & AD[*:1]>=25 & AC>=\$nac & DP>=10" -g8 -G10 -Ob -o bcftools/fpg.filt.bcf
	       
               #bcftools filter --threads ${task.cpus} -s 'LowQual' -e  'QUAL/DP<2.0 || FS>60 || MQ<40 || DP<10' -g8 -G10 -Ob | \\
	       #bcftools filter --threads ${task.cpus} -s 'LowQual' -i  'QUAL>=50 && AD[*:1]>=25 && AC>=3 && DP>=10' -g8 -G10 -Ob -o bcftools/fpg.filt.bcf
          
               #| bcftools filter --threads ${task.cpus} -s 'LowQual' -e  'QUAL/DP<2.0 || FS>60 || MQ<40 || DP<10' -g8 -G10 -Ob \\
    	  
	        #// 'QUAL/DP>2.0 && FS<=60 && MQ>=40 && DP>=10 && AD>0.8 && QUAL>=30'
    		#bcftools view -v 'snps' --threads ${task.cpus} $bcf | \\
       		#	bcftools filter -s 'LowQual' -i 'QUAL/DP>2.0 && FS<=60 && MQ>=30 && DP>=10' -g8 -G10 -Ob -o bcftools/fpg.filt.bcf
                
                #// bcftools view --threads ${task.cpus} $bcf -Ob -o bcftools/fpg.filt.bcf   
    	  #// && MQ>=30 && F_MISSING<0.9

   	  bcftools index bcftools/fpg.filt.bcf

    """
   }
   else {

    """
    #!/usr/bin/env bash

    if ! [[ -d bcftools ]]; then mkdir bcftools; fi


    #bcftools view -v 'snps' --threads ${task.cpus} $bcf | \\
     #  bcftools filter -s 'LowQual' -i 'FS<=60 && MQ>=40 && DP>=10 && QUAL>=30' -g8 -G10 -Ob -o bcftools/fpg.filt.bcf

    bcftools view --threads ${task.cpus} $bcf -Ob -o bcftools/fpg.filt.bcf
    bcftools index bcftools/fpg.filt.bcf

    """

  }
}




process FILTERVCF {
 
   cpus 12
   executor 'slurm'

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
     path("bcftools/fpg.sampleList.txt"), emit: sample_list


   script:

    """
    #!/usr/bin/env bash
   
    if ! [[ -d bcftools ]]; then mkdir bcftools; fi

    #bcftools index $bcf
 
    bcftools view --threads ${task.cpus} -i 'FILTER="PASS"' $bcf -Ob -o bcftools/fpg.filt.norm.pass.bcf
    # bcftools view --threads ${task.cpus} $bcf -Ob -o bcftools/fpg.filt.norm.pass.bcf
    bcftools index bcftools/fpg.filt.norm.pass.bcf

    # create vcf file n bgzip it
    bcftools view --threads ${task.cpus} -i 'FILTER="PASS"' $bcf -Ov -o bcftools/fpg.filt.norm.pass.vcf
    bcftools index -t bcftools/fpg.filt.norm.pass.vcf
    

    #bcftools view --threads ${task.cpus} -Ov -o bcftools/fpg.filt.norm.pass.vcf $bcf
    
    # get count of samples in bcf file
    bcftools query --list-samples bcftools/fpg.filt.norm.pass.bcf > bcftools/fpg.sampleList.txt


    """
}



process VCFSNPS2FASTA {

   conda "$params.cacheDir/fpgVcf2FastaEnv"
   publishDir "$params.bcftl", mode: 'copy'


   input:
     path(vcf)
     path(fpg)
     //path(idx)


   output:
     path("SNPfasta/fpg_snp_aln.fa"), emit: snp_msa
     //path("SNPfasta/fpg.filt.norm.pass.vcf"), emit: vcf_filt
     path("SNPfasta/vcf_infile"), emit: infile


   script:

    """
    #!/usr/bin/env bash
    
    if ! [[ -d SNPfasta ]]; then mkdir SNPfasta; fi   

    
    readlink -f $vcf > SNPfasta/vcf_infile
    nsamples=\$(cat "$fpg" | wc -l)
    amb_samples=`awk "BEGIN {print (\$nsamples/100)*10}" | awk '{ print int(\$1) }'`  #// `echo \$(( \$nsamples/100*10 ))`   #// `awk "BEGIN {print (\$nsamples/100)*10}"`

    python $projectDir/templates/broad-fungalgroup/scripts/SNPs/vcfSnpsToFasta.py --max_amb_samples \$amb_samples SNPfasta/vcf_infile > SNPfasta/fpg_snp_aln.fa
    
    """

}





process VCF2PHYLIP {

   conda "$params.cacheDir/fpgVcfTools"
   publishDir "$params.bcftl", mode: 'copy'


   input:
     path(vcf)
     //path(idx)


   output:
     path("vcf2phylip/vcfSNPs.min1.fasta"), emit: snp_aln
     path("vcf2phylip/vcfSNPs.min1.phy"), emit: phy
     path("vcf2phylip/vcfSNPs.min1.used_sites.tsv"), emit: sites
     path("vcf2phylip/vcfSNPs.min1.fold.fasta"), emit: fold_aln

   script:

    """
    #!/usr/bin/env bash

    if ! [[ -d vcf2phylip ]]; then mkdir vcf2phylip; fi
  
     mainPath=`echo "$projectDir"`

    \$mainPath/templates/vcf2phylip.py -i $vcf -f -w -m 1 --output-folder vcf2phylip --output-prefix vcfSNPs

    if [[ -f "vcf2phylip/vcfSNPs.min1.fasta" ]]; then 
	seqtk seq -Cl60 vcf2phylip/vcfSNPs.min1.fasta > vcf2phylip/vcfSNPs.min1.fold.fasta
	
    else
	echo ">no variants to process" > vcf2phylip/vcfSNPs.min1.fasta
        echo ">no variants to process" > vcf2phylip/vcfSNPs.min1.phy
        echo ">no variants to process" > vcf2phylip/vcfSNPs.min1.used_sites.tsv
	echo ">no variants to process" > vcf2phylip/vcfSNPs.min1.fold.fasta
    fi

    """
}




workflow BCFTOOLS {
  take:
    fa
    bam

  main:
     //CALLVARIANTSgrp(fa,bam)
     // //INDEXBCF(CALLVARIANTS.out.vcfs)
     CALLVARIANTS(fa,bam)
     //CALLVARIANTS(fa,bam.
     //			 map { file -> def key = file.name.toString().tokenize('_').get(0).replace('[','')
     //        			return tuple(key, file)} 
     //)
     FILTERSAMPLE(CALLVARIANTS.out.bcf,CALLVARIANTS.out.bcf_idx)
     BCFMERGE(FILTERSAMPLE.out.smpl_filt.collect(),FILTERSAMPLE.out.smpl_idx.collect())
     REHEADERVCF(BCFMERGE.out.mge,BCFMERGE.out.mge_idx)
     //REHEADERVCF(CALLVARIANTSgrp.out.vcfs,CALLVARIANTSgrp.out.idx)
     BCFNORM(REHEADERVCF.out.reh_bcf,REHEADERVCF.out.reh_idx,fa)
     SOFTFILTERVCF(BCFNORM.out.bcf_norm,BCFNORM.out.idx_norm)
     //BCFNORM(SOFTFILTERVCF.out.bcf_filt,SOFTFILTERVCF.out.idx_filt,fa) 
     FILTERVCF(SOFTFILTERVCF.out.bcf_filt,SOFTFILTERVCF.out.idx_filt)
     VCFSNPS2FASTA(FILTERVCF.out.vcf_pass,FILTERVCF.out.sample_list) //,FILTERVCF.out.tbi_pass)
     VCF2PHYLIP(FILTERVCF.out.vcf_pass) //,FILTERVCF.out.tbi_pass)

  emit:
    //bcf_raw = CALLVARIANTSgrp.out.vcfs
    //bcf_idx = CALLVARIANTSgrp.out.idx
    bcf_ind = CALLVARIANTS.out.bcf
    bcf_rawIdx = CALLVARIANTS.out.bcf_idx
    smpl_fbcf = FILTERSAMPLE.out.smpl_filt
    smpl_fidx = FILTERSAMPLE.out.smpl_idx
    bcf_mge = BCFMERGE.out.mge
    bcf_idx = BCFMERGE.out.mge_idx
    bcf_reh = REHEADERVCF.out.reh_bcf
    bcf_rehIdx = REHEADERVCF.out.reh_idx
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

