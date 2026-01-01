#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GETREPEATS {
   conda "$params.cacheDir/fpgMaskRepeats"
  // publishDir "$params.refMasked", mode: 'copy'
  
  label 'process_medium'

  input:
    file(ref)

  output:
    path('out.delta'), emit: delta
    path('out.coords'), emit: coord

  script:


  """
  #!/usr/bin/env bash
  
  #if ! [[ -d bwa ]]; then mkdir rpts; fi
  
  nucmer --coords $ref $ref 

  """
}




process REPEATSBED {
   conda "$params.cacheDir/fpgMaskRepeats"
  // publishDir "$params.refMasked", mode: 'copy'

  label 'process_medium'

  input:
    file(repeats)

  output:
    path('repeats.bed'), emit: rpts_bed

  script:


  """
  #!/usr/bin/env bash

  #if ! [[ -d bwa ]]; then mkdir rpts; fi


  show-coords -r -T -H $repeats \
  | awk '{if (\$1 != \$3 && \$2 != \$4) print \$0}' \
  | awk '{print \$8\"\\t\"\$1"\\t"\$2}' \
  | sort -k1,1V -k2,3n  > repeats.bed


  """
}





process MASKREF {
   conda "$params.cacheDir/fpgMaskRepeats"
   label 'process_medium'

   publishDir "$params.refMasked", mode: 'copy'


  input:
    file(ref)
    file(bed)

  output:
    path('*masked.fasta'), emit: masked_fa

  script:


  """
  #!/usr/bin/env bash

  if ! [[ -d masked ]]; then mkdir masked; fi

  maskFastaFromBed -fi $ref -bed $bed -fo ${ref.simpleName}_masked.fasta

  """
}




process INDEXREF {
  conda "$params.cacheDir/fpgAlign"
  label 'process_medium'

  publishDir "$params.refIndex", mode: 'copy'

  input:
  file(ref)

  output:
   path(bwa) , emit: indx
   val true, emit: prs

  script:
  """
  #!/usr/bin/env bash
  
  if ! [[ -d bwa ]]; then mkdir bwa; fi

  bwa index -p bwa/${ref.simpleName} $ref

  """
}


/*
workflow {
  GETREPEATS(file("$params.refseq"))  
  REPEATSBED(GETREPEATS.out.delta)  
  MASKREF(file("$params.refseq"),REPEATSBED.out.rpts_bed)
  INDEXREF(MASKREF.out.masked_fa) 
}
*/


