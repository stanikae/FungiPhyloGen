#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process RUNSNPEFF {
    tag "snpeff"
    label 'process_high'
    
    conda "$params.cacheDir/fpgCallVariants"
    publishDir "$params.ann", mode: 'copy'

    input:
    path(vcf)
    file(ref)
    file(gbk)

    output:
    path("snpeff/snpeff_ann.vcf"), emit: ann
    // path("snpeff/snpEff_summary.html"), emit: summary, optional: true

    script:
    """
    #!/usr/bin/env bash

    # 1. Create directory structure
    if ! [[ -d snpeff/ref/genomes/ref ]]; then mkdir -p snpeff/ref/genomes/ref; fi

    # 2. Setup Variables
    # FIX: Use \$(...) syntax instead of backticks to avoid syntax errors
    CONDA_BASE="$params.cacheDir/fpgCallVariants"
    
    # We use 'find' to locate the config and jar safely
    configPath=\$(find \${CONDA_BASE}/share -name "snpEff.config" | head -n 1)
    binDir=\$(dirname \$(find \${CONDA_BASE}/share -name "snpEff.jar" | head -n 1))

    # 3. Stage Reference Files
    cp "$ref" snpeff/ref/genomes/ref/ref.fa
    cp "$gbk" snpeff/ref/genomes/ref/genes.gbk
    cp \$configPath snpeff/ref/snpEff.config

    # --- HANDLE VCF COMPRESSION ---
    # Copy input to a .gz file
    cp "$vcf" snpeff/ref/fpg.vcf.gz
    
    # Decompress to plain text for SnpEff
    gunzip -f snpeff/ref/fpg.vcf.gz
    
    # 4. Configure SnpEff
    echo -e "ref.genome : FPG Reference" >> snpeff/ref/snpEff.config

    # 5. Run SnpEff
    cd ./snpeff/ref

    if [[ -f genomes/ref/genes.gbk ]]; then
        # Build database
        java -jar \$binDir/snpEff.jar build -genbank -dataDir genomes -c snpEff.config ref
        
        # Annotate (reading plain text fpg.vcf)
        java -jar \$binDir/snpEff.jar ann ref \
            -noLog -nodownload -onlyProtein \
            -dataDir genomes \
            -c snpEff.config \
            fpg.vcf > ../snpeff_ann.vcf
    fi
    """
}

workflow WFANNOTATESNP {
    take:
    ch_vcf
    ch_ref
    ch_gbk

    main:
    RUNSNPEFF(ch_vcf, ch_ref, ch_gbk)

    emit:
    eff = RUNSNPEFF.out.ann
}
