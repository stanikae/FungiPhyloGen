#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --- MODULE: VARIANT CALLING ---

process CALLVARIANTS {
    tag "${meta.id}"
    publishDir "$params.bcftl/raw_per_sample", mode: 'symlink'

    conda "$params.cacheDir/fpgCallVariants"
    label 'process_high'

    input:
    tuple val(meta), path(bam)
    path ref

    output:
    tuple val(meta), path("${meta.id}.raw.bcf"), emit: bcf
    tuple val(meta), path("${meta.id}.raw.bcf.csi"), emit: bcf_idx

    script:
    """
    # FIX: Added 'ADF,ADR' to the annotation tags (-a).
    # This provides the raw forward/reverse counts needed to calculate FS later.
    
    bcftools mpileup \
        --threads ${task.cpus} \
        -a AD,DP,SP,ADF,ADR \
        -Q 30 \
        -f $ref \
        $bam \
        -Ou | \
    bcftools call \
        --threads ${task.cpus} \
        -a GQ \
        --ploidy $params.ploidy \
        -m \
        -Ob \
        -o ${meta.id}.raw.bcf

    bcftools index ${meta.id}.raw.bcf
    """
}

process FILTERSAMPLE {
    tag "${meta.id}"
    publishDir "$params.bcftl/filtered_per_sample", mode: 'symlink'

    conda "$params.cacheDir/fpgCallVariants"
    label 'process_high'


    input:
    tuple val(meta), path(bcf), path(bcf_idx)

    output:
    path("${meta.id}.filtered.bcf"), emit: smpl_filt
    path("${meta.id}.filtered.bcf.csi"), emit: smpl_idx

    script:
    """
    # Debug print
    echo "Processing input BCF: $bcf"
    
    bcftools filter \
        --threads ${task.cpus} \
        -i 'TYPE="snp" && QUAL>=20 && INFO/DP>=5' \
        -Ob \
        -o ${meta.id}.filtered.bcf \
        $bcf

    bcftools index ${meta.id}.filtered.bcf
    """
}

process BCFMERGE {
    tag "merge"
    publishDir "$params.bcftl/merged", mode: 'symlink'

    conda "$params.cacheDir/fpgCallVariants"
    label 'process_high'

    input:
    path bcf
    path idx

    output:
    path "merged.bcf", emit: mge
    path "merged.bcf.csi", emit: mge_idx

    script:
    """
    # Per-sample inputs are sparse SNP-only BCFs. The -0 option represents
    # absent SNP records as reference genotypes when constructing the
    # population SNP matrix.
    #
    # Sparse record absence cannot distinguish confident reference calls
    # from insufficient coverage. Therefore, downstream profiles use AC,
    # AN and AF for population representation and do not use F_MISSING.

    bcftools merge \
        --threads ${task.cpus} \
        -0 \
        ${bcf} \
        -Ob \
        -o merged.bcf

    bcftools index --threads ${task.cpus} merged.bcf
    """
}

process REHEADERVCF {
    tag "reheader"
    publishDir "$params.bcftl/reheadered", mode: 'symlink'

    conda "$params.cacheDir/fpgCallVariants"
    label 'process_high'


    input:
    path bcf
    path idx

    output:
    path "reheadered.bcf", emit: reh_bcf
    path "reheadered.bcf.csi", emit: reh_idx

    script:
    """
    bcftools query -l "${bcf}" | awk '{ old=\$1; sub(/\\.filtered\$/, "", \$1); print old "\t" \$1 }' > reheader.txt

    bcftools reheader \
        -s reheader.txt \
        -o reheadered.bcf \
        "${bcf}"

    bcftools index reheadered.bcf
    """
}

process BCFNORM {
    tag "normalize"
    publishDir "$params.bcftl/normalized", mode: 'copy'

    conda "$params.cacheDir/fpgCallVariants"
    label 'process_high'

    input:
    path bcf
    path idx
    path ref

    output:
    path "normalized.bcf", emit: bcf_norm
    path "normalized.bcf.csi", emit: idx_norm
    path "norm.log", emit: log_norm

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    bcftools norm \
        -m -any \
        --threads ${task.cpus} \
        -f "$ref" \
        "$bcf" \
        -Ou \
        2> norm.log \
    | bcftools +fill-tags \
        -Ob \
        -o normalized.bcf \
        -- \
        -t AC,AN,AF

    bcftools index \
        --threads ${task.cpus} \
        normalized.bcf
    """
}


process SOFTFILTERVCF {
    tag "soft_filter"
    publishDir "$params.bcftl/soft_filtered", mode: 'copy'

    conda "$params.cacheDir/fpgCallVariants"
    label 'process_high'

    input:
    path bcf
    path idx

    output:
    path "soft_filtered.bcf", emit: bcf_filt
    path "soft_filtered.bcf.csi", emit: idx_filt

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Apply the selected haploid population-SNP profile to the normalized
    # multi-sample BCF. AC, AN, and AF have been materialized by fill-tags.
    #
    # Effective exclusion expression:
    # ${params.active_filter}

    bcftools filter \
        --threads ${task.cpus} \
        -s FAIL \
        -e '${params.active_filter}' \
        -Ob \
        -o soft_filtered.bcf \
        "$bcf"

    bcftools index \
        --threads ${task.cpus} \
        soft_filtered.bcf
    """
}


process FILTERVCF {
    tag "extract_pass"
    publishDir "$params.bcftl/final_vcf", mode: 'copy'

    conda "$params.cacheDir/fpgCallVariants"
    label 'process_high'

    input:
    path bcf
    path idx

    output:
    path "final.pass.bcf", emit: bcf_pass, optional: true
    path "final.pass.bcf.csi", emit: idx_pass, optional: true
    path "final.pass.vcf.gz", emit: vcf_pass, optional: true
    path "final.pass.vcf.gz.tbi", emit: tbi_pass, optional: true

    script:
    """
    #!/usr/bin/env bash
    
    bcftools view --threads ${task.cpus} -f 'PASS' -Ob -o tmp.pass.bcf $bcf

    N_VARIANTS=\$(bcftools view -H tmp.pass.bcf | wc -l) 

    if [ "\$N_VARIANTS" -gt 0 ]; then
        mv tmp.pass.bcf final.pass.bcf
        bcftools index final.pass.bcf

        bcftools view -Oz -o final.pass.vcf.gz final.pass.bcf
        bcftools index -t final.pass.vcf.gz

    fi
    """    
}


process VCFSNPS2FASTA {
    tag "vcf_to_fasta"
    publishDir "$params.bcftl/msa", mode: 'copy'

    conda "$params.cacheDir/fpgVcf2FastaEnv"
    label 'process_medium'

    input:
    path vcf

    output:
    path "fpg_snp_aln.fa", emit: snp_msa
    path "fpg_snp_sites.tsv", emit: sites
    path "fpg_snp_alignment_qc.tsv", emit: alignment_qc
    path "fpg_snp_distances.tsv", emit: distances
    path "fpg_snp_comparable_sites.tsv", emit: comparable_sites
    path "fpg_snp_manifest.json", emit: manifest

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    GIT_COMMIT=\$(git -C "$projectDir" rev-parse HEAD 2>/dev/null || printf 'unknown')

    python "$projectDir/scripts/fpg_vcf_to_snp_alignment.py" \
        --vcf "$vcf" \
        --output-prefix fpg_snp \
        --pipeline-version-file "$projectDir/VERSION" \
        --git-commit "\$GIT_COMMIT"
    """
}


process VCF2PHYLIP {
    tag "vcf_to_phylip"
    publishDir "$params.bcftl/phylip", mode: 'copy'

    conda "$params.cacheDir/fpgVcfTools"
    label 'process_medium'

    input:
    path vcf

    output:
    path "vcf2phylip/vcfSNPs.min1.fasta", emit: snp_aln
    path "vcf2phylip/vcfSNPs.min1.phy", emit: phy
    path "vcf2phylip/vcfSNPs.min1.used_sites.tsv", emit: sites
    path "vcf2phylip/vcfSNPs.min1.fold.fasta", emit: fold_aln

    script:
    """
    #!/usr/bin/env bash
    mkdir vcf2phylip
    vcf2phylip.py -i $vcf -f -w -m 1 --output-folder vcf2phylip --output-prefix vcfSNPs

    if [[ -f "vcf2phylip/vcfSNPs.min1.fasta" ]]; then
        seqtk seq -Cl60 vcf2phylip/vcfSNPs.min1.fasta > vcf2phylip/vcfSNPs.min1.fold.fasta
    else
        echo ">no_variants_found" > vcf2phylip/vcfSNPs.min1.fasta
        touch vcf2phylip/vcfSNPs.min1.phy
        touch vcf2phylip/vcfSNPs.min1.used_sites.tsv
        echo ">no_variants_found" > vcf2phylip/vcfSNPs.min1.fold.fasta
    fi
    """
}

workflow BCFTOOLS {
    take:
    ch_bam 
    ch_ref 

    main:
    CALLVARIANTS(ch_bam, ch_ref)
    FILTERSAMPLE(CALLVARIANTS.out.bcf.join(CALLVARIANTS.out.bcf_idx))

    smpl_filt_ch = FILTERSAMPLE.out.smpl_filt.collect()
    smpl_idx_ch  = FILTERSAMPLE.out.smpl_idx.collect()
    BCFMERGE(smpl_filt_ch, smpl_idx_ch)

    REHEADERVCF(BCFMERGE.out.mge, BCFMERGE.out.mge_idx)
    BCFNORM(REHEADERVCF.out.reh_bcf, REHEADERVCF.out.reh_idx, ch_ref)
    SOFTFILTERVCF(BCFNORM.out.bcf_norm, BCFNORM.out.idx_norm)
    FILTERVCF(SOFTFILTERVCF.out.bcf_filt, SOFTFILTERVCF.out.idx_filt)
    VCFSNPS2FASTA(FILTERVCF.out.vcf_pass)
    VCF2PHYLIP(FILTERVCF.out.vcf_pass)

    emit:
    pass_vcf             = FILTERVCF.out.vcf_pass
    msa_snp              = VCFSNPS2FASTA.out.snp_msa
    snp_sites            = VCFSNPS2FASTA.out.sites
    snp_alignment_qc     = VCFSNPS2FASTA.out.alignment_qc
    snp_distances        = VCFSNPS2FASTA.out.distances
    snp_comparable_sites = VCFSNPS2FASTA.out.comparable_sites
    snp_manifest         = VCFSNPS2FASTA.out.manifest
    aln_snp              = VCF2PHYLIP.out.snp_aln
    aln_fold             = VCF2PHYLIP.out.fold_aln
}
