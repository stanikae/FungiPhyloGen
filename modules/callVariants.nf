#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --- MODULE: VARIANT CALLING ---

process CALLVARIANTS {
    tag "${meta.id}"
    publishDir "$params.bcftl/raw_per_sample", mode: 'copy'

    conda "$params.cacheDir/fpgCallVariants"
    cpus 8
    executor 'slurm'

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
    publishDir "$params.bcftl/filtered_per_sample", mode: 'copy'

    conda "$params.cacheDir/fpgCallVariants"
    cpus 8
    executor 'slurm'

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
    publishDir "$params.bcftl/merged", mode: 'copy'

    conda "$params.cacheDir/fpgCallVariants"
    cpus 12
    executor 'slurm'

    input:
    path bcf
    path idx

    output:
    path "merged.bcf", emit: mge
    path "merged.bcf.csi", emit: mge_idx

    script:
    """
    # List of files passed directly (no -l flag)
    bcftools merge \
        --threads ${task.cpus} \
        -0 \
        ${bcf} \
        -Ob \
        -o merged.bcf

    bcftools index merged.bcf
    """
}

process REHEADERVCF {
    tag "reheader"
    publishDir "$params.bcftl/reheadered", mode: 'copy'

    conda "$params.cacheDir/fpgCallVariants"
    cpus 8
    executor 'slurm'

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
    cpus 12
    executor 'slurm'

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
    bcftools norm \
        -m -any \
        --threads ${task.cpus} \
        -f $ref \
        $bcf \
        -Ob \
        -o normalized.bcf 2> norm.log

    bcftools index normalized.bcf
    """
}



process SOFTFILTERVCF {
    tag "soft_filter"
    publishDir "$params.bcftl/soft_filtered", mode: 'copy'

    conda "$params.cacheDir/fpgCallVariants"
    cpus 12
    executor 'slurm'

    input:
    path bcf
    path idx

    output:
    path "soft_filtered.bcf", emit: bcf_filt
    path "soft_filtered.bcf.csi", emit: idx_filt

    script:
    filter_expression = params.filters[params.genus] ?: params.filters['default']
    """
    # FIX: Removed 'fill-tags' entirely.
    # BCFtools handles strand bias by lowering the QUAL score, 
    # so we rely on QUAL < 30 to catch bias issues.
    
    bcftools filter \
        --threads ${task.cpus} \
        -s 'FAIL' \
        -e '${filter_expression}' \
        -Ob \
        -o soft_filtered.bcf \
        $bcf

    bcftools index soft_filtered.bcf
    """
}



process FILTERVCF {
    tag "extract_pass"
    publishDir "$params.bcftl/final_vcf", mode: 'copy'

    conda "$params.cacheDir/fpgCallVariants"
    cpus 12
    executor 'slurm'

    input:
    path bcf
    path idx

    output:
    path "final.pass.bcf", emit: bcf_pass, optional: true
    path "final.pass.bcf.csi", emit: idx_pass, optional: true
    path "final.pass.vcf.gz", emit: vcf_pass, optional: true
    path "final.pass.vcf.gz.tbi", emit: tbi_pass, optional: true
    path "sample_list.txt", emit: sample_list, optional: true

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

        bcftools query -l final.pass.bcf > sample_list.txt
    fi
    """    
}


process VCFSNPS2FASTA {
    tag "vcf_to_fasta"
    publishDir "$params.bcftl/msa", mode: 'copy'

    conda "$params.cacheDir/fpgVcf2FastaEnv"
    cpus 2
    executor 'slurm'

    input:
    path vcf
    path samples

    output:
    path "fpg_snp_aln.fa", emit: snp_msa
    path "vcf_infile", emit: infile

    script:
    """
    # --- FIX START ---
    # Decompress the VCF because the Python script cannot read .gz
    # We pipe to a new file named 'input.vcf'
    gunzip -c $vcf > input.vcf
    
    # Point the python script to the uncompressed file
    readlink -f input.vcf > vcf_infile
    # --- FIX END ---

    nsamples=\$(cat "$samples" | wc -l)
    amb_samples=\$(awk "BEGIN {print (\$nsamples/100)*10}" | awk '{ print int(\$1) }')

    python $projectDir/templates/broad-fungalgroup/scripts/SNPs/vcfSnpsToFasta.py \
        --max_amb_samples \$amb_samples \
        vcf_infile > fpg_snp_aln.fa
    """
}





process VCF2PHYLIP {
    tag "vcf_to_phylip"
    publishDir "$params.bcftl/phylip", mode: 'copy'

    conda "$params.cacheDir/fpgVcfTools"
    cpus 2
    executor 'slurm'

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
    VCFSNPS2FASTA(FILTERVCF.out.vcf_pass, FILTERVCF.out.sample_list)
    VCF2PHYLIP(FILTERVCF.out.vcf_pass)

    emit:
    pass_vcf  = FILTERVCF.out.vcf_pass
    msa_snp   = VCFSNPS2FASTA.out.snp_msa
    aln_snp   = VCF2PHYLIP.out.snp_aln
    aln_fold  = VCF2PHYLIP.out.fold_aln
}
