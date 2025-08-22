#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --- MODULE: VARIANT CALLING ---

process CALLVARIANTS {
    tag "$sample_id"
    publishDir "$params.bcftl/raw_per_sample", mode: 'copy'

    conda "$params.cacheDir/fpgCallVariants"
    cpus 8
    executor 'slurm'

    input:
    tuple val(sample_id), path(bam)
    path ref

    output:
    tuple val(sample_id), path("${sample_id}.raw.bcf"), emit: bcf
    tuple val(sample_id), path("${sample_id}.raw.bcf.csi"), emit: bcf_idx

    script:
    """
    bcftools mpileup \\
        --threads ${task.cpus} \\
        -a AD,DP,SP \\
        --output-tags MQ,SOR \\
        -Q 30 \\
        -f $ref \\
        $bam \\
        -Ou | \\
    bcftools call \\
        --threads ${task.cpus} \\
        -a GQ \\
        --ploidy $params.ploidy \\
        -m \\
        -Ob \\
        -o ${sample_id}.raw.bcf

    bcftools index ${sample_id}.raw.bcf
    """
}

process FILTERSAMPLE {
    tag "$sample_id"
    publishDir "$params.bcftl/filtered_per_sample", mode: 'copy'

    conda "$params.cacheDir/fpgCallVariants"
    cpus 8
    executor 'slurm'

    input:
    tuple val(sample_id), path(bcf)
    path bcf_idx // Dummy input to ensure index is staged

    output:
    path("${sample_id}.filtered.bcf"), emit: smpl_filt
    path("${sample_id}.filtered.bcf.csi"), emit: smpl_idx

    script:
    // Apply a consistent, minimal hard filter to remove obvious junk before merging.
    """
    bcftools filter \\
        --threads ${task.cpus} \\
        -i 'TYPE="snp" && QUAL>=20 && DP>=5' \\
        -Ob \\
        -o ${sample_id}.filtered.bcf \\
        $bcf

    bcftools index ${sample_id}.filtered.bcf
    """
}

process BCFMERGE {
    tag "merge_${bcf.size()}_samples"
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
    bcftools merge \\
        --threads ${task.cpus} \\
        -0 \\
        -l ${bcf.join(' ')} \\
        -Ob \\
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
    # Create a mapping file to rename samples (e.g., remove suffix)
    bcftools query -l "$bcf" | sed -E 's/\\.filtered$//' | awk '{print \$1 "\t" \$1}' > reheader.txt

    bcftools reheader \\
        -s reheader.txt \\
        -o reheadered.bcf \\
        "$bcf"

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
    bcftools norm \\
        -m -any \\
        --threads ${task.cpus} \\
        -f $ref \\
        $bcf \\
        -Ob \\
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
    // Dynamically select the filter expression from the config file
    filter_expression = params.filters[params.genus] ?: params.filters['default']
    """
    # Step 1: Add the QD (Quality by Depth) tag, which is crucial for good filtering.
    bcftools +fill-tags $bcf -Ob -o tmp.bcf -- -t QD
    bcftools index tmp.bcf

    # Step 2: Apply a single, comprehensive soft filter.
    # Failing variants are marked with 'FAIL' but not removed yet.
    bcftools filter \\
        --threads ${task.cpus} \\
        -s 'FAIL' \\
        -e '${filter_expression}' \\
        -Ob \\
        -o soft_filtered.bcf \\
        tmp.bcf

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
    
    # Create a temporary VCF with only PASSing variants
    bcftools view --threads ${task.cpus} -f 'PASS' -Ob -o tmp.pass.bcf $bcf


    # Count the number of variant records (excluding the header)
    N_VARIANTS=\$(bcftools view -H tmp.pass.bcf | wc -l) 



    # Only create final output files if variants exist
    if [ "\$N_VARIANTS" -gt 0 ]; then
        # Rename temp file to final output
        mv tmp.pass.bcf final.pass.bcf
        bcftools index final.pass.bcf

        # Create the compressed VCF and its index
        bcftools view -Oz -o final.pass.vcf.gz final.pass.bcf
        bcftools index -t final.pass.vcf.gz

        # Get the final list of samples
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
    readlink -f $vcf > vcf_infile
    nsamples=\$(cat "$samples" | wc -l)
    amb_samples=\$(awk "BEGIN {print (\$nsamples/100)*10}" | awk '{ print int(\$1) }')

    python $projectDir/templates/broad-fungalgroup/scripts/SNPs/vcfSnpsToFasta.py --max_amb_samples \$amb_samples vcf_infile > fpg_snp_aln.fa
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
    ch_bam // Expects channel: [ [id], bam ]
    ch_ref // Expects channel: [ ref ]

    main:
    // 1. Call variants per sample
    CALLVARIANTS(ch_bam, ch_ref)

    // 2. Initial hard filter on each sample's BCF
    FILTERSAMPLE(CALLVARIANTS.out.bcf.join(CALLVARIANTS.out.bcf_idx))

    // 3. Collect and merge all sample BCFs
    smpl_filt_ch = FILTERSAMPLE.out.smpl_filt.collect()
    smpl_idx_ch  = FILTERSAMPLE.out.smpl_idx.collect()
    BCFMERGE(smpl_filt_ch, smpl_idx_ch)

    // 4. Reheader to clean sample names
    REHEADERVCF(BCFMERGE.out.mge, BCFMERGE.out.mge_idx)

    // 5. Normalize variants (split multiallelic sites)
    BCFNORM(REHEADERVCF.out.reh_bcf, REHEADERVCF.out.reh_idx, ch_ref)

    // 6. Apply comprehensive soft filters
    SOFTFILTERVCF(BCFNORM.out.bcf_norm, BCFNORM.out.idx_norm)

    // 7. Extract PASSing variants into final BCF/VCF files
    FILTERVCF(SOFTFILTERVCF.out.bcf_filt, SOFTFILTERVCF.out.idx_filt)

    // 8. Convert final VCF to other formats
    VCFSNPS2FASTA(FILTERVCF.out.vcf_pass, FILTERVCF.out.sample_list)
    VCF2PHYLIP(FILTERVCF.out.vcf_pass)

    emit:
    pass_vcf  = FILTERVCF.out.vcf_pass
    msa_snp   = VCFSNPS2FASTA.out.snp_msa
    aln_snp   = VCF2PHYLIP.out.snp_aln
    aln_fold  = VCF2PHYLIP.out.fold_aln
}
