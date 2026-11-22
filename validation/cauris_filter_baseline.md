# C. auris filter-validation baseline

## Repository baseline

- Production commit: 38d4380509bbb7c0f41266917623fd5dcfea293a
- Validation branch: fix/cauris-filter-validation
- Historical development revision: 78bda17ecb
- Validation date: 2026-11-20

## Validated ERG11 feature

- Reference assembly: B8441 V3
- ERG11 locus: B9J08_03698
- Previous V2 locus: B9J08_001448
- Validated site: CM076440.1:1475923 A>G
- Expected consequence: c.374T>C, p.Val125Ala

## Historical run finding

The ERG11 variant was present through merging, reheadering, normalization and soft filtering, but was labelled LowQD and removed during PASS-only selection.

- QUAL: 284.59
- Cohort INFO/DP: 9125
- QUAL/INFO/DP: 0.0311233
- QD cutoff: 0.035
- Final PASS VCF: variant absent
- SnpEff and targeted report: variant unavailable

## Current production behaviour

The production workflow resolves and applies filtering as follows:

    params.filter_profile
    -> params.filters[params.filter_profile]
    -> params.active_filter
    -> bcftools filter -e
    -> FAIL assignment
    -> bcftools view -f PASS

Complementary tests confirmed that bcftools correctly evaluates INFO/AF and F_MISSING from genotype data.

Both current cauris_small and cauris_batch profiles classify the validated ERG11 site as FAIL because both retain the cohort-level QUAL/INFO/DP threshold.

## Required end-to-end validation

Before changing the production profiles, run independent end-to-end validation for:

- 20 samples using cauris_small
- 30 samples as the proposed profile boundary
- 82 samples using cauris_batch

Each run must begin from aligned BAMs or per-sample raw BCFs and continue through calling, merging, normalization, production filtering, annotation, targeted reporting, SNP distances and phylogenetic reconstruction.
