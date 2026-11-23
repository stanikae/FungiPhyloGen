# VCF-to-SNP alignment compatibility audit

## Scope

This audit evaluates the Broad-derived `vcfSnpsToFasta.py` converter against
the normalized, multi-sample, haploid `final.pass.vcf.gz` produced by
FungiPhyloGen's bcftools workflow.

## Confirmed baseline behaviour

The controlled fixture contains:

- an explicit haploid reference genotype;
- an explicit haploid alternate genotype;
- an explicit missing genotype;
- a multiallelic SNP with haploid allele indices 1 and 2;
- a heterozygous genotype in a haploid callset;
- an invariant reference-only site.

The existing converter:

- represents explicit missing genotypes as `N`;
- represents heterozygous genotypes as `N`;
- decodes haploid multiallelic allele indices 1 and 2;
- excludes invariant reference-only sites;
- produces equal-length sequences.

## Confirmed defect

`--max_amb_samples 0` produces the same alignment as
`--max_amb_samples 1`.

The zero threshold is not applied because the script tests the integer option
by truth value. Zero is false, so the default ambiguity allowance remains in
effect.

For FungiPhyloGen cohorts containing fewer than 10 samples, the workflow
calculates a 10% ambiguity allowance of zero. The current implementation
therefore does not enforce the intended ambiguity threshold for those
cohorts.

## Replacement requirements

The FungiPhyloGen-specific converter must:

1. consume `final.pass.vcf.gz` directly;
2. preserve all sites from the final PASS VCF without applying an additional ambiguity threshold;
3. preserve explicit `GT=0` as REF;
4. convert explicit missing genotypes to `N`;
5. handle unexpected heterozygous calls explicitly;
6. validate haploid allele indices;
7. reject duplicate sample identifiers;
8. reject malformed and duplicate sites;
9. emit deterministic, equal-length SNP-only sequences;
10. report retained and excluded sites;
11. report ambiguity per sample;
12. produce pairwise comparable-site counts alongside SNP differences.

## Production nine-isolate input audit

The validated nine-isolate `final.pass.vcf.gz` contains:

- 22,750 PASS SNP records;
- 103,524 haploid alternate genotypes (`GT=1`);
- 101,226 haploid reference genotypes (`GT=0`);
- no explicit missing genotypes;
- no diploid or phased genotypes;
- no invalid allele indices;
- 641 duplicated `CHROM:POS` keys.

All 641 duplicated keys contain two normalized biallelic records representing
one biological multiallelic nucleotide position. The normalized VCF therefore
contains 22,109 unique genomic positions.

A controlled split-site regression fixture confirmed that the current
Broad-derived converter correctly reconstructs normalized records sharing the
same chromosome and position as one biological alignment column.

The current converter produced the expected canonical FASTA and pairwise SNP
distance matrix. Split multiallelic positions are therefore not double-counted
by the current implementation.

Any replacement converter must preserve this validated reconstruction
behaviour and emit one alignment column per unique genomic position.
