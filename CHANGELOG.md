# Changelog

All notable changes to FungiPhyloGen are documented here.

## [Unreleased]

### Fixed

- Revised haploid population-SNP filtering to remove the cohort-depth-dependent `QUAL/INFO/DP` proxy, which was not equivalent to GATK QD after independently called sample BCFs were merged.
- Materialized cohort `AC`, `AN`, and `AF` annotations after normalization using `bcftools +fill-tags`, making population allele-frequency filtering explicit and auditable.
- Added `cauris` as the canonical *Candida auris* population-SNP profile while retaining `cauris_small` and `cauris_batch` as backward-compatible aliases.
- Removed profile-level `F_MISSING` filtering pending correction and validation of missing-genotype handling with `bcftools merge -0`.
- Retained the established species-specific population allele-frequency thresholds for haploid SNP selection.
- Corrected the `cparapsilosis` SNP-filtering profile by replacing the ambiguous unqualified `DP<10` expression with `INFO/DP<10`, removing the depth-dependent `(QUAL/MAX(AD[:1]))<=2.0` criterion, and adding an explicit `QUAL<30` exclusion threshold. This prevents `bcftools filter` failures when both `INFO/DP` and `FORMAT/DP` are defined and avoids rejecting well-supported variants solely because of high alternate-read depth.

### Changed

- Large workflow intermediates, including trimmed reads, alignment files, per-sample BCFs, merged BCFs, and reheadered BCFs, are now published as symbolic links to avoid unnecessary duplication of Nextflow work files.
- Corrected `samtools sort` and `samtools index` argument ordering for compatibility with current Samtools releases.
- Made SNP-distance output cleanup compatible with both BSD `sed` on macOS and GNU `sed` on Linux.
