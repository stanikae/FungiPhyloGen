# Changelog

All notable changes to FungiPhyloGen are documented here.

## [Unreleased]

## [1.2.0] - 2026-11-23

### Added

- Added a caller-neutral Python 3.10 converter that constructs the canonical SNP alignment directly from `final.pass.vcf.gz` without applying secondary biological filtering.
- Added `fpg_snp_sites.tsv` for alignment-column provenance and split-multiallelic reconstruction metadata.
- Added `fpg_snp_alignment_qc.tsv` with alignment dimensions, source-record totals, reconstructed split positions, and per-taxon ambiguity counts.
- Added `fpg_snp_distances.tsv` and `fpg_snp_comparable_sites.tsv` to report raw pairwise SNP differences together with their comparable-site denominators.
- Added regression tests for split multiallelic reconstruction, explicit missing-genotype retention, comparable-site accounting, and rejection of non-haploid genotypes.


### Fixed

- Revised haploid population-SNP filtering to remove the cohort-depth-dependent `QUAL/INFO/DP` proxy, which was not equivalent to GATK QD after independently called sample BCFs were merged.
- Materialized cohort `AC`, `AN`, and `AF` annotations after normalization using `bcftools +fill-tags`, making population allele-frequency filtering explicit and auditable.
- Added `cauris` as the canonical *Candida auris* population-SNP profile while retaining `cauris_small` and `cauris_batch` as backward-compatible aliases.
- Removed profile-level `F_MISSING` filtering pending correction and validation of missing-genotype handling with `bcftools merge -0`.
- Retained the established species-specific population allele-frequency thresholds for haploid SNP selection.
- Corrected the `cparapsilosis` SNP-filtering profile by replacing the ambiguous unqualified `DP<10` expression with `INFO/DP<10`, removing the depth-dependent `(QUAL/MAX(AD[:1]))<=2.0` criterion, and adding an explicit `QUAL<30` exclusion threshold. This prevents `bcftools filter` failures when both `INFO/DP` and `FORMAT/DP` are defined and avoids rejecting well-supported variants solely because of high alternate-read depth.
- Validated the revised `cparapsilosis` profile using a targeted five-isolate audit and a clean nine-isolate end-to-end workflow run; soft-filter classifications reconciled exactly with the final PASS BCF in both cohorts.

### Changed

- Replaced the production invocation of the Broad-derived `vcfSnpsToFasta.py` workflow with direct processing of the normalized multi-sample PASS VCF.
- Removed downstream ambiguous-site filtering, temporary VCF decompression, `vcf_infile`, and `sample_list.txt`; all biological filtering remains authoritative in the upstream `bcftools` workflow.
- Preserved one canonical SNP alignment as the input to IQ-TREE and RapidNJ, while validating the generated SNP-distance matrix against `snp-dists` 0.8.2.

- Large workflow intermediates, including trimmed reads, alignment files, per-sample BCFs, merged BCFs, and reheadered BCFs, are now published as symbolic links to avoid unnecessary duplication of Nextflow work files.
- Corrected `samtools sort` and `samtools index` argument ordering for compatibility with current Samtools releases.
- Made SNP-distance output cleanup compatible with both BSD `sed` on macOS and GNU `sed` on Linux.

### Validation

- Reconstructed 22,750 PASS SNP records from a nine-isolate *Candida parapsilosis* integration dataset into 22,109 genomic alignment positions, including 641 split multiallelic positions and zero ambiguous bases.
- Confirmed exact agreement across all 100 matrix cells between the new converter, independent `snp-dists`, and the previous production alignment's SNP-distance matrix.
