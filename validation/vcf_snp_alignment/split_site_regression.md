# Split-site regression test

## Purpose

Verify that normalized biallelic records sharing the same chromosome,
position and reference allele are reconstructed as one biological nucleotide
position.

## Fixture

The fixture contains:

- one ordinary SNP at `chr1:100`;
- two normalized records at `chr1:200` representing `A>C` and `A>T`.

At `chr1:200`, samples S1-S3 carry C and sample S4 carries T.

## Required behaviour

The converter must:

1. group both `chr1:200` records into one genomic position;
2. reconstruct one haploid nucleotide per sample;
3. emit one alignment column for `chr1:200`;
4. preserve the alternate alleles C and T;
5. produce a two-column alignment overall;
6. count one difference between C and T at `chr1:200`.

## Confirmed current behaviour

The Broad-derived converter produced the expected canonical alignment:

- reference: `GA`
- S1: `GC`
- S2: `AC`
- S3: `AC`
- S4: `GT`

The current and canonical SNP-distance matrices were identical.

The split-site reconstruction behaviour therefore passes and must be
preserved by any replacement implementation.

## Alignment filtering policy

The replacement converter will not implement an ambiguous-sample threshold.
All biological variant filtering is completed upstream by the bcftools
workflow. Explicit missing genotypes will be represented as `N`, retained in
the alignment and reported through alignment and pairwise-comparability QC.
