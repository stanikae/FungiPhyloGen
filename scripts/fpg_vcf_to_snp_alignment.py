#!/usr/bin/env python3
"""Construct a canonical haploid SNP alignment from a FungiPhyloGen PASS VCF.

The input is expected to be the normalized, multi-sample, SNP-only
final.pass.vcf.gz produced by FungiPhyloGen.

This program does not apply biological variant-quality filters. It translates
the authoritative PASS VCF into:

- a canonical SNP FASTA alignment;
- a retained-site manifest;
- alignment QC;
- a raw SNP-difference matrix;
- a pairwise comparable-site matrix.

Normalized records sharing CHROM, POS and REF are reconstructed as one
biological nucleotide position.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from collections import Counter
from pathlib import Path
from typing import Iterator, TextIO

CANONICAL_BASES = frozenset("ACGT")
MISSING_GENOTYPES = frozenset({".", "./.", ".|."})
CONVERTER_VERSION = "1.2.0"


class AlignmentError(RuntimeError):
    """Raised when the input violates the FungiPhyloGen VCF contract."""


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Construct a canonical haploid SNP alignment and pairwise SNP "
            "distance matrices from a FungiPhyloGen final PASS VCF."
        )
    )
    parser.add_argument(
        "--vcf", required=True, type=Path,
        help="Input final.pass.vcf.gz or uncompressed VCF."
    )
    parser.add_argument(
        "--output-prefix", required=True, type=Path,
        help="Output prefix, for example results/fpg_snp."
    )
    parser.add_argument(
        "--line-width", type=int, default=80,
        help="FASTA line width. Default: 80."
    )
    parser.add_argument(
        "--pipeline-version-file", type=Path,
        help="File containing the FungiPhyloGen version, normally VERSION."
    )
    parser.add_argument(
        "--git-commit", default="unknown",
        help="FungiPhyloGen Git commit tested by the workflow."
    )
    return parser.parse_args()


def open_text(path: Path) -> TextIO:
    if not path.is_file():
        raise AlignmentError(f"Input VCF does not exist: {path}")
    if path.stat().st_size == 0:
        raise AlignmentError(f"Input VCF is empty: {path}")
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def parse_header(path: Path) -> tuple[list[str], list[str]]:
    samples: list[str] = []
    contigs: list[str] = []
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("##contig=<ID="):
                token = line.split("##contig=<ID=", 1)[1]
                contig = token.split(",", 1)[0].split(">", 1)[0]
                contigs.append(contig)
            elif line.startswith("#CHROM"):
                fields = line.split("\t")
                if len(fields) < 10:
                    raise AlignmentError("The input VCF contains no sample columns.")
                samples = fields[9:]
                break
    if not samples:
        raise AlignmentError("The VCF #CHROM header was not found.")
    duplicates = sorted(name for name, count in Counter(samples).items() if count > 1)
    if duplicates:
        raise AlignmentError("Duplicate sample identifiers: " + ", ".join(duplicates))
    if "reference" in samples:
        raise AlignmentError("Sample name 'reference' is reserved.")
    return samples, contigs


def iter_records(path: Path, sample_count: int) -> Iterator[dict]:
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            expected = 9 + sample_count
            if len(fields) < expected:
                raise AlignmentError(
                    f"Line {line_number}: expected at least {expected} fields, found {len(fields)}."
                )
            chrom, pos_text, _rid, ref, alt_text, qual, filt, info, fmt = fields[:9]
            sample_fields = fields[9:expected]
            try:
                pos = int(pos_text)
            except ValueError as exc:
                raise AlignmentError(
                    f"Line {line_number}: invalid POS value {pos_text!r}."
                ) from exc
            if filt != "PASS":
                raise AlignmentError(
                    f"{chrom}:{pos}: expected FILTER=PASS, found {filt!r}."
                )
            ref = ref.upper()
            alts = [allele.upper() for allele in alt_text.split(",")]
            if len(ref) != 1 or ref not in CANONICAL_BASES:
                raise AlignmentError(
                    f"{chrom}:{pos}: REF must be one of A,C,G,T; found {ref!r}."
                )
            if not alts or any(len(a) != 1 or a not in CANONICAL_BASES for a in alts):
                raise AlignmentError(
                    f"{chrom}:{pos}: every ALT must be one of A,C,G,T; found {alt_text!r}."
                )
            if ref in alts:
                raise AlignmentError(f"{chrom}:{pos}: REF also occurs in ALT.")
            if len(set(alts)) != len(alts):
                raise AlignmentError(f"{chrom}:{pos}: duplicate ALT alleles were found.")
            format_keys = fmt.split(":")
            if "GT" not in format_keys:
                raise AlignmentError(f"{chrom}:{pos}: FORMAT does not contain GT.")
            gt_index = format_keys.index("GT")
            genotypes: list[str] = []
            for sample_number, sample_field in enumerate(sample_fields, start=1):
                values = sample_field.split(":")
                if gt_index >= len(values):
                    raise AlignmentError(
                        f"{chrom}:{pos}: sample column {sample_number} has no GT value."
                    )
                genotypes.append(values[gt_index])
            yield {
                "chrom": chrom, "pos": pos, "ref": ref, "alts": alts,
                "qual": qual, "info": info, "genotypes": genotypes,
                "line_number": line_number,
            }


def decode_haploid_record(record: dict, samples: list[str]) -> list[str]:
    alleles = [record["ref"], *record["alts"]]
    decoded: list[str] = []
    for sample_name, genotype in zip(samples, record["genotypes"]):
        if genotype in MISSING_GENOTYPES:
            decoded.append("N")
            continue
        if "/" in genotype or "|" in genotype:
            raise AlignmentError(
                f"{record['chrom']}:{record['pos']}: sample {sample_name!r} "
                f"has non-haploid genotype {genotype!r}."
            )
        try:
            allele_index = int(genotype)
        except ValueError as exc:
            raise AlignmentError(
                f"{record['chrom']}:{record['pos']}: sample {sample_name!r} "
                f"has invalid genotype {genotype!r}."
            ) from exc
        if allele_index < 0 or allele_index >= len(alleles):
            raise AlignmentError(
                f"{record['chrom']}:{record['pos']}: sample {sample_name!r} "
                f"has allele index {allele_index}, but {len(alleles)} alleles are defined."
            )
        decoded.append(alleles[allele_index])
    return decoded


def reconstruct_position(records: list[dict], samples: list[str]) -> tuple[str, list[str], list[str]]:
    first = records[0]
    chrom, pos, ref = first["chrom"], first["pos"], first["ref"]
    seen_alts: set[str] = set()
    observed_alts: list[str] = []
    sample_alleles = [ref] * len(samples)
    sample_missing = [False] * len(samples)
    sample_alt_seen = [False] * len(samples)

    for record in records:
        if record["chrom"] != chrom or record["pos"] != pos:
            raise AlignmentError("Internal position-grouping error.")
        if record["ref"] != ref:
            raise AlignmentError(
                f"{chrom}:{pos}: split records have conflicting REF alleles."
            )
        for alt in record["alts"]:
            if alt in seen_alts:
                raise AlignmentError(
                    f"{chrom}:{pos}: ALT {alt!r} occurs in more than one split record."
                )
            seen_alts.add(alt)
            observed_alts.append(alt)

        decoded = decode_haploid_record(record, samples)
        for index, allele in enumerate(decoded):
            if allele == "N":
                sample_missing[index] = True
            elif allele != ref:
                if sample_alt_seen[index] and sample_alleles[index] != allele:
                    raise AlignmentError(
                        f"{chrom}:{pos}: sample {samples[index]!r} carries conflicting "
                        f"alternate alleles {sample_alleles[index]!r} and {allele!r}."
                    )
                sample_alleles[index] = allele
                sample_alt_seen[index] = True

    for index, is_missing in enumerate(sample_missing):
        if is_missing:
            sample_alleles[index] = "N"

    if not any(allele not in {ref, "N"} for allele in sample_alleles):
        raise AlignmentError(
            f"{chrom}:{pos}: PASS record group is invariant after reconstruction."
        )
    return ref, sample_alleles, observed_alts


def construct_alignment(path: Path, samples: list[str]) -> tuple[list[str], dict[str, list[str]], list[dict]]:
    sequences: dict[str, list[str]] = {
        "reference": [], **{sample: [] for sample in samples}
    }
    retained_sites: list[dict] = []
    current_key: tuple[str, int] | None = None
    current_records: list[dict] = []
    seen_keys: set[tuple[str, int]] = set()

    def flush_group() -> None:
        if not current_records:
            return
        ref, sample_alleles, observed_alts = reconstruct_position(current_records, samples)
        first = current_records[0]
        sequences["reference"].append(ref)
        for sample, allele in zip(samples, sample_alleles):
            sequences[sample].append(allele)
        retained_sites.append({
            "alignment_column": len(retained_sites) + 1,
            "chrom": first["chrom"], "pos": first["pos"], "ref": ref,
            "observed_alt_alleles": ",".join(observed_alts),
            "source_record_count": len(current_records),
            "ambiguous_sample_count": sum(a == "N" for a in sample_alleles),
        })

    for record in iter_records(path, len(samples)):
        key = (record["chrom"], record["pos"])
        if current_key is None:
            current_key = key
            seen_keys.add(key)
        if key != current_key:
            flush_group()
            if key in seen_keys:
                raise AlignmentError(
                    f"{record['chrom']}:{record['pos']}: records for the same position "
                    "are not contiguous; VCF is not correctly sorted."
                )
            seen_keys.add(key)
            current_key = key
            current_records = []
        current_records.append(record)
    flush_group()

    if not retained_sites:
        raise AlignmentError("No variable PASS SNP positions were found.")
    expected_length = len(retained_sites)
    for taxon, sequence in sequences.items():
        if len(sequence) != expected_length:
            raise AlignmentError(
                f"Internal alignment error: {taxon!r} has length {len(sequence)}, "
                f"expected {expected_length}."
            )
    return ["reference", *samples], sequences, retained_sites


def write_fasta(path: Path, taxa: list[str], sequences: dict[str, list[str]], line_width: int) -> None:
    if line_width < 1:
        raise AlignmentError("--line-width must be at least 1.")
    with path.open("w", encoding="utf-8") as handle:
        for taxon in taxa:
            sequence = "".join(sequences[taxon])
            handle.write(f">{taxon}\n")
            for start in range(0, len(sequence), line_width):
                handle.write(sequence[start:start + line_width] + "\n")


def write_sites(path: Path, retained_sites: list[dict]) -> None:
    columns = [
        "alignment_column", "chrom", "pos", "ref", "observed_alt_alleles",
        "source_record_count", "ambiguous_sample_count",
    ]
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for site in retained_sites:
            handle.write("\t".join(str(site[column]) for column in columns) + "\n")


def pairwise_matrices(taxa: list[str], sequences: dict[str, list[str]]) -> tuple[list[list[int]], list[list[int]]]:
    strings = {taxon: "".join(sequences[taxon]) for taxon in taxa}
    size = len(taxa)
    differences = [[0] * size for _ in range(size)]
    comparable = [[0] * size for _ in range(size)]
    for left_index in range(size):
        left = strings[taxa[left_index]]
        for right_index in range(left_index, size):
            right = strings[taxa[right_index]]
            difference_count = 0
            comparable_count = 0
            for left_base, right_base in zip(left, right):
                if "N" in (left_base, right_base):
                    continue
                comparable_count += 1
                if left_base != right_base:
                    difference_count += 1
            differences[left_index][right_index] = difference_count
            differences[right_index][left_index] = difference_count
            comparable[left_index][right_index] = comparable_count
            comparable[right_index][left_index] = comparable_count
    return differences, comparable


def write_matrix(path: Path, taxa: list[str], matrix: list[list[int]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("Taxon\t" + "\t".join(taxa) + "\n")
        for taxon, row in zip(taxa, matrix):
            handle.write(taxon + "\t" + "\t".join(map(str, row)) + "\n")


def write_qc(path: Path, input_vcf: Path, taxa: list[str], sequences: dict[str, list[str]], retained_sites: list[dict]) -> None:
    unique_positions = len(retained_sites)
    source_records = sum(site["source_record_count"] for site in retained_sites)
    split_positions = sum(site["source_record_count"] > 1 for site in retained_sites)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"input_vcf\t{input_vcf}\n")
        handle.write(f"sample_count\t{len(taxa) - 1}\n")
        handle.write(f"alignment_taxa\t{len(taxa)}\n")
        handle.write(f"alignment_columns\t{unique_positions}\n")
        handle.write(f"source_vcf_records\t{source_records}\n")
        handle.write(f"split_multiallelic_positions\t{split_positions}\n")
        for taxon in taxa:
            ambiguous = sum(base == "N" for base in sequences[taxon])
            handle.write(f"ambiguous_bases[{taxon}]\t{ambiguous}\n")



def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_pipeline_version(path: Path | None) -> str:
    if path is None:
        return "unknown"
    if not path.is_file():
        raise AlignmentError(f"Pipeline version file does not exist: {path}")
    version = path.read_text(encoding="utf-8").strip()
    if not version:
        raise AlignmentError(f"Pipeline version file is empty: {path}")
    return version


def write_manifest(
    path: Path,
    args: argparse.Namespace,
    pipeline_version: str,
    taxa: list[str],
    retained_sites: list[dict],
    output_paths: list[Path],
) -> None:
    source_records = sum(site["source_record_count"] for site in retained_sites)
    split_positions = sum(site["source_record_count"] > 1 for site in retained_sites)
    manifest = {
        "schema_version": "1.0",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "fungiphylogen": {
            "version": pipeline_version,
            "git_commit": args.git_commit,
        },
        "converter": {
            "name": "fpg_vcf_to_snp_alignment.py",
            "version": CONVERTER_VERSION,
        },
        "runtime": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "platform": platform.platform(),
        },
        "command": [sys.executable, *sys.argv],
        "input": {
            "path": str(args.vcf),
            "size_bytes": args.vcf.stat().st_size,
            "sha256": sha256_file(args.vcf),
        },
        "parameters": {
            "line_width": args.line_width,
            "secondary_biological_filtering": False,
            "missing_genotype_character": "N",
        },
        "alignment": {
            "sample_count": len(taxa) - 1,
            "taxon_count_including_reference": len(taxa),
            "columns": len(retained_sites),
            "source_vcf_records": source_records,
            "split_multiallelic_positions": split_positions,
        },
        "outputs": [
            {
                "path": str(output),
                "size_bytes": output.stat().st_size,
                "sha256": sha256_file(output),
            }
            for output in output_paths
        ],
    }
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

def main() -> int:
    args = parse_arguments()
    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    samples, _contigs = parse_header(args.vcf)
    taxa, sequences, retained_sites = construct_alignment(args.vcf, samples)

    fasta_path = Path(f"{prefix}_aln.fa")
    sites_path = Path(f"{prefix}_sites.tsv")
    qc_path = Path(f"{prefix}_alignment_qc.tsv")
    distances_path = Path(f"{prefix}_distances.tsv")
    comparable_path = Path(f"{prefix}_comparable_sites.tsv")
    manifest_path = Path(f"{prefix}_manifest.json")
    pipeline_version = read_pipeline_version(args.pipeline_version_file)

    write_fasta(fasta_path, taxa, sequences, args.line_width)
    write_sites(sites_path, retained_sites)
    write_qc(qc_path, args.vcf, taxa, sequences, retained_sites)
    differences, comparable = pairwise_matrices(taxa, sequences)
    write_matrix(distances_path, taxa, differences)
    write_matrix(comparable_path, taxa, comparable)
    write_manifest(
        manifest_path,
        args,
        pipeline_version,
        taxa,
        retained_sites,
        [fasta_path, sites_path, qc_path, distances_path, comparable_path],
    )

    print(f"Alignment: {fasta_path}")
    print(f"Sites: {sites_path}")
    print(f"QC: {qc_path}")
    print(f"SNP distances: {distances_path}")
    print(f"Comparable sites: {comparable_path}")
    print(f"Manifest: {manifest_path}")
    print(f"Samples: {len(samples)}")
    print(f"Alignment columns: {len(retained_sites)}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AlignmentError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
