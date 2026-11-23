#!/usr/bin/env python3
"""Regression tests for the FungiPhyloGen VCF-to-SNP alignment converter."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONVERTER = ROOT / "scripts" / "fpg_vcf_to_snp_alignment.py"
FIXTURES = ROOT / "validation" / "vcf_snp_alignment" / "fixtures"
EXPECTED = ROOT / "validation" / "vcf_snp_alignment" / "expected"


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def read_matrix(path: Path) -> tuple[list[str], dict[str, dict[str, int]]]:
    with path.open(encoding="utf-8") as handle:
        rows = [
            line.rstrip("\n").split("\t")
            for line in handle
            if line.strip()
        ]
    taxa = rows[0][1:]
    matrix = {
        row[0]: dict(zip(taxa, map(int, row[1:])))
        for row in rows[1:]
    }
    return taxa, matrix


class ConverterTests(unittest.TestCase):
    maxDiff = None

    def run_converter(
        self,
        fixture: str,
        prefix: Path,
        expect_success: bool = True,
        extra_args: list[str] | None = None,
    ) -> subprocess.CompletedProcess[str]:
        command = [
            sys.executable,
            str(CONVERTER),
            "--vcf",
            str(FIXTURES / fixture),
            "--output-prefix",
            str(prefix),
        ]
        if extra_args:
            command.extend(extra_args)

        result = subprocess.run(
            command,
            text=True,
            capture_output=True,
            check=False,
        )

        if expect_success and result.returncode != 0:
            self.fail(
                f"converter failed\nSTDOUT:\n{result.stdout}"
                f"\nSTDERR:\n{result.stderr}"
            )
        if not expect_success and result.returncode == 0:
            self.fail("converter unexpectedly succeeded")
        return result

    def test_split_site_reconstruction_and_distances(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "fpg_snp"
            self.run_converter("split_site_test.vcf", prefix)

            self.assertEqual(
                read_text(Path(f"{prefix}_aln.fa")),
                read_text(EXPECTED / "split_site_expected.fa"),
            )
            self.assertEqual(
                read_text(Path(f"{prefix}_distances.tsv")),
                read_text(EXPECTED / "split_site_expected_distances.tsv"),
            )
            self.assertEqual(
                read_text(Path(f"{prefix}_comparable_sites.tsv")),
                read_text(
                    EXPECTED /
                    "split_site_expected_comparable_sites.tsv"
                ),
            )
            sites = read_text(Path(f"{prefix}_sites.tsv"))
            self.assertIn("2\tchr1\t200\tA\tC,T\t2\t0\n", sites)

    def test_missing_call_is_retained_as_n(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "fpg_snp"
            self.run_converter("missing_call_test.vcf", prefix)

            fasta = read_text(Path(f"{prefix}_aln.fa"))
            self.assertIn(">S3\nNTG\n", fasta)

            qc = read_text(Path(f"{prefix}_alignment_qc.tsv"))
            self.assertIn("alignment_columns\t3\n", qc)
            self.assertIn("ambiguous_bases[S3]\t1\n", qc)

            taxa, matrix = read_matrix(
                Path(f"{prefix}_comparable_sites.tsv")
            )
            self.assertEqual(
                taxa,
                ["reference", "S1", "S2", "S3", "S4"],
            )
            self.assertEqual(matrix["reference"]["S3"], 2)
            self.assertEqual(matrix["S2"]["S3"], 2)
            self.assertEqual(matrix["reference"]["S2"], 3)

    def test_non_haploid_genotype_fails_without_partial_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "fpg_snp"
            result = self.run_converter(
                "haploid_test.vcf",
                prefix,
                expect_success=False,
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("non-haploid genotype '0/1'", result.stderr)
            self.assertFalse(Path(f"{prefix}_aln.fa").exists())
            self.assertFalse(Path(f"{prefix}_distances.tsv").exists())

    def test_no_secondary_filter_options(self):
        result = subprocess.run(
            [sys.executable, str(CONVERTER), "--help"],
            text=True,
            capture_output=True,
            check=True,
        )
        for forbidden in (
            "max_amb",
            "QUAL",
            "MQ",
            "FS",
            "INFO/DP",
            "INFO/AF",
            "INFO/AC",
        ):
            self.assertNotIn(forbidden, result.stdout)

    def test_manifest_records_version_commit_and_checksums(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "fpg_snp"
            self.run_converter(
                "split_site_test.vcf",
                prefix,
                extra_args=[
                    "--pipeline-version-file",
                    str(ROOT / "VERSION"),
                    "--git-commit",
                    "test-commit",
                ],
            )

            manifest_path = Path(f"{prefix}_manifest.json")
            self.assertTrue(manifest_path.is_file())
            manifest = json.loads(
                manifest_path.read_text(encoding="utf-8")
            )

            self.assertEqual(
                manifest["fungiphylogen"]["version"],
                "1.2.0",
            )
            self.assertEqual(
                manifest["fungiphylogen"]["git_commit"],
                "test-commit",
            )
            self.assertEqual(
                manifest["converter"]["version"],
                "1.2.0",
            )
            self.assertFalse(
                manifest["parameters"][
                    "secondary_biological_filtering"
                ]
            )
            self.assertEqual(manifest["alignment"]["columns"], 2)
            self.assertEqual(
                manifest["alignment"]["source_vcf_records"],
                3,
            )
            self.assertEqual(
                manifest["alignment"][
                    "split_multiallelic_positions"
                ],
                1,
            )
            self.assertEqual(len(manifest["input"]["sha256"]), 64)
            self.assertEqual(len(manifest["outputs"]), 5)

            for output in manifest["outputs"]:
                self.assertEqual(len(output["sha256"]), 64)


if __name__ == "__main__":
    unittest.main(verbosity=2)
