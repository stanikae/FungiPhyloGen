#!/usr/bin/env python3
"""Compare two tab-delimited SNP distance matrices by taxon name."""

from pathlib import Path
import argparse
import sys


def read_matrix(path: Path):
    rows = [
        line.rstrip("\n").split("\t")
        for line in path.open(encoding="utf-8")
        if line.strip()
    ]
    if len(rows) < 2:
        raise ValueError(f"Matrix is empty or incomplete: {path}")

    taxa = rows[0][1:]
    matrix = {}
    for row in rows[1:]:
        if len(row) != len(taxa) + 1:
            raise ValueError(
                f"Malformed row in {path}: expected {len(taxa)+1} columns, "
                f"found {len(row)}"
            )
        matrix[row[0]] = dict(zip(taxa, map(int, row[1:])))

    if set(matrix) != set(taxa):
        raise ValueError(
            f"Row and column taxa differ in {path}: "
            f"rows-only={sorted(set(matrix)-set(taxa))}, "
            f"columns-only={sorted(set(taxa)-set(matrix))}"
        )
    return taxa, matrix


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("left", type=Path)
    parser.add_argument("right", type=Path)
    parser.add_argument("--left-label", default="left")
    parser.add_argument("--right-label", default="right")
    args = parser.parse_args()

    left_taxa, left = read_matrix(args.left)
    right_taxa, right = read_matrix(args.right)

    left_set = set(left_taxa)
    right_set = set(right_taxa)
    print(f"{args.left_label} taxa: {len(left_set)}")
    print(f"{args.right_label} taxa: {len(right_set)}")
    print(f"{args.left_label}-only taxa: {sorted(left_set-right_set)}")
    print(f"{args.right_label}-only taxa: {sorted(right_set-left_set)}")

    if left_set != right_set:
        print("FAIL: taxon sets differ.", file=sys.stderr)
        return 1

    mismatches = []
    for a in left_taxa:
        for b in left_taxa:
            lv = left[a][b]
            rv = right[a][b]
            if lv != rv:
                mismatches.append((a, b, lv, rv))

    print(f"Matrix cells: {len(left_taxa) ** 2}")
    print(f"Distance mismatches: {len(mismatches)}")
    for a, b, lv, rv in mismatches[:20]:
        print(f"{a}\t{b}\t{args.left_label}={lv}\t{args.right_label}={rv}")

    if mismatches:
        print("DIFFERENCE: matrices are not equivalent.", file=sys.stderr)
        return 1

    print("PASS: matrices are exactly equivalent by taxon name.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
