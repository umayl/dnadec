#!/usr/bin/env python3
"""Evaluate planted synthetic circuits; this is not a biological accuracy benchmark."""

import json
import re
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROMOTER = "TTGACAAAAAAAAAAAAAAAAAATATAAT"
RBS = "AGGAGGAAAAAAA"
TERMINATOR = "GCGCGCGCGAAAGCGCGCGCTTTTTT"
OPERATOR = "GATATATATATC"
START_STOP = re.compile(r"@span:\s*(\d+)\.\.(\d+)\s+\(([+-])\)")
SITE = re.compile(
    r"CANDIDATE_OPERATOR_SITE gene_\d+ motif=([ACGT]+) "
    r"span=(\d+)\.\.(\d+) strand=([+-])"
)
MASK_64 = (1 << 64) - 1


def random_dna(length, seed):
    bases = "ACGT"
    sequence = []
    state = seed
    for _ in range(length):
        state = (state * 6364136223846793005 + 1442695040888963407) & MASK_64
        sequence.append(bases[(state >> 33) % 4])
    return "".join(sequence)


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def make_positive(index):
    seed = 20261125 + index * 7919
    prefix = random_dna(24, seed)
    repeat_count = 1 + index % 3
    operator_block = OPERATOR * repeat_count
    codon = ("GCT", "GCC", "GAA")[index % 3]
    body = codon * (40 + index % 41)
    oriented = prefix + PROMOTER + operator_block + RBS + "ATG" + body + "TAA" + TERMINATOR
    start = len(prefix) + len(PROMOTER) + len(operator_block) + len(RBS)
    end = start + 3 + len(body) + 3
    expected_sites = [
        (OPERATOR, len(prefix) + len(PROMOTER) + offset, "+")
        for offset in range(0, len(operator_block), len(OPERATOR))
    ]

    if index % 2:
        sequence = reverse_complement(oriented)
        expected_gene = (len(oriented) - end, len(oriented) - start, "-")
        expected_sites = [
            (motif, len(oriented) - position - len(motif), "-")
            for motif, position, _ in expected_sites
        ]
    else:
        sequence = oriented
        expected_gene = (start, end, "+")

    return sequence, expected_gene, expected_sites


def make_negative(index):
    seed = 911_000 + index * 1013
    prefix = "C" * 24
    body = "GCT" * 55
    if index % 4 == 0:
        return random_dna(600, seed)
    if index % 4 == 1:
        return prefix + PROMOTER + "C" * 30 + "ATG" + body + "TAA" + TERMINATOR
    if index % 4 == 2:
        return prefix + PROMOTER + RBS + "ATG" + body + "C" * 60
    return prefix + PROMOTER + RBS + "ATG" + "GCT" * 16 + "TAA" + TERMINATOR


def f1(precision, recall):
    return 2 * precision * recall / (precision + recall) if precision + recall else 0.0


def main():
    subprocess.run(
        ["cargo", "build", "--release", "--bin", "dnadec"],
        cwd=ROOT,
        check=True,
    )
    binary = ROOT / "target" / "release" / "dnadec.exe"
    if not binary.exists():
        binary = ROOT / "target" / "release" / "dnadec"
    if not binary.exists():
        raise RuntimeError(f"Release binary not found: {binary}")

    positives = [make_positive(index) for index in range(40)]
    negatives = [make_negative(index) for index in range(40)]
    predicted_gene_count = 0
    exact_gene_matches = 0
    predicted_sites = set()
    expected_sites = set()

    with tempfile.TemporaryDirectory(prefix="dnadec-synthetic-") as directory:
        temp = Path(directory)
        cases = [
            (f"positive_{index}", sequence, expected_gene, sites)
            for index, (sequence, expected_gene, sites) in enumerate(positives)
        ]
        cases.extend(
            (f"negative_{index}", sequence, None, [])
            for index, sequence in enumerate(negatives)
        )

        for name, sequence, expected_gene, sites in cases:
            fasta = temp / f"{name}.fasta"
            fasta.write_text(f">{name}\n{sequence}\n", encoding="ascii")
            result = subprocess.run(
                [str(binary), str(fasta), "--non-strict"],
                cwd=ROOT,
                text=True,
                capture_output=True,
                timeout=60,
                check=False,
            )
            if result.returncode != 0:
                raise RuntimeError(
                    f"dnadec failed on {name}:\n{result.stderr[-4000:]}"
                )

            genes = {
                (int(start), int(end), strand)
                for start, end, strand in START_STOP.findall(result.stdout)
            }
            predicted_gene_count += len(genes)
            if expected_gene is not None and expected_gene in genes:
                exact_gene_matches += 1
            predicted_sites.update(
                (name, motif, int(start), strand)
                for motif, start, _, strand in SITE.findall(result.stdout)
            )
            expected_sites.update(
                (name, motif, start, strand)
                for motif, start, strand in sites
            )

    false_positive_genes = predicted_gene_count - exact_gene_matches
    gene_precision = exact_gene_matches / predicted_gene_count if predicted_gene_count else 0.0
    gene_recall = exact_gene_matches / len(positives)
    matched_sites = len(predicted_sites & expected_sites)
    site_precision = matched_sites / len(predicted_sites) if predicted_sites else 0.0
    site_recall = matched_sites / len(expected_sites) if expected_sites else 0.0
    print(
        json.dumps(
            {
                "scope": "held-out deterministic synthetic constructs; not biological validation",
                "positive_constructs": len(positives),
                "negative_constructs": len(negatives),
                "predicted_genes": predicted_gene_count,
                "exact_gene_matches": exact_gene_matches,
                "false_positive_genes": false_positive_genes,
                "missed_genes": len(positives) - exact_gene_matches,
                "gene_precision": gene_precision,
                "gene_recall": gene_recall,
                "gene_f1": f1(gene_precision, gene_recall),
                "predicted_operator_sites": len(predicted_sites),
                "expected_operator_sites": len(expected_sites),
                "matched_operator_sites": matched_sites,
                "missed_operator_sites": len(expected_sites) - matched_sites,
                "operator_site_precision": site_precision,
                "operator_site_recall": site_recall,
                "operator_site_f1": f1(site_precision, site_recall),
                "temporary_files_deleted": not temp.exists(),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
