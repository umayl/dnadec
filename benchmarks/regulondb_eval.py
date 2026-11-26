#!/usr/bin/env python3
"""Evaluate dnadec against official, noncommercial RegulonDB exports.

The downloaded reference data is held only in a temporary directory and removed
when the run ends. RegulonDB terms permit academic/noncommercial use; do not
redistribute the downloaded records.
"""

import bisect
import hashlib
import json
import math
import re
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from urllib.parse import urlsplit
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[1]
REGULONDB_HOST = "regulondb.ccg.unam.mx"
GRAPHQL_URL = f"https://{REGULONDB_HOST}/graphql"
GENOME_URL = f"https://{REGULONDB_HOST}/media/raw/e_coli_k12.fna"
MAX_DOWNLOAD_BYTES = 20 * 1024 * 1024
START_CODONS = {"ATG", "GTG", "TTG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}
SPAN_RE = re.compile(r"@span:\s*(\d+)\.\.(\d+)\s+\(([+-])\)")
REGULATORY_EDGE_RE = re.compile(r"REGULATORY_EDGE gene_(\d+) -\| gene_(\d+)")
OPERATOR_SITE_RE = re.compile(
    r"CANDIDATE_OPERATOR_SITE gene_\d+ motif=[ACGT]+ "
    r"span=(\d+)\.\.(\d+) strand=[+-]"
)
RELEASE_RE = re.compile(r"RegulonDB Release:\s*([^\r\n]+)")


def fetch_bytes(url, *, body=None, max_bytes=MAX_DOWNLOAD_BYTES):
    headers = {"User-Agent": "dnadec-benchmark/1.0"}
    if body is not None:
        headers["Content-Type"] = "application/json"
    request = Request(url, data=body, headers=headers)
    with urlopen(request, timeout=60) as response:
        if urlsplit(response.geturl()).hostname != REGULONDB_HOST:
            raise RuntimeError(f"Refusing unexpected redirect to {response.geturl()}")
        data = response.read(max_bytes + 1)
        if len(data) > max_bytes:
            raise RuntimeError(f"Download exceeded {max_bytes} bytes: {url}")
        if not 200 <= response.status < 300:
            raise RuntimeError(f"HTTP {response.status} from {url}")
        return data


def fetch_export(name):
    query = (
        "{ getDataOfFile(fileName: "
        + json.dumps(name)
        + ") { fileName version creationDate citation content } }"
    )
    request_body = json.dumps({"query": query}).encode("utf-8")
    response = json.loads(fetch_bytes(GRAPHQL_URL, body=request_body))
    if response.get("errors"):
        raise RuntimeError(f"RegulonDB GraphQL error for {name}: {response['errors']}")
    result = response.get("data", {}).get("getDataOfFile")
    if not result or result.get("fileName") != name or not result.get("content"):
        raise RuntimeError(f"RegulonDB did not return the expected {name} export")
    return result


def parse_tsv(path, expected_columns):
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    if not lines:
        raise RuntimeError(f"Empty reference export: {path.name}")
    header = lines[0].split("\t")
    if len(header) != expected_columns:
        raise RuntimeError(
            f"Unexpected {path.name} schema: expected {expected_columns} columns, "
            f"got {len(header)}"
        )
    rows = []
    for line_number, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != expected_columns:
            raise RuntimeError(
                f"Malformed {path.name} row {line_number}: "
                f"expected {expected_columns} columns, got {len(fields)}"
            )
        rows.append(fields)
    return rows


def parse_reference_genes(path):
    rows = parse_tsv(path, 12)
    genes = []
    gene_names = []
    sequence_length_mismatches = 0
    missing_coordinates = 0
    for row in rows:
        if not row[2].isdigit() or not row[3].isdigit():
            missing_coordinates += 1
            continue
        left, right = int(row[2]), int(row[3])
        strand = "+" if row[4] == "forward" else "-" if row[4] == "reverse" else None
        start_codon, stop_codon = row[7].upper(), row[8].upper()
        if strand is None or start_codon not in START_CODONS or stop_codon not in STOP_CODONS:
            continue
        if right < left:
            raise RuntimeError(f"Invalid RegulonDB gene coordinates for {row[1]}")
        sequence = row[9].upper()
        if len(sequence) != right - left + 1:
            sequence_length_mismatches += 1
        genes.append((left - 1, right, strand))
        gene_names.append(row[1])
    return genes, gene_names, len(rows), sequence_length_mismatches, missing_coordinates


def parse_reference_sites(path):
    rows = parse_tsv(path, 27)
    all_sites = set()
    confirmed_or_strong_sites = set()
    missing_coordinates = 0
    descending_coordinates = 0
    for row in rows:
        left, right = row[6].strip(), row[7].strip()
        if not left.isdigit() or not right.isdigit():
            missing_coordinates += 1
            continue
        left, right = int(left), int(right)
        if left > right:
            descending_coordinates += 1
        start, end = min(left, right) - 1, max(left, right)
        if end <= start:
            raise RuntimeError(f"Invalid RegulonDB binding-site coordinates: {left}..{right}")
        interval = (start, end)
        all_sites.add(interval)
        if row[19].strip() in {"C", "S"}:
            confirmed_or_strong_sites.add(interval)
    return (
        sorted(all_sites),
        sorted(confirmed_or_strong_sites),
        len(rows),
        missing_coordinates,
        descending_coordinates,
    )


def parse_predictions(output):
    predictions = []
    for start, end, strand in SPAN_RE.findall(output):
        start, end = int(start), int(end)
        if end > start:
            predictions.append((start, end, strand))
    return predictions


def parse_operator_site_predictions(output):
    sites = {
        (int(start), int(end))
        for start, end in OPERATOR_SITE_RE.findall(output)
        if int(end) > int(start)
    }
    return sorted(sites)


def candidate_annotation_counts(output):
    return {
        "candidate_site_annotation_groups": output.count(
            "[CANDIDATE_OPERATOR_SITES:"
        ),
        "candidate_site_to_gene_annotations": output.count("operator_site("),
        "candidate_operator_site_coordinates": len(
            parse_operator_site_predictions(output)
        ),
        "unterminated_readthrough_calls": output.count(
            "[WARN: UNTERMINATED_READTHROUGH]"
        ),
    }


def parse_regulatory_edges(output):
    return [
        (int(regulator), int(target), "-")
        for regulator, target in REGULATORY_EDGE_RE.findall(output)
    ]


def interval_iou(left, right):
    overlap = max(0, min(left[1], right[1]) - max(left[0], right[0]))
    union = (left[1] - left[0]) + (right[1] - right[0]) - overlap
    return overlap / union if union else 0.0


def maximum_matching(adjacency, right_count):
    """Hopcroft-Karp maximum-cardinality matching on a sparse interval graph."""
    sys.setrecursionlimit(max(sys.getrecursionlimit(), len(adjacency) + 100))
    left_match = [-1] * len(adjacency)
    right_match = [-1] * right_count
    distance = [0] * len(adjacency)
    infinity = len(adjacency) + right_count + 1

    def breadth_first_search():
        queue = []
        head = 0
        shortest = infinity
        for left, matched in enumerate(left_match):
            if matched == -1:
                distance[left] = 0
                queue.append(left)
            else:
                distance[left] = infinity
        while head < len(queue):
            left = queue[head]
            head += 1
            if distance[left] >= shortest:
                continue
            for right in adjacency[left]:
                next_left = right_match[right]
                if next_left == -1:
                    shortest = distance[left] + 1
                elif distance[next_left] == infinity:
                    distance[next_left] = distance[left] + 1
                    queue.append(next_left)
        return shortest

    def depth_first_search(left, shortest):
        for right in adjacency[left]:
            next_left = right_match[right]
            if next_left == -1 and distance[left] + 1 == shortest:
                left_match[left] = right
                right_match[right] = left
                return True
            if next_left != -1 and (
                distance[next_left] == distance[left] + 1
                and depth_first_search(next_left, shortest)
            ):
                left_match[left] = right
                right_match[right] = left
                return True
        distance[left] = infinity
        return False

    matched = 0
    while True:
        shortest = breadth_first_search()
        if shortest == infinity:
            return [
                (left, right)
                for left, right in enumerate(left_match)
                if right != -1
            ]
        for left, right in enumerate(left_match):
            if right == -1 and depth_first_search(left, shortest):
                matched += 1


def gene_match_pairs(predictions, reference, threshold):
    by_strand = {}
    for gold_index, gold in enumerate(reference):
        by_strand.setdefault(gold[2], []).append((gold_index, gold))
    gold_starts = {}
    max_gold_length = max((end - start for start, end, _ in reference), default=0)
    for strand, intervals in by_strand.items():
        intervals.sort(key=lambda item: item[1][0])
        gold_starts[strand] = [item[1][0] for item in intervals]

    adjacency = [[] for _ in predictions]
    for pred_index, pred in enumerate(predictions):
        intervals = by_strand.get(pred[2], [])
        starts = gold_starts.get(pred[2], [])
        low = bisect.bisect_left(starts, pred[0] - max_gold_length)
        high = bisect.bisect_left(starts, pred[1])
        for local_index in range(low, high):
            gold_index, gold = intervals[local_index]
            if gold[1] > pred[0] and interval_iou(pred, gold) >= threshold:
                adjacency[pred_index].append(gold_index)

    return maximum_matching(adjacency, len(reference))


def gene_metrics(predictions, reference, threshold):
    true_positives = len(gene_match_pairs(predictions, reference, threshold))
    precision = true_positives / len(predictions) if predictions else None
    recall = true_positives / len(reference) if reference else None
    f1 = (
        2 * precision * recall / (precision + recall)
        if precision is not None and recall is not None and precision + recall
        else 0.0
        if precision is not None and recall is not None
        else None
    )
    return {
        "iou_threshold": threshold,
        "matched": true_positives,
        "predicted": len(predictions),
        "reference": len(reference),
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


def operator_site_match_pairs(predictions, reference, tolerance):
    reference_centers = sorted(
        ((start + end, index) for index, (start, end) in enumerate(reference))
    )
    center_values = [center for center, _ in reference_centers]
    adjacency = [[] for _ in predictions]
    doubled_tolerance = 2 * tolerance
    for prediction_index, (start, end) in enumerate(predictions):
        center = start + end
        low = bisect.bisect_left(
            center_values, center - doubled_tolerance
        )
        high = bisect.bisect_right(
            center_values, center + doubled_tolerance
        )
        adjacency[prediction_index] = [
            reference_centers[index][1] for index in range(low, high)
        ]
    return maximum_matching(adjacency, len(reference))


def operator_site_metrics(predictions, reference):
    metrics = {}
    for tolerance in (0, 5, 10):
        true_positives = len(
            operator_site_match_pairs(predictions, reference, tolerance)
        )
        precision = true_positives / len(predictions) if predictions else None
        recall = true_positives / len(reference) if reference else None
        f1 = (
            2 * precision * recall / (precision + recall)
            if precision is not None
            and recall is not None
            and precision + recall
            else 0.0
            if precision is not None and recall is not None
            else None
        )
        metrics[f"center_within_{tolerance}_bp"] = {
            "matched": true_positives,
            "predicted": len(predictions),
            "reference": len(reference),
            "precision": precision,
            "recall": recall,
            "f1": f1,
        }
    return metrics


def reference_regulatory_edges(path):
    rows = parse_tsv(path, 7)
    edge_confidence = {}
    confirmed_signed_edges = set()
    for row in rows:
        regulator_gene, target_gene, effect, confidence = (
            row[2].strip(),
            row[4].strip(),
            row[5].strip(),
            row[6].strip(),
        )
        if regulator_gene and target_gene and effect in {"+", "-", "-+"}:
            edge_confidence.setdefault((regulator_gene, target_gene, effect), set()).add(
                confidence
            )
            if effect in {"+", "-"} and confidence in {"C", "S"}:
                confirmed_signed_edges.add((regulator_gene, target_gene, effect))
    confirmed = sum("C" in confidence for confidence in edge_confidence.values())
    confirmed_or_strong = sum(
        bool(confidence & {"C", "S"}) for confidence in edge_confidence.values()
    )
    confirmed_or_strong_regulators = {
        regulator
        for (regulator, _, _), confidence in edge_confidence.items()
        if confidence & {"C", "S"}
    }
    stats = {
        "source_rows": len(rows),
        "unique_signed_gene_edges_all_confidences": len(edge_confidence),
        "confirmed_edges": confirmed,
        "confirmed_or_strong_edges": confirmed_or_strong,
        "confirmed_or_strong_unambiguous_signed_edges": len(confirmed_signed_edges),
        "confirmed_or_strong_regulators": len(confirmed_or_strong_regulators),
    }
    return stats, confirmed_signed_edges


def signed_edge_metrics(predicted_edges, gene_names_by_prediction, reference_edges):
    mapped_predictions = {
        (
            gene_names_by_prediction[regulator],
            gene_names_by_prediction[target],
            effect,
        )
        for regulator, target, effect in predicted_edges
        if regulator in gene_names_by_prediction and target in gene_names_by_prediction
    }
    true_positives = len(mapped_predictions & reference_edges)
    precision = true_positives / len(predicted_edges) if predicted_edges else None
    mapped_precision = (
        true_positives / len(mapped_predictions) if mapped_predictions else None
    )
    recall = true_positives / len(reference_edges) if reference_edges else None
    f1 = (
        2 * precision * recall / (precision + recall)
        if precision is not None and recall is not None and precision + recall
        else 0.0
        if precision is not None and recall is not None
        else None
    )
    return {
        "matched_signed_edges": true_positives,
        "predicted_signed_edges": len(predicted_edges),
        "mapped_predicted_signed_edges": len(mapped_predictions),
        "reference_confirmed_or_strong_unambiguous_signed_edges": len(reference_edges),
        "precision": precision,
        "precision_among_gene_mapped_predictions": mapped_precision,
        "recall": recall,
        "f1": f1,
        "unmapped_prediction_edges": len(predicted_edges) - len(mapped_predictions),
        "precision_note": (
            "Overall precision counts gene-unmapped predictions as unmatched; "
            "mapped-only precision is undefined when no predicted edge maps to genes."
        ),
    }


def benchmark_scaling(binary, sequence, temp_path):
    sizes = [100_000, 250_000, 500_000, 1_000_000]
    points = []
    for size in sizes:
        fasta_path = temp_path / f"scaling_{size}.fasta"
        fasta_path.write_text(
            f">NC_000913.3_prefix_{size}\n{sequence[:size]}\n",
            encoding="ascii",
        )
        result, elapsed = run_dnadec(binary, fasta_path, True)
        if result.returncode != 0:
            raise RuntimeError(
                f"dnadec failed on {size}-bp scaling input: "
                + "\n".join(result.stderr.splitlines()[-20:])
            )
        points.append({"bases": size, "runtime_seconds": elapsed})

    log_bases = [math.log(point["bases"]) for point in points]
    log_times = [math.log(point["runtime_seconds"]) for point in points]
    mean_x = statistics.mean(log_bases)
    mean_y = statistics.mean(log_times)
    denominator = sum((x - mean_x) ** 2 for x in log_bases)
    exponent = (
        sum((x - mean_x) * (y - mean_y) for x, y in zip(log_bases, log_times))
        / denominator
    )
    return {
        "method": "single-run non-strict timing on prefixes of the same reference genome",
        "points": points,
        "log_log_runtime_exponent": exponent,
        "expected_linear_exponent": 1.0,
    }


def run_dnadec(binary, fasta_path, non_strict):
    args = [str(binary), str(fasta_path)]
    if non_strict:
        args.append("--non-strict")
    started = time.perf_counter()
    result = subprocess.run(
        args,
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=300,
        check=False,
    )
    elapsed = time.perf_counter() - started
    return result, elapsed


def sha256(data):
    return hashlib.sha256(data).hexdigest()


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

    temp = tempfile.TemporaryDirectory(prefix="dnadec-regulondb-")
    temp_path = Path(temp.name)
    try:
        genome_bytes = fetch_bytes(GENOME_URL)
        genome_path = temp_path / "e_coli_k12.fna"
        genome_path.write_bytes(genome_bytes)
        fasta_header = genome_bytes.decode("ascii").splitlines()[0]
        if not fasta_header.startswith(">NC_000913.3"):
            raise RuntimeError(f"Unexpected reference assembly: {fasta_header}")
        genome_sequence = "".join(
            line.strip()
            for line in genome_bytes.decode("ascii").splitlines()
            if line and not line.startswith(">")
        ).upper()
        if not genome_sequence or set(genome_sequence) - set("ACGT"):
            raise RuntimeError("Reference FASTA is empty or contains unexpected symbols")

        query = """
        query BenchmarkReferences {
          genes: getDataOfFile(fileName: "Gene_sequence") {
            fileName version creationDate citation content
          }
          network: getDataOfFile(fileName: "NetworkRegulatorGene") {
            fileName version creationDate citation content
          }
          sites: getDataOfFile(fileName: "TF-RISet") {
            fileName version creationDate citation content
          }
        }
        """
        request_body = json.dumps({"query": query}).encode("utf-8")
        response = json.loads(fetch_bytes(GRAPHQL_URL, body=request_body))
        if response.get("errors"):
            raise RuntimeError(f"RegulonDB GraphQL error: {response['errors']}")
        exports = response.get("data", {})
        gene_export = exports.get("genes")
        network_export = exports.get("network")
        site_export = exports.get("sites")
        if not gene_export or not network_export or not site_export:
            raise RuntimeError("RegulonDB did not return all requested reference exports")
        if gene_export["fileName"] != "Gene_sequence":
            raise RuntimeError("Unexpected gene-sequence export")
        if network_export["fileName"] != "NetworkRegulatorGene":
            raise RuntimeError("Unexpected regulatory-network export")
        if site_export["fileName"] != "TF-RISet":
            raise RuntimeError("Unexpected TF binding-site export")
        if len(
            {
                gene_export["version"],
                network_export["version"],
                site_export["version"],
            }
        ) != 1:
            raise RuntimeError("Reference export versions do not match")

        gene_path = temp_path / "Gene_sequence.tsv"
        network_path = temp_path / "NetworkRegulatorGene.tsv"
        site_path = temp_path / "TF-RISet.tsv"
        gene_bytes = gene_export["content"].encode("utf-8")
        network_bytes = network_export["content"].encode("utf-8")
        site_bytes = site_export["content"].encode("utf-8")
        gene_path.write_bytes(gene_bytes)
        network_path.write_bytes(network_bytes)
        site_path.write_bytes(site_bytes)

        (
            reference_genes,
            reference_gene_names,
            annotated_rows,
            sequence_length_mismatches,
            missing_coordinate_rows,
        ) = parse_reference_genes(gene_path)
        regulatory_reference, reference_signed_edges = reference_regulatory_edges(
            network_path
        )
        (
            reference_sites,
            confirmed_or_strong_reference_sites,
            site_annotation_rows,
            missing_site_coordinate_rows,
            descending_site_coordinate_rows,
        ) = parse_reference_sites(site_path)
        genome_length = len(genome_sequence)
        if any(end > genome_length for _, end, _ in reference_genes):
            raise RuntimeError("RegulonDB gene coordinates exceed the reference genome")
        if any(end > genome_length for _, end in reference_sites):
            raise RuntimeError("RegulonDB binding-site coordinates exceed the reference genome")

        strict_result, strict_seconds = run_dnadec(binary, genome_path, False)
        non_strict_runs = [
            run_dnadec(binary, genome_path, True)
            for _ in range(4)
        ]
        non_strict_result, _ = non_strict_runs[-1]
        if non_strict_result.returncode != 0:
            raise RuntimeError(
                "dnadec failed in --non-strict mode: "
                + "\n".join(non_strict_result.stderr.splitlines()[-20:])
            )
        predictions = parse_predictions(non_strict_result.stdout)
        candidate_annotations = candidate_annotation_counts(non_strict_result.stdout)
        candidate_sites = parse_operator_site_predictions(non_strict_result.stdout)
        site_metrics = operator_site_metrics(
            candidate_sites, confirmed_or_strong_reference_sites
        )
        run_times = [seconds for _, seconds in non_strict_runs[1:]]
        predicted_edges = parse_regulatory_edges(non_strict_result.stdout)
        gene_matches = gene_match_pairs(predictions, reference_genes, 0.9)
        gene_names_by_prediction = {
            prediction_index: reference_gene_names[reference_index]
            for prediction_index, reference_index in gene_matches
        }
        logic_metrics = signed_edge_metrics(
            predicted_edges,
            gene_names_by_prediction,
            reference_signed_edges,
        )
        scaling_metrics = benchmark_scaling(binary, genome_sequence, temp_path)

        release_match = RELEASE_RE.search(gene_export["citation"])
        report = {
            "source": "official RegulonDB GraphQL exports and reference FASTA",
            "regulondb_release": release_match.group(1).strip() if release_match else "unknown",
            "export_creation_date": gene_export["creationDate"],
            "reference_assembly": "NC_000913.3",
            "genome_length_bp": genome_length,
            "gene_export_sha256": sha256(gene_bytes),
            "network_export_sha256": sha256(network_bytes),
            "site_export_sha256": sha256(site_bytes),
            "genome_fasta_sha256": sha256(genome_bytes),
            "gene_annotation": {
                "records": annotated_rows,
                "protein_coding_orf_records": len(reference_genes),
                "rows_missing_coordinates": missing_coordinate_rows,
                "sequence_coordinate_length_mismatches": sequence_length_mismatches,
            },
            "gene_detection_non_strict": {
                "median_runtime_seconds": statistics.median(run_times),
                "runs_including_warmup": len(non_strict_runs),
                "predicted_orfs": len(predictions),
                "iou_1_00": gene_metrics(predictions, reference_genes, 1.0),
                "iou_0_90": gene_metrics(predictions, reference_genes, 0.9),
                "iou_0_50": gene_metrics(predictions, reference_genes, 0.5),
            },
            "strict_mode": {
                "exit_code": strict_result.returncode,
                "runtime_seconds": strict_seconds,
                "predicted_orfs": len(parse_predictions(strict_result.stdout)),
                "stderr_lines": len(strict_result.stderr.splitlines()),
            },
            "regulatory_edges": {
                **regulatory_reference,
                **logic_metrics,
                "inference_rule": (
                    "Edges are inferred only between adjacent, same-strand, "
                    "non-coupled, confident genes with shared candidate motifs "
                    "or a local operator pair containing a 3x tandem site."
                ),
            },
            "candidate_binding_sites": {
                "reference_rows": site_annotation_rows,
                "reference_unique_site_coordinates": len(reference_sites),
                "confirmed_or_strong_unique_site_coordinates": len(
                    confirmed_or_strong_reference_sites
                ),
                "rows_missing_coordinates": missing_site_coordinate_rows,
                "rows_with_descending_coordinates": descending_site_coordinate_rows,
                "candidate_unique_site_coordinates": len(candidate_sites),
                "against_confirmed_or_strong_reference": site_metrics,
            },
            "complexity_scaling": scaling_metrics,
            "abstentions": {
                **candidate_annotations,
                "candidate_edge_accuracy_scored": True,
                "reason": (
                    "Unlinked operator candidates remain unresolved; inferred local edges "
                    "are scored against confirmed/strong signed RegulonDB gene edges."
                ),
            },
            "data_retained": False,
            "data_terms": "RegulonDB exports are for academic/noncommercial use; no data is redistributed.",
        }
    finally:
        temp.cleanup()

    report["temporary_downloads_deleted"] = not temp_path.exists()
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
