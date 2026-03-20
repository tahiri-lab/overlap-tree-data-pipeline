#!/usr/bin/env python3
"""
Validate the illustrative supertree results.

Expected folder layout:
supertree_illustration/
├── input_datasets/
│   ├── multiset_1.txt
│   ├── ...
│   └── multiset_100.txt
├── supertrees/
│   ├── supertrees_sfit.txt
│   ├── supertrees_dfit.txt
│   ├── supertrees_nj.txt
│   ├── supertrees_mrplus.txt
│   └── supertrees_scs.txt
└── supertree_validation.py

Outputs:
- supertree_validation_results_summary.csv
- supertree_validation_results_detailed.csv
"""

from __future__ import annotations

import argparse
import csv
import statistics
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set

from Bio import Phylo


METHOD_FILES = {
    "Split Fit": "supertrees/supertrees_sfit.txt",
    "Most Similar": "supertrees/supertrees_dfit.txt",
    "Average NJ": "supertrees/supertrees_nj.txt",
    "Majority-Rule": "supertrees/supertrees_mrplus.txt",
    "Spectral Clustering": "supertrees/supertrees_scs.txt",
}


@dataclass
class InputDatasetInfo:
    dataset_id: int
    path: Path
    parsed_ok: bool
    tree_count: int
    expected_taxa: Set[str]
    expected_taxa_count: int
    duplicate_leaf_labels_found: bool
    empty_leaf_labels_found: bool
    error_message: str = ""


@dataclass
class OutputAuditRow:
    method: str
    dataset_id: int
    output_tree_present: bool
    parsed_ok: bool
    output_taxa_count: int
    expected_taxa_count: int
    all_expected_taxa_covered: bool
    exact_taxon_set_match: bool
    duplicate_leaf_labels_found: bool
    empty_leaf_labels_found: bool
    analysis_ready: bool
    success: bool
    error_message: str = ""


def default_base_dir() -> Path:
    if "__file__" in globals():
        return Path(__file__).resolve().parent
    return Path.cwd()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate the illustrative supertree outputs against the 100 input datasets."
    )
    parser.add_argument(
        "base_dir",
        nargs="?",
        type=Path,
        default=default_base_dir(),
        help=(
            "Path to supertree_illustration/. "
            "Default: the folder containing this script (or the current working directory)."
        ),
    )
    parser.add_argument(
        "--expected-datasets",
        type=int,
        default=100,
        help="Expected number of multiset files and output supertrees per method. Default: 100",
    )
    parser.add_argument(
        "--expected-input-trees",
        type=int,
        default=30,
        help="Expected number of input trees per multiset. Default: 30",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Directory for output files. Default: base_dir",
    )
    return parser.parse_args()


def split_newick_trees(text: str) -> List[str]:
    """
    Split a text blob containing one or more Newick trees terminated by ';'.
    """
    trees: List[str] = []
    current: List[str] = []
    for ch in text:
        current.append(ch)
        if ch == ";":
            tree = "".join(current).strip()
            if tree:
                trees.append(tree)
            current = []
    tail = "".join(current).strip()
    if tail:
        raise ValueError("Trailing non-empty text found after the last ';'")
    return trees


def parse_tree_taxa(newick: str) -> tuple[Set[str], int, bool, bool]:
    """
    Parse a Newick tree with Biopython and return:
    - set of terminal taxon labels
    - terminal count
    - duplicate label flag
    - empty label flag
    """
    tree = Phylo.read(StringIO(newick), "newick")
    labels: List[str] = []
    empty_found = False

    for terminal in tree.get_terminals():
        label = terminal.name
        if label is None:
            empty_found = True
            label = ""
        label = label.strip()
        if not label:
            empty_found = True
        labels.append(label)

    label_set = set(labels)
    duplicate_found = len(labels) != len(label_set)
    if "" in label_set:
       
        label_set.remove("")

    return label_set, len(labels), duplicate_found, empty_found


def load_input_datasets(
    input_dir: Path,
    expected_datasets: int,
    expected_input_trees: int,
) -> List[InputDatasetInfo]:
    infos: List[InputDatasetInfo] = []

    for dataset_id in range(1, expected_datasets + 1):
        path = input_dir / f"multiset_{dataset_id}.txt"
        if not path.exists():
            infos.append(
                InputDatasetInfo(
                    dataset_id=dataset_id,
                    path=path,
                    parsed_ok=False,
                    tree_count=0,
                    expected_taxa=set(),
                    expected_taxa_count=0,
                    duplicate_leaf_labels_found=False,
                    empty_leaf_labels_found=False,
                    error_message="Input dataset file not found",
                )
            )
            continue

        try:
            trees = split_newick_trees(path.read_text(encoding="utf-8"))
            expected_taxa: Set[str] = set()
            duplicate_found = False
            empty_found = False

            for newick in trees:
                taxa, _, dup, empty = parse_tree_taxa(newick)
                expected_taxa |= taxa
                duplicate_found = duplicate_found or dup
                empty_found = empty_found or empty

            infos.append(
                InputDatasetInfo(
                    dataset_id=dataset_id,
                    path=path,
                    parsed_ok=len(trees) == expected_input_trees and not empty_found,
                    tree_count=len(trees),
                    expected_taxa=expected_taxa,
                    expected_taxa_count=len(expected_taxa),
                    duplicate_leaf_labels_found=duplicate_found,
                    empty_leaf_labels_found=empty_found,
                    error_message=(
                        ""
                        if len(trees) == expected_input_trees and not empty_found
                        else f"Expected {expected_input_trees} input trees, found {len(trees)}"
                    ),
                )
            )
        except Exception as exc:
            infos.append(
                InputDatasetInfo(
                    dataset_id=dataset_id,
                    path=path,
                    parsed_ok=False,
                    tree_count=0,
                    expected_taxa=set(),
                    expected_taxa_count=0,
                    duplicate_leaf_labels_found=False,
                    empty_leaf_labels_found=False,
                    error_message=str(exc),
                )
            )

    return infos


def audit_method_outputs(
    method_name: str,
    method_path: Path,
    inputs: Sequence[InputDatasetInfo],
) -> List[OutputAuditRow]:
    rows: List[OutputAuditRow] = []

    if not method_path.exists():
        for info in inputs:
            rows.append(
                OutputAuditRow(
                    method=method_name,
                    dataset_id=info.dataset_id,
                    output_tree_present=False,
                    parsed_ok=False,
                    output_taxa_count=0,
                    expected_taxa_count=info.expected_taxa_count,
                    all_expected_taxa_covered=False,
                    exact_taxon_set_match=False,
                    duplicate_leaf_labels_found=False,
                    empty_leaf_labels_found=False,
                    analysis_ready=False,
                    success=False,
                    error_message=f"Method file not found: {method_path}",
                )
            )
        return rows

    try:
        output_trees = split_newick_trees(method_path.read_text(encoding="utf-8"))
    except Exception as exc:
        for info in inputs:
            rows.append(
                OutputAuditRow(
                    method=method_name,
                    dataset_id=info.dataset_id,
                    output_tree_present=False,
                    parsed_ok=False,
                    output_taxa_count=0,
                    expected_taxa_count=info.expected_taxa_count,
                    all_expected_taxa_covered=False,
                    exact_taxon_set_match=False,
                    duplicate_leaf_labels_found=False,
                    empty_leaf_labels_found=False,
                    analysis_ready=False,
                    success=False,
                    error_message=f"Failed to read method file: {exc}",
                )
            )
        return rows

    for idx, info in enumerate(inputs):
        if idx >= len(output_trees):
            rows.append(
                OutputAuditRow(
                    method=method_name,
                    dataset_id=info.dataset_id,
                    output_tree_present=False,
                    parsed_ok=False,
                    output_taxa_count=0,
                    expected_taxa_count=info.expected_taxa_count,
                    all_expected_taxa_covered=False,
                    exact_taxon_set_match=False,
                    duplicate_leaf_labels_found=False,
                    empty_leaf_labels_found=False,
                    analysis_ready=False,
                    success=False,
                    error_message="No output tree present for this dataset index",
                )
            )
            continue

        newick = output_trees[idx]

        try:
            output_taxa, output_taxa_count, duplicate_found, empty_found = parse_tree_taxa(newick)

            # "All expected taxa covered":
            # every expected taxon from the corresponding input set must be present.
            covered = info.expected_taxa.issubset(output_taxa)

            # Exact match is useful for auditing.
            exact_match = output_taxa == info.expected_taxa

            # "Analysis-ready" operationalization:
            # valid parse + non-empty taxon labels + no duplicate leaf labels
            analysis_ready = not duplicate_found and not empty_found

            # "Success rate" is the proportion of datasets for which
            # a supertree was generated and processed without error.
            success = analysis_ready

            rows.append(
                OutputAuditRow(
                    method=method_name,
                    dataset_id=info.dataset_id,
                    output_tree_present=True,
                    parsed_ok=True,
                    output_taxa_count=output_taxa_count,
                    expected_taxa_count=info.expected_taxa_count,
                    all_expected_taxa_covered=covered,
                    exact_taxon_set_match=exact_match,
                    duplicate_leaf_labels_found=duplicate_found,
                    empty_leaf_labels_found=empty_found,
                    analysis_ready=analysis_ready,
                    success=success,
                    error_message="",
                )
            )
        except Exception as exc:
            rows.append(
                OutputAuditRow(
                    method=method_name,
                    dataset_id=info.dataset_id,
                    output_tree_present=True,
                    parsed_ok=False,
                    output_taxa_count=0,
                    expected_taxa_count=info.expected_taxa_count,
                    all_expected_taxa_covered=False,
                    exact_taxon_set_match=False,
                    duplicate_leaf_labels_found=False,
                    empty_leaf_labels_found=False,
                    analysis_ready=False,
                    success=False,
                    error_message=str(exc),
                )
            )

    return rows


def summarize_method(rows: Sequence[OutputAuditRow]) -> Dict[str, str]:
    successful = [r for r in rows if r.success]
    parsed = [r for r in rows if r.parsed_ok]

    success_rate = 100.0 * len(successful) / len(rows) if rows else 0.0
    avg_taxa = statistics.mean(r.output_taxa_count for r in parsed) if parsed else 0.0
    all_expected_yes = all(r.all_expected_taxa_covered for r in rows) if rows else False
    analysis_ready_yes = all(r.analysis_ready for r in rows) if rows else False

    return {
        "Method": rows[0].method if rows else "",
        "Success rate": f"{success_rate:.0f}%",
        "Average number of taxa in output trees": f"{avg_taxa:.1f}",
        "All expected taxa covered": "Yes" if all_expected_yes else "No",
        "Analysis-ready": "Yes" if analysis_ready_yes else "No",
    }


def write_detailed_csv(rows: Sequence[OutputAuditRow], path: Path) -> None:
    fieldnames = [
        "method",
        "dataset_id",
        "output_tree_present",
        "parsed_ok",
        "output_taxa_count",
        "expected_taxa_count",
        "all_expected_taxa_covered",
        "exact_taxon_set_match",
        "duplicate_leaf_labels_found",
        "empty_leaf_labels_found",
        "analysis_ready",
        "success",
        "error_message",
    ]
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(
                {
                    "method": r.method,
                    "dataset_id": r.dataset_id,
                    "output_tree_present": r.output_tree_present,
                    "parsed_ok": r.parsed_ok,
                    "output_taxa_count": r.output_taxa_count,
                    "expected_taxa_count": r.expected_taxa_count,
                    "all_expected_taxa_covered": r.all_expected_taxa_covered,
                    "exact_taxon_set_match": r.exact_taxon_set_match,
                    "duplicate_leaf_labels_found": r.duplicate_leaf_labels_found,
                    "empty_leaf_labels_found": r.empty_leaf_labels_found,
                    "analysis_ready": r.analysis_ready,
                    "success": r.success,
                    "error_message": r.error_message,
                }
            )


def write_summary_csv(summary_rows: Sequence[Dict[str, str]], path: Path) -> None:
    fieldnames = [
        "Method",
        "Success rate",
        "Average number of taxa in output trees",
        "All expected taxa covered",
        "Analysis-ready",
    ]
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)


def write_latex_table(summary_rows: Sequence[Dict[str, str]], path: Path) -> None:
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\\begin{tabular}{lcccc}\n")
        fh.write("\\toprule\n")
        fh.write(
            "\\textbf{Method} & "
            "\\textbf{\\shortstack{Success\\\\ rate}} & "
            "\\textbf{\\shortstack{Average number of\\\\ taxa in output trees}} & "
            "\\textbf{\\shortstack{All expected\\\\ taxa covered}} & "
            "\\textbf{\\shortstack{Analysis-\\\\ ready}} \\\\\n"
        )
        fh.write("\\midrule\n")
        for row in summary_rows:
            fh.write(
                f"{row['Method']} & "
                f"{row['Success rate']} & "
                f"{row['Average number of taxa in output trees']} & "
                f"{row['All expected taxa covered']} & "
                f"{row['Analysis-ready']} \\\\\n"
            )
        fh.write("\\bottomrule\n")
        fh.write("\\end{tabular}\n")


def print_summary(summary_rows: Sequence[Dict[str, str]]) -> None:
    headers = [
        "Method",
        "Success rate",
        "Average number of taxa in output trees",
        "All expected taxa covered",
        "Analysis-ready",
    ]
    widths = {h: len(h) for h in headers}
    for row in summary_rows:
        for h in headers:
            widths[h] = max(widths[h], len(row[h]))

    def fmt(row: Dict[str, str]) -> str:
        return " | ".join(row[h].ljust(widths[h]) for h in headers)

    print(fmt({h: h for h in headers}))
    print("-+-".join("-" * widths[h] for h in headers))
    for row in summary_rows:
        print(fmt(row))


def main() -> int:
    args = parse_args()
    base_dir = args.base_dir.expanduser().resolve()
    outdir = (args.outdir.expanduser().resolve() if args.outdir else base_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    input_dir = base_dir / "input_datasets"
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input folder not found: {input_dir}")

    inputs = load_input_datasets(
        input_dir=input_dir,
        expected_datasets=args.expected_datasets,
        expected_input_trees=args.expected_input_trees,
    )

    all_rows: List[OutputAuditRow] = []
    summary_rows: List[Dict[str, str]] = []

    for method_name, relpath in METHOD_FILES.items():
        rows = audit_method_outputs(method_name, base_dir / relpath, inputs)
        all_rows.extend(rows)
        summary_rows.append(summarize_method(rows))

    write_detailed_csv(all_rows, outdir / "supertree_validation_results_detailed.csv")
    write_summary_csv(summary_rows, outdir / "supertree_validation_results_summary.csv")
    #write_latex_table(summary_rows, outdir / "supertree_results_table.tex")

    print_summary(summary_rows)
    print(f"\nWrote outputs to: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())