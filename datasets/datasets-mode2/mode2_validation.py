#!/usr/bin/env python3
"""
Audit Mode 2 benchmark datasets.

This script computes, for every dataset:
- number of input trees
- input-tree size min/max (number of taxa per tree)
- pairwise Jaccard overlap min/max across all unordered tree pairs
- minimum number of shared taxa across all unordered tree pairs
- number of taxa shared by all input trees (the observed anchor core)
- overlap-constraint violations
- minimum-shared-taxa violations
- full coverage against the corresponding reference tree
- duplicate leaf label and branch-length sanity checks

The script writes:
- one CSV with all datasets across all groups
- one CSV with group-level aggregate summaries
- one CSV with selected rows (default: datasets 1, 5, 10, 15, 20)
- a LaTeX table per species group for the selected rows
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence

LEAF_LABEL_RE = re.compile(r"(?<=[,(])([^():;,\[\]'\"\s]+)(?=\s*[:),;])")
BRANCH_LENGTH_RE = re.compile(r":\s*(-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)")
MULTISET_RE = re.compile(r"multiset_(\d+)\.txt$")

DEFAULT_SELECTED = (1, 5, 10, 15, 20)


@dataclass
class DatasetAudit:
    species_group: str
    dataset_id: int
    n_input_trees: int
    n_reference_trees_in_file: int
    reference_taxa_count: int
    union_taxa_count: int
    size_min: int
    size_max: int
    overlap_min: float
    overlap_max: float
    min_shared_taxa: int
    anchors_shared_by_all_trees: int
    total_pairs: int
    overlap_violating_pairs: int
    min_shared_violating_pairs: int
    duplicate_leaf_labels_found: bool
    non_positive_branch_lengths_found: bool
    all_input_trees_subset_of_reference: bool
    full_coverage_ok: bool
    reference_match_ok: bool


@dataclass
class GroupSummary:
    species_group: str
    datasets_checked: int
    total_input_trees: int
    total_pairs: int
    overlap_violating_pairs_total: int
    min_shared_violating_pairs_total: int
    overlap_min_across_group: float
    overlap_max_across_group: float
    min_shared_taxa_across_group: int
    anchors_shared_by_all_trees_min: int
    anchors_shared_by_all_trees_max: int
    full_coverage_failures: int
    subset_failures: int
    duplicate_label_failures: int
    non_positive_branch_length_failures: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit the Mode 2 datasets in the overlap treedata pipeline."
    )
    parser.add_argument(
        "input_path",
        type=Path,
        help=(
            "Path to either the repository root or directly to the "
            "datasets/datasets-mode2 directory."
        ),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("mode2_audit_outputs"),
        help="Directory where CSV and LaTeX outputs will be written.",
    )
    parser.add_argument(
        "--jaccard-min",
        type=float,
        default=0.30,
        help="Configured lower Jaccard bound. Default: 0.30",
    )
    parser.add_argument(
        "--jaccard-max",
        type=float,
        default=0.70,
        help="Configured upper Jaccard bound. Default: 0.70",
    )
    parser.add_argument(
        "--required-shared",
        type=int,
        default=10,
        help=(
            "Minimum required number of shared taxa per tree pair. "
            "For the released datasets this is 10 because 10 anchors were used."
        ),
    )
    parser.add_argument(
        "--selected-datasets",
        type=int,
        nargs="+",
        default=list(DEFAULT_SELECTED),
        help="Dataset IDs to export into the paper-style selected-rows table.",
    )
    return parser.parse_args()


def resolve_mode2_root(path: Path) -> Path:
    path = path.expanduser().resolve()
    if path.name == "datasets-mode2" and path.is_dir():
        return path
    candidate = path / "datasets" / "datasets-mode2"
    if candidate.is_dir():
        return candidate
    raise FileNotFoundError(
        f"Could not find datasets/datasets-mode2 under {path}. "
        "Pass either the repository root or the datasets-mode2 directory itself."
    )


def extract_newick_trees(text: str) -> List[str]:
    """Extract Newick trees terminated by ';' from a text blob."""
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
        raise ValueError("Input contains trailing non-empty text after the last ';'.")
    return trees


def extract_leaf_labels(newick: str) -> List[str]:
    return LEAF_LABEL_RE.findall(newick)


def extract_branch_lengths(newick: str) -> List[float]:
    return [float(x) for x in BRANCH_LENGTH_RE.findall(newick)]


def dataset_number(path: Path) -> int:
    match = MULTISET_RE.search(path.name)
    if not match:
        raise ValueError(f"Could not parse dataset number from {path.name}")
    return int(match.group(1))


def format_range_int(low: int, high: int) -> str:
    return f"{low}--{high}"


def format_range_float(low: float, high: float, digits: int = 3) -> str:
    return f"{low:.{digits}f}--{high:.{digits}f}"


def discover_species_dirs(mode2_root: Path) -> List[Path]:
    species_dirs = [
        p for p in mode2_root.iterdir() if p.is_dir() and p.name.endswith("_datasets")
    ]
    species_dirs.sort(key=lambda p: p.name)
    return species_dirs


def species_name_from_dir(species_dir: Path) -> str:
    name = species_dir.name
    if not name.endswith("_datasets"):
        raise ValueError(f"Unexpected species directory name: {name}")
    return name[: -len("_datasets")]


def compute_dataset_audit(
    species_group: str,
    dataset_id: int,
    multiset_file: Path,
    reference_trees: Sequence[str],
    jaccard_min: float,
    jaccard_max: float,
    required_shared: int,
) -> DatasetAudit:
    multiset_text = multiset_file.read_text(encoding="utf-8")
    input_trees = extract_newick_trees(multiset_text)
    if dataset_id < 1 or dataset_id > len(reference_trees):
        raise IndexError(
            f"Dataset {dataset_id} for {species_group} has no matching reference tree. "
            f"Reference file contains {len(reference_trees)} trees."
        )
    reference_tree = reference_trees[dataset_id - 1]

    tree_taxa_sets: List[set[str]] = []
    sizes: List[int] = []
    duplicate_leaf_labels_found = False
    non_positive_branch_lengths_found = False

    for tree in input_trees:
        labels = extract_leaf_labels(tree)
        label_set = set(labels)
        tree_taxa_sets.append(label_set)
        sizes.append(len(label_set))
        if len(labels) != len(label_set):
            duplicate_leaf_labels_found = True
        branch_lengths = extract_branch_lengths(tree)
        if any(x <= 0 for x in branch_lengths):
            non_positive_branch_lengths_found = True

    ref_labels = extract_leaf_labels(reference_tree)
    ref_set = set(ref_labels)
    reference_match_ok = len(ref_labels) == len(ref_set)

    if not tree_taxa_sets:
        raise ValueError(f"No trees found in {multiset_file}")

    intersection_all = set.intersection(*tree_taxa_sets)
    union_all = set.union(*tree_taxa_sets)
    subset_ok = all(taxa.issubset(ref_set) for taxa in tree_taxa_sets)
    full_coverage_ok = union_all == ref_set

    overlaps: List[float] = []
    shared_counts: List[int] = []
    overlap_violating_pairs = 0
    min_shared_violating_pairs = 0

    for a, b in itertools.combinations(tree_taxa_sets, 2):
        inter = len(a & b)
        union = len(a | b)
        jaccard = inter / union if union else math.nan
        overlaps.append(jaccard)
        shared_counts.append(inter)
        if not (jaccard_min <= jaccard <= jaccard_max):
            overlap_violating_pairs += 1
        if inter < required_shared:
            min_shared_violating_pairs += 1

    total_pairs = len(overlaps)
    if total_pairs == 0:
        overlap_min = math.nan
        overlap_max = math.nan
        min_shared = 0
    else:
        overlap_min = min(overlaps)
        overlap_max = max(overlaps)
        min_shared = min(shared_counts)

    return DatasetAudit(
        species_group=species_group,
        dataset_id=dataset_id,
        n_input_trees=len(input_trees),
        n_reference_trees_in_file=len(reference_trees),
        reference_taxa_count=len(ref_set),
        union_taxa_count=len(union_all),
        size_min=min(sizes),
        size_max=max(sizes),
        overlap_min=overlap_min,
        overlap_max=overlap_max,
        min_shared_taxa=min_shared,
        anchors_shared_by_all_trees=len(intersection_all),
        total_pairs=total_pairs,
        overlap_violating_pairs=overlap_violating_pairs,
        min_shared_violating_pairs=min_shared_violating_pairs,
        duplicate_leaf_labels_found=duplicate_leaf_labels_found,
        non_positive_branch_lengths_found=non_positive_branch_lengths_found,
        all_input_trees_subset_of_reference=subset_ok,
        full_coverage_ok=full_coverage_ok,
        reference_match_ok=reference_match_ok,
    )


def build_group_summary(audits: Sequence[DatasetAudit]) -> GroupSummary:
    if not audits:
        raise ValueError("Cannot build a group summary from an empty audit list")
    return GroupSummary(
        species_group=audits[0].species_group,
        datasets_checked=len(audits),
        total_input_trees=sum(a.n_input_trees for a in audits),
        total_pairs=sum(a.total_pairs for a in audits),
        overlap_violating_pairs_total=sum(a.overlap_violating_pairs for a in audits),
        min_shared_violating_pairs_total=sum(a.min_shared_violating_pairs for a in audits),
        overlap_min_across_group=min(a.overlap_min for a in audits),
        overlap_max_across_group=max(a.overlap_max for a in audits),
        min_shared_taxa_across_group=min(a.min_shared_taxa for a in audits),
        anchors_shared_by_all_trees_min=min(a.anchors_shared_by_all_trees for a in audits),
        anchors_shared_by_all_trees_max=max(a.anchors_shared_by_all_trees for a in audits),
        full_coverage_failures=sum(not a.full_coverage_ok for a in audits),
        subset_failures=sum(not a.all_input_trees_subset_of_reference for a in audits),
        duplicate_label_failures=sum(a.duplicate_leaf_labels_found for a in audits),
        non_positive_branch_length_failures=sum(
            a.non_positive_branch_lengths_found for a in audits
        ),
    )


def write_all_datasets_csv(audits: Sequence[DatasetAudit], outpath: Path) -> None:
    fieldnames = [
        "species_group",
        "dataset_id",
        "n_input_trees",
        "n_reference_trees_in_file",
        "reference_taxa_count",
        "union_taxa_count",
        "size_min",
        "size_max",
        "overlap_min",
        "overlap_max",
        "min_shared_taxa",
        "anchors_shared_by_all_trees",
        "total_pairs",
        "overlap_violating_pairs",
        "min_shared_violating_pairs",
        "duplicate_leaf_labels_found",
        "non_positive_branch_lengths_found",
        "all_input_trees_subset_of_reference",
        "full_coverage_ok",
        "reference_match_ok",
    ]
    with outpath.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for a in audits:
            writer.writerow(
                {
                    "species_group": a.species_group,
                    "dataset_id": a.dataset_id,
                    "n_input_trees": a.n_input_trees,
                    "n_reference_trees_in_file": a.n_reference_trees_in_file,
                    "reference_taxa_count": a.reference_taxa_count,
                    "union_taxa_count": a.union_taxa_count,
                    "size_min": a.size_min,
                    "size_max": a.size_max,
                    "overlap_min": f"{a.overlap_min:.6f}",
                    "overlap_max": f"{a.overlap_max:.6f}",
                    "min_shared_taxa": a.min_shared_taxa,
                    "anchors_shared_by_all_trees": a.anchors_shared_by_all_trees,
                    "total_pairs": a.total_pairs,
                    "overlap_violating_pairs": a.overlap_violating_pairs,
                    "min_shared_violating_pairs": a.min_shared_violating_pairs,
                    "duplicate_leaf_labels_found": a.duplicate_leaf_labels_found,
                    "non_positive_branch_lengths_found": a.non_positive_branch_lengths_found,
                    "all_input_trees_subset_of_reference": a.all_input_trees_subset_of_reference,
                    "full_coverage_ok": a.full_coverage_ok,
                    "reference_match_ok": a.reference_match_ok,
                }
            )


def write_group_summary_csv(summaries: Sequence[GroupSummary], outpath: Path) -> None:
    fieldnames = [
        "species_group",
        "datasets_checked",
        "total_input_trees",
        "total_pairs",
        "overlap_violating_pairs_total",
        "min_shared_violating_pairs_total",
        "overlap_min_across_group",
        "overlap_max_across_group",
        "min_shared_taxa_across_group",
        "anchors_shared_by_all_trees_min",
        "anchors_shared_by_all_trees_max",
        "full_coverage_failures",
        "subset_failures",
        "duplicate_label_failures",
        "non_positive_branch_length_failures",
    ]
    with outpath.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for s in summaries:
            writer.writerow(
                {
                    "species_group": s.species_group,
                    "datasets_checked": s.datasets_checked,
                    "total_input_trees": s.total_input_trees,
                    "total_pairs": s.total_pairs,
                    "overlap_violating_pairs_total": s.overlap_violating_pairs_total,
                    "min_shared_violating_pairs_total": s.min_shared_violating_pairs_total,
                    "overlap_min_across_group": f"{s.overlap_min_across_group:.6f}",
                    "overlap_max_across_group": f"{s.overlap_max_across_group:.6f}",
                    "min_shared_taxa_across_group": s.min_shared_taxa_across_group,
                    "anchors_shared_by_all_trees_min": s.anchors_shared_by_all_trees_min,
                    "anchors_shared_by_all_trees_max": s.anchors_shared_by_all_trees_max,
                    "full_coverage_failures": s.full_coverage_failures,
                    "subset_failures": s.subset_failures,
                    "duplicate_label_failures": s.duplicate_label_failures,
                    "non_positive_branch_length_failures": s.non_positive_branch_length_failures,
                }
            )


def write_selected_rows_csv(
    audits: Sequence[DatasetAudit], selected_ids: Sequence[int], outpath: Path
) -> None:
    selected = [a for a in audits if a.dataset_id in set(selected_ids)]
    selected.sort(key=lambda a: (a.species_group, a.dataset_id))
    fieldnames = [
        "species_group",
        "dataset_label",
        "size_min_max",
        "overlap_min_max",
        "min_shared_taxa",
        "anchors",
    ]
    with outpath.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for a in selected:
            writer.writerow(
                {
                    "species_group": a.species_group,
                    "dataset_label": f"Dataset {a.dataset_id}",
                    "size_min_max": format_range_int(a.size_min, a.size_max),
                    "overlap_min_max": format_range_float(a.overlap_min, a.overlap_max),
                    "min_shared_taxa": a.min_shared_taxa,
                    "anchors": a.anchors_shared_by_all_trees,
                }
            )


def write_latex_tables(
    audits: Sequence[DatasetAudit], selected_ids: Sequence[int], outdir: Path
) -> None:
    by_group: dict[str, List[DatasetAudit]] = {}
    selected_set = set(selected_ids)
    for audit in audits:
        if audit.dataset_id in selected_set:
            by_group.setdefault(audit.species_group, []).append(audit)

    for group, rows in by_group.items():
        rows.sort(key=lambda a: a.dataset_id)
        path = outdir / f"{group}_selected_validation_table.tex"
        with path.open("w", encoding="utf-8") as fh:
            fh.write("\\begin{tabular}{lcccc}\n")
            fh.write("\\hline\n")
            fh.write("Datasets & Size min-max & Overlap min-max & Min shared taxa & Anchors \\\\\n")
            fh.write("\\hline\n")
            for a in rows:
                fh.write(
                    f"Dataset {a.dataset_id} & "
                    f"{format_range_int(a.size_min, a.size_max)} & "
                    f"{format_range_float(a.overlap_min, a.overlap_max)} & "
                    f"{a.min_shared_taxa} & "
                    f"{a.anchors_shared_by_all_trees} \\\\\n"
                )
            fh.write("\\hline\n")
            fh.write("\\end{tabular}\n")


def print_console_summary(
    group_summaries: Sequence[GroupSummary],
    all_audits: Sequence[DatasetAudit],
) -> None:
    print("\nMode 2 dataset audit completed.\n")
    for s in group_summaries:
        print(
            f"[{s.species_group}] datasets={s.datasets_checked}, "
            f"input_trees={s.total_input_trees}, pairs={s.total_pairs}, "
            f"overlap violations={s.overlap_violating_pairs_total}, "
            f"min-shared violations={s.min_shared_violating_pairs_total}, "
            f"observed overlap range={s.overlap_min_across_group:.3f}-{s.overlap_max_across_group:.3f}, "
            f"anchors shared by all trees={s.anchors_shared_by_all_trees_min}-{s.anchors_shared_by_all_trees_max}"
        )
    print()

    key_groups = sorted({a.species_group for a in all_audits})
    for group in key_groups:
        rows = [a for a in all_audits if a.species_group == group and a.dataset_id in DEFAULT_SELECTED]
        if not rows:
            continue
        rows.sort(key=lambda a: a.dataset_id)
        print(group)
        print("dataset\tsize_min-max\toverlap_min-max\tmin_shared\tanchors")
        for a in rows:
            print(
                f"{a.dataset_id}\t"
                f"{format_range_int(a.size_min, a.size_max)}\t"
                f"{format_range_float(a.overlap_min, a.overlap_max)}\t"
                f"{a.min_shared_taxa}\t"
                f"{a.anchors_shared_by_all_trees}"
            )
        print()


def main() -> int:
    args = parse_args()
    mode2_root = resolve_mode2_root(args.input_path)
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    species_dirs = discover_species_dirs(mode2_root)
    if not species_dirs:
        raise FileNotFoundError(f"No '*_datasets' directories found in {mode2_root}")

    all_audits: List[DatasetAudit] = []
    group_summaries: List[GroupSummary] = []

    for species_dir in species_dirs:
        species_group = species_name_from_dir(species_dir)
        ref_files = sorted(species_dir.glob("*reference_trees20.txt"))
        if len(ref_files) != 1:
            raise FileNotFoundError(
                f"Expected exactly one '*reference_trees20.txt' in {species_dir}, found {len(ref_files)}"
            )
        reference_trees = extract_newick_trees(ref_files[0].read_text(encoding="utf-8"))

        multiset_files = sorted(species_dir.glob("multiset_*.txt"), key=dataset_number)
        if not multiset_files:
            raise FileNotFoundError(f"No multiset_*.txt files found in {species_dir}")

        species_audits: List[DatasetAudit] = []
        for multiset_file in multiset_files:
            dataset_id = dataset_number(multiset_file)
            audit = compute_dataset_audit(
                species_group=species_group,
                dataset_id=dataset_id,
                multiset_file=multiset_file,
                reference_trees=reference_trees,
                jaccard_min=args.jaccard_min,
                jaccard_max=args.jaccard_max,
                required_shared=args.required_shared,
            )
            species_audits.append(audit)
            all_audits.append(audit)

        species_audits.sort(key=lambda a: a.dataset_id)
        group_summaries.append(build_group_summary(species_audits))

    all_audits.sort(key=lambda a: (a.species_group, a.dataset_id))
    group_summaries.sort(key=lambda s: s.species_group)

    write_all_datasets_csv(all_audits, outdir / "mode2_audit_all_datasets.csv")
    write_group_summary_csv(group_summaries, outdir / "mode2_audit_group_summaries.csv")
    write_selected_rows_csv(
        all_audits,
        args.selected_datasets,
        outdir / "mode2_audit_selected_rows.csv",
    )
    write_latex_tables(all_audits, args.selected_datasets, outdir)
    print_console_summary(group_summaries, all_audits)
    print(f"Wrote outputs to: {outdir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
