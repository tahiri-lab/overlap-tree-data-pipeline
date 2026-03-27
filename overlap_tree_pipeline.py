#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Overlap tree dataset pipeline with two modes.

Mode "service"
    Builds 10 overlapping taxa subsets from a selected base set of species and downloads
    pruned tree samples via the VertLife PhyloSubsets web service (Selenium).

Mode "reference"
    Uses one or more locally available broad-coverage VertLife trees as starting trees,
    prunes them to base taxa lists to create reference trees, then prunes each reference
    tree to multiple partially intersecting taxa subsets to create overlapping input trees.
    Optional topology and branch-length perturbations can be applied after pruning.

Usage:

  # Mode "service" (download pruned tree sets from VertLife PhyloSubsets)
  python overlap_dataset_pipeline.py service amphibians 120 550 name@example.com --seed 7

  # Mode "reference" (local full tree(s) -> reference trees + overlapping input trees)
  python overlap_dataset_pipeline.py reference \
      --group amphibians \
      --full_trees amphibians_full_trees20.txt \
      --outdir ./benchmarks \
      --base_sizes 50:145:5 \
      --n_input_trees 30 \
      --pairwise_overlap_range 0.30 0.70 \
      --anchor_taxa_count 10 \
      --nni_moves 3 \
      --length_scale_range 1.003 1.009 \
      --seed 7

Dependencies:
  - Python 3.9+
  - pandas
  - ete3
  - selenium + Chrome/Chromedriver for mode "service"

Author: AK
Last updated: March 2026
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

# Local modules are expected to sit alongside this script.
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))


def _require_ete3():
    try:
        from ete3 import Tree
        return Tree
    except Exception as e:
        raise RuntimeError(
            "ete3 is required. Install with: pip install ete3."
        ) from e


def _require_reference_modules():
    # Import heavy dependencies lazily so 'service' mode and '--help' work without ete3.
    Tree = _require_ete3()
    import gen_pruned_trees as gp
    import make_base_species_lists_phylo_stratified as strat
    return Tree, gp, strat


# Helpers

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def read_wide_species_csv(path: str, group: str) -> List[str]:
    df = pd.read_csv(path)
    col_map = {
        "amphibians": "Amphibians",
        "birds": "Birds",
        "mammals": "Mammals",
        "sharks": "Sharks",
        "squamates": "Squamates",
    }
    col = col_map.get(group.lower())
    if col is None or col not in df.columns:
        raise ValueError(f"Group '{group}' not found in {path}. Expected column: {col}")
    return df[col].dropna().astype(str).tolist()


def normalize_spaces_to_underscores(name: str) -> str:
    name = str(name).strip()
    name = " ".join(name.split())
    return name.replace(" ", "_")


def split_newick_file(path: str) -> List[str]:
    """
    Load one or more Newick trees from a file.
    Splits on semicolons and returns list of strings each ending with ';'.
    """
    with open(path, "r", encoding="utf-8") as f:
        content = f.read()
    parts = [p.strip() for p in re.split(r";\s*", content) if p.strip()]
    return [p + ";" for p in parts]


def write_newick_list(path: Path, newicks: Sequence[str]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for n in newicks:
            n = n.strip()
            if not n:
                continue
            if not n.endswith(";"):
                n += ";"
            f.write(n + "\n")


def prune_newick(newick_str: str, keep_taxa: Sequence[str]) -> Optional[str]:
    """
    Prune an ETE3 tree to keep_taxa while preserving branch lengths.
    Returns a Newick string or None if fewer than 2 taxa remain.
    """
    Tree = _require_ete3()
    t = Tree(newick_str, format=1)
    present = {lf.name for lf in t.iter_leaves()}
    keep = [x for x in keep_taxa if x in present]
    if len(keep) < 2:
        return None
    t.prune(keep, preserve_branch_length=True)
    # format=5 keeps internal supports + lengths (if present) and avoids quoting leaves
    return t.write(format=5).strip()


def compute_leafset(newick_str: str) -> set:
    Tree = _require_ete3()
    t = Tree(newick_str, format=1)
    return {lf.name for lf in t.iter_leaves() if lf.name}


def jaccard(a: set, b: set) -> float:
    u = len(a | b)
    return 0.0 if u == 0 else len(a & b) / u


def audit_multiset_lines(
    newicks: List[str],
    lo: float,
    hi: float,
    min_shared: int,
    anchors: Optional[List[str]] = None
) -> Dict:
    anchors_set = set(anchors or [])
    leafsets = [compute_leafset(n) for n in newicks]
    n = len(leafsets)

    sizes = [len(s) for s in leafsets]
    size_min = min(sizes) if sizes else 0
    size_max = max(sizes) if sizes else 0

    anchors_ok = None
    if anchors_set:
        anchors_ok = sum(1 for s in leafsets if anchors_set.issubset(s))

    overlaps = []
    inter_sizes = []
    out_of_range = 0
    below_min_shared = 0

    for i in range(n):
        for j in range(i + 1, n):
            A, B = leafsets[i], leafsets[j]
            inter = len(A & B)
            inter_sizes.append(inter)
            if inter < min_shared:
                below_min_shared += 1
            ov = jaccard(A, B)
            overlaps.append(ov)
            if ov < lo or ov > hi:
                out_of_range += 1

    total_pairs = n * (n - 1) // 2
    return {
        "n_trees": n,
        "sizes": {"min": size_min, "max": size_max},
        "overlap": {"min": min(overlaps) if overlaps else 0.0, "max": max(overlaps) if overlaps else 0.0},
        "intersection": {"min": min(inter_sizes) if inter_sizes else 0, "max": max(inter_sizes) if inter_sizes else 0},
        "out_of_range_pairs": out_of_range,
        "below_min_shared_pairs": below_min_shared,
        "trees_with_all_anchors": anchors_ok,
        "total_pairs": total_pairs,
    }


# Mode "service"

def run_mode_service(args: argparse.Namespace) -> None:
    """
    Delegate to the existing Mode A script for downloading pruned tree sets
    and assembling a combined dataset.
    """
    modeA = SCRIPT_DIR / "overlap_tree_pipeline_modeA.py"
    if not modeA.exists():
        raise FileNotFoundError(
            f"Missing {modeA}. Place overlap_tree_pipeline_modeA.py alongside this script."
        )

    cmd = [
        sys.executable,
        str(modeA),
        args.species_group,
        str(args.n),
        str(args.number_of_trees),
        args.email,
        "--selection_mode",
        args.selection_mode,
        "--stratify_by",
        args.stratify_by,
    ]

    if args.species_list_file:
        cmd += ["--species_list_file", args.species_list_file]
    if args.max_per_stratum is not None:
        cmd += ["--max_per_stratum", str(args.max_per_stratum)]
    if args.seed is not None:
        cmd += ["--seed", str(args.seed)]
    if args.prepare_only:
        cmd += ["--prepare_only"]
    if args.use_existing_nexus:
        cmd += ["--use_existing_nexus", args.use_existing_nexus]

    print("\n[MODE service] Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)



# Mode "reference" (reference tree + controlled noise)

@dataclass
class ReferenceConfig:
    group: str
    outdir: Path
    all_species_csv: Path

    # Base lists
    base_sizes: List[int]
    seed: int
    k_strata: int
    reuse_fraction: float
    max_per_species: int

    # Reference trees
    full_trees_file: Path
    base_lists_xlsx: Optional[Path]
    base_sheet: Optional[str]
    base_lists_csv: Optional[Path]

    # Overlap set generation
    n_input_trees: int
    pairwise_overlap_range: Tuple[float, float]
    min_shared_leaves_per_pair: int
    enforce_full_coverage: bool
    anchor_taxa_count: int
    prune_mode: str
    clade_selection_bias: float
    leaf_prune_blockiness: float
    contract_degree2: bool

    # Target sizes
    target_size_frac: Tuple[float, float]

    # Noise
    topology_noise: str
    nni_moves: int
    protect_anchors_in_noise: bool
    length_scaling: str
    length_scale_range: Tuple[float, float]
    renormalize_root_height: str


def parse_sizes(spec: str) -> List[int]:
    spec = spec.strip()
    if ":" in spec:
        a, b, c = map(int, spec.split(":"))
        return list(range(a, b + 1, c))
    return [int(x) for x in spec.split(",") if x.strip()]


def load_base_lists_from_xlsx(path: Path, sheet: str) -> List[List[str]]:
    xl = pd.ExcelFile(path)
    if sheet not in xl.sheet_names:
        matches = [s for s in xl.sheet_names if s.lower() == sheet.lower()]
        if matches:
            sheet = matches[0]
        else:
            raise ValueError(f"Sheet '{sheet}' not found in {path}. Sheets: {xl.sheet_names}")
    df = xl.parse(sheet)

    base_cols = [
        c for c in df.columns
        if re.match(r"^base[_\s-]?\d+$", str(c).strip(), flags=re.IGNORECASE)
    ]
    if not base_cols:
        raise ValueError("No Base_* columns found in the provided workbook sheet.")

    def _idx(col):
        m = re.search(r"(\d+)$", str(col).strip())
        return int(m.group(1)) if m else 0

    base_cols.sort(key=_idx)

    bases = []
    for c in base_cols:
        raw = df[c].dropna().astype(str).tolist()
        names = [normalize_spaces_to_underscores(x) for x in raw if str(x).strip()]
        names = list(dict.fromkeys(names))
        bases.append(names)
    return bases


def load_base_lists_from_csv(path: Path) -> List[List[str]]:
    df = pd.read_csv(path)
    base_cols = [
        c for c in df.columns
        if re.match(r"^base[_\s-]?\d+$", str(c).strip(), flags=re.IGNORECASE)
    ]
    if not base_cols:
        raise ValueError("No Base_* columns found in the provided CSV.")

    def _idx(col):
        m = re.search(r"(\d+)$", str(col).strip())
        return int(m.group(1)) if m else 0

    base_cols.sort(key=_idx)

    bases = []
    for c in base_cols:
        raw = df[c].dropna().astype(str).tolist()
        names = [normalize_spaces_to_underscores(x) for x in raw if str(x).strip()]
        names = list(dict.fromkeys(names))
        bases.append(names)
    return bases


def generate_phylo_stratified_bases(
    full_tree_newick: str,
    universe: List[str],
    sizes: List[int],
    seed: int,
    k_strata: int,
    reuse_fraction: float,
    max_per_species: int
) -> List[List[str]]:

    Tree, _gp, strat = _require_reference_modules()

    # Convert universe to tree naming convention (underscores)
    uni_norm = [normalize_spaces_to_underscores(x) for x in universe]
    uni_norm = list(dict.fromkeys(uni_norm))

    full_tree = Tree(full_tree_newick, format=1)
    tree_leaves = set(full_tree.get_leaf_names())

    # Restrict to taxa present in the tree
    uni_norm = [x for x in uni_norm if x in tree_leaves]
    if len(uni_norm) < 2:
        raise ValueError("Universe intersection with full tree leaves has <2 taxa.")

    strata_nodes = strat.make_strata(full_tree, k_strata)
    strata_leaves = strat.strata_leaf_sets(strata_nodes)

    bases, _used_count, _strata_counts_df = strat.build_phylo_stratified_bases(
        universe=uni_norm,
        strata_leaves=strata_leaves,
        sizes=sizes,
        seed=seed,
        max_per_species=max_per_species,
        reuse_fraction=reuse_fraction,
    )
    return bases


def write_bases_csv(path: Path, bases: List[List[str]]) -> None:
    max_len = max(len(b) for b in bases) if bases else 0
    cols = {}
    for i, b in enumerate(bases, start=1):
        cols[f"Base_{i:02d}"] = b + [""] * (max_len - len(b))
    pd.DataFrame(cols).to_csv(path, index=False)


def run_mode_reference(cfg: ReferenceConfig) -> None:
    Tree, gp, strat = _require_reference_modules()
    ensure_dir(cfg.outdir)

    group_dir = cfg.outdir / cfg.group
    ensure_dir(group_dir)

    paths = {
        "bases_dir": group_dir / "base_taxa_lists",
        "ref_dir": group_dir / "reference_trees",
        "multisets_dir": group_dir / "input_multisets",
        "meta_dir": group_dir / "metadata",
        "logs_dir": group_dir / "logs",
    }
    for p in paths.values():
        ensure_dir(p)

    # Load species universe from CSV
    universe = read_wide_species_csv(str(cfg.all_species_csv), cfg.group)
    universe = [normalize_spaces_to_underscores(x) for x in universe]

    # Load full trees (one or many)
    full_newicks = split_newick_file(str(cfg.full_trees_file))
    if not full_newicks:
        raise ValueError(f"No trees found in {cfg.full_trees_file}")

    # Base taxa lists
    bases: List[List[str]]
    bases_source = "generated"
    if cfg.base_lists_xlsx:
        sheet = cfg.base_sheet or cfg.group
        bases = load_base_lists_from_xlsx(cfg.base_lists_xlsx, sheet)
        bases_source = f"xlsx:{cfg.base_lists_xlsx.name}:{sheet}"
    elif cfg.base_lists_csv:
        bases = load_base_lists_from_csv(cfg.base_lists_csv)
        bases_source = f"csv:{cfg.base_lists_csv.name}"
    else:
        # generate from the FIRST full tree (strata construction uses one tree)
        bases = generate_phylo_stratified_bases(
            full_tree_newick=full_newicks[0],
            universe=universe,
            sizes=cfg.base_sizes,
            seed=cfg.seed,
            k_strata=cfg.k_strata,
            reuse_fraction=cfg.reuse_fraction,
            max_per_species=cfg.max_per_species,
        )

    # Write bases to disk (CSV)
    bases_csv = paths["bases_dir"] / f"{cfg.group}_base_taxa_lists.csv"
    write_bases_csv(bases_csv, bases)
    
    # Align base lists with full trees
    n_bases = len(bases)
    n_full = len(full_newicks)

    if n_full == 1 and n_bases >= 1:
        n_ref = n_bases
        bases_used = bases
        full_used = [full_newicks[0]] * n_ref
        print(f"\n=== {cfg.group.upper()} ===  ({n_ref} reference trees)")
        if n_bases > 1:
            print(f"[INFO] Using a single starting tree for all {n_bases} base taxa lists.")
    elif n_bases == 1 and n_full >= 1:
        n_ref = n_full
        bases_used = [bases[0]] * n_ref
        full_used = full_newicks
        print(f"\n=== {cfg.group.upper()} ===  ({n_ref} reference trees)")
        if n_full > 1:
            print(f"[INFO] Using one base taxa list for all {n_full} starting trees.")
    elif n_bases == n_full:
        n_ref = n_bases
        bases_used = bases
        full_used = full_newicks
        print(f"\n=== {cfg.group.upper()} ===  ({n_ref} reference trees)")
    else:
        n_ref = min(n_bases, n_full)
        bases_used = bases[:n_ref]
        full_used = full_newicks[:n_ref]
        print(f"\n=== {cfg.group.upper()} ===  ({n_ref} reference trees)")
        print(f"[WARN] Bases ({n_bases}) and starting trees ({n_full}) differ. Processing first {n_ref}.")

    # Prune starting trees to bases -> reference trees
    reference_newicks: List[str] = []
    skipped_ref = 0
    for i in range(n_ref):
        ref = prune_newick(full_used[i], bases_used[i])
        if ref is None:
            skipped_ref += 1
            reference_newicks.append("")
        else:
            reference_newicks.append(ref)

    # Write reference trees as individual files + a combined file
    combined_ref_path = paths["ref_dir"] / f"{cfg.group}_reference_trees{n_ref}.txt"
    write_newick_list(combined_ref_path, [x for x in reference_newicks if x])

    for i, nwk in enumerate(reference_newicks, start=1):
        if not nwk:
            continue
        (paths["ref_dir"] / f"reference_{i:02d}.nwk").write_text(
            nwk + ("\n" if not nwk.endswith("\n") else ""), encoding="utf-8"
        )

    if skipped_ref:
        print(f"[WARN] Skipped {skipped_ref} reference trees due to <2 taxa after intersection.")

    # Prepare generator base config (gen_pruned_trees)
    base_cfg = {
        "random_seed": int(cfg.seed),
        "n_trees": int(cfg.n_input_trees),

        "overlap_metric": "jaccard",
        "pairwise_overlap_range": [float(cfg.pairwise_overlap_range[0]), float(cfg.pairwise_overlap_range[1])],
        "min_shared_leaves_per_pair": int(cfg.min_shared_leaves_per_pair),
        "per_leaf_min_coverage": 1,
        "enforce_full_coverage": bool(cfg.enforce_full_coverage),
        "min_leaves_per_tree": 2,  # overwritten per reference tree
        "prune_mode": str(cfg.prune_mode),
        "clade_selection_bias": float(cfg.clade_selection_bias),
        "clade_size_pref": "by_size",
        "leaf_prune_blockiness": float(cfg.leaf_prune_blockiness),
        "contract_degree2": bool(cfg.contract_degree2),

        "anchor_taxa_count": int(cfg.anchor_taxa_count),
        "preserve_paths_between_anchors": False,

        "topology_noise": str(cfg.topology_noise),
        "nni_moves": int(cfg.nni_moves),
        "swap_fraction": 0.0,
        "protect_anchors_in_noise": bool(cfg.protect_anchors_in_noise),

        "length_scaling": str(cfg.length_scaling),
        "length_scale_params": {"low": float(cfg.length_scale_range[0]), "high": float(cfg.length_scale_range[1])},
        "renormalize_root_height": str(cfg.renormalize_root_height),

        # placeholders expected by module
        "base_tree_path": "__inline__",
        "out_dir": "__inline__",
    }
    base_cfg = gp.resolve_defaults(base_cfg)

    # For each reference tree, generate one multiset file
    for i in range(n_ref):
        ref_nwk = reference_newicks[i]
        if not ref_nwk:
            continue

        ref_tree = Tree(ref_nwk, format=1)

        # Ensure positive branch lengths
        for n in ref_tree.traverse():
            try:
                n.dist = max(float(getattr(n, "dist", 0.0)), 1e-8)
            except Exception:
                n.dist = 1e-8

        leaves = gp.get_leaf_names(ref_tree)
        base_n = len(leaves)
        if base_n < 2:
            continue

        # Per-reference-tree seed
        seed_here = int(cfg.seed) + 1000 * (i + 1)
        cfg_i = dict(base_cfg)
        cfg_i["random_seed"] = seed_here
        random.seed(seed_here)

        # Target size range as fractions of the reference tree taxa count
        lo_frac, hi_frac = cfg.target_size_frac
        lo = max(int(cfg.min_shared_leaves_per_pair) + int(cfg.anchor_taxa_count), int(lo_frac * base_n))
        hi = max(lo, int(hi_frac * base_n))
        cfg_i["target_tree_size_range"] = [lo, hi]
        cfg_i["min_leaves_per_tree"] = lo

        # Generate leaf sets satisfying overlap constraints
        leaf_sets, notes = gp.generate_sets_of_leaves(ref_tree, cfg_i, leaves)

        anchors = set(cfg_i.get("anchors_explicit", []))
        base_h = gp.root_height(ref_tree)

        K = int(cfg_i["n_trees"])
        lines: List[str] = []
        per_tree_meta: List[dict] = []

        for k in range(K):
            keep = leaf_sets[k]
            pruned = gp.prune_to_leaves(ref_tree, keep, contract_degree2=bool(cfg_i.get("contract_degree2", True)))
            gp.apply_topology_noise(pruned, cfg_i, anchors)
            gp.apply_length_scaling(pruned, cfg_i, anchors, base_h)

            nwk = pruned.write(format=1).strip()
            if not nwk.endswith(";"):
                nwk += ";"
            lines.append(nwk)

            per_tree_meta.append({
                "tree_index": k + 1,
                "n_taxa": len(keep),
            })

        out_path = paths["multisets_dir"] / f"multiset_{i+1}.txt"
        write_newick_list(out_path, lines)

        # Write dataset-level metadata
        meta = {
            "group": cfg.group,
            "dataset_index": i + 1,
            "bases_source": bases_source,
            "base_taxa_list_size": len(bases_used[i]),
            "base_taxa_list": bases_used[i],
            "reference_tree_taxa_count": base_n,
            "n_input_trees": K,
            "subset_size_range": [lo, hi],
            "pairwise_overlap_range": list(cfg.pairwise_overlap_range),
            "min_shared_leaves_per_pair": cfg.min_shared_leaves_per_pair,
            "anchor_taxa_count": cfg.anchor_taxa_count,
            "anchors_explicit": sorted(list(anchors)),
            "prune_mode": cfg.prune_mode,
            "topology_noise": cfg.topology_noise,
            "nni_moves": cfg.nni_moves,
            "length_scaling": cfg.length_scaling,
            "length_scale_range": list(cfg.length_scale_range),
            "seed_base": cfg.seed,
            "seed_dataset": seed_here,
            "notes": notes,
            "input_trees": per_tree_meta,
            "paths": {
                "reference_tree": str((paths["ref_dir"] / f"reference_{i+1:02d}.nwk").resolve()),
                "multiset_file": str(out_path.resolve()),
                "base_taxa_lists_csv": str(bases_csv.resolve()),
            },
        }
        meta_path = paths["meta_dir"] / f"dataset_{i+1:02d}.json"
        meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")

        # Optional audit
        audit = audit_multiset_lines(
            lines,
            lo=float(cfg.pairwise_overlap_range[0]),
            hi=float(cfg.pairwise_overlap_range[1]),
            min_shared=int(cfg.min_shared_leaves_per_pair),
            anchors=sorted(list(anchors)) if anchors else None
        )
        audit_path = paths["logs_dir"] / f"audit_{i+1:02d}.json"
        audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")

        if (i + 1) % 5 == 0 or (i + 1) == n_ref:
            print(
                f"  wrote {out_path.name}  [{K} trees]  "
                f"(subset size range: {lo}-{hi}, anchors: {len(anchors)})"
            )

    print(f"[DONE] Wrote {n_ref} multiset files under {paths['multisets_dir']}")
    print(f"[DONE] Reference trees under {paths['ref_dir']}")
    print(f"[DONE] Metadata under {paths['meta_dir']}")



# CLI

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Unified overlap dataset pipeline with two modes: service and reference."
    )
    sub = ap.add_subparsers(dest="mode", required=True)

    # service
    apA = sub.add_parser(
        "service",
        help="VertLife PhyloSubsets mode (download pruned tree sets)."
    )
    apA.add_argument("species_group", choices=["amphibians", "birds", "mammals", "sharks", "squamates"])
    apA.add_argument("n", type=int)
    apA.add_argument("number_of_trees", type=int)
    apA.add_argument("email", type=str)
    apA.add_argument("--selection_mode", choices=["user_list", "random", "stratified"], default="stratified")
    apA.add_argument("--species_list_file", type=str, default=None)
    apA.add_argument("--stratify_by", choices=["genus"], default="genus")
    apA.add_argument("--max_per_stratum", type=int, default=None)
    apA.add_argument("--seed", type=int, default=None)
    apA.add_argument(
        "--prepare_only",
        action="store_true",
        help="Generate selected_species.csv and <group>_overlapping_subsets.csv, then stop before VertLife submission."
    )
    apA.add_argument(
        "--use_existing_nexus",
        type=str,
        default=None,
        help="Skip VertLife submission and use an existing folder of Nexus files to finish the dataset."
    )

    # Reference
    apB = sub.add_parser(
        "reference",
        help="Reference-tree benchmarking mode (local full trees -> reference + overlapping inputs)."
    )
    apB.add_argument("--group", required=True, choices=["amphibians", "birds", "mammals", "sharks", "squamates"])
    apB.add_argument(
        "--full_trees",
        required=True,
        help="Path to a file containing one or more full-tree Newicks (semicolon-separated)."
    )

    apB.add_argument("--outdir", default="./benchmark_datasets", help="Output root directory.")
    apB.add_argument("--all_species_csv", default="all_species_lists.csv", help="Wide CSV with group columns.")

    # base list options
    apB.add_argument("--base_sizes", default="50:145:5", help="Base taxa list sizes as start:end:step or comma list.")
    apB.add_argument("--seed", type=int, default=7)

    apB.add_argument("--k_strata", type=int, default=12, help="Phylogenetic strata count for base list generation.")
    apB.add_argument("--reuse_fraction", type=float, default=0.15)
    apB.add_argument("--max_per_species", type=int, default=6)

    apB.add_argument("--base_lists_xlsx", default=None, help="Optional workbook with Base_01.. columns.")
    apB.add_argument("--base_sheet", default=None, help="Sheet name to read from base_lists_xlsx. Default is the group name.")
    apB.add_argument("--base_lists_csv", default=None, help="Optional CSV with Base_01.. columns.")

    # Overlap set generation
    apB.add_argument("--n_input_trees", type=int, default=30)
    apB.add_argument("--pairwise_overlap_range", type=float, nargs=2, default=[0.30, 0.70])
    apB.add_argument("--min_shared", type=int, default=2)
    apB.add_argument("--enforce_full_coverage", action="store_true", default=True)
    apB.add_argument("--no_full_coverage", action="store_false", dest="enforce_full_coverage")

    apB.add_argument("--anchor_taxa_count", type=int, default=10)
    apB.add_argument("--prune_mode", choices=["leaves", "clades", "mixed"], default="mixed")
    apB.add_argument("--clade_selection_bias", type=float, default=0.6)
    apB.add_argument("--leaf_prune_blockiness", type=float, default=0.3)
    apB.add_argument("--no_contract_degree2", action="store_true", help="Disable degree-2 contraction after pruning.")

    # Subset size fractions relative to reference tree
    apB.add_argument(
        "--target_size_frac",
        type=float,
        nargs=2,
        default=[0.50, 0.75],
        help="Low and high fractions of taxa retained per input tree."
    )

    # Noise controls
    apB.add_argument("--topology_noise", choices=["none", "swap_labels", "nni"], default="nni")
    apB.add_argument("--nni_moves", type=int, default=3)
    apB.add_argument("--no_protect_anchors", action="store_true", help="Allow anchors to move during topology noise.")

    apB.add_argument("--length_scaling", choices=["none", "global", "global_uniform"], default="global_uniform")
    apB.add_argument("--length_scale_range", type=float, nargs=2, default=[1.003, 1.009])
    apB.add_argument("--renormalize_root_height", choices=["none", "to_base"], default="none")

    return ap


def main() -> None:
    ap = build_parser()
    args = ap.parse_args()

    if args.mode == "service":
        run_mode_service(args)
        return

    if args.mode == "reference":
        cfg = ReferenceConfig(
            group=args.group.lower(),
            outdir=Path(args.outdir),
            all_species_csv=Path(args.all_species_csv),

            base_sizes=parse_sizes(args.base_sizes),
            seed=int(args.seed),
            k_strata=int(args.k_strata),
            reuse_fraction=float(args.reuse_fraction),
            max_per_species=int(args.max_per_species),

            full_trees_file=Path(args.full_trees),
            base_lists_xlsx=Path(args.base_lists_xlsx) if args.base_lists_xlsx else None,
            base_sheet=args.base_sheet,
            base_lists_csv=Path(args.base_lists_csv) if args.base_lists_csv else None,

            n_input_trees=int(args.n_input_trees),
            pairwise_overlap_range=(float(args.pairwise_overlap_range[0]), float(args.pairwise_overlap_range[1])),
            min_shared_leaves_per_pair=int(args.min_shared),
            enforce_full_coverage=bool(args.enforce_full_coverage),
            anchor_taxa_count=int(args.anchor_taxa_count),
            prune_mode=str(args.prune_mode),
            clade_selection_bias=float(args.clade_selection_bias),
            leaf_prune_blockiness=float(args.leaf_prune_blockiness),
            contract_degree2=not bool(args.no_contract_degree2),

            target_size_frac=(float(args.target_size_frac[0]), float(args.target_size_frac[1])),

            topology_noise=str(args.topology_noise),
            nni_moves=int(args.nni_moves),
            protect_anchors_in_noise=not bool(args.no_protect_anchors),
            length_scaling=str(args.length_scaling),
            length_scale_range=(float(args.length_scale_range[0]), float(args.length_scale_range[1])),
            renormalize_root_height=str(args.renormalize_root_height),
        )
        run_mode_reference(cfg)
        return

    raise ValueError(f"Unknown mode: {args.mode}")

if __name__ == "__main__":
    main()
