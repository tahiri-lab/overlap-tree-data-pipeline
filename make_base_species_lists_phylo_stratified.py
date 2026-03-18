#!/usr/bin/env python3
"""
Phylogenetically stratified base-taxa sampler

Steps:
1) Read one full Newick tree and get its leaf set (taxa universe).
2) Create K non-overlapping phylogenetic strata (each stratum is a clade) using
   a greedy top-down partition:
   - Start with the full tree as 1 bin.
   - Repeatedly split the currently largest bin into its two child clades
     until there are K bins.
   - If the tree has polytomies, they are first resolved to binary
     (adds 0-length structure) so each split increases bins by exactly 1.
3) For each desired base size N (e.g., 50,55,...,145):
   - Allocate N picks across the K strata to be roughly equal and guarantee
     at least 1 pick per stratum when N >= K.
   - Optionally seed a fraction of taxa from a previous base to maintain
     non-zero overlaps (reuse_fraction).
   - Fill remaining slots by weighted sampling that favors underused taxa across
     all bases, while respecting a soft reuse cap (max_per_species).
   - If any stratum runs out of eligible taxa, redistribute the deficit to
     other strata with remaining candidates.

Outputs:
- An Excel workbook with:
  - a base-list sheet (default prefix: "Bases")
  - a "<prefix>_Jaccard" sheet with pairwise overlaps between bases
  - a "<prefix>_Reuse" sheet with reuse-count summary stats
  - a "<prefix>_Strata" sheet with the K strata (leaf lists per column)
  - a "<prefix>_StrataCounts" sheet with per-base counts by stratum

Dependencies:
- ete3 (for Tree)
- pandas, numpy
- xlsxwriter (used by pandas ExcelWriter engine)

"""

import argparse
from collections import defaultdict
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from ete3 import Tree


# Utilities

def parse_sizes(spec: str) -> List[int]:
    """Parse sizes like '50:145:5' -> [50,55,...,145]. Or comma list: '50,60,75'."""
    spec = spec.strip()
    if ":" in spec:
        start, end, step = map(int, spec.split(":"))
        return list(range(start, end + 1, step))
    return [int(x) for x in spec.split(",") if x.strip()]


def load_first_newick(path: str) -> str:
    """Load the first Newick tree from a file that may contain multiple trees."""
    with open(path, "r", encoding="utf-8") as f:
        txt = f.read().strip()
    parts = [p.strip() for p in txt.split(";") if p.strip()]
    if not parts:
        raise ValueError(f"No Newick trees found in {path}")
    return parts[0] + ";"


def normalize_name_to_tree(name: str) -> str:
    """Convert 'Genus species' -> 'Genus_species' and normalize whitespace."""
    name = str(name).strip()
    name = " ".join(name.split())
    return name.replace(" ", "_")


def display_name_from_tree(name: str) -> str:
    """Convert 'Genus_species' -> 'Genus species' for nicer Excel display."""
    return str(name).replace("_", " ").strip()


def make_excel_safe_sheet_name(name: str, max_len: int = 31) -> str:
    """
    Improve an Excel sheet name:
    - remove invalid characters []:*?/\\
    - trim whitespace
    - limit length to 31 characters
    """
    invalid = set('[]:*?/\\')
    cleaned = "".join(ch for ch in str(name) if ch not in invalid).strip()
    if not cleaned:
        cleaned = "Sheet"
    return cleaned[:max_len]


def sheet_name(prefix: str, suffix: str = "") -> str:
    """Compose a safe Excel sheet name from a prefix and optional suffix."""
    base = prefix if not suffix else f"{prefix}_{suffix}"
    return make_excel_safe_sheet_name(base)


def weighted_choice_without_replacement(
    candidates: List[str], weights: np.ndarray, k: int, rng: np.random.Generator
) -> List[str]:
    """Sample k distinct items from candidates according to weights (non-negative)."""
    if k <= 0 or not candidates:
        return []
    w = np.asarray(weights, dtype=float)
    w[w < 0] = 0
    k = min(k, len(candidates))
    if w.sum() <= 0:
        idx = rng.choice(len(candidates), size=k, replace=False)
        return [candidates[i] for i in idx]
    p = w / w.sum()
    idx = rng.choice(len(candidates), size=k, replace=False, p=p)
    return [candidates[i] for i in idx]


def jaccard(a: List[str], b: List[str]) -> float:
    sa, sb = set(a), set(b)
    return 0.0 if not sa and not sb else len(sa & sb) / len(sa | sb)


def jaccard_matrix(bases: List[List[str]]) -> pd.DataFrame:
    n = len(bases)
    M = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            M[i, j] = 1.0 if i == j else jaccard(bases[i], bases[j])
    idx = [f"Base_{i+1:02d}" for i in range(n)]
    return pd.DataFrame(M, index=idx, columns=idx)


def to_wide_df(bases: List[List[str]], display: bool = True) -> pd.DataFrame:
    """Bases -> columns Base_01..Base_N, padded with blanks."""
    max_len = max(len(b) for b in bases) if bases else 0
    data = {}
    for i, b in enumerate(bases):
        col = [display_name_from_tree(x) for x in b] if display else list(b)
        data[f"Base_{i+1:02d}"] = col + [""] * (max_len - len(col))
    return pd.DataFrame(data)


# Strata construction

def make_strata(tree: Tree, k: int) -> List[Tree]:
    """
    Greedy top-down partition into k clade-strata.
    Assumes binary splits (polytomies are resolved beforehand).
    """
    if k < 1:
        raise ValueError("k must be >= 1")

    t = tree.copy(method="deepcopy")
    t.resolve_polytomy(recursive=True)

    bins: List[Tree] = [t]
    while len(bins) < k:
        splittable = [n for n in bins if (not n.is_leaf()) and len(n.children) >= 2 and len(n) > 1]
        if not splittable:
            break
        node = max(splittable, key=lambda n: len(n))
        bins.remove(node)
        bins.extend(node.children)

    if len(bins) < k:
        raise ValueError(
            f"Could only create {len(bins)} strata (requested {k}). "
            f"Try a smaller --k-strata value."
        )

    return bins[:k]


def strata_leaf_sets(strata: List[Tree]) -> List[List[str]]:
    """Return leaf-name lists for each stratum."""
    return [s.get_leaf_names() for s in strata]


def build_stratum_index(strata_leaves: List[List[str]]) -> Dict[str, int]:
    """Map each species name to its stratum index."""
    idx = {}
    for i, leaves in enumerate(strata_leaves):
        for sp in leaves:
            idx[sp] = i
    return idx


# Sampling

def allocate_across_strata_equal(
    size: int, k: int, strata_sizes: List[int], ensure_min1: bool = True
) -> List[int]:
    """
    Allocate 'size' picks across k strata:
    - roughly equal
    - optionally ensure >=1 per stratum if size>=k
    - remainder goes first to larger strata
    """
    if k <= 0:
        return []

    if ensure_min1 and size >= k:
        alloc = [1] * k
        remaining = size - k
    else:
        alloc = [0] * k
        remaining = size

    if remaining <= 0:
        return alloc

    base = remaining // k
    rem = remaining % k
    alloc = [a + base for a in alloc]

    order = sorted(range(k), key=lambda i: strata_sizes[i], reverse=True)
    for i in order[:rem]:
        alloc[i] += 1
    return alloc


def pick_from_stratum(
    stratum_species: List[str],
    need: int,
    chosen: set,
    used_count: Dict[str, int],
    cap: int,
    rng: np.random.Generator,
) -> List[str]:
    """Pick taxa from one stratum, prioritizing underused taxa and respecting cap when possible."""
    if need <= 0:
        return []

    eligible = [s for s in stratum_species if s not in chosen and used_count[s] < cap]
    if len(eligible) < need:
        eligible = [s for s in stratum_species if s not in chosen]

    usage = np.array([used_count[s] for s in eligible], dtype=float)
    weights = 1.0 / (1.0 + usage)
    return weighted_choice_without_replacement(eligible, weights, need, rng)


def build_phylo_stratified_bases(
    universe: List[str],
    strata_leaves: List[List[str]],
    sizes: List[int],
    seed: int = 0,
    max_per_species: int = 6,
    reuse_fraction: float = 0.15,
) -> Tuple[List[List[str]], Dict[str, int], pd.DataFrame]:
    """
    Build bases of varying sizes using phylogenetic strata, controlled overlap,
    and underuse-weighted sampling with a soft reuse cap.
    """
    rng = np.random.default_rng(seed)

    seen = set()
    uni = []
    for s in universe:
        s = s.strip()
        if s and s not in seen:
            seen.add(s)
            uni.append(s)

    strata = []
    for leaves in strata_leaves:
        leafset = set(leaves)
        strata.append([s for s in uni if s in leafset])

    k = len(strata)
    strata_sizes = [len(x) for x in strata]
    sp_to_stratum = build_stratum_index(strata)

    used_count = defaultdict(int)
    bases: List[List[str]] = []
    counts_rows = []

    for base_i, size in enumerate(sizes, start=1):
        cap = max_per_species + (1 if base_i > 0.75 * len(sizes) else 0)

        chosen: List[str] = []
        chosen_set = set()

        if base_i > 1 and bases and reuse_fraction > 0:
            prev = bases[int(rng.integers(0, len(bases)))]
            seed_cnt = max(0, min(size - 1, int(round(reuse_fraction * min(size, len(prev))))))
            if seed_cnt > 0:
                seed_taxa = rng.choice(prev, size=seed_cnt, replace=False).tolist()
                seed_taxa = [s for s in seed_taxa if s in seen]
                chosen.extend(seed_taxa)
                chosen_set.update(seed_taxa)

        alloc = allocate_across_strata_equal(size=size, k=k, strata_sizes=strata_sizes, ensure_min1=True)

        seeded_counts = [0] * k
        for s in chosen:
            si = sp_to_stratum.get(s, None)
            if si is not None:
                seeded_counts[si] += 1
        need_per = [max(0, alloc[i] - seeded_counts[i]) for i in range(k)]

        deficit = 0
        for i in range(k):
            take = pick_from_stratum(strata[i], need_per[i], chosen_set, used_count, cap, rng)
            chosen.extend(take)
            chosen_set.update(take)
            deficit += max(0, need_per[i] - len(take))

        if deficit > 0:
            pool = [s for s in uni if s not in chosen_set and used_count[s] < cap]
            usage = np.array([used_count[s] for s in pool], dtype=float)
            weights = 1.0 / (1.0 + usage)
            extra = weighted_choice_without_replacement(pool, weights, deficit, rng)
            chosen.extend(extra)
            chosen_set.update(extra)

        if len(chosen) < size:
            pool = [s for s in uni if s not in chosen_set]
            extra = weighted_choice_without_replacement(pool, np.ones(len(pool)), size - len(chosen), rng)
            chosen.extend(extra)
            chosen_set.update(extra)

        chosen = chosen[:size]

        bases.append(chosen)
        for s in chosen:
            used_count[s] += 1

        row = {"Base": f"Base_{base_i:02d}", "Size": len(chosen)}
        comp = [0] * k
        for s in chosen:
            si = sp_to_stratum.get(s, None)
            if si is not None:
                comp[si] += 1
        for i in range(k):
            row[f"Stratum_{i+1:02d}"] = comp[i]
        counts_rows.append(row)

    strata_counts_df = pd.DataFrame(counts_rows)
    return bases, used_count, strata_counts_df


# Main

def strata_to_wide_df(strata_leaves: List[List[str]]) -> pd.DataFrame:
    """Write strata as columns Stratum_01..Stratum_K."""
    max_len = max(len(x) for x in strata_leaves) if strata_leaves else 0
    data = {}
    for i, leaves in enumerate(strata_leaves, start=1):
        col = [display_name_from_tree(x) for x in leaves]
        data[f"Stratum_{i:02d}"] = col + [""] * (max_len - len(col))
    return pd.DataFrame(data)


def read_wide_species_csv_group(path: str, group_col: str) -> List[str]:
    """
    Read one group column from a wide CSV (one column per group; rows are species names).
    Returns names normalized to the tree naming convention (underscores).
    """
    df = pd.read_csv(path)
    if group_col not in df.columns:
        raise ValueError(f"Column '{group_col}' not found in {path}. Columns: {list(df.columns)}")
    col = df[group_col].dropna().astype(str).tolist()
    out = []
    seen = set()
    for s in col:
        nm = normalize_name_to_tree(s)
        if nm and nm not in seen:
            seen.add(nm)
            out.append(nm)
    return out


# Backward-compatible alias for older imports
read_wide_species_csv_birds = read_wide_species_csv_group


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--tree",
        required=True,
        help="Path to a full Newick tree file (single tree recommended)",
    )
    ap.add_argument("--out", required=True, help="Output Excel path")
    ap.add_argument("--sizes", default="50:145:5", help="Sizes as start:end:step or comma list")
    ap.add_argument("--k-strata", type=int, default=12, help="Number of phylogenetic strata (clade bins)")
    ap.add_argument("--seed", type=int, default=20251104, help="RNG seed")

    ap.add_argument(
        "--reuse-fraction",
        type=float,
        default=0.15,
        help="Fraction of taxa to reuse from a previous base",
    )
    ap.add_argument(
        "--max-per-species",
        type=int,
        default=6,
        help="Soft per-species reuse cap across all bases",
    )

    ap.add_argument(
        "--csv",
        default=None,
        help="Optional wide CSV with one column per group to define allowed taxa",
    )
    ap.add_argument(
        "--group-col",
        default=None,
        help="Column name in the wide CSV to use when --csv is provided",
    )
    ap.add_argument(
        "--sheet-prefix",
        default=None,
        help="Optional prefix for Excel sheet names (defaults to the group column or 'Bases')",
    )

    # Hidden backward-compatible alias
    ap.add_argument("--birds-col", dest="group_col", help=argparse.SUPPRESS)

    args = ap.parse_args()

    sizes = parse_sizes(args.sizes)

    newick = load_first_newick(args.tree)
    full_tree = Tree(newick, format=1)

    tree_leaves = full_tree.get_leaf_names()
    tree_leaf_set = set(tree_leaves)

    if args.csv:
        if not args.group_col:
            raise ValueError(
                "When --csv is used, --group-col must be provided "
                "(for example: --group-col Amphibians)."
            )
        csv_group = read_wide_species_csv_group(args.csv, group_col=args.group_col)
        universe = [s for s in csv_group if s in tree_leaf_set]
        if len(universe) < 2:
            raise ValueError(
                f"After intersecting CSV column '{args.group_col}' with tree leaves, <2 taxa remain."
            )
    else:
        universe = tree_leaves[:]

    strata_nodes = make_strata(full_tree, args.k_strata)
    strata_leaves = strata_leaf_sets(strata_nodes)

    bases, used_count, strata_counts_df = build_phylo_stratified_bases(
        universe=universe,
        strata_leaves=strata_leaves,
        sizes=sizes,
        seed=args.seed,
        max_per_species=args.max_per_species,
        reuse_fraction=args.reuse_fraction,
    )

    J = jaccard_matrix(bases)
    reuse_summary = pd.Series(dict(used_count)).describe().to_frame("reuse_count")

    if args.sheet_prefix:
        prefix = args.sheet_prefix
    elif args.group_col:
        prefix = args.group_col
    else:
        prefix = "Bases"

    base_sheet = sheet_name(prefix)
    jaccard_sheet = sheet_name(prefix, "Jaccard")
    reuse_sheet = sheet_name(prefix, "Reuse")
    strata_sheet = sheet_name(prefix, "Strata")
    strata_counts_sheet = sheet_name(prefix, "StrataCounts")

    with pd.ExcelWriter(args.out, engine="xlsxwriter") as writer:
        to_wide_df(bases, display=True).to_excel(writer, sheet_name=base_sheet, index=False)
        J.to_excel(writer, sheet_name=jaccard_sheet)
        reuse_summary.to_excel(writer, sheet_name=reuse_sheet)
        strata_to_wide_df(strata_leaves).to_excel(writer, sheet_name=strata_sheet, index=False)
        strata_counts_df.to_excel(writer, sheet_name=strata_counts_sheet, index=False)

    M = J.values
    upper = M[np.triu_indices_from(M, k=1)]
    prefix_display = prefix
    print(
        f"[DONE] Wrote {len(bases)} bases to {args.out}\n"
        f"Group/sheet prefix: {prefix_display}\n"
        f"Sizes: {list(map(len, bases))}\n"
        f"Jaccard min={upper.min():.3f}, median={np.median(upper):.3f}, max={upper.max():.3f}\n"
        f"Strata (K) = {len(strata_leaves)}"
    )


if __name__ == "__main__":
    main()