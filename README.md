# Overlap dataset pipeline

## Overview

This project provides a reproducible pipeline for constructing phylogenetic tree datasets with partial taxon overlap.

It supports two complementary dataset-construction modes:

- **Mode 1 (`service`)** builds overlapping taxon subsets and retrieves pruned phylogenetic tree samples from the VertLife PhyloSubsets web service.
- **Mode 2 (`reference`)** builds benchmark datasets from one or more locally available broad-coverage phylogenetic trees by pruning them to create reference trees and then pruning again to create overlapping input trees, with optional controlled topology and branch-length perturbation.

The repository is intended for constructing empirical overlap-controlled datasets for tasks such as supertree construction, tree comparison, clustering, and related benchmarking workflows.

The repository contains:

- **`overlap_tree_pipeline.py`**  
  Unified script for both modes

- **`overlap_tree_pipeline_modeA.py`**  
  Mode 1 implementation for generating 10 overlapping subsets and retrieving VertLife tree samples

- **`gen_pruned_trees.py`**  
  Core Mode 2 tree-pruning and benchmark-generation logic

- **`make_base_species_lists_phylo_stratified.py`**  
  Helper for generating phylogeny-stratified base taxa lists for Mode 2

- **`all_species_lists.csv`**  
  Species list file used by the pipeline

- **`datasets-mode1/`**  
  Datasets constructed using mode 1

- **`datasets-mode2/`**  
  Datasets constructed using mode 2

The unified script runs the appropriate workflow depending on the selected mode and passes parameters through to the underlying implementation.

---

## How the pipeline works

### Mode 1: Overlap tree datasets

Mode 1 constructs a base set of species, generates exactly 10 overlapping taxon subsets, submits these subsets to the VertLife PhyloSubsets service, downloads the resulting pruned Nexus files, samples trees from each subset, and combines the sampled trees into one Newick dataset file. The script enforces equal tree sampling across the 10 subsets, which is why the requested total number of trees must be divisible by 10.

### Mode 2: reference-tree benchmark datasets

Mode 2 starts from one or more local full-tree Newick strings, creates reference trees by pruning to selected base taxa lists, generates overlapping input-tree multisets under user-defined overlap constraints, optionally applies topology and branch-length perturbations, and writes dataset metadata and audit logs alongside the generated trees. If explicit base lists are not supplied, they are generated from the first full tree using a phylogeny-stratified procedure.

---

## Installation

### Python

Python 3.9+ is recommended. 

### Required packages

Install the main dependencies with:

```bash
pip install pandas requests biopython selenium ete3 openpyxl
```

### Browser requirements for Mode 1

Mode 1 uses Selenium to submit jobs through the VertLife web form. It requires:

* Chrome or Chromium
* ChromeDriver available on `PATH`
* A ChromeDriver version matching the installed browser major version 

---

## Input data

### Species list file

The pipeline expects `all_species_lists.csv` in wide format, with one column per group:

* `Amphibians`
* `Birds`
* `Mammals`
* `Sharks`
* `Squamates`

Each column contains one species name per row.

Mode 1 expects this file in the current working directory. Mode 2 uses the same file by default, but the path can be overridden with `--all_species_csv`.

### Full-tree file for Mode 2

Mode 2 requires a file containing one or more Newick trees. Multiple trees may be stored in one file, separated by semicolons. The pipeline splits the file on semicolons and treats each Newick string as a candidate starting tree.

Example:

```text
(full_tree_1_newick...);
(full_tree_2_newick...);
...
```

### Taxon naming convention

Mode 2 normalizes taxa names by converting spaces to underscores when reading species lists. For best results, the leaf labels in the full-tree file should use underscore-separated names as well. If the full trees use spaces, either convert the leaf labels before running the pipeline or adjust the normalization logic in the script.

---

## How to run the pipeline

### General command structure

```bash
python overlap_tree_pipeline.py <mode> [options]
```

Available modes:

* `service`
* `reference`

Help can be shown with:

```bash
python overlap_tree_pipeline.py --help
python overlap_tree_pipeline.py service --help
python overlap_tree_pipeline.py reference --help
```

---

## Mode 1: building overlap datasets from VertLife (`service`)

### Purpose

This mode constructs 10 overlapping taxon subsets and retrieves pruned tree samples for each subset from VertLife. It is suitable for building empirical overlap-controlled datasets with branch lengths preserved in the downloaded trees.

### Example

```bash
python overlap_tree_pipeline.py service amphibians 120 550 name@example.com --seed 7
```

### Main positional arguments

* `species_group`
  One of `amphibians`, `birds`, `mammals`, `sharks`, `squamates`

* `n`
  Size of the base species set; ignored when `--selection_mode user_list`

* `number_of_trees`
  Total number of output trees in the final combined dataset

* `email`
  Email address required by VertLife

### Main optional arguments

* `--selection_mode {user_list,random,stratified}`
  Default: `stratified`

* `--species_list_file PATH`
  Used with `user_list`

* `--stratify_by genus`
  Currently only genus is supported

* `--max_per_stratum INT`
  Optional per-genus cap in stratified sampling

* `--seed INT`
  Random seed for reproducible species selection and tree sampling

### Constraints

* `number_of_trees` must be divisible by 10
* `number_of_trees` must be between 10 and 1000 inclusive 

### Mode 1 outputs

Mode 1 writes to the current working directory:

* `selected_species.csv`
* `<group>_overlapping_subsets.csv`
* `<group>_nexus/`
* `overlapping_dataset_<group>.txt`

---

## Mode 2: building reference-tree benchmark datasets (`reference`)

### Purpose

This mode constructs benchmark datasets that include:

* a fixed reference tree
* a multiset of overlapping input trees
* metadata describing the dataset construction
* audit logs for validation and reproducibility

It is intended for benchmarking methods against a known target under controlled overlap and optional perturbation (noise).

### Example

```bash
python overlap_tree_pipeline.py reference \
  --group amphibians \
  --full_trees full_trees/amphibians_full_trees20.txt \
  --outdir datasets-mode2 \
  --base_sizes 50:145:5 \
  --n_input_trees 30 \
  --pairwise_overlap_range 0.30 0.70 \
  --anchor_taxa_count 10 \
  --nni_moves 3 \
  --length_scale_range 1.003 1.009 \
  --seed 7
```

### Core inputs

* `--group`
  One of `amphibians`, `birds`, `mammals`, `sharks`, `squamates`

* `--full_trees PATH`
  File containing one or more semicolon-separated Newick trees

* `--outdir PATH`
  Output root directory
  Default: `./benchmark_datasets`

* `--all_species_csv PATH`
  Species-list CSV
  Default: `all_species_lists.csv`

### Base taxa list configuration

Generated base taxa lists can be configured with:

* `--base_sizes SPEC`
  Format `start:end:step` or comma-separated list
  Default: `50:145:5`

* `--seed INT`
  Default: `7`

Instead of generating them, explicit base lists may be supplied with:

* `--base_lists_xlsx PATH`
* `--base_sheet NAME`
* `--base_lists_csv PATH`

If explicit base lists are not supplied, the pipeline generates them from the first full tree using the phylogeny-stratified helper script.

### Phylogeny-stratified base generation

When base lists are generated automatically, the following controls are available:

* `--k_strata INT`
  Default: `12`

* `--reuse_fraction FLOAT`
  Default: `0.15`

* `--max_per_species INT`
  Default: `6`

### Overlap-set generation

The overlapping input-tree sets are controlled through:

* `--n_input_trees INT`
  Number of input trees per dataset
  Default: `30`

* `--pairwise_overlap_range LO HI`
  Jaccard overlap bounds
  Default: `0.30 0.70`

* `--min_shared INT`
  Minimum shared taxa per tree pair
  Default: `2`

* `--enforce_full_coverage` / `--no_full_coverage`
  Full coverage is enabled by default

* `--anchor_taxa_count INT`
  Number of anchor taxa included in every input tree
  Default: `10` 

### Pruning behavior

Tree pruning can be configured with:

* `--prune_mode {leaves,clades,mixed}`
  Default: `mixed`

* `--clade_selection_bias FLOAT`
  Default: `0.6`

* `--leaf_prune_blockiness FLOAT`
  Default: `0.3`

* `--no_contract_degree2`
  Disables contraction of degree-2 nodes after pruning

### Input-tree size range

* `--target_size_frac LO HI`
  Fraction of reference-tree taxa retained per input tree
  Default: `0.50 0.75` 

### Noise controls

Topology and branch-length perturbation can be configured with:

* `--topology_noise {none,swap_labels,nni}`
  Default: `nni`

* `--nni_moves INT`
  Default: `3`

* `--no_protect_anchors`
  Allows anchor taxa to move during topology perturbation

* `--length_scaling {none,global,global_uniform}`
  Default: `global_uniform`

* `--length_scale_range LO HI`
  Default: `1.003 1.009`

* `--renormalize_root_height {none,to_base}`
  Default: `none` 

### Mode 2 outputs

When `--outdir datasets-mode2` is used, the output structure is:

```text
datasets-mode2/
  amphibians/
    base_taxa_lists/
      amphibians_base_taxa_lists.csv
    reference_trees/
      reference_01.nwk
      reference_02.nwk
      ...
      amphibians_reference_trees20.txt
    input_multisets/
      multiset_1.txt
      multiset_2.txt
      ...
    metadata/
      dataset_01.json
      dataset_02.json
      ...
    logs/
      audit_01.json
      audit_02.json
      ...
```

Main output files:

* `reference_XX.nwk`
  Reference tree for dataset `XX`

* `<group>_reference_trees<N>.txt`
  All successfully created reference trees, one per line

* `multiset_X.txt`
  Overlapping input trees for dataset `X`, one Newick tree per line

* `dataset_XX.json`
  Dataset metadata including parameters, seeds, taxa counts, and file paths

* `audit_XX.json`
  Audit statistics including overlap extrema, shared-taxa extrema, total pair counts, and anchor-coverage summaries

If a full tree and a base taxa list intersect in fewer than 2 taxa, that dataset is skipped. 

---

## Configuration summary

The most important configuration settings in this project are:

### For Mode 1

* species group
* base-set size
* total number of sampled trees
* species selection mode
* seed

### For Mode 2

* taxonomic group
* full-tree source file
* base taxa list specification
* number of input trees
* overlap bounds
* minimum shared taxa
* full-coverage requirement
* anchor taxa count
* pruning mode
* topology perturbation settings
* branch-length perturbation settings
* seed

These settings determine both dataset composition and reproducibility metadata written by the pipeline.

---

## Reproducibility

The pipeline is designed for reproducible dataset generation.

* In Mode 1, `--seed` controls base-species selection and random tree sampling from the downloaded Nexus files.
* In Mode 2, `--seed` controls generated base taxa lists and overlap-set generation.
* In Mode 2, each dataset additionally receives a deterministic dataset-specific seed derived from the main seed and dataset index.

For reproducible use, record:

* the full command line
* the main seed
* the generated metadata file for the dataset, such as `metadata/dataset_XX.json`

---

## Validation and auditing

Mode 2 writes an audit file for each generated dataset:

```text
logs/audit_XX.json
```

The audit records include:

* number of trees written
* minimum and maximum taxa per tree
* minimum and maximum pairwise Jaccard overlap
* minimum and maximum pairwise intersection size
* number of pairs outside the requested overlap bounds
* number of pairs below the minimum shared-taxa threshold
* anchor coverage counts when anchors are explicit 

These outputs support automated checking of the generated benchmark datasets.

---

## Troubleshooting

### Mode 1

**ChromeDriver errors**

* Confirm Chrome or Chromium is installed
* Confirm ChromeDriver is installed and on `PATH`
* Confirm the ChromeDriver major version matches the installed browser major version 

**VertLife delays or download failures**

Mode 1 depends on an external web service. Temporary failures may occur because of service load or network issues.

**Errors with `number_of_trees`**

Mode 1 requires `number_of_trees` to be divisible by 10 and to lie between 10 and 1000 inclusive. 

**Missing `all_species_lists.csv`**

Mode 1 expects this file in the current working directory.

### Mode 2

**`ete3 is required`**

Install ETE3 with:

```bash
pip install ete3
```

**No trees found in `--full_trees`**

Confirm the file contains valid Newick strings and at least one semicolon.

**Too few taxa after pruning**

If the overlap between a base taxa list and the full tree is too small, the corresponding dataset may be skipped. Matching taxon naming conventions between the species lists and the full-tree leaf labels is important.

---

## Citation

If VertLife data are used, cite the relevant [VertLife paper(s)](https://vertlife.org/data/) and the corresponding group-level phylogeny sources. The Mode 1 script also prints the relevant data-source citations after a run.

---

## License

This project is licensed under the [MIT License](LICENSE).
