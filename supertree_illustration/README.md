# Illustrative application to supertree construction

This pilot supertree illustartion utilized overlapping mammal datasets (constructed by mode 1). The experimental design involved:

* Datasets: 100 distinct input sets.
* Taxon subsets: Each set comprised 10 overlap-controlled taxon subsets, with target overlap levels ranging from 10% to 90%.
* Tree composition: 30 phylogenetic input trees per set (3 trees per subset).
* Scale: Each tree included 30 taxa sampled from a total of 55 unique species.

Content:

- **`input_datasets/`**

  Contains 100 input sets of phylogenetic trees used for the supertree reconstruction demonstration. Each file (`multiset_X.txt`) consists of 30 trees with overlapping taxa, constructed using the proposed pipeline (mode 1).

- **`supertrees/`**

  This folder includes output supertrees produced by five different methods:
  - `supertrees_sfit.txt` — Split Fit algorithm
  - `supertrees_dfit.txt` — Most Similar Supertree
  - `supertrees_nj.txt` — Average NJ
  - `supertrees_mrplus.txt` — Majority-Rule
  - `supertrees_scs.txt` — Spectral Clustering
 
 - **`supertree_validation.ipynb`**
   
   Contains validation results.
