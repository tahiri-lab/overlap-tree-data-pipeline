# Pipeline algorithm pseudocodes: construction of phylogenetic tree datasets with partial taxon overlap

This document provides high-level pseudocodes of the implemented pipeline, organized by construction mode and processing stage. It complements the runnable code, and it is not an alternative implementation or a specification of behavior that the code does not enforce.

**Implementation files:**

- `overlap_tree_pipeline.py` — unified entry point and Mode 2 part.
- `overlap_tree_pipeline_modeA.py` — Mode 1 selection, subset construction, VertLife retrieval, and dataset assembly.
- `make_base_species_lists_phylo_stratified.py` — optional phylogeny-stratified taxon list generation for Mode 2.
- `gen_pruned_trees.py` — Mode 2 subset selection, coverage repair, pruning, and perturbation.

All input trees are assumed to be available as Newick or Nexus data, as appropriate. `SET(T)` means the set of terminal taxon labels in a tree `T`; `PRUNE(T,S)` denotes phylogenetic pruning retaining taxa `S` (with degree-2 contraction when configured). Set operations use their standard mathematical meanings.

## Common definitions

For two taxon sets `A` and `B`, the **Jaccard similarity** is

\[
J(A,B)=\frac{|A\cap B|}{|A\cup B|}.
\]

Mode 1 specifies *target* Jaccard values between each generated taxon subset and one base subset. Mode 2 tests *interval constraints* for each new subset against all previously accepted subsets. Neither mode is described here as achieving exact targets unconditionally: Mode 1 has rounding/wrapping effects, and Mode 2 may relax or bypass candidate constraints and subsequently repair coverage.

---

## Algorithm 1 — Unified pipeline dispatcher

**Input:** Mode `service` or `reference`; mode-specific configuration.  
**Output:** Mode 1 overlapping-tree dataset, or Mode 2 reference trees and input-tree multisets.

```text
PROCEDURE CONSTRUCT_DATASETS(mode, configuration)
    IF mode = "service" THEN
        CALL MODE1_CONSTRUCT_SUBSETS(configuration)
        IF configuration.prepare_only THEN
            RETURN generated species and subset lists
        END IF
        CALL MODE1_RETRIEVE_AND_ASSEMBLE(configuration)
    ELSE IF mode = "reference" THEN
        CALL MODE2_PREPARE_REFERENCES_AND_GENERATE(configuration)
    ELSE
        RAISE invalid-mode error
    END IF
END PROCEDURE
```

## Mode 1 — Empirical-source trees with overlap targets relative to a base subset

### Algorithm 2 — Species selection and construction of ten overlapping subsets

**Input:** A taxonomic group, a species universe or user-provided list, requested base-list size `n`, species selection mode, nine target overlaps `P = (0.9, 0.8, ..., 0.1)`.  
**Output:** Ordered base list `B` and ten equal-size taxon subsets `S1, ..., S10`.

```text
PROCEDURE MODE1_CONSTRUCT_SUBSETS(configuration)
    B <- SELECT_BASE_SPECIES(
             group species universe,
             mode in {user_list, random, genus_stratified},
             requested n, optional seed and per-genus cap)
    n <- LENGTH(B)
    P <- [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]

    # Find the first k whose successive window offsets cover B.
    k <- NONE
    FOR candidate_k <- 1 TO 1000 DO
        required_span <- candidate_k
        previous_offset <- 0
        FOR p in P DO
            nominal_shared <- ROUND(2*candidate_k*p/(1+p))
            offset <- candidate_k - nominal_shared
            required_span <- required_span + CEIL(offset - previous_offset)
            previous_offset <- offset
        END FOR
        IF required_span >= n THEN
            k <- candidate_k
            BREAK
        END IF
    END FOR
    IF k = NONE THEN
        REPORT inability to construct the subsets
        RETURN
    END IF

    S1 <- first k entries of B
    FOR i <- 2 TO 10 DO
        p <- P[i-1]
        nominal_shared <- ROUND(2*k*p/(1+p))
        offset <- k - nominal_shared
        Si <- k consecutive entries of B beginning at offset,
              wrapping to the start of B if necessary
    END FOR

    WRITE B to selected_species.csv
    WRITE S1,...,S10 to <group>_overlapping_subsets.csv
    RETURN B, S1,...,S10
END PROCEDURE
```

`ROUND` denotes Python's `round()` applied to the floating-point expression: halfway cases use ties-to-even, while floating-point evaluation can affect values very close to a halfway value. In the current ordered target sequence, the offset sequence is nondecreasing and the implemented cumulative-size calculation simplifies to `k + max(offset)`. The construction targets `J(S1,Si) = p` but does not independently specify the overlap of `Si` and `Sj` for arbitrary `i,j`.

For base-list length `n` and subset size `k`, feasibility requires `|S1 ∩ Si| >= max(0, 2k-n)`. In the implemented circular-window construction, the realized common-taxon count is `max(nominal_shared, 2k-n)`; hence realized Jaccard values may differ from their targets.


### Algorithm 3 — VertLife retrieval, local sampling, and export

**Input:** Ten taxon subsets; total requested number `N` of trees divisible by ten; VertLife access or a local folder of previously downloaded Nexus files; optional random seed for local sampling.  
**Output:** Combined Newick dataset and retained subset lists.

```text
PROCEDURE MODE1_RETRIEVE_AND_ASSEMBLE(configuration)
    REQUIRE N between 10 and 1000, and N mod 10 = 0
    t <- N / 10

    IF an existing Nexus folder was supplied THEN
        nexus_files <- prepare ten existing Nexus inputs
    ELSE
        FOR each of the ten subsets DO
            SUBMIT subset and group to VertLife PhyloSubsets
            DOWNLOAD the resulting tree file or archive
        END FOR
        nexus_files <- prepare the ten downloaded Nexus inputs
    END IF

    combined <- empty sequence
    FOR each Nexus file in the prepared list DO
        trees <- PARSE_NEXUS(file)
        IF LENGTH(trees) < t THEN
            WARN and use all available trees from that file
        ELSE
            trees <- SAMPLE_WITHOUT_REPLACEMENT(trees, t)
        END IF
        CONVERT selected trees to Newick (retain branch lengths)
        APPEND the converted Newick strings to combined
    END FOR

    WRITE combined, one Newick tree per line,
          to overlapping_dataset_<group>.txt
END PROCEDURE
```

Reproducibility boundary: The local seed can control local species selection and tree subsampling, but not the initial tree sampling performed remotely by VertLife. To reproduce an exact Mode 1 release, retain the original Nexus inputs. The current code warns and samples fewer than `t` trees when an input file is too short; it does not silently create missing trees.

---

## Mode 2 — Reference-based benchmark collections with pairwise overlap intervals

### Algorithm 4 — Reference-tree preparation and dataset-level orchestration

**Input:** One or more broad-coverage source Newick trees, taxonomic universe, base-taxon lists or base-list generation settings, pairwise overlap interval `[lo,hi]`, input-tree count `K`, perturbation settings.  
**Output:** One reference tree and one overlapping input-tree multiset for each processed dataset.

```text
PROCEDURE MODE2_PREPARE_REFERENCES_AND_GENERATE(configuration)
    full_trees <- READ_NEWICK_TREES(configuration.full_trees)
    IF explicit base-taxon lists were supplied THEN
        bases <- READ_BASE_LISTS()
    ELSE
        bases <- PHYLOGENY_STRATIFIED_BASE_LISTS(
                     first full tree, species universe,
                     requested sizes, seed, strata and reuse settings)
    END IF

    # Alignment rule used by the wrapper:
    IF only one full tree THEN
        pair that full tree with every base list
    ELSE IF only one base list THEN
        pair that list with every full tree
    ELSE IF counts match THEN
        pair corresponding entries
    ELSE
        WARN and pair only the first min(counts) entries
    END IF

    FOR each paired (full_tree, base_list), indexed by i DO
        R <- PRUNE(full_tree, base_list)
        IF fewer than two taxa remain THEN
            SKIP this reference
            CONTINUE
        END IF
        SAVE R as the dataset's reference tree
        SET_RANDOM_SEED(base_seed + 1000*i)   # i starts at 1
        configure the integer target-size interval from
            the reference-taxon fraction bounds

        leaf_sets, notes <- MODE2_GENERATE_INPUT_TAXON_SETS(R, configuration)
        CALL MODE2_PRUNE_PERTURB_EXPORT(R, leaf_sets, notes, configuration)
    END FOR
END PROCEDURE
```

When base lists are generated automatically, the helper partitions the first source phylogeny into clade-based strata and samples taxa approximately across strata, using reuse/underuse weighting and a soft reuse cap. Supplied base lists bypass that generation step. The retained reference taxa are the intersection between the specified base list and taxa actually present in the source tree.


### Algorithm 5 — Candidate subset generation, overlap checks, relaxation, and fallback

**Input:** Reference tree `R`; requested input-tree count `K`; target size range; pairwise overlap interval `[lo,hi]`; minimum pairwise shared taxa `q`; number of common anchor taxa `a`; selection settings; retry limit `M` (default 200); optional relaxation width `delta` (default 0.05).  
**Output:** Input taxon sets `S1,...,SK`, common anchor set `A`, and pre-repair selection notes.

```text
PROCEDURE MODE2_GENERATE_INPUT_TAXON_SETS(R, configuration)
    U <- SET(R)
    A <- UNIFORM_RANDOM_SAMPLE(U, min(a, |U|))
    target_sizes <- SAMPLE_INTEGER_TARGET_SIZES(K, configured size interval)
    accepted_sets <- empty sequence
    notes <- empty sequence

    FOR i <- 1 TO K DO
        target_size <- max(target_sizes[i], |A| + q)
        accepted <- FALSE

        # First pass: test against the originally requested bounds.
        FOR attempt <- 1 TO M DO
            C, method <- GENERATE_CANDIDATE(R, target_size,
                                           prune_mode, clade_bias, blockiness)
            C <- FORCE_ANCHORS_AND_TRIM_NONANCHORS(C, A, target_size)
            IF VALID_AGAINST_ALL(C, accepted_sets, lo, hi, q) THEN
                APPEND C to accepted_sets
                RECORD method and relaxed_bounds = NONE
                accepted <- TRUE
                BREAK
            END IF
        END FOR

        # Despite its name, strict_constraints does not guarantee strict failure.
        IF NOT accepted AND strict_constraints = TRUE THEN
            relaxed_lo <- max(0, lo - delta)
            relaxed_hi <- min(1, hi + delta)
            FOR attempt <- 1 TO M DO
                C, method <- GENERATE_CANDIDATE(...)
                C <- FORCE_ANCHORS_AND_TRIM_NONANCHORS(C, A, target_size)
                IF VALID_AGAINST_ALL(C, accepted_sets,
                                     relaxed_lo, relaxed_hi, q) THEN
                    APPEND C to accepted_sets
                    RECORD method and relaxed_bounds = [relaxed_lo,relaxed_hi]
                    accepted <- TRUE
                    BREAK
                END IF
            END FOR
        END IF

        IF NOT accepted THEN
            C, method <- GENERATE_CANDIDATE(...)
            C <- FORCE_ANCHORS_AND_TRIM_NONANCHORS(C, A, target_size)
            APPEND C to accepted_sets WITHOUT an overlap check
            RECORD method and relaxed_bounds = "failed_all"
        END IF
    END FOR

    IF full_coverage_enabled OR per_leaf_min_coverage > 0 THEN
        CALL MODE2_REPAIR_COVERAGE(U, accepted_sets, A, configuration)
    END IF
    RETURN accepted_sets, notes, A
END PROCEDURE
```

`GENERATE_CANDIDATE` selects taxa either by random leaves (optionally biased toward one clade), by a greedy union of clades followed by random filling, or by a configured mixture. `VALID_AGAINST_ALL` requires both `|C ∩ S| >= q` and `lo <= J(C,S) <= hi` for every previously accepted `S`. The code also supports a different overlap metric (`base_fraction`) through direct generator configuration; the unified Mode 2 workflow uses Jaccard.

Important: The widened pass is attempted only when `strict_constraints=True`, and an unconditional fallback exists regardless of that flag. The fallback does not raise a constraint-satisfaction exception. `relaxed_bounds` notes describe *candidate acceptance before coverage repair*, not a verified final constraint status.


### Algorithm 6 — Coverage repair after candidate selection

**Input:** Reference taxon set `U`, accepted input taxon sets `S1,...,SK`, anchor set `A`, requested minimum per-taxon coverage `c`, original overlap constraints.  
**Output:** Modified input taxon sets (best effort to retain overlap constraints).

```text
PROCEDURE MODE2_REPAIR_COVERAGE(U, sets, A, configuration)
    missing <- SHUFFLE(U minus UNION(sets))
    FOR each taxon x in missing DO
        i <- index of a smallest set in sets
        trial <- sets[i] UNION {x}
        IF trial satisfies constraints against all other sets THEN
            sets[i] <- trial
        ELSE
            removable <- taxa of sets[i] not in A,
                         ordered by descending frequency across sets
            swapped <- FALSE
            FOR each y in removable DO
                trial <- (sets[i] minus {y}) UNION {x}
                IF trial satisfies constraints against all other sets THEN
                    sets[i] <- trial
                    swapped <- TRUE
                    BREAK
                END IF
            END FOR
            IF NOT swapped THEN
                sets[i] <- sets[i] UNION {x}  # unconditional last resort
            END IF
        END IF
    END FOR

    IF c > 0 THEN
        FOR each reference taxon x whose current coverage is below c DO
            FOR eligible sets not yet containing x, following the
                implementation's set-ordering rule DO
                trial <- set UNION {x}
                IF trial satisfies constraints against all other sets THEN
                    accept the addition
                    IF coverage(x) >= c THEN BREAK END IF
                END IF
            END FOR
        END FOR
    END IF
    RETURN sets
END PROCEDURE
```

The union-coverage stage can introduce an out-of-bounds pair if no constraint-preserving repair is found. The per-taxon minimum-coverage stage is best effort and has no equivalent unconditional fallback. The generator calls this repair whenever either `enforce_full_coverage` is true or `per_leaf_min_coverage > 0`; within the repair function, union coverage is attempted first.


### Algorithm 7 — Rooted NNI perturbation with optional anchor protection

**Input:** An input tree `T` after pruning, requested number `m` of NNI moves, common anchor set `A`, anchor-protection flag.  
**Output:** Perturbed tree `T`.

```text
PROCEDURE MODE2_NNI_PERTURB(T, m, A, protect_anchors)
    FOR step <- 1 TO m DO
        eligible <- empty sequence
        FOR each nonroot node v of T DO
            u <- parent(v)
            IF u has exactly two children AND v has exactly two children THEN
                a <- sibling(v) under u
                b, c <- children(v)
                IF protect_anchors AND
                   (a or b or c contains any taxon in A) THEN
                    CONTINUE
                END IF
                APPEND (u,v,a,b,c) to eligible
            END IF
        END FOR

        IF eligible is empty THEN
            RETURN T  # fewer than m successful moves
        END IF

        (u,v,a,b,c) <- UNIFORM_RANDOM_CHOICE(eligible)
        z <- UNIFORM_RANDOM_CHOICE({b,c})
        old_a_length <- length of edge (u,a)
        old_z_length <- length of edge (v,z)
        DETACH a and z
        ATTACH z under u, and a under v
        ASSIGN z the former length of (u,a)
        ASSIGN a the former length of (v,z)
    END FOR
    RETURN T
END PROCEDURE
```

The locally eligible configuration is a rooted binary parent-child pair, with a swap between the sibling of the child and one grandchild. Eligibility is recomputed after each move. No successful-move counter is returned or stored by this implementation, and final Robinson–Foulds differences cannot reconstruct the historical number of executed NNIs. The alternative `swap_labels` mode swaps selected leaf names and pendant lengths instead of applying NNI.


### Algorithm 8 — Final tree creation, branch-length modification, and post-generation audit

**Input:** Reference tree `R`; selected taxon sets `S1,...,SK`; anchors and perturbation settings.  
**Output:** Per-dataset input multiset, reference tree, metadata, and final-overlap audit.

```text
PROCEDURE MODE2_PRUNE_PERTURB_EXPORT(R, sets, notes, configuration)
    H_ref <- maximum root-to-tip distance of R
    newicks <- empty sequence

    FOR each taxon set Si DO
        T <- PRUNE(R, Si)  # contract degree-2 nodes when requested

        IF topology_noise = "nni" THEN
            T <- MODE2_NNI_PERTURB(T, requested_moves, anchors, protection)
        ELSE IF topology_noise = "swap_labels" THEN
            SWAP_SELECTED_LEAF_LABELS_AND_PENDANT_LENGTHS(T)
        END IF

        IF branch_length_scaling = "global_uniform" THEN
            f <- UNIFORM_RANDOM(low, high)
            MULTIPLY_EVERY_EDGE_LENGTH(T, f)
        ELSE IF branch_length_scaling selects another supported mode THEN
            APPLY_THE_CONFIGURED_LENGTH_TRANSFORMATION(T)
        END IF

        IF branch_length_scaling != "none" AND
           renormalize_root_height = "to_base" THEN
            H <- maximum root-to-tip distance of T
            IF H > 0 THEN
                MULTIPLY_EVERY_EDGE_LENGTH(T, H_ref / H)
            END IF
        END IF

        APPEND NEWICK(T) to newicks
    END FOR

    WRITE reference R and combined input-tree multiset
    WRITE metadata: parameters, base list, anchor taxa,
          selected-set construction type, and pre-repair relaxation notes
    audit <- CHECK_FINAL_PAIRWISE_OVERLAPS(newicks, original [lo,hi], q)
    audit <- audit plus final input-tree size range and anchor membership
    WRITE audit to the dataset audit JSON
END PROCEDURE
```

This part computes realized pairwise Jaccard overlap minima/maxima, numbers of out-of-range and insufficient-shared-taxon pairs, subset-size extrema, and the count of trees containing all anchors. It does not enforce a strict failure when violations are found, or itself compare the union of input taxa with the full reference taxon set. Full taxon coverage, rooted RF, and branch-length deviations can be checked in separate post-hoc validation analyses.

Important branch-length detail: `global_uniform` uses one factor per tree. If subsequent root-height normalization rescales all edges to match the full reference's root height, the initial uniform global factor cancels algebraically (apart from numerical precision); this combination does not preserve an independently measurable global scaling effect.