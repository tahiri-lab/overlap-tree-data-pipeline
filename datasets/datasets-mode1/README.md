## Datasets of phylogenetic trees defined on different but overlapping sets of taxa (Mode 1)

This folder contains datasets of phylogenetic trees defined on different but overlapping sets of taxa constructed by a proposed data pipeline in mode 1.

### Data description :deciduous_tree:

* The biological data utilized in this study was obtained from [vertlife.org](https://vertlife.org/phylosubsets/) [[1]](#ref1), which offers a straightforward method for acquiring tree distributions with specified subsets of taxa. The tool initially prunes a comprehensive dataset to a smaller subset and then samples trees from the selected pseudoposterior distribution.
* The following four groups are provided: :frog: **amphibians** [[2]](#ref2), :eagle: **birds** [[3]](#ref3), :monkey: **mammals** [[1]](#ref1), and :shark: **sharks** [[4]](#ref4).
* The number of species included in each group varies (overall). In particular, there are 7239 species of amphibians, 9993 species of birds, 5911 species of mammals, and 1192 species of sharks (see [all_species_lists.csv](https://github.com/tahiri-lab/overlap-tree-data-pipeline/blob/main/all_species_lists.csv).
* The final dataset for each group is comprised of 10 subsets that overlap by species, with the degree of overlap ranging from 10% to 90% (see the following Figure).

![Levels of overlap among subsets for 4 groups of species](https://github.com/tahiri-lab/KNCL/blob/main/data/images/overlaps_subsets.png "Levels of overlap among subsets for 4 groups of species")

The biological data consists of the following four files:

1. :frog: [**`amphibians_trees.txt`**](https://github.com/tahiri-lab/overlap-tree-data-pipeline/blob/main/datasets/datasets-mode1/amphibians_trees.txt):
   - Contains 550 phylogenetic trees (Newick) of *Amphibians* with overlap levels ranging from 10% to 100%.
   - Total number of unique species: 120.
   - Average level of overlap: 58.71%.
   - Number of unique pairs of trees: 150975.
   - Number of unique pairs of trees with 100% overlap: 14850 or 9.84%.
  
2. :eagle: [**`birds_trees.txt`**](https://github.com/tahiri-lab/overlap-tree-data-pipeline/blob/main/datasets/datasets-mode1/birds_trees.txt):
   - Contains 600 phylogenetic trees (Newick) of *Birds* with overlap levels ranging from 10% to 100%.
   - Total number of unique species: 135.
   - Average level of overlap: 58.78%.
   - Number of unique pairs of trees: 179700.
   - Number of unique pairs of trees with 100% overlap: 17700 or 9.85%.

3. :monkey: [**`mammals_trees.txt`**](https://github.com/tahiri-lab/overlap-tree-data-pipeline/blob/main/datasets/datasets-mode1/mammals_trees.txt):
   - Contains 500 phylogenetic trees (Newick) of *Mammals* with overlap levels ranging from 10% to 100%.
   - Total number of unique species: 105.
   - Average level of overlap: 58.88%.
   - Number of unique pairs of trees: 124750.
   - Number of unique pairs of trees with 100% overlap: 12250 or 9.82%.

4. :shark: [**`sharks_trees.txt`**](https://github.com/tahiri-lab/overlap-tree-data-pipeline/blob/main/datasets/datasets-mode1/sharks_trees.txt):
   - Contains 450 phylogenetic trees (Newick) of *Sharks* with overlap levels ranging from 10% to 100%.
   - Total number of unique species: 95.
   - Average level of overlap: 58.81%.
   - Number of unique pairs of trees: 101025.
   - Number of unique pairs of trees with 100% overlap: 9900 or 9.80%.

The distributions of overlap levels of biological data for all possible pairs of trees (without repeats) are presented below.

![The distributions for the levels of overlap in the biological data](https://github.com/tahiri-lab/KNCL/blob/main/data/images/overlap_levels_comparison.png "The distributions for the levels of overlap in the biological data")

### References

1. <a id="ref1"></a> Upham, N. S., J. A. Esselstyn, and W. Jetz. 2019. Inferring the mammal tree: species-level sets of phylogenies for questions in ecology, evolution, and conservation. *PLOS Biology*. [https://doi.org/10.1371/journal.pbio.3000494](https://doi.org/10.1371/journal.pbio.3000494)

2. <a id="ref2"></a> Jetz, W., and R. A. Pyron. 2018. The interplay of past diversification and evolutionary isolation with present imperilment across the amphibian tree of life. *Nature Ecology & Evolution*, 1. [https://www.nature.com/articles/s41559-018-0515-5](https://www.nature.com/articles/s41559-018-0515-5)

3. <a id="ref3"></a> Jetz, W., G. H. Thomas, J. B. Joy, K. Hartmann, and A. O. Mooers. 2012. The global diversity of birds in space and time. *Nature*, 491:444–448. [http://www.nature.com/nature/journal/v491/n7424/abs/nature11631.html](http://www.nature.com/nature/journal/v491/n7424/abs/nature11631.html)

4. <a id="ref4"></a> Stein, R. W., Mull, C. G., Kuhn, T. S., Aschliman, N. C., Davidson, L. N. K., Joy, J. B., Smith, G. J., Dulvy, N. K., & Mooers, A. O. 2018. Global priorities for conserving the evolutionary history of sharks, rays, and chimaeras. *Nature Ecology & Evolution*, 2:288–298. [http://dx.doi.org/10.1038/s41559-017-0448-4](http://dx.doi.org/10.1038/s41559-017-0448-4)
