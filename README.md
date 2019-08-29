# REDassignment
Benny Lin, Connor Morgan-Lang, Steven J. Hallam

## Overview:
Script for assigning relative evolutionary distances (RED) to taxonomic trees and generating linear, logistic, and random forest classifier models from the derived RED values versus taxonomic ranks. The script was implemented into the TreeSAPP python package and evaluated on classification performance.

TreeSAPP is a functional and taxonomic annotation tool that leverages phylogenetic placements and distances to provide accurate taxonomic assignments (https://github.com/hallamlab/TreeSAPP). It currently uses PD to recommend taxonomic rank, rather than relying on LCA alone. The rank recommendendation algorithm in TreeSAPP is currently computationally burdensome: iterative clade exclusion analysis is performed where at each iteration sequences that represent a taxon are removed from the tree, then inserted and their PDs and taxonomic rank they represent are recorded for downstream regression. By replacing PD with RED values, TreeSAPP's rank recommendation process becomes more efficient. 

## Running REDassigment
Note: compatible only with Python 3

To list all required inputs run `py REDAssignment.py`.

To decorate a tree's nodes with RED values and generate a model for those values:
```
py REDAssignment.py -i taxonomic_ids_of_tree.txt -t newick_tree_file_of_interest.txt -r Strings From Lineages Not Wanted -l 2 -m linear
```
