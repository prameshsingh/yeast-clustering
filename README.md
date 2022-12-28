# yeast-clustering

This repository contains codes for construction of **S. cerevisiae** (yeast) gene co-expression networks and community detection in them.


To use the submodules run the following commands:

```
git submodule init
git submodule update
```

## Network Inference 

The gene co-expression network is inferred using the following three measures:

- Pearson correlation (R): 
We use `r` function `cor()` with `method = "pearson"` 

- Mutual Information (MI):
Use the submodule `mutual-information-relatedness` for computing *MI*.
See `README.md` in `mutual-information-relatedness` for more informaiton.


- Relateness (f):
Use the submodule `mutual-information-relatedness` for computing *f*.
See `README.md` in `mutual-information-relatedness` for more informaiton.

## Community detection (Clustering)

To detect communities in the inferred gene co-expression networks, we use the following four methods:

- Markov Cluster Algorithm (MCL):
For MCL, we use `r` package `mcl()`. See https://cran.r-project.org/web/packages/MCL/index.html for more info.

- Modularity (Q) and Generalize Modularity Density ($`Q_g`$):
To find communities using *Q* or *Qg*, use an implementation of the RenEEL algorithm from the submodule `generalized-modularity-density`. *Q* is a special case of *Qg* when the exponent \chi is 0. `See README.md` in `generalized-modularity-density` for more info. 

- Excess Modularity Density (Qx): 
For community detection by maximizing Qx, use the submodule `excess-modularity-density`. 

## Evaluation

To evaluate the modules by Average Adujested Rand Index (AARI) run `src/evaluation.R`. 
All clustering results (one column for each method) are contained in the file `data/all_partitions.xlsx` and the leaf GO-terms associations are in `data/gt_all_solution.xlsx`. 

Required libraries: 
- `openxlsx`
- `igraph`
- `ggplot2` (if plotting).

