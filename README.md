# Introduction
we proposed PRMeth, a method to deconvolve tumor mixtures using partially available DNA methylation data. By adopting an iteratively optimized non-negative matrix factorization framework, PRMeth took DNA methylation profiles of a portion of the cell types in the tissue mixtures (including blood and solid tumors) as input to estimate the proportions of all cell types as well as the methylation profiles of unknown cell types simultaneously.
# How to use
## Install and load these packages
''
    install.packages("matrixStats")
    install.packages("quadprog")
    library(matrixStats)
    library(quadprog)
''

![](https://github.com/hedingqin/PRMeth/blob/main/PRMeth/picture/optimalK.png)
