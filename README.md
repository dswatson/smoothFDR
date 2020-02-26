FDR Smoothing for Genomic Data
================

# Overview

False discovery rate (FDR) smoothing is an empirical Bayes method for
exploiting spatial structure in large multiple-testing problems. FDR
smoothing automatically finds spatially localized regions of significant
test statistics. It then relaxes the threshold of statistical
significance within these regions and tightens it elsewhere, ensuring
overall FDR control at a given level. This results in increased power
and cleaner spatial separation of signals from noise. The approach
requires solving a nonstandard high-dimensional optimization problem,
for which an efficient augmented-Lagrangian algorithm is implemented.
For details, see [Tansey et al.’s (2018)](https://bit.ly/2msEdWH) paper.

# Installation

To install the package, run the following in R:

``` r
# Install devtools if you have not already
install.packages('devtools')

# Then install directly from GitHub
devtools::install_github('dswatson/smoothFDR')
```

# Example

``` r
# Load library
library(smoothFDR)

# Set seed
set.seed(123)

# Import DNA methylation data, distributed with the package
data('DNAm')

# Run FDR smoothing
res <- smoothFDR(DNAm, probe = 'cpg', parallel = FALSE)

# How many significant CpG sites at 5% FDR according to Benjamini-Hochberg?
sum(res$BH_q.value <= 0.05)
```

    ## [1] 0

``` r
# How many according to the smooth FDR estimate?
sum(res$q.value <= 0.05)
```

    ## [1] 255
