
<!-- README.md is generated from README.Rmd. Please edit that file -->

# momaz

<!-- badges: start -->

<!-- badges: end -->

**Consensus gene ranking from multi-omic association studies using
Fisher, SSz, and cwMCD**

## Installation

You can install the development version of momaz from
[GitHub](https://github.com/doraaobodo/momaz) with:

``` r
# install.packages("remotes")
remotes::install_github("doraaobodo/momaz")
```

## Description

The package implements and compares three statistical approaches for
consolidating significance rankings from multi-omic association studies
into consensus gene rankings:

- **Fisher’s method**  
- **Sum of squared z-statistics (SSz)**  
- **Cellwise Minimum Covariance Determinant (cellMCD)**

The package includes functions to compute diagnostics, evaluate method
performance in simulations, and generate visualizations such as
bivariate scatterplots, uniform QQ plots, and ranked heatmaps.
