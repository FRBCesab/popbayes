
<!-- README.md is generated from README.Rmd. Please edit that file -->

# popbayes <img src="man/figures/hexsticker.png" height="120" align="right"/>

<!-- badges: start -->

[![R CMD
Check](https://github.com/frbcesab/popbayes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/frbcesab/popbayes/actions/workflows/R-CMD-check.yaml)
[![Website
deployment](https://github.com/frbcesab/popbayes/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/frbcesab/popbayes/actions/workflows/pkgdown.yaml)
[![License: GPL (>=
2)](https://img.shields.io/badge/License-GPL%20%28%3E%3D%202%29-blue.svg)](https://choosealicense.com/licenses/gpl-2.0/)
[![LifeCycle](man/figures/lifecycle/lifecycle-stable.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Project Status:
Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

The goal of the R package `popbayes` is to estimate population trends
from counts series

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
## Install 'remotes' package (if not already installed) ----
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

## Install dev version of 'popbayes' from GitHub ----
remotes::install_github("frbcesab/popbayes", build_vignettes = TRUE)
```

Then you can attach the package `popbayes`:

``` r
library("popbayes")
```

## Overview

![](man/figures/popbayes-diagram.png)

## Get started

Please read the
[Vignette](https://frbcesab.github.io/popbayes/articles/popbayes.html)

## Citation

Please cite this package as:

> Casajus N. & Pradel R. (2021) popbayes: Bayesian model to estimate
> populations trend. R package version 0.1.

You can also run:

``` r
citation("popbayes")

## A BibTeX entry for LaTeX users is:
## 
## @Manual{,
##   title  = {{popbayes}: {B}ayesian model to estimate populations trend,
##   author = {{Casajus N.}, and {Pradel R.}},
##   year   = {2021},
##   note   = {R package version 0.1},
##   url    = {https://github.com/frbcesab/popbayes},
## }
```
