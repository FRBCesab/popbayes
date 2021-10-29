
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

The goal of the R package `popbayes` is to infer trends of one or
several populations over time from series of counts. It does so by
accounting for count precision (provided or inferred based on expert
knowledge, e.g. for guesstimates), smoothing the population rate of
increase over time, and accounting for the maximum demographic potential
of species. Inference is carried out in a Bayesian framework.

## Installation

**Before using this package, users need to install the freeware
[JAGS](https://mcmc-jags.sourceforge.io/).**

You can install the development version from
[GitHub](https://github.com/) with:

``` r
## Install 'remotes' package (if not already installed) ----
#  install.packages("remotes")

## Install development version of 'popbayes' ----
remotes::install_github("frbcesab/popbayes", build_vignettes = TRUE)
```

## Overview

![](vignettes/docs/popbayes-diagram.png)

## Get started

Please read the [Get
started](https://frbcesab.github.io/popbayes/articles/popbayes.html)
vignette.

## Citation

Please cite this package as:

> Casajus N. & Pradel R. (2021) popbayes: Bayesian model to estimate
> population trends from counts series. R package version 1.0. URL:
> <https://frbcesab.github.io/popbayes/>.

You can also run:

``` r
citation("popbayes")

## A BibTeX entry for LaTeX users is:
## 
## @Manual{,
##   title  = {{popbayes}: {B}ayesian model to estimate population trends from counts series,
##   author = {{Casajus N.}, and {Pradel R.}},
##   year   = {2021},
##   note   = {R package version 1.0},
##   url    = {https://frbcesab.github.io/popbayes/},
## }
```
