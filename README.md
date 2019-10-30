# popbayes

[![Build Status](https://travis-ci.org/FRBCesab/popbayes.svg?branch=master)](https://travis-ci.org/FRBCesab/popbayes) [![](https://img.shields.io/badge/licence-GPLv3-8f10cb.svg)](http://www.gnu.org/licenses/gpl.html)

Overview
--------

In construction

Installation
--------

To install the package `popbayes` from GitHub, first install the package [`devtools``](http://cran.r-project.org/web/packages/devtools/index.html) from the CRAN.

```r
### Install the < devtools > package
install.packages("devtools", dependencies = TRUE)
```

Then install the `popbayes` package:

```r
### Install the < popbayes > package from GitHub
devtools::install_github("frbcesab/popbayes", build_vignettes = TRUE)

### Load the < popbayes > package
library(popbayes)
```

Getting started
--------

```r
### Browse the < popbayes > package vignette
vignette(topic = "popbayes")
```
