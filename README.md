[![Build Status](https://travis-ci.org/mskilab/JaBbA.svg?branch=master)](https://travis-ci.org/mskilab/JaBbA)
[![codecov.io](https://img.shields.io/codecov/c/github/mskilab/gUtils.svg)](https://codecov.io/github/mskilab/JaBbA?branch=master)

# JaBbA 

Inferring balanced cancer genome graphs using mixed integer programming analysis
of read depth and junction patterns in WGS data. 

Installation
------------

1. Install dependent packages and latest Bioconductor (if you haven't already)

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
```

2. Install devtools from CRAN (if you don't have it already)

```{r}
install.packages('devtools')
```

3. Install dependent mskilab dependencies

```{r}
devtools::install_github('mskilab/gUtils')
devtools::install_github('mskilab/gGnome')
```

3. Install JaBbA

```{r}
devtools::install_github('mskilab/JaBbA)
```
