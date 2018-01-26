[![Build Status](https://travis-ci.org/mskilab/JaBbA.svg?branch=master)](https://travis-ci.org/mskilab/JaBbA)
[![codecov.io](https://img.shields.io/codecov/c/github/mskilab/gUtils.svg)](https://codecov.io/github/mskilab/JaBbA?branch=master)

# JaBbA 

Inferring balanced integer cancer genome graphs using mixed integer programming analysis
of read depth and junction patterns in WGS data. . 

Installation
------------

1. Install dependent packages and latest Bioconductor (if you haven't already)

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
biocLite("DNAcopy")
```

2. Install devtools from CRAN (if you don't have it already)

```{r}
install.packages('devtools')  
install.packages('slam')  
install.packages('Rcplex')  
## Requires IBM ILOG CPLEX: install the community version, Download Rcplex source,
## note the include dir, and build Rcplex from source if auto install fails
## R CMD INSTALL --configure-args="--with-cplex-dir=/Applications/CPLEX_Studio_Community128/cplex" Rcplex_0.3-3.tar.gz
```

3. Install dependent mskilab dependencies

```{r}
devtools::install_github('mskilab/gUtils')
devtools::install_github('gTrack')
devtools::install_github('mskilab/gGnome')
```

3. Install JaBbA

```{r}
devtools::install_github('mskilab/JaBbA)
```

