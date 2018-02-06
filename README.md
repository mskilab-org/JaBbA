[![Build Status](https://travis-ci.org/mskilab/JaBbA.svg?branch=master)](https://travis-ci.org/mskilab/JaBbA)
[![codecov.io](https://img.shields.io/codecov/c/github/mskilab/JaBbA.svg)](https://codecov.io/github/mskilab/JaBbA?branch=master)

# JaBbA 

Inferring balanced cancer genome graphs using mixed integer programming analysis
of read depth and junction patterns in WGS data. 

Installation
------------
1. Install IBM ILOG
   [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio).
   The software is proprietary, but can be obtained for free under [IBM's academic
   initiative](https://www.ibm.com/products/ilog-cplex-optimization-studio/pricing).

2. Install dependent packages and latest Bioconductor (if you haven't already)

```{r}
install.packages('Rcplex')
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
```

3. Install devtools from CRAN (if you don't have it already)

```{r}
install.packages('devtools')
```

4. Install dependent mskilab R packages

```{r}
devtools::install_github('mskilab/gUtils')
devtools::install_github('mskilab/gGnome')
```

5. Install JaBbA

```{r}
devtools::install_github('mskilab/JaBbA)
```
