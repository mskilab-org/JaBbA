[![Build Status](https://travis-ci.org/mskilab/JaBbA.svg?branch=master)](https://travis-ci.org/mskilab/JaBbA)
[![codecov.io](https://img.shields.io/codecov/c/github/mskilab/JaBbA.svg)](https://codecov.io/github/mskilab/JaBbA?branch=master)

# JaBbA 

Inferring balanced cancer genome graphs with mixed-integer programming analysis
of read depth and junction patterns in WGS data. 

Installation
------------
1. Install IBM ILOG
   [CPLEX Studio](https://www.ibm.com/products/ilog-cplex-optimization-studio).
   The software is proprietary, but can be obtained for free under [IBM's academic
   initiative](https://www.ibm.com/products/ilog-cplex-optimization-studio/pricing).

2. Set CPLEX_DIR variable (e.g. in your shell or .bash_profile) to your CPLEX
   Studio installation

```{sh}
export CPLEX_DIR=/path/to/your/copy/of/CPLEX_Studio/
```

3. Install dependent packages and latest Bioconductor (if you haven't already)

```{r}
install.packages(DNAcopy')
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
```

4. Install devtools from CRAN (if you don't have it already)

```{r}
install.packages('devtools')
```

5. Install dependent mskilab R packages

```{r}
devtools::install_github('mskilab/gUtils')
devtools::install_github('mskilab/gGnome')
devtools::install_github('mskilab/Ppurple')
```

6. Install JaBbA

```{r}
devtools::install_github('mskilab/JaBbA)
```
