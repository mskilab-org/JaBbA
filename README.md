[![Build Status](https://travis-ci.org/mskilab/JaBbA.svg?branch=master)](https://travis-ci.org/mskilab/JaBbA)
[![codecov.io](https://img.shields.io/codecov/c/github/mskilab/JaBbA.svg)](https://codecov.io/github/mskilab/JaBbA?branch=master)

# JaBbA
## (Junction Balance Analysis)

Inferring balanced cancer genome graphs with mixed-integer programming analysis
of read depth and junction patterns in WGS data. 
 
Installation (R package)
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


Installation (command line)
------------

7. Pull JaBbA git and add pulled directory to PATH

```{bash}
$ git clone git@github:mskilab/JaBbA
$ export PATH="$PATH:$PWD/JaBbA"
```

8. test run jba executable on provided data
```{bash}

$ jba JaBbA/inst/extdata/junctions.vcf JaBbA/inst/extdata/coverage.txt -w -u -v
 _____         ___    _      _____ 
(___  )       (  _`\ ( )    (  _  )
    | |   _ _ | (_) )| |_   | (_) |
 _  | | /'_` )|  _ <'| '_`\ |  _  |
( )_| |( (_| || (_) )| |_) )| | | |
`\___/'`\__,_)(____/'(_,__/'(_) (_)

(Junction     Balance     Analysis)

JaBbA 2018-02-13 21:29:50: Located junction file JaBbA/inst/extdata/junctions.vcf
JaBbA 2018-02-13 21:29:50: Located coverage file JaBbA/inst/extdata/coverage.txt
JaBbA 2018-02-13 21:29:50: Loading packages ...
JaBbA 2018-02-13 21:30:00: Starting analysis in ./jbaout
JaBbA 2018-02-13 21:32:13: Done .. job output in: ./jbaout

```

9. command line documentation

```{bash}
 _____         ___    _      _____ 
(___  )       (  _`\ ( )    (  _  )
    | |   _ _ | (_) )| |_   | (_) |
 _  | | /'_` )|  _ <'| '_`\ |  _  |
( )_| |( (_| || (_) )| |_) )| | | |
`\___/'`\__,_)(____/'(_,__/'(_) (_)

(Junction     Balance     Analysis)

Usage: jba [options] JUNCTIONS COVERAGE
 	JUNCTIONS can be BND style vcf, bedpe, rds of GrangesList
 	COVERAGE is a .wig, .bw, .bedgraph, .bed., .rds of a granges, or .txt file that is coercible to a granges, use --field=FIELD argument if using specific field of a multi-column table')


Options:
	-s SEG, --seg=SEG
		Path to .rds file of GRanges object or .bed file or  .txt / .csv file of intervals corresponding to initial segmentation (optional, will use CBS of coverage to compute if not provided)

	-f FIELD, --field=FIELD
		Name of meta data field or column of coverage file to use for coverage signal from coverage file, may be required if coverage file has several fields

	-t TFIELD, --tfield=TFIELD
		Name of meta data field of ra GRanges or data frame that specifies tiers of junctions, where tier 1 is forced to be included, tier 2 is optional, and tier 3 junctions are only used in when --iterate is set

	-i HETS, --hets=HETS
		Path to tab or commadelimited hets file output of het counts with columns /fields $seqnames, $start, $end, $alt, $ref

	-o OUTDIR, --outdir=OUTDIR
		Directory to dump output into (default JaBbA)

	-k SLACK, --slack=SLACK
		Slack penalty to apply per loose end copy

	-z SUBSAMPLE, --subsample=SUBSAMPLE
		Numeric value between 0 and 1 specifying whether to subsample coverage for intra segment variance estimation

	-l TILIM, --tilim=TILIM
		Time limit for JaBbA MIP

	-p PLOIDY, --ploidy=PLOIDY
		Ploidy guess

	-q PURITY, --purity=PURITY
		Purity guess

	-c CORES, --cores=CORES
		Number of cores for JaBBa MIP

	-m ITERATE, --iterate=ITERATE
		How many times to iterate through tiers

	-r WINDOW, --window=WINDOW
		Integer window (in bp) to dumpster dive and rescue junctions near loose ends (default 10,000)

	-e EDGENUDGE, --edgenudge=EDGENUDGE
		Edge nudge for optimization, to be multiplied by edge specific confidence score if provided

	-b NSEG, --nseg=NSEG
		Path to .rds file of GRanges object of intervals corresponding or .txt / .tsv / .csv with GRanges fields (seqnames, start, end) specifying copy number for normal tissue, needs to have $cn field

	-x, --strict
		restricting input junctions to only the subset overlapping seg

	-u, --gurobi
		flag to use gurobi (gurobi R package must be installed) instead of default CPLEX (cplex must be installed prior to library installation)

	-w, --overwrite
		Flag whether to overwrite previous directory

	-v, --verbose
		verbose output

	-y, --nudgebalanced
		Manually nudge balanced junctions into the model.

	-h, --help
		Show this help message and exit



```
