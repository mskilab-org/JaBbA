[![Build Status](https://travis-ci.org/mskilab/JaBbA.svg?branch=master)](https://travis-ci.org/mskilab/JaBbA)
[![codecov.io](https://img.shields.io/codecov/c/github/mskilab/JaBbA.svg)](https://codecov.io/github/mskilab/JaBbA?branch=master)

# JaBbA (Junction Balance Analysis)

Inferring balanced cancer genome graphs with mixed-integer programming analysis
of read depth and junction patterns in WGS data. 
 
Installation (R package)
------------
1. Install IBM ILOG
   [CPLEX Studio](https://www.ibm.com/products/ilog-cplex-optimization-studio).
   The software is proprietary, but can be obtained for free under [IBM's academic
   initiative](https://www.ibm.com/products/ilog-cplex-optimization-studio/pricing).

2. Set `CPLEX_DIR` variable to your CPLEX Studio installation directory

```{sh}
export CPLEX_DIR=/path/to/your/copy/of/CPLEX_Studio/
```
**NOTE: if `CPLEX_DIR` is set correctly then `$CPLEX_DIR/cplex/include` and `$CPLEX_DIR/cplex/lib` should both exist.**

3. Install JaBbA 

```{r}
devtools::install_github('mskilab/JaBbA')
```


Installation (jba executable)
------------

4. (after installing R package) Pull JaBbA git and add pulled directory to PATH

```{bash}
$ git clone git@github:mskilab/JaBbA
$ export PATH="$PATH:$PWD/JaBbA"
```

5. test run jba executable on provided data
```{bash}

$ jba JaBbA/inst/extdata/junctions.vcf JaBbA/inst/extdata/coverage.txt 

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
JaBbA 2018-02-13 21:30:00: Starting analysis in ./jba_out
JaBbA 2018-02-13 21:32:13: Done .. job output in: ./jba_out

```


Usage (jba executable)
------------

```{bash}
Usage: jba [options] JUNCTIONS COVERAGE
	JUNCTIONS can be BND style vcf, bedpe, rds of GrangesList
 	COVERAGE is a .wig, .bw, .bedgraph, .bed., .rds of a granges, or .tsv  .csv /.txt  file that is coercible to a GRanges
	use --field=FIELD argument so specify which column to use if specific meta field of a multi-column table


Options:
	-s SEG, --seg=SEG
		Path to .rds file of GRanges object of intervals corresponding to initial segmentation (required)

	-a ABU, --abu=ABU
		Path to .rds, file of GRanges object of fine scale genomic coverage / abundance as tiled intervals (100 - 5000 bp) along genome (required)

	-r RA, --ra=RA
		Path to rearrangement file, which can be VCF breakend format, dRanger tab delim output,  or an rds of GRangesList of signed locus pairs pointing AWAY from junction (required)

	-j CFIELD, --cfield=CFIELD
		Name of meta data field of ra GRanges or data frame that specifies a junction confidence score, this will be used as a score to nudge (ie reward) each copy of that junction in the optimization

	-i TFIELD, --tfield=TFIELD
		Name of meta data field of ra GRanges or data frame that specifies tiers of junctions, where tier 1 is forced to be included

	-b NSEG, --nseg=NSEG
		Path to .rds file of GRanges object of intervals corresponding to normal tissue copy number, needs to have $cn field

	-d HETS, --hets=HETS
		Path to tab delimited hets file output of pileup with fields seqnames, start, end, alt.count.t, ref.count.t, alt.count.n, ref.count.n

	-l LIBDIR, --libdir=LIBDIR
		Directory containing karyoMIP.R file (eg default GIT.HOME/isva)

	-o OUTDIR, --outdir=OUTDIR
		Directory to dump output into (default jba_out)

	-n NAME, --name=NAME
		Sample / Individual name

	-f FIELD, --field=FIELD
		Name of meta data field or column to use for abundance / coverage signal from abundance / coverage soignal file

	-k SLACK, --slack=SLACK
		Slack penalty to apply per loose end copy

	-z SUBSAMPLE, --subsample=SUBSAMPLE
		Numeric value between 0 and 1 specifying whether to subsample coverage for intra segment variance estimation

	-t TILIM, --tilim=TILIM
		Time limit for JaBbA MIP

	-p PLOIDY, --ploidy=PLOIDY
		Ploidy guess

	-q PURITY, --purity=PURITY
		Purity guess

	-c CORES, --cores=CORES
		Number of cores for JaBBa MIP

	-m ITERATE, --iterate=ITERATE
		How many times to iterate through tiers

	-w WINDOW, --window=WINDOW
		Window to dumpster dive for junctions around loose ends

	-e EDGENUDGE, --edgenudge=EDGENUDGE
		Edge nudge for optimization, to be multiplied by edge specific confidence score if provided

	--ppmethod=PPMETHOD
		choose from sequenza, ppurple, or ppgrid to estimate purity ploidy

	--indel
		if TRUE will force the small isolated junctions in tier 2 have non-zero copy numbers

	--allin
		if TRUE will put all tiers in the first round of iteration

	--boolean
		if TRUE will use Boolean loose end penalty

	--epgap=EPGAP
		threshold for calling convergence

	-x STRICT, --strict=STRICT
		if TRUE will restrict input junctions to only the subset overlapping seg

	-u GUROBI, --gurobi=GUROBI
		if TRUE will use gurobi (gurobi R package must be installed) and if FALSE will use CPLEX (cplex must be installed prior to library installation)

	-v, --verbose
		verbose output

	-y, --nudgebalanced
		Manually nudge balanced junctions into the model.

	-h, --help
		Show this help message and exit

```

Output (R package)
------------
1. `jabba.simple[.rds|.png|.cnv.vcf|.gg.rds]`

   Main results, the optimized and simplified rearrangement graph. The four formats are R list object, PNG image of the graph generated by gTrack, VCF of the copy number variations, and
   gGraph object constructed with gGnome. In the list output, field "segstats" is the `GRanges` object of the nodes (including loose ends), field "adj" is the adjacency matrix, field "edges"
   is the edge table, field "gtrack" is the gTrack object used to generate the plot in the PNG file.

2. `karyograph.rds.ppfit.png`

   This plot illustrates the distribution of the raw segmental mean of the coverage signal, with red dashed vertical lines indicating the grid of integer copy number states. When the grid align well
   with the peaks in the underlying histogram, it indicates the purity/ploidy estimation is relatively successful.

3. `jabba.seg.txt`

   SEG format file of the final segmental copy numbers, compatible with IGV.

4. `opt.report.rds`

   This file contains an R `data.table` object of the convergence statistics of all the sub-problems (identified by "cl" column). The column "convergence" indicates the state of the final solution:
   - 1: converged quickly, within short time limit (input `tilim/10`) to the stringent epgap (input `epgap/1000`)
   - 2: converged roughly, within short time limit to the relaxed epgap (input `epgap`)
   - 3: converged after a second round, within long time limit to the relaxed epgap
   - 4: hardly converged after a second round, even after long time limit still above the relaxed epgap

   For detailed explanation of `tilim` and `epgap` please read our manuscript and CPLEX help doc.

Attributions
------------
> Marcin Imielinski - Assistant Professor, Weill Cornell Medicine
> Core Member, New York Genome Center.

> Xiaotong Yao - Graduate Research Assistant, Weill Cornell Medicine, New York
> Genome Center.


Funding sources
------------

<img
src="https://static1.squarespace.com/static/562537a8e4b0bbf0e0b819f1/5ad81984575d1f7d69517350/5ad819f02b6a28750f79597c/1524111879079/DDCF.jpeg?format=1500w"
height="150" class ="center"> <img
src="https://static1.squarespace.com/static/562537a8e4b0bbf0e0b819f1/5ad81984575d1f7d69517350/5ad819b8aa4a996c2d584594/1524111841815/BWF.png?format=500w"
height="150" class ="center">




```
