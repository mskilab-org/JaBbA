[![Build Status](https://travis-ci.org/mskilab/JaBbA.svg?branch=master)](https://travis-ci.org/mskilab/JaBbA)
[![codecov.io](https://img.shields.io/codecov/c/github/mskilab/JaBbA.svg)](https://codecov.io/github/mskilab/JaBbA?branch=master)

# JaBbA (Junction Balance Analysis)
```
 _____         ___    _      _____ 
(___  )       (  _`\ ( )    (  _  )
    | |   _ _ | (_) )| |_   | (_) |
 _  | | /'_` )|  _ <'| '_`\ |  _  |
( )_| |( (_| || (_) )| |_) )| | | |
`\___/'`\__,_)(____/'(_,__/'(_) (_)

(Junction     Balance     Analysis)

```

JaBbA builds a genome graph based on junctions and read depth from whole genome sequencing, inferring optimal copy numbers for both vertices (DNA segments) and edges (bonds between segments). It can be used for discovering various patterns of structural variations.

If you use JaBbA in your work, please cite: [Distinct Classes of Complex Structural Variation Uncovered across Thousands of Cancer Genome Graphs](https://doi.org/10.1016/j.cell.2020.08.006)

## Table of contents
- [JaBbA (Junction Balance Analysis)](#jabba-junction-balance-analysis)
	- [Table of contents](#table-of-contents)
	- [Installation](#installation)
	- [Usage](#usage)
	- [Output](#output)
	- [Best practice](#best-practice)
	- [FAQ](#faq)
	- [The design of JaBbA](#the-design-of-jabba)
	- [Attributions](#attributions)
	- [Funding sources](#funding-sources)
	- [Fun fact](#fun-fact)

## Installation
1. Install IBM ILOG
   [CPLEX Studio](https://www.ibm.com/products/ilog-cplex-optimization-studio).
   The software is proprietary, but can be obtained for free under [IBM's academic
   initiative](https://www.ibm.com/products/ilog-cplex-optimization-studio/pricing).

2. Set `CPLEX_DIR` variable to your CPLEX Studio installation directory

```{bash}
export CPLEX_DIR=/path/to/your/copy/of/CPLEX_Studio/
```
**NOTE: if `CPLEX_DIR` is set correctly then `$CPLEX_DIR/cplex/include` and `$CPLEX_DIR/cplex/lib` should both exist.**

3. Install JaBbA 

```{r}
devtools::install_github('mskilab/JaBbA')
```

4. For convenience, add `jba` executable to your `PATH`

```{bash}
$ JABBA_PATH=$(Rscript -e 'cat(paste0(installed.packages()["JaBbA", "LibPath"], "/JaBbA/extdata/"))')
$ export PATH=${PATH}:${JABBA_PATH}
$ jba ## to see usage
```

5. Test run jba executable on provided toy data
```{bash}
$ jba ${JABBA_PATH}/junctions.vcf ${JABBA_PATH}/coverage.txt 
```

## Usage
```
 _____         ___    _      _____ 
(___  )       (  _`\ ( )    (  _  )
    | |   _ _ | (_) )| |_   | (_) |
 _  | | /'_` )|  _ <'| '_`\ |  _  |
( )_| |( (_| || (_) )| |_) )| | | |
`\___/'`\__,_)(____/'(_,__/'(_) (_)

(Junction     Balance     Analysis)

Usage: jba JUNCTIONS COVERAGE [options]
	JUNCTIONS can be BND style vcf, bedpe, rds of GrangesList
 	COVERAGE is a .wig, .bw, .bedgraph, .bed., .rds of a granges, or .tsv  .csv /.txt  file that is coercible to a GRanges
	use --field=FIELD argument so specify which column to use if specific meta field of a multi-column table


Options:
	--j.supp=J.SUPP
		supplement junctions in the same format as junctions

	--blacklist.junctions=BLACKLIST.JUNCTIONS
		rearrangement junctions to be excluded from consideration

	--whitelist.junctions=WHITELIST.JUNCTIONS
		rearrangement junctions to be forced to be incorporated

	--geno
		whether the junction has genotype information

	--indel=INDEL
		character of the decision to 'exclude' or 'include' small(< min.nbins * coverage bin width) isolated INDEL-like events into the model. Default NULL, do nothing.

	--cfield=CFIELD
		junction confidence meta data field in ra

	--tfield=TFIELD
		tier confidence meta data field in ra. tiers are 1 = must use, 2 = may use, 3 = use only in iteration>1 if near loose end. Default 'tier'.

	--iterate=ITERATE
		the number of extra re-iterations allowed, to rescue lower confidence junctions that are near loose end. Default 0. This requires junctions to be tiered via a metadata field tfield.

	--window=WINDOW
		window size in bp within which to look for lower confidence junctions. Default 1000.

	--nudgebalanced=NUDGEBALANCED
		whether to attempt to add a small incentive for chains of quasi-reciprocal junctions.

	--edgenudge=EDGENUDGE
		numeric hyper-parameter of how much to nudge or reward aberrant junction incorporation. Default 0.1 (should be several orders of magnitude lower than average 1/sd on individual segments), a nonzero value encourages incorporation of perfectly balanced rearrangements which would be equivalently optimal with 0 copies or more copies.

	--strict
		if used will only include junctions that exactly overlap segs

	--allin
		if TRUE will put all tiers in the first round of iteration

	--field=FIELD
		name of the metadata column of coverage that contains the data. Default 'ratio' (coverage ratio between tumor and normal). If using dryclean, it is 'foreground'.

	-s SEG, --seg=SEG
		Path to .rds file of GRanges object of intervals corresponding to initial segmentation (required)

	--maxna=MAXNA
		Any node with more NA than this fraction will be ignored

	--blacklist.coverage=BLACKLIST.COVERAGE
		Path to .rds, BED, TXT, containing the blacklisted regions of the reference genome

	--nseg=NSEG
		Path to .rds file of GRanges object of intervals corresponding to normal tissue copy number, needs to have $cn field

	--hets=HETS
		Path to tab delimited hets file output of pileup with fields seqnames, start, end, alt.count.t, ref.count.t, alt.count.n, ref.count.n

	--ploidy=PLOIDY
		Ploidy guess, can be a length 2 range

	--purity=PURITY
		Purity guess, can be a length 2 range

	--ppmethod=PPMETHOD
		select from 'ppgrid', 'ppurple', and 'sequenza' to infer purity and ploidy if not both given. Default, 'sequenza'.

	--cnsignif=CNSIGNIF
		alpha value for CBS

	--slack=SLACK
		Slack penalty to apply per loose end

	--linear
		if TRUE will use L1 loose end penalty

	-t TILIM, --tilim=TILIM
		Time limit for JaBbA MIP

	--epgap=EPGAP
		threshold for calling convergence

	-o OUTDIR, --outdir=OUTDIR
		Directory to dump output into

	-n NAME, --name=NAME
		Sample / Individual name

	--cores=CORES
		Number of cores for JaBBa MIP

	-v, --verbose
		verbose output

	-h, --help
		Show this help message and exitsage: jba JUNCTIONS COVERAGE [options]
        JUNCTIONS can be BND style vcf, bedpe, rds of GrangesList
        COVERAGE is a .wig, .bw, .bedgraph, .bed., .rds of a granges, or .tsv  .csv /.txt  file that is coercible to a GRanges
        use --field=FIELD argument so specify which column to use if specific meta field of a multi-column table
```

## Output

1. `jabba.simple.gg.rds` 

	The out put genome graph with integer copy numbers. For more information of how to analyze it please refer to [gGnome tutorial](http://mskilab.com/gGnome/tutorial.html#applications)

2. `karyograph.rds.ppfit.png`

   This plot illustrates the distribution of the raw segmental mean of the coverage signal, with red dashed vertical lines indicating the grid of integer copy number states. When the grid align well with the peaks in the underlying histogram, it indicates the purity/ploidy estimation is relatively successful.

3. `jabba.seg.txt`

   SEG format file of the final segmental copy numbers, compatible with IGV/ABSOLUTE/GISTIC and many more.

4. `opt.report.rds`

   This file contains an R `data.table` object of the convergence statistics of all the sub-problems (identified by "cl" column). The column "convergence" indicates the state of the final solution:
   - 1= converged quickly, within short time limit (input `tilim/10`) to the stringent epgap (input `epgap/1000`)
   - 2= converged roughly, within short time limit to the relaxed epgap (input `epgap`)
   - 3= converged after a second round, within long time limit to the relaxed epgap
   - 4= hardly converged after a second round, even after long time limit still above the relaxed epgap

   For detailed explanation of `tilim` and `epgap` please read our manuscript and CPLEX help doc.

## Best practice
We are working on a best practice pipeline setting to take you from BAMs to good quality reconstructed genome graphs. For the time being please refer to FAQ for practical guidances.

## FAQ

1. "CPLEX Error 1016: Community Edition. Problem size limits exceeded."
	
	You are using a free trial version of CPLEX, please contact [IBM's academic initiative](https://www.ibm.com/products/ilog-cplex-optimization-studio/pricing) for a full license for academic use. We will work on supporting Gurobi in the future.

2. How to prepare the genome-wide coverage input?
	
	We used [*fragCounter*](https://github.com/mskilab/fragCounter) at
    200bp resolution on hg19 reference genome in our
    [paper](https://doi.org/10.1016/j.cell.2020.08.006) (details in
    STAR Methods), which summarized the ratio between the numbers of
    reads mapped to a bin in tumor versus normal sample, and then
    corrected for GC content and mappability. Recently we've adopted
    [*dryclean*](https://github.com/mskilab/dryclean) ([Deshpande et
    al., BioRxiv, 2019](https://doi.org/10.1101/847681)) to denoise
    coverage data using robust PCA with a panel of normal (PON). The
    default bin width is 1 kbp and the best practice protocol with
    dryclean is in preparation (expected Dec 2020).
	
	Please make sure that the `field` argument to JaBbA is set to the
    column name of the coverage data in your coverage input file.

3. How to prepare the junction input?
	
	We used [SvABA](https://github.com/walaj/svaba) ([Wala et al.,
    Genome Research, 2018](https://doi.org/10.1101/gr.221028.117)),
    and we support all junctions callers whose output either conforms
    to BND style in
    [VCF4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf)
    specifications or
    [BEDPE](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format)
    format for junctions. Some popular junction callers that do not
    conform to either are supported too, namely
    [Delly](https://github.com/dellytools/delly),
    [Lumpy](https://github.com/arq5x/lumpy-sv),
    [Novobreak](https://sourceforge.net/projects/novobreak/). Besides,
    the input format can also be a RDS file containing *GRangesList*
    of junctions.
	

	The *tfield* argument indicates the metadata column name in the
    junction input reflecting confidence in the call. Take advantage
    of this to customize your own consensus junction merging, e.g. set
    junctions called by all callers to tier 1 (highest confidence,
    JaBbA must use these junctions), called by at least 2 callers but
    not all to tier 2 (normal confidence, JaBbA will decide whether to
    use based on coverage change), and the junctions supported by only
    one caller to tier 3 (low confidence, JaBbA will only search for
    plausible candidates from them if close to a loose end, must set
    *iteration*>0).
	
	
4. How to prepare the segmentation input?

	Optional, as JaBbA will infer internally with
    [*CBS*](https://www.bioconductor.org/packages/release/bioc/html/DNAcopy.html),
    but if needed you can also provide tumor and/or normal
    segmentation through *seg* and *nseg* arguments.
	

5. How to prepare the purity and ploidy input?
   
	Optional, as JaBbA will infer internally with one of *ppgrid*,
    [*sequenza*](https://cran.r-project.org/web/packages/sequenza/index.html),
    and [*Ppurple*](https://github.com/mskilab/ppurple) of your choice
    through *ppmethod* argument, depending on the availability of
    necessary input files. Purity and ploidy estimation is arguably
    the most influential hyper-parameter of a JaBbA run, as it
    dictates the relationship between coverage data to copy number
    space, so if you want to make high quality graphs, start with
    better estiamtions. It is a very hard problem disguised as an easy
    one. Other external tools that helps: ABSOLUTE, ASCAT-wgs, TITAN.
	

6. How to choose slack penalty?
   
   The default used in the paper is 100, which in theory correspond to
   a prior belief that 1 in 100 breakends have loose end, or
   unexplained copy number change without consistent junction. In
   practice though, there are many source of noise in the coverage
   data that could mix with the desired signal of copy number
   change. After running JaBbA, plot the output copy number alongside
   the input coverage, If your output segmentation looks too "rigid",
   i.e. missing obvious clean copy number change points, you might
   want to consider dropping the slack penalty, as it indicates that
   JaBbA is so reluctant to add a loose end that it ignores the true
   signal from the coverage. For dryclean coverage input, which is at
   1kbp resolution and denoised, we currently recommend trying slack
   penalty around 20.
   

## The design of JaBbA
The key to understanding what JaBbA is doing lies in its name,
"balancing junctions". When analyzing SVs, junctions and segmental
copy numbers have been treated separately in most if not all large
scale WGS analyses, but what's not being addressed directly is that
these are just two measurable features of the same DNA sequence. The
structure of DNA tells us that it is a string, every segment in it are
just joined in tandem, hence every copy of a segment should have
exactly one upstream neighbor and one downstream neighbor. When adding
all copies a segment has, the simple rule that cannot be broken is
there must be same number of copies of up/downstream neighbors. If we
treat segments as vertices and the 3'-5' phosphodiester bond between
segments as edges, we get the *junction balance constraints* that
couple them together.


Of course there would be copy number change points where we can't find
a matching junction, and we fill them with loose ends. To make the
copy number more correct without a junction we have to use these loose
ends like placeholders so the junction balance constraint is still
met. Then the construction of the objective function is clear: we want
the segment copy number estiamte to be as close to the data's center
as possible (minimize residual sum of square) while limit the places
where we had to use loose ends for better fit (minimize the number of
loose ends).


## Attributions

> Marcin Imielinski - Assistant Professor, Weill Cornell Medicine
> Core Member, New York Genome Center.

> Xiaotong Yao - Graduate Research Assistant, Weill Cornell Medicine, New York
> Genome Center.

## Funding sources

<img
src="https://static1.squarespace.com/static/562537a8e4b0bbf0e0b819f1/5ad81984575d1f7d69517350/5ad819f02b6a28750f79597c/1524111879079/DDCF.jpeg?format=1500w"
height="150" class ="center"> <img
src="https://static1.squarespace.com/static/562537a8e4b0bbf0e0b819f1/5ad81984575d1f7d69517350/5ad819b8aa4a996c2d584594/1524111841815/BWF.png?format=500w"
height="150" class ="center">

## Fun fact
Why the name? Keith Bradnam started a parody award called [JABBA](http://www.acgt.me/blog/2020/3/31/three-cheers-for-jabba-awards) and we thought we'd won. However, jokes aside, did you notice the palindrome there? That's what would happen to a genome if it were to undergo BFB cycles.
