[![Build Status](https://travis-ci.org/mskilab/JaBbA.svg?branch=master)](https://travis-ci.org/mskilab/JaBbA)
[![codecov.io](https://img.shields.io/codecov/c/github/mskilab/JaBbA.svg)](https://codecov.io/github/mskilab/JaBbA?branch=master)

# JaBbA (Junction Balance Analysis)

JaBbA builds a genome graph based on junctions and fits integer copy numbers for both vertices (DNA segments) and edges (bonds between segments) from whole genome sequencing data. The two required inputs are junctions and binned read depth coverage.

## Citation

If you use JaBbA in your work, please cite: 
 
## Installation
------------
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

        --rescue.window=RESCUE.WINDOW
                window size in bp within which to look for lower confidence junctions. Default 1000.

        --rescue.all
                Attempt to rescue all loose ends regardless of internal quality filter

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
                Show this help message and exit
```
## Output (R package)
------------
1. `jabba.simple.gg.rds` 

	The out put genome graph with integer copy numbers. For more information of how to analyze it please refer to [gGnome tutorial](http://mskilab.com/gGnome/tutorial.html#applications)

2. `jabba.simple[.rds|.png|.cnv.vcf]`

   Main results, the optimized and simplified rearrangement graph. The four formats are R list object, PNG image of the graph generated by gTrack, VCF of the copy number variations, and gGraph object constructed with gGnome ([learn more about gGnome here](http://mskilab.com/gGnome/tutorial.html#introduction)). In the list output, field "segstats" is the `GRanges` object of the nodes (including loose ends), field "adj" is the adjacency matrix, field "edges" is the edge table, field "gtrack" is the gTrack object used to generate the plot in the PNG file.

3. `karyograph.rds.ppfit.png`

   This plot illustrates the distribution of the raw segmental mean of the coverage signal, with red dashed vertical lines indicating the grid of integer copy number states. When the grid align well with the peaks in the underlying histogram, it indicates the purity/ploidy estimation is relatively successful.

4. `jabba.seg.txt`

   SEG format file of the final segmental copy numbers, compatible with IGV/ABSOLUTE/GISTIC and many more.

5. `opt.report.rds`

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
