#############################################################################
## Marcin Imielinski
## The Broad Institute of MIT and Harvard / Cancer program.
## marcin@broadinstitute.org
##
## Weill-Cornell Medical College
## mai9037@med.cornell.edu
##
## New York Genome Center
## mimielinski@nygenome.org
##
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.

## You should have received a copy of the GNU Lesser General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

#' @import gTrack
#' @import GenomicRanges
#' @import igraph
#' @import Matrix
#' @import parallel
#' @import data.table
#' @import Rcplex
#' @import DNAcopy
#' @import gUtils

##
##
## to install igraph, may need "use GCC-4.8"
## installing RCplex need to run:
##
## install.packages('Rcplex', configure.args = c("--with-cplex-include=/broad/software/free/Linux/redhat_5_x86_64/pkgs/cplex_12.2/cplex/include/", "--with-cplex-lib='-L/broad/software/free/Linux/redhat_5_x86_64/pkgs/cplex_12.2/cplex/lib/x86-64_sles10_4.1/static_pic/ -lcplex -lm -lpthread'", "--with-cplex-cflags='-m64 -fPIC -I/usr/include'"))
##
## install.packages('Rcplex', configure.args = c("--with-cplex-include=/broad/software/free/Linux/redhat_5_x86_64/pkgs/cplex_12.6/cplex/include/", "--with-cplex-lib='-L/broad/software/free/Linux/redhat_5_x86_64/pkgs/cplex_12.6/cplex/lib/x86-64_linux/static_pic/ -lcplex -lm -lpthread'", "--with-cplex-cflags='-m64 -fPIC -I/usr/include'"))
##
## install.packages('Rcplex', configure.args = c("--with-cplex-include=/nfs/home/mimielinski/Software/CPLEX/CPLEX_Studio/cplex/include/", "--with-cplex-lib='-L/nfs/home/mimielinski/Software/CPLEX/CPLEX_Studio/cplex/lib/x86-64_linux/static_pic/ -lcplex -lm -lpthread'", "--with-cplex-cflags='-m64 -fPIC -I/usr/include'"))
## 
##
## to install igraph, need a "use GCC-4.8"
##
##

#' @name JaBbA
#' @title JaBbA
#' @details
#' Module to run jbaMIP + preprocessing from text file or rds input and dump files out to text.
#' 
#' Generates the following files in the output directory:
#' 
#' karyograph.rds --- file of unpopulated karyograph as an RDS file of a list object storing the output of karyograph
#' 
#' 
#' jabba.rds --- file storing JaBbA object
#' 
#' jabba.simple.rds --- file storing JaBbA object simplified so that segments containing all unpopulated aberrant junctions are merged
#' 
#' jabba.raw.rds --- storing raw jbaMIP solution, this may be useful for debugging and QC
#' 
#' jabba.png, jabba.simple.png --- gTrack images of the above reconstructions
#' 
#' jabba.seg.txt --- tsv file with jabba.simple solution segments
#' 
#' jabba.seg.rds --- GRanges rds with jabba.simple solution segments
#' 
#' jabba.adj.txt --- tsv file with edges (i.e. node pairs) of adjacency matrix populated with inferred copy numbers and node ids indexing segments in jabba.seg.txt
#' 
#' jabba.vcf, jabba.simple.vcf --- BND-style vcf output of junctions in JaBbA output populated with rearrangement and interval copy numbers
#' 
#' jabba.cnv.vcf, jabba.simple.cnv.vcf --- cfopy number style VCF showing jabba copy number output
#' 
#' 
#' @param ra  GRangesList of junctions  (i.e. bp pairs with strands oriented AWAY from break) OR path to junction VCF file (BND format), dRanger txt file or rds of GRangesList 
#' @param coverage  GRanges of coverage OR path to cov file, rds of GRanges or .wig / .bed file of (normalized, GC corrected) fragment density
#' @param field  field of coverage GRanges to use as fragment density signal (only relevant if coverage is GRanges rds file)
#' @param seg  optional path to existing segmentation, if missing then will segment abu using DNACopy with standard settings 
#' @param cfield  character, junction confidence meta data field in ra
#' @param tfield  character, tier confidence meta data field in ra
#' @param outdir  out directory to dump into, default ./
#' @param nseg  optional path to normal seg file with $cn meta data field
#' @param hets  optional path to hets.file which is tab delimited text file with fields seqnames, start, end, alt.count.t, ref.count.t, alt.count.n, ref.count.n
#' @param name  prefix for sample name to be output to seg file
#' @param cores  number of cores to use (default 1)
#' @param nseg  path to data.frame or GRanges rds of normal seg file with coordinates and $cn data field specifying germline integer copy number
#' @param subsample  numeric between 0 and 1 specifying how much to sub-sample high confidence coverage data
#' @param tilim  timeout for jbaMIP computation (default 1200 seconds)
#' @param edgenudge  numeric hyper-parameter of how much to nudge or reward aberrant junction incorporation, default 0.1 (should be several orders of magnitude lower than average 1/sd on individual segments), a nonzero value encourages incorporation of perfectly balanced rearrangements which would be equivalently optimal with 0 copies or more copies.
#' @param slack.penalty  penalty to put on every loose.end copy, should be calibrated with respect to 1/(k*sd)^2 for each segment, i.e. that we are comfortable with junction balance constraints introducing k copy number deviation from a segments MLE copy number assignment (the assignment in the absence of junction balance constraints)
#' @param overwrite  flag whether to overwrite existing output directory contents or just continue with existing files.
#' @export
#' @import DNAcopy
############################################
JaBbA = function(
    ra, # path to junction VCF file, dRanger txt file or rds of GRangesList of junctions (with strands oriented pointing AWAY from breakpoint)
    coverage, # path to cov file, rds of GRanges
    seg = NULL, # path to seg file, rds of GRanges
    cfield = NULL, # character, junction confidence meta data field in ra
    tfield = NULL, # character, tier confidence meta data field in ra
    outdir = './', # out directory to dump into
    nseg = NULL, # path to normal seg file with $cn meta data field
    hets = NULL, # path to hets.file which is tab delimited text file with fields seqnames, start, end, alt.count.t, ref.count.t, alt.count.n, ref.count.n
    name = 'tumor', ## prefix for sample name to be output to seg file
    cores = 4, # default 1 
    field = 'ratio', ## character, meta data field to use from coverage object to indicate numeric coveragendance, coverage,
    subsample = NULL, ## numeric scalar between 0 and 1, how much to subsample coverage per segment 
    tilim = 1200, ## timeout for MIP portion
    edgenudge = 0.1, ## hyper-parameter of how much to "nudge" or reward edge use, will be combined with cfield information if provided
    slack.penalty = 1e3, ## nll penalty for each loose end copy
    overwrite = F ## whether to overwrite existing output
)
{
    kag.file = paste(outdir, 'karyograph.rds', sep = '/')
    junctions.txt.file = paste(outdir, 'junctions.txt', sep = '/')
    junctions.rds.file = paste(outdir, 'junctions.rds', sep = '/')
    jabba.raw.rds.file = paste(outdir, 'jabba.raw.rds', sep = '/')
    jabba.rds.file = paste(outdir, 'jabba.rds', sep = '/')
    jabba.vcf.file = paste(outdir, 'jabba.vcf', sep = '/')
    jabba.cnv.vcf.file = paste(outdir, 'jabba.cnv.vcf', sep = '/')
    jabba.simple.rds.file = paste(outdir, 'jabba.simple.rds', sep = '/')
    jabba.simple.vcf.file = paste(outdir, 'jabba.simple.vcf', sep = '/')
    jabba.simple.cnv.vcf.file = paste(outdir, 'jabba.simple.cnv.vcf', sep = '/')
    jabba.png.file = paste(outdir, 'jabba.png', sep = '/')
    jabba.simple.png.file = paste(outdir, 'jabba.simple.png', sep = '/')
    seg.tab.file = paste(outdir, 'jabba.seg.txt', sep = '/')
    seg.gr.file = paste(outdir, 'jabba.seg.rds', sep = '/')
    seg.adj.file = paste(outdir, 'jabba.adj.txt', sep = '/')
    nozzle.file = paste(outdir, 'nozzle', sep = '/')


    message('Reading in coverage')
    if (is.character(coverage))
    {
        if (grepl('\\.rds$', coverage))
        {
            cov = readRDS(coverage)
        }
        else
        {
            message('Importing seg from UCSC format')
            cov = import.ucsc(coverage)
            field = 'score';
            cov = gr.fix(cov)
        }
    }
    else
        cov = coverage

    
    if (!(field %in% names(values(cov))))
    {
        new.field = names(values(cov))[1]
        warning(paste0('Field ', field, ' not found in coverage GRanges metadata so using ', new.field, ' instead'))
        field = new.field
    }
    
    if (is.null(seg))
    {
        message('No segmentation provided via seg variable, so performing segmentation using CBS')
        vals = values(cov)[, field]
        ix = which(!is.na(vals))
        cna = CNA(log(vals[ix]), as.character(seqnames(cov))[ix], start(cov)[ix], data.type = 'logratio')
        cat('finished making cna\n')
        seg = segment(smooth.CNA(cna), alpha = 1e-5, verbose = T)
        cat('finished segmenting\n')
        seg = seg2gr(print(seg), seqlengths(cov)) ## remove seqlengths that have not been segmented
        seg = gr.fix(seg, seqlengths(cov), drop = T)
        cat(length(seg), ' segments produced\n')
        names(seg) = NULL
    }
    else
    {
        if (is.character(seg))
        {
            if (grepl('\\.rds$', seg))
            {
                cov = readRDS(seg)
            }
            else
            {
                message('Importing seg from UCSC format')
                cov = import.ucsc(seg)
                field = 'score';
            }
        }
        else
            cov = seg
    }
       
    
    if (!is.null(hets))
        if (!file.exists(hets))
            {
                warning(sprintf('hets file "%s" not found, ignoring hets\n'))
                hets = NULL
            }

    if (!is.character(ra))
        {
            if (overwrite | !file.exists(kag.file))
                karyograph_stub(seg, coverage, ra = ra, out.file = kag.file, nseg.file = nseg, field = field, subsample = subsample, het.file = hets)
            else
                warning("Skipping over karyograph creation because file already exists and overwrite = FALSE")
        }
    else
    {
        if (grepl('rds$', ra) | grepl('vcf$', ra) | grepl('vcf\\.gz$', ra))
        {
            if (overwrite | !file.exists(kag.file))
                karyograph_stub(seg, coverage, ra.file = ra, out.file = kag.file, nseg.file = nseg, field = field, subsample = subsample, het.file = hets)
            else
                warning("Skipping over karyograph creation because file already exists and overwrite = FALSE")
        } else  {
            if (overwrite | !file.exists(kag.file))
                karyograph_stub(seg, coverage, dranger.file = gsub('all.mat$', 'somatic.txt', ra), nseg.file = nseg, out.file = kag.file, field = field, subsample = subsample, het.file = hets)
            else
                warning("Skipping over karyograph creation because file already exists and overwrite = FALSE")
        }
    }
        
    kag = readRDS(kag.file)
    
    if (!is.null(cfield))
        {  
            if (cfield %in% names(values(kag$junctions)))
                {
                    val = values(kag$junctions)[, cfield]
                    val[is.na(val)] = 0      
                    edgenudge = val * edgenudge
                }
        }
    
    ab.force = NULL
    if (!is.null(tfield))
        {
            if (tfield %in% names(values(kag$junctions)))
                {
                    ab.force = which(gsub('tier', '', as.character(values(kag$junctions)[, tfield]))=='1')
                    cat('Found tier field enforcing >=1 CN at ', length(ab.force), 'junctions\n')
                }
        }

    gc()
    
    if (overwrite | !file.exists(jabba.raw.rds.file))
        ramip_stub(kag.file, jabba.raw.rds.file, mc.cores = cores, tilim = tilim, edge.nudge = edgenudge, ab.force = ab.force, mem = 5, slack.prior = 1/slack.penalty)
    
    kag = readRDS(kag.file)
    jab = readRDS(jabba.raw.rds.file)

    if (!file.exists(jabba.rds.file))  
        jabd = JaBbA.digest(jab, kag)
    else
        jabd = readRDS(jabba.rds.file)
    
    jabd$purity = jab$purity
    jabd$ploidy = jab$ploidy

    if (!file.exists(jabba.simple.rds.file))
        jabd.simple = JaBbA.digest(jab, kag, keep.all = F) ## simplified
    else
        jabd.simple = readRDS(jabba.simple.rds.file)

    jabd.simple$purity = jab$purity
    jabd.simple$ploidy = jab$ploidy
    junctions = kag$junctions
    values(junctions)$cn = jab$adj[rbind(kag$ab.edges[, 1:2, 1])]
    jabd.simple$junctions = jabd$junctions = jab$junctions = junctions
    
    jab$ab.edges = kag$ab.edges
    seg.out = cbind(sample = name, as.data.frame(jabd$segstats))
    names(seg.out)[1:4] = c('track.name', 'chrom', 'start', 'end')
    seg.out$seg.id = 1:nrow(seg.out)
    cols = c('track.name', 'chrom', 'start', 'end', 'cn', 'seg.id')
    seg.out = seg.out[, c(cols, setdiff(names(seg.out), cols))]
    write.tab(seg.out, seg.tab.file)
    jabd$segstats$seg.id = 1:length(jabd$segstats)

    cat('Checking for hets\n')
    if (!is.null(hets))
        if (file.exists(hets))
            tryCatch(
                {
                    cat('Loading hets\n')
                    hetsd = fread(hets)
                    hets.gr = seg2gr(hetsd[pmin(ref.frac.n, 1-ref.frac.n) > 0.2 & (ref.count.n + alt.count.n)>20, ])
                    cat('Computing alleles for jabd\n')
                    jabd = c(jabd, jabba.alleles(jabd, hets.gr, verbose = TRUE, uncoupled=TRUE)[c('asegstats', 'aadj', 'atd')])
                    cat('Computing alleles for jabd simple\n')
                    jabd.simple = c(jabd.simple, jabba.alleles(jabd.simple, hets.gr, verbose = TRUE, uncoupled=TRUE)[c('asegstats', 'aadj', 'atd')])
                    cat('Done computing alleles\n')
                },
                error = function(e) print("Jabba allelic generation failed"))
    
    saveRDS(jabd$segstats, seg.gr.file)
    saveRDS(jab, jabba.raw.rds.file)
    saveRDS(jabd, jabba.rds.file)
    saveRDS(jabd.simple, jabba.simple.rds.file)

    tryCatch(
        {
            jabba2vcf(jabd, jabba.vcf.file);
            jabba2vcf(jabd, jabba.cnv.vcf.file, cnv = TRUE)
            jabba2vcf(jabd.simple, jabba.simple.vcf.file)
            jabba2vcf(jabd.simple, jabba.simple.cnv.vcf.file, cnv = TRUE)            
        }, error = function(e) print("Jabba VCF generation failed"))

            
    seg.adj = cbind(data.frame(sample = rep(name, nrow(jabd$edges))), jabd$edges[, c('from', 'to', 'cn', 'type')])
    write.tab(seg.adj, seg.adj.file)

    values(kag$junctions)$cn.jabba = jab$adj[rbind(jab$ab.edges[, 1:2, 1])]

    if (length(kag$junctions)>0)
        {
            tmp = grl.pivot(kag$junctions)
            names(tmp[[1]]) = 1:length(kag$junctions)
            names(tmp[[2]]) = 1:length(kag$junctions)
            ra1 = as.data.frame(tmp[[1]])
            ra2 = as.data.frame(tmp[[2]])
            names(ra1) = paste('bp1_', names(ra1), sep = '')
            names(ra2) = paste('bp2_', names(ra2), sep = '')
            junc.txt = as.data.frame(values(kag$junctions))
            write.tab(cbind(ra1, ra2, junc.txt), junctions.txt.file)
        }
    else
        writeLines(c("\t"), junctions.txt.file)

    saveRDS(kag$junctions, junctions.rds.file)
    
    tmp.cov = sample(cov, pmin(length(cov), 5e5))
    tmp.cov = gr.fix(tmp.cov, jabd$segstats)
    
    y1 = pmax(5, max(jabd$segstats$cn)*1.1)
    jabd$td$y1 = y1
    jabd.simple$td$y1 = y1
    
    td.cov = gTrack(tmp.cov, y.field = field, col = alpha('black', 0.2), name = 'Cov', y1 = (y1 + jab$gamma)/jab$beta)

    cat('Generating figures\n')

    if (!file.exists(jabba.png.file))
        {
            if (is.character(tryCatch(png(jabba.png.file, width = 2000, height = 1000), error = function(e) 'bla')))
                pdf(gsub('png$', 'pdf', jabba.png.file), width = 10, height = 10)

            if (is.null(jabd$atd))
                plot(c(td.cov, jabd$td))
            else
                plot(c(jabd$atd, td.cov, jabd$td))
            
            dev.off()
        }
    
    if (!file.exists(jabba.simple.png.file))
        {
            if (is.character(tryCatch(png(jabba.simple.png.file, width = 2000, height = 1000), error = function(e) 'bla')))
                pdf(gsub("png$", "pdf", jabba.simple.png.file), width = 20, height = 10)
            
            if (is.null(jabd.simple$atd))
                plot(c(td.cov, jabd.simple$td))
            else
                plot(c(jabd.simple$atd, td.cov, jabd.simple$td))
            
            dev.off()      
        }

    return('done')
}


#################
# ra_tier
#
# Classify full set of dRanger rearrangements into "tiers" of confidence
#
# (1) Tier 1 BPresult>0 and somatic_score>min.score1
# (2) Tier 2 BPresult=0 and somatic_score>min.score1
# (3) Tier 3 min.score2<=somatic_score<=min.score2 & tumreads>min.reads
#
##################
ra_tier = function(dra, min.score1 = 10, min.score2 = 4, min.treads1 = 10, min.treads2 = 3, max.nreads = Inf)
  {
    if (is(dra, 'GRangesList'))
      dra = values(dra)         
    dra$BPresult[is.na(dra$BPresult)] = -1
    dra$T_SWreads[is.na(dra$T_SWreads)] = 0
    dra$N_SWreads[is.na(dra$N_SWreads)] = 0
    tier = rep(3, nrow(dra))
    tier[dra$BPresult==1 & dra$somatic_score>=min.score1 & (dra$tumreads >= min.treads1 | dra$T_SWreads > min.treads1 | dra$T_BWAreads > min.treads1) & dra$normreads==0 & dra$N_SWreads < 0] = 1
    tier[dra$BPresult>=0 & tier!=1 & dra$normreads == 0 & 
         ((dra$tumreads + dra$T_SWreads) >= min.treads2 |
          (dra$somatic_score >= min.score2 | (dra$tumreads >= min.treads2))
         )] = 2
           
    return(tier)
  }


##########
#' karyograph
#'
#' @details
#' builds graph from rearrangement breakpoints +/- copy number endpoints
#' used for downstream jbaMIP and karyoMIP functions
#' 
#' Input bpp is a GRangesList of signed locus pairs describing aberrant adjacencies.
#' The convention is as follows: Each locus in the input breakpoint pair points to the direction that
#' is being joined by the adjacencies i.e.
#' (-) bp points to "left" or preceding segment 
#' (+) bp points to the  "right" or the following segment
#' 
#' eg imagine a|bp1|b
#'            c|bp2|d
#'  "+" bp point to the right (eg b or d), "-" bp point to the left (a or c)
#' 
#' Input "tile" is a set of intervals whose endpoints are also used to partition the genome prior to the building of the
#' karyograph.
#' 
#' Output karyograph connects signed genomic intervals (in a signed tiling of the reference genome) with "aberrant" and "reference" edges.
#' Reference edges connect intervals that are adjacent in the reference genome, and aberrant edges are inferred (upstream
#' of this) through cancer genome paired end analysis.
#' Note that every node, edge, and path in this karyograph has a "reciprocal path" 
#' 
#' @param junctions GRangesList of junctions, where each item is a length GRanges of signed locations
#' @param tile GRanges optional existing tiling of the genome (eg a copy number segmentation) from which additional segments will be created
#' @return
#'  a list with the following fields
#' $tile = GRanges of length 2*n tiling of the genome corresponding to union of rearrangement breakpoints and copy number endpoints
#' $G = igraph object representing karyograph, here are the edge and vertex features 
#'      vertex features: $chrom, $start, $end, $width, $strand, $size, $shape, $border.width, $label, $chrom.ord, $y, $col, $weight
#'      edge features:$ $bp.id, $weight, $from, $to, $col, $type, $line.style, $arrow.shape, $width, 
#'      important:  $type specifies which edges are "aberrant" and "reference", $bp.id specifies which input rearrangement (item in junctions)
#'      a given aberrant edge came from (and is NA for reference edges)
#' $adj = 2n x 2n adjacency matrix whose nonzero entries ij show the edge.id in $G
#' $ab.adj = 2n x 2n binary matrix specifying aberrant edges
#' $ab.edges = length(junctions) x {'from', 'to'} x {'+', '-'} mapping junction id's (indices into input junctions lists) to source and sink vertices,
#'             in both orientations
#' @export
############################################
karyograph = function(junctions, ## this is a grl of breakpoint pairs (eg output of ra_breaks(dranger.df) where dranger is df of dranger output)
    tile = NULL, # pre-existing set of intervals on top of which to build a graph (eg endpoints from a copy number based segmentation)
    label.edges = FALSE        
  )
  {
      require(gplots)
      require(igraph)
      require(Matrix)
      require(data.table)
      require(RColorBrewer)

    if (length(junctions)>0)
      {
        bp.p = grl.pivot(junctions)        
        bp1 = gr.end(gr.fix(bp.p[[1]]), 1, ignore.strand = F)
        bp2 = gr.start(gr.fix(bp.p[[2]]), 1, ignore.strand = F)


        #' mimielinski Sunday, Aug 06, 2017 06:46:15 PM
        #' fix added to handle strange [0 0] junctions outputted 
        #' by Snowman ... which failed to match to any tile
        #' todo: may want to also handle junctions that
        #' fall off the other side of the chromosome
        end(bp1) = pmax(1, end(bp1))
        start(bp1) = pmax(1, start(bp1))        
        end(bp1) = pmax(1, end(bp1))
        start(bp1) = pmax(1, start(bp1))
    
        if (any(as.logical(strand(bp1) == '*') | as.logical(strand(bp2) == '*')))
          stop('bp1 and bp2 must be signed intervals (i.e. either + or -)')
        
        if (length(bp1) != length(bp2))
          stop('bp1 and bp2 inputs must have identical lengths')
        
                                        #    if (sum(width(reduce(bp1))) != sum(width(bp1)) | sum(width(reduce(bp2))) != sum(width(bp2)))
                                        #      stop('bp1 or bp2 cannot have duplicates / overlaps (with respect to location AND strand)')
        
        values(bp1)$bp.id = 1:length(bp1);    
        values(bp2)$bp.id = 1:length(bp1)+length(bp1);
    
        pgrid = sgn1 = c('-'=-1, '+'=1)[as.character(strand(bp1))]
        sgn2 = c('-'=-1, '+'=1)[as.character(strand(bp2))]
        
        ### HACK HACK to force seqlengths to play with each other if malformedo
        tmp.sl = seqlengths(grbind(bp1, bp2))
        tmp.sl.og = tmp.sl
                                        #        tmp.sl = gr2dt(grbind(bp1, bp2))[, max(end, na.rm = TRUE), keyby = seqnames][, sl := pmax(V1+2, tmp.sl[as.character(seqnames)], na.rm = TRUE)][, structure(sl, names = as.character(seqnames))]
        tmp.sl = gr2dt(grbind(bp1, bp2))[, max(end+1, na.rm = TRUE), keyby = seqnames][names(tmp.sl.og), structure(pmax(V1, tmp.sl.og, na.rm = TRUE), names = names(tmp.sl.og))]
        bp1 = gr.fix(bp1, tmp.sl)
        bp2 = gr.fix(bp2, tmp.sl)
    # first we tile the genome around the combined breakpoints
      }
    else
      {
        if (is.null(tile))
            {
                tile = si2gr(junctions)
                if (length(tile)==0)
                    {
                        warning('Empty input given, producing empty output')
                        return(NULL)
                    }
                A = sparseMatrix(1,1, x = 0, dims = rep(length(tile), 2))
                return(
                    list(tile = tile, adj = A,
                         G = graph.adjacency(A), ab.adj = A != 0, ab.edges = NULL, junctions = junctions))
            }
        
        junctions = GRangesList()
        bp1 = bp2 = GRanges()
      }

    if (!is.null(tile))
      {
        ## find disjoint union of tile and join with gaps
        tile = gr.fix(tile)
        tile = gr.fix(tile, bp1)
        bp1 = gr.fix(bp1, tile) ## argh argh argh .. more pain avoiding hacks
        strand(tile) = '+'
        tile = disjoin(tile)
        tile = sort(c(tile, gaps(tile)))

        ## make sure seqlevels / seqinfo are identical
        if (!identical(sort(seqlevels(tile)), seqlevels(junctions)))
          {
            tile = gr.fix(tile, junctions)
            junctions = gr.fix(junctions, tile)
          }

        if(length(junctions)>0)
            {
                tbp = setdiff(gr.stripstrand(gr.trim(tile, 1)), gr.stripstrand(grbind(bp1, bp2)))
                bp1 = gr.fix(bp1, tbp)
                bp2 = gr.fix(bp2, tbp) ## seqlengths pain
                tbp = gr.fix(tbp, bp1)
            }
        else
          tbp = gr.stripstrand(gr.trim(tile, 1))

        tbp = tbp[start(tbp)!=1]
        
        if (length(tbp)>0)
            tbp$seg.bp = TRUE
      }
    else
        tbp = NULL;
      
      if (length(junctions)>0)
          if (length(tbp)>0)
              g = gaps(gr.stripstrand(sort(c(bp1[, c()], bp2[, c()], tbp[, c()]))))
          else
              g = gaps(gr.stripstrand(sort(c(bp1[, c()], bp2[, c()]))))
    else
      g = gaps(gr.stripstrand(sort(tbp)));

    g = g[as.logical( strand(g)=='*' )];
    strand(g) = '+';
    
    values(g)$bp.id = NA
    values(g)$seg.bp = NA

      ## combine tiles and find disjoint set
      tile.og = tile
    tile = grbind(bp1, bp2, g, tbp);
    tile = disjoin(gr.stripstrand(tile[order(gr.stripstrand(tile))]))
    strand(tile) = '+'
    tile = gr.fix(tile);
    tile$is.tel = start(tile)==1 | end(tile) == seqlengths(tile)[as.character(seqnames(tile))]
      values(tile)$tile.id = 1:length(tile);      

    # find "breakpoint" i.e. bp associated intervals, i.e. width 1 intervals that end with a bp1 or bp2 location   
    junc.bp = grbind(bp1, bp2)
    junc.bpix = numeric()
    if (length(junc.bp)>0)
      junc.bpix = which(paste(seqnames(tile), end(tile)) %in% paste(seqnames(junc.bp), start(junc.bp)))

      ## make sure all seqlenths are compatible (so effing annoying)
      tile = gr.fix(tile, bp1)
      tile = gr.fix(tile, bp2)
      bp1 = gr.fix(bp1, tile)
      bp2 = gr.fix(bp2, tile)
      
    
    ## also keep track of tbp associatd bp.ix
    all.bp = grbind(bp1, bp2, tbp)
    all.bpix = numeric()

    if (length(all.bp)>0)
      all.bpix = which(paste(seqnames(tile), end(tile)) %in% paste(seqnames(all.bp), start(all.bp)))

    ## now to build the graph, we would like to fuse all the bp associated intervals with their previous interval
    ## UNLESS they are preceded by another bp associated interval
    ##
    if (length(all.bpix>0))
      {
        to.fuse = all.bpix[which(all.bpix>1 & !((all.bpix-1) %in% all.bpix))]
        end(tile)[to.fuse-1] = end(tile)[to.fuse-1]+1
        tile = tile[-to.fuse]
      }
    
    if (length(junc.bpix)>0)
      {        
        ## we have a partition of genomic segments flanked by tile endpoints and/or ra junctions
        ##
        ## Input junction syntax is interpreted as follows:
        ## a- b+ junctions connect seg ending with position a to seg starting with b+1
        ## a- b- junctions connect seg ending with position a to seg ending with position b (on neg strand)
        ## a+ b+ junctions connect seg starting with position a+1 (on negative strand) to seg starting with position b+1
        ## a+ b- junctions connect seg starting with position a+1 (on negative strand) to seg ending with position b (on neg strand)
        
        # collect all pairwise adjacencies implied by breakpoints
        # eg imagine a|bp1|b
        #            c|bp2|d
        # "+" bp point to the right (eg b or d), "-" bp point to the left (a or c)
          
          ab.pairs = cbind(
              ifelse(as.logical(strand(bp1)=='+'), gr.match(GenomicRanges::shift(gr.start(bp1), 1), gr.start(tile)), 
                     gr.match(gr.start(bp1), gr.end(tile))),
              ifelse(as.logical(strand(bp2)=='+'), gr.match(GenomicRanges::shift(gr.start(bp2), 1), gr.start(tile)),
                     gr.match(gr.start(bp2), gr.end(tile)))                    
              )

          
        ## ab.pairs = cbind(
        ##   ifelse(as.logical(strand(bp1)=='+'), match(paste(seqnames(bp1), start(bp1)+1), paste(seqnames(tile), start(tile))),
        ##          match(paste(seqnames(bp1), start(bp1)), paste(seqnames(tile), end(tile)))),
        ##   ifelse(as.logical(strand(bp2)=='+'), match(paste(seqnames(bp2), start(bp2)+1), paste(seqnames(tile), start(tile))),
        ##          match(paste(seqnames(bp2), start(bp2)), paste(seqnames(tile), end(tile))))          
        ##   )
        ab.pairs.bpid = bp1$bp.id    
        pp = (sgn1*sgn2)>0 & sgn1>0;
        mm = (sgn1*sgn2)>0 & sgn1<0;
        mp = sgn1>0 & sgn2<0
        ab.pairs[pp,1] = -ab.pairs[pp,1] # ++ breakpoints --> (-b)d adjacency
        ab.pairs[mm,2] = -ab.pairs[mm,2] # -- breakpoints --> a(-c) adjacency
        ab.pairs[mp, ] = -ab.pairs[mp, ] # +- breakpoints --> (-b)(-c) adjacency

        # clean up adj pairs
        # remove any that have crossed a chromosome boundary from their breakpoint
        # this will occur in cases of badly formed breakpoint input (eg breakpoints that point outward
        # from their telomeres)
        edge.id = rep(1:nrow(ab.pairs), 2)
        ab.pairs = rbind(ab.pairs, cbind(-ab.pairs[,2], -ab.pairs[,1]));
        ab.pairs.bpid = c(ab.pairs.bpid, ab.pairs.bpid)
    
        # build "aberrant" adjacency matrix representing directed graph of edges connecting
        # <signed> nodes.
        # note: indices of matrix represent edge labels
          adj.ab = Matrix( 0,
                          nrow = 2*length(tile),
                          ncol = 2*length(tile),
                          dimnames = rep( list( as.character(c(1:length(tile), -(1:length(tile))))), 2))
          tmp.ix = cbind(match(as.character(ab.pairs[,1]), rownames(adj.ab)),
          match(as.character(ab.pairs[,2]), colnames(adj.ab)))        
        adj.ab[tmp.ix[!duplicated(tmp.ix), , drop = F]] = ab.pairs.bpid[!duplicated(tmp.ix)]
      }
    else
      {
        ab.pairs.bpid = edge.id = c()
        ab.pairs = matrix(nrow = 0, ncol = 2);        
        adj.ab = Matrix(FALSE, nrow = 2*length(tile), ncol = 2*length(tile),
          dimnames = rep(list(as.character(c(1:length(tile), -(1:length(tile))))), 2))
      }
    
    # build reference adjacency matrix (representing consecutive segments on the reference genome)
    # note: indices of matrix represent edge labels
    seg.ix = 1:length(tile)
    ref.pairs = cbind(seg.ix[1:(length(seg.ix)-1)], seg.ix[2:(length(seg.ix))])
    # ref.pairs = ref.pairs[ref.pairs[,1]>0 & ref.pairs[,2]!=length(tile), ]
    ref.pairs = ref.pairs[which(as.character(seqnames(tile[ref.pairs[,1]])) == as.character(seqnames(tile[ref.pairs[,2]]))), ]

    if (nrow(ref.pairs)>0)
      {
        edge.id = c(edge.id, max(edge.id) + rep(1:nrow(ref.pairs), 2))        
        ref.pairs = rbind(ref.pairs, cbind(-ref.pairs[,2], -ref.pairs[,1])) # reverse ref pairs        
        adj.ref = Matrix(0, nrow = 2*length(tile), ncol = 2*length(tile),
          dimnames = rep(list(as.character(c(1:length(tile), -(1:length(tile))))), 2))
        adj.ref[cbind(match(as.character(ref.pairs[,1]), rownames(adj.ref)),
                      match(as.character(ref.pairs[,2]), colnames(adj.ref)))] = nrow(ab.pairs)+1:nrow(ref.pairs)
      }
    else
      {
        adj.ref = Matrix(FALSE, nrow = 2*length(tile), ncol = 2*length(tile),
          dimnames = rep(list(as.character(c(1:length(tile), -(1:length(tile))))), 2))
      }

    ## current tile is partition of genome only in positive orientation + dummy intervals for breakpoints
    ## output tile is forward partition and followed by reverse partition
    ## (this is what is currently referenced by adj.ref and adj.ab)
    ## TODO: clean up this part
    tmp.nm = as.character(c(1:length(tile), -(1:length(tile))))
    tile = c(tile, gr.flipstrand(tile))
    names(tile) = tmp.nm

    ## apply ix to adj.ref and adj.ab, and create "adj" which has union of reference and aberrant junctions
    ## and adj.source which remembers whether edge ij was reference (value = 1) or aberrant (value = 2)
    adj.source = sign(adj.ref)+2*sign(adj.ab)
    adj = sign(adj.ref)+sign(adj.ab)
    tryres <- try( edges <- which(adj!=0, arr.ind=T), silent=T ) ## num edge x 2 matrix of vertex pairs  
    if ( class( tryres ) == "try-error" ) {
          relib( "Matrix" ); relib( "gUtils" )
          edges <- which(adj!=0, arr.ind=T)
    }
    adj[edges] = 1:nrow(edges) ## re number edges across edge set 
    rownames(adj) = colnames(adj) = 1:nrow(adj)
    G = graph.adjacency(adj ,weighted = 'edge.ix') ## edge.ix will allow us to match up edges in the adj matrix with edges in the igraph
    node.ind = abs(as.numeric(V(G)$name))

    ## add vertex features including formatting to igraph
    V(G)$chrom = as.character(seqnames(tile))[node.ind]
    V(G)$start = start(tile)[node.ind]
    V(G)$end = end(tile)[node.ind]
    V(G)$width = width(tile)[node.ind]
    V(G)$strand = sign(as.numeric(V(G)$name))
    V(G)$size = 5;
    V(G)$shape= c('rectangle', 'crectangle')[1 + as.numeric(V(G)$strand<0)];
    V(G)$border.width = c(1, 2)[1 + as.numeric(V(G)$strand=='-')] ;
    V(G)$label = paste(V(G)$chrom, ':', round(V(G)$start/1e6,0), '-', round(V(G)$end/1e6,0), sep = '')
    V(G)$label[V(G)$strand<0] = paste(V(G)$chrom, ':', round(V(G)$end/1e6,0), '-', round(V(G)$start/1e6,0), sep = '')[V(G)$strand<0]
    col.map = structure(brewer.master(length(seqlevels(tile))), names = seqlevels(tile))
    V(G)$chrom.ord = levapply(as.numeric(V(G)$start), list(V(G)$chrom), 'rank')
    V(G)$y = V(G)$chrom.ord*30
    V(G)$x = chr2num(V(G)$chrom)*300 + 100*rep(c(0,1), each = length(tile)/2)
    V(G)$col = col.map[V(G)$chrom]

    ## add edge features including formatting to igraph
    E(G)$weight = 1
    E(G)$from = edges[E(G)$edge.ix, 1]
    E(G)$to = edges[E(G)$edge.ix, 2]
    E(G)$col = c(col2hex('gray20'), col2hex('red'))[adj.source[edges[E(G)$edge.ix, ]]]
    E(G)$type = c('reference', 'aberrant', 'aberrant')[adj.source[edges[E(G)$edge.ix, ]]]
    E(G)$line.style = 'SEPARATE_ARROW'
    E(G)$arrow.shape = 'ARROW'
    E(G)$width = 1
    ab.ix = E(G)$type=='aberrant'  ## keep track of bp.id leading to edge
    E(G)$bp.id = NA;
    if (length(ab.pairs.bpid)>0)          
      E(G)$bp.id[ab.ix] = ab.pairs.bpid[adj.ab[cbind(E(G)$from[ab.ix], E(G)$to[ab.ix])]]
    E(G)$eid = NA; ## what is edge ID??? how is different from edge.ix?
    E(G)$eid[ab.ix] = edge.id[adj.ab[cbind(E(G)$from[ab.ix], E(G)$to[ab.ix])]] 
    E(G)$eid[!ab.ix] = edge.id[adj.ref[cbind(E(G)$from[!ab.ix], E(G)$to[!ab.ix])]]
    values(tile) = values(tile)[, c('tile.id', 'is.tel')]
    tile$ab.source = 1:length(tile) %in% E(G)$from[ab.ix]
    tile$ab.target = 1:length(tile) %in% E(G)$to[ab.ix]

    # important: map input ra to aberrant graph edges, i.e. ab.edges matrix with $from $to and $edge.ix columns
    # and one row for each aberrant edge    
    ab.edges = array(NA, dim = c(length(junctions), 3, 2), dimnames = list(NULL, c('from', 'to', 'edge.ix'), c('+', '-')))
    dupped = duplicated(ab.pairs.bpid)
    ab.edges[,1:2,1] = cbind(match(ab.pairs[!dupped,1], names(tile)), match(ab.pairs[!dupped,2], names(tile)))
    ab.edges[,1:2,2] = cbind(match(ab.pairs[dupped,1], names(tile)), match(ab.pairs[dupped,2], names(tile)))
    ab.edges[,3, 1] = match(paste(ab.edges[,1,1], '|', ab.edges[,2,1]), paste(E(G)$from, '|', E(G)$to)) ## must be easier way to perform this taks
    ab.edges[,3, 2] = match(paste(ab.edges[,1,1], '|', ab.edges[,2,1]), paste(E(G)$from, '|', E(G)$to))

    if (label.edges & nrow(ab.edges)>0)
        {
            ix = c(ab.edges[,1,1], ab.edges[,2,1], ab.edges[,1,2], ab.edges[,2,2])
            tile$edges.out = tile$edges.in = ''            
            tile$edges.in[ix]= sapply(ix,
                function(x) {ix = which(adj[,x]!=0); paste(ix, '->', sep = '', collapse = ',')})
            tile$edges.out[ix] = sapply(ix,
                function(x) {ix = which(adj[x, ]!=0); paste('->', ix,  sep = '', collapse = ',')})
        }
    
    return(list(tile = tile, adj = adj, G = G, ab.adj = adj.ab != 0, ab.edges = ab.edges, junctions = junctions))
  }


#' @name karyotrack
#' @title karyo track
#' @details
#' 
#' Takes karyograph and outputs gTrack +/- highlighting of one or more paths defined as GRanges or GRangesList (for multiple paths)
#' Edges will only be highlighted when the exact interval pair corresponding to the edge is included in
#' the graph
#'
#' @param kag  output of karyograph
#' @param paths GRanges or GRangesList
#' @return gTrack of karyograph with particular nodes / edges colored with specified colors
#' @export
############################################
karyotrack = function(kag, paths = NULL, col = 'red', pad = 0)
    {
        if (length(paths)==0)
            paths = NULL
                
        edge.ix = which(kag$adj!=0, arr.ind = T) ## collect aberrant edges
                
        ## convert to "simplified form"
        edges = data.frame(from = edge.ix[,1], to = edge.ix[,2])

        estr = paste(edge.ix[,1], edge.ix[,2])
        abestr = paste(kag$ab.edges[,1,1:2], kag$ab.edges[,2,1:2])

        if (nrow(edges)>0)
            {
                edges$type = 'reference'
                if (any(ix <- estr %in% abestr))
                    edges$type[ix] = 'aberrant'
                
                edges$col = ifelse(edges$type == 'aberrant', alpha('purple', 0.4), alpha('gray10', 0.4))                
                edges$h = 1                
                edges$lwd = ifelse(edges$type == 'aberrant', 2, 1)
                edges$lty = 1
                edges$cex.arrow = 0
                edges$v = 1
                edges$not.flat = edges$type == 'aberrant'
                edges$v[edges$type == 'aberrant'] = 2
                edges$h[edges$type == 'aberrant'] = 2
                edges$dangle.w = 0.5
            }       
        
        pos.ix = which( as.logical( strand(kag$tile)=='+') )
        kag$tile$tile.id = match(gr.stripstrand(kag$tile), gr.stripstrand(kag$tile[pos.ix]))
        ss = kag$tile
        ss$col = alpha('gray', 0.2)
        ss$border = alpha('black', 0.5)
        ss$ywid = 0.8
        ss$bin = suppressWarnings(disjointBins(ss+1+pad, ignore.strand = FALSE))

        out = gTrack()
        if (!is.null(paths))
            {                
                if (is(paths, 'GRanges'))
                    paths = split(paths, 1)
                
                if (is.null(values(paths)$col))
                    values(paths)$col = col

                paths.u = grl.unlist(paths)[, c('col', 'grl.ix')] %**% ss[, c()]

                ss$col[paths.u$subject.id] = paths.u$col

                edges = as.data.table(edges)

                edges[, end.from := ifelse(as.logical(strand(ss)[from]=='+'), end(ss)[from], start(ss)[from])]
                edges[, start.to := ifelse(as.logical(strand(ss)[to]=='+'), start(ss)[to], end(ss)[to])]
                
                paths.u = gr2dt(paths.u)[, ix := 1:length(seqnames)]
                paths.up = paths.u[ , list(ix.from = (ix)[-(length(ix))], ix.to = (ix)[-1], color = col[1]), by = grl.ix]
                paths.up[, from := paths.u$subject.id[ix.from]][, to := paths.u$subject.id[ix.to]]
                paths.up[, strand.from := paths.u$strand[ix.from]][, strand.to := paths.u$strand[ix.to]]
                paths.up[ , end.from := ifelse(strand.from == '+', paths.u$end[ix.from], paths.u$start[ix.from])]
                paths.up[ , start.to := ifelse(strand.to == '+', paths.u$start[ix.to], paths.u$end[ix.to])]
                setkeyv(edges, c('from', 'to', 'end.from', 'start.to'))
                setkeyv(paths.up, c('from', 'to', 'end.from', 'start.to'))

                edges[, col := alpha(col, 0.2)];
                edges[, path.match := 0]
                                        #      edges[paths.up, col := color] ## update colors if there is a match on from to, end and start
                edges[paths.up, lwd := lwd*3] ## update colors if there is a match on from to, end and start
                edges[paths.up, col := alpha(col, 0.8)] ## update colors if there is a match on from to, end and start
                edges[paths.up, path.match := 1] ## update colors if there is a match on from to, end and start
                edges = edges[order(path.match), ]
                out = gTrack(paths, draw.paths = TRUE, name = 'Paths')
            }                            

        out = c(gTrack(ss, y.field = 'bin', edges = edges, name ='Karyograph', angle = 0, labels.suppress = TRUE, yaxis = FALSE), out)
                                                                                    
        return(out)
    }


##############################
#' karyoseg
#' 
#' performs preliminary segmentation on preliminary karyograph using CBS on coverage in tiles of coverage samples on the genome
#' only segments inter-junction reference intervals above min.width (which is usually defined as 20 times the coverage width)
#' returns the new karyograph with the addition partitions (and reference interval pairs)
#' 
###############################
karyoseg = function(kag, cov)
  {
  }


##############################
#' segstats is a step in the JaBbA pipeline
#'
#' @details
#' computes posterior mean's and sd's for a target tiling GRanges of segments (target)
#' target must be a non-overlapping gapless strandless or two stranded tiling of the genome (eg output of gr.tile)
#' if two stranded, then every stranded interval must have a mirror image interval included 
#'
#' given a GRanges of signals using value field "field" of signal GRanges
#' assuming that the signals inside each interval in "target" are independent samples from a
#' normal distribution of unknown mean and variance
#' 
#' will also compute means ands std deviations on a gamma posterior of the
#' the "high" and "low" alleles for a granges "asignal" representing allelic signal ref count and tot count
#' across a set of locations (i.e. modeling the counts as a poisson random variable)
#' Fields are specified by afields
#' 
#' outputs target GRanges with "$mean" and "$sd" fields populated
#'
#' @param target GRanges of segments on which segstats will be computed
#' @param signal GRanges of coverage from which samples will be taken
#' @param field  field of "signal" GRanges from which coverage signal will be pulled
#' @param asignal optional GRanges corresponding width 1 bialellic snp allele counts across the genome,
#' @param afields length 2 character vector meta data fields of asignal GRanges that will be used to get allele counts (default is ref.count, alt.count)
#' @param subsample number between 0 and 1 with which to subsample per segment for coverage (useful for superdense coverage eg 50 bases to avoid correlations between samples due to read overlap)
#' @param mc.cores number of cores to run on (default 1)
#' @export
############################################
segstats = function(target,
  signal = NULL, 
  field = 'signal',
  asignal = NULL, ## granges corresponding to width 1 snp allele counts across the genome
  afields = c('ref.count', 'alt.count'), ## length 2 character vector specifying the two allele fields of asignal
  prior_weight = 1,
  prior_mean = NA, # if NA will compute prior "empirically"
  prior_alpha = NA, # priors for inverse gamma for variance inference
  prior_beta = NA,
  max.chunk = 1e8,
  max.slice = 2e4,
  na.thresh = 0.2,    
  subsample = NULL, ## number between 0 and 1 to subsample per segment for coverage (useful for dense coverage)
  mc.cores = 1, 
  nsamp_prior = 1e3, ## number of data samples to estimate alpha / beta prior value
  ksamp_prior = 100  ## size of data samples to estimate alpha / beta prior values
  )
{
    if (!is.null(asignal))
      {
        if (all(afields %in% names(values(asignal))))
          {
            asignal$low.count = pmin(values(asignal)[, afields[1]], values(asignal)[, afields[2]])
            asignal$high.count = pmax(values(asignal)[, afields[1]], values(asignal)[, afields[2]])
            asignal$ix = gr.match(asignal, target, max.slice = max.slice)

            aprior_alpha = 1 #mean(c(asignal$low.count, signal$high.count), na.rm = T)
            aprior_beta = 1
            
            .postalpha = function(x)
              aprior_alpha + sum(x, na.rm = T)

            .postbeta = function(x)
              aprior_beta + sum(!is.na(x))
            
            asignal.df = as.data.frame(asignal)
            alpha_high = vaggregate(high.count ~ ix, asignal.df, .postalpha)[as.character(1:length(target))]
            beta_high = vaggregate(high.count ~ ix, asignal.df, .postbeta)[as.character(1:length(target))]
            alpha_low = vaggregate(low.count ~ ix, asignal.df, .postalpha)[as.character(1:length(target))]
            beta_low = vaggregate(low.count ~ ix, asignal.df, .postbeta)[as.character(1:length(target))]
            
            target$mean_high = alpha_high / beta_high
            target$sd_high = sqrt(alpha_high / (beta_high)^2)
            target$mean_low = alpha_low / beta_low
            target$sd_low = sqrt(alpha_low / (beta_low)^2)                   
          }
        else
          stop('One or more of the afields ', paste(afields, collapse = ', '), ' not found as meta data columns of asignal')
      }

    if (!is.null(signal))
      {

        if (!(field %in% names(values(signal))))
          stop('Field not found in signal GRanges')

        utarget = unique(gr.stripstrand(target))
                
        map = gr.tile.map(utarget, signal, verbose = T, mc.cores = mc.cores)               
        val = values(signal)[, field]
        val[is.infinite(val)] = NA
        vall = lapply(map, function(x) val[x])

        ## do subsampling before we dedup
##         if (!is.null(subsample))
##           {
##             subsample = pmax(0, pmin(1, subsample))
##             vall = lapply(vall, function(x) if (length(x)==0) NA else if (all(is.na(x))) NA else sample(x, ceiling(length(x)*subsample)))
##           }
        
        vall = vall[match(gr.stripstrand(target), utarget)]       
        
        if (is.na(prior_mean))
          {            
            #          prior_mean = mean(vaggregate(signal ~ query.id, data = as.data.frame(r), mean))
            #            tmp = vaggregate(signal ~ query.id, data = as.data.frame(r), mean)
            tmp = sapply(vall, mean, na.rm = T)
            tmp.ix = !is.na(tmp)
            prior_mean = sum(as.numeric(width(target)*tmp)[tmp.ix])/sum(as.numeric(width(target))[tmp.ix])
          }
        
        if (is.na(prior_alpha) | is.na(prior_beta))
          {
            ## guessing alpha and beta priors for inverse gamma (for variance estimation)
            ## we want to bias this estimate by looking at a bunch of local "variance" samples
            ## from signal .. i.e. a hundred adjacent markers
#            tmp.sig = r$signal[!is.na(r$signal)]
            tmp.sig = val[!is.na(val)]

            if (length(tmp.sig)<ksamp_prior)
              stop('Something wrong, number of samples smaller than data sample size')
            
            vars = sapply(1:nsamp_prior, function(x) var(sample(tmp.sig, ksamp_prior)))            
            E_var = mean(vars, na.rm = T)
            var_var = var(vars, na.rm = T)
            prior_alpha = E_var^2/var_var + 2
            prior_beta = E_var * (E_var^2/var_var + 1)
          }
        
        .postmean = function(x)
          {
            x_bar = mean(x, na.rm = T)
            n = sum(!is.na(x))
            if (is.na(x_bar))
              x_bar = 0
            return((prior_mean*prior_weight + n*x_bar)/(n + prior_weight))
          }
        
        .postsd = function(x, subsample = NULL)
          {
            ## TODO: fix prior weight here .. 
            var_bar = var(x, na.rm = T)
            n = sum(!is.na(x))
            if (!is.null(subsample)) ## apply penalty to simulate subsampling
              {
#                print('subsampling')
                n = ceiling(n*subsample)
              }
            post_alpha = prior_alpha + n/2 - 1/2
            post_beta = prior_beta + 1/2*var_bar*n
              
            sigma2_mode = post_beta/(post_alpha + 1) ## only using modal value of igamma posterior
            if (is.na(post_beta))
              sigma2_mode = prior_beta/(prior_alpha+1)
            return(sqrt(sigma2_mode/(n + prior_weight)))
          }
        
        mu = sapply(vall, .postmean)
        sd = sapply(vall, .postsd, subsample)
        pc.na = sapply(vall, function(x) sum(is.na(x))/length(x))
       
        mu[which(pc.na>na.thresh)] = NA
        sd[which(pc.na>na.thresh)] = NA
        
        target$mean = prior_mean;
        target$sd = sqrt(prior_beta / (prior_alpha + 1));

        ix = !is.na(mu)
        target$mean[ix] = mu[ix]
        target$sd[ix] = sd[ix]
    }

    return(target) 
}

####################################################################
#' jbaMIP
#'
#' @details
#' primary "heavy" lifting task of JaBbA.  Sets up optimization problem given an input graph and segstats input
#' and sends to CPLEX via RCplex
#'
#' combines edge-conservation constraints from karyograph (n x n adjacency matrix connecting n genomic intervals) with segment abundance data
#' (segstats - length n GRanges object with mean and sd fields corresponding to posterior means and sd's on the
#' the relative "concentration" of each interval) to infer
#' 
#' (1) interval and edge absolute copy numbers on the karyograph
#' (2) purity and ploidy
#' (3) slack edges (if any needed)
#' 
#' basically solves ABSOLUTE problem (fitting integer grid to continuous segment intensities)
#' while enforcing edge-conservation constraints. 
#' 
#' Most important optional parameters include
#' (1) cn.sd (expected deviation of absolute copy number from ploidy)
#' (2) ploidy.min and ploidy.max --> useful for probing alternate solutions, but can be generously set
#' (3) adj.lb --> enforces minimal edge absolute copy number, eg to force aberrant adjacency use
#' (4) edge.slack - logical variable to determine whether or not to allow penalized relaxation of edge conservation constraints
#' (5) nsolutions - number of alternate solutions
#'
#' @param adj n x n adjacency matrix interpreted as binary (this is the $adj output of karyograph)
#' @param segstats n x 1 GRanges object with "mean" and "sd" value fields
#' @param beta numeric guess for beta (i.e. from ppgrid)
#' @param gamma numeric guess for gamma (i.e. from ppgrid)
#' @param slack.prior 1/slack.prior = penalty for each additional copy number of each slack edge, the higher slack.prior the more slack we allow in the reconstruction, should be intuitively calibrated to the expected "incompleteness" of the reconstruction, 1/slack.prior should be calibrated with respect to 1/(k*sd)^2 for each segment, so that we are comfortable with junction balance constraints introducing k copy number deviation from a segments MLE copy number assignment (the assignment in the absence of junction balance constraints)
#' @param field.ncn this field takes into account normal copy number in relative to absolute conversion
#' @param adj.lb nxn matrix of lower bounds on particular copy numbers - this is used to force certain junctions into the graph
#' @param adj.nudge nxn adjacency matrix of "nudge" rewards on individual junctions, NOTE: maximum value in this matrix
#' 
#' @return
#' output is a Rcplex solution or list of Rcplex solution with additional fields, each Rcplex solution is a list and the additional fields
#' added by jbaMIP are
#' $adj input n x n adjacency matrix populated with integer copy numbers
#' $segstats input segstats vector populated with meta data fields $cn, $ecn.in, $ecn.out, $edges.out, $eslack.in, $eslack.out
#' $purity purity value associated with relative-absolute affine copy number conversion for this solution
#' $ploidy purity value associated with relative-absolute affine copy number conversion for this solution
#' $
#' 
#' 
#' Additional fields for qc / technical debugging:
#' $nll.cn negative log likelihood corresponding to the CN fit in this solution
#' $nll.opt negative log likelihood correpsonding to the MLE CN without junction constraints
#' $residual = value of residual between copy solution and MLE fit without junction constraints
#' $beta beta value associated with relative-absolute affine copy number conversion for this solution
#' $gamma gamma value associated with relative-absolute affine copy number conversion for this solution
#' $gap.cn total gap between MLE fit without junction constraints and JaBbA fit
#' $ploidy.constraints input ploidy constraints
#' $beta.constraints input beta constrinats
#' $cn.prior input cn.prior
#' $slack.prior input slack.prior
#' @export
############################################
jbaMIP = function(
  adj, # binary n x n adjacency matrix ($adj output of karyograph)
  segstats, # n x 1 GRanges object with "mean" and "sd" value fields

  ########### optional args
  beta = NA, # beta guess (optional)  
  gamma = NA, # gamma guess (optional)
  field.ncn = 'ncn', # will use this field to take into account normal copy number in transformation of relative to integer copy number
  tilim = 20, mipemphasis = 0, epgap = 0.01, # MIP params
  ploidy.min = 0.1, # ploidy bounds (can be generous)
  ploidy.max = 20,
#  purity.guess = NA,
#  ploidy.guess = NA,
  beta.guess = beta, 
  beta.min = beta.guess,
  beta.max = beta.guess,
  gamma.guess = NA, 
  gamma.min = gamma.guess,
  gamma.max = gamma.guess,
  cn.sd = 1, # sd of cn prior (mean is ploidy) - i.e. sd of local copy number from ploidy
  cn.prior = cn.sd,
  partition = T, ## whether to partition the problem into MIP subproblems depending on the relationships of the segment standard deviation and the value of the slack.prior (only works if gamma.guess, beta.guess are specified)
  purity.prior.mean = NA,
  purity.prior.sd = 0.3,
  purity.prior.strength = 1,
  purity.prior = c(purity.prior.mean, purity.prior.sd),
  cn.fix = rep(NA, length(segstats)), ## vector of NA's and (integer) values to which to "fix" copy states, only non NA's are incorporated
  cn.lb = cn.fix,  
  cn.ub = cn.fix,
  loose.ends = c(), ## integer vector specifies indices of "loose ends", slack won't be penalized at these vertices
  adj.lb = 0*adj, # lower bounds for adjacency matrix
  adj.nudge = 0*adj, # linear objective function coefficients for edges (only which(adj!=0) components considered)
  na.node.nudge = TRUE, 
  ecn.out.ub = rep(NA, length(segstats)), ## upper bound for cn of edges leaving nodes
  ecn.in.ub = rep(NA, length(segstats)),  ## upper bound for cn of edges entering nodes  
  gurobi = F, # otherwise will use cplex
  nsolutions = 1,
  verbose = F,
  debug = F,
  mc.cores = 1, ## only matters if partition = T
  ignore.edge = FALSE, ignore.cons = TRUE, edge.slack = TRUE, slack.prior = 1, 
  ... # passed to optimizer
  )
{
  require(Rcplex)
  if (length(segstats) != nrow(adj))
    stop('length(segstats) !=  nrow(adj)')

  if (is.null(adj.lb))
      adj.lb = 0*adj
  
#  message('Gunes!!! We are enforcing ', sum(adj.lb!=0), ' lower bound constraints!!!!!')
  
  #####
  # wrapper that calls jbaMIP recursively on subgraphs after "fixing"
  #####
  if (partition & !is.na(gamma.guess) & !is.na(beta.guess))
    {
      require(igraph)
      
#      m = segstats$mean*beta.guess - gamma.guess 
      m = rel2abs(segstats, gamma = gamma.guess, beta = beta.guess, field = 'mean', field.ncn = field.ncn)
      cnmle = round(m) ## MLE estimate for CN

      residual.min = ((m-cnmle)/(segstats$sd))^2
      residual.other = apply(cbind((m-cnmle-1)/segstats$sd, (m-cnmle+1)/segstats$sd)^2, 1, min)
      residual.diff = residual.other - residual.min ## penalty for moving to closest adjacent copy state
      
      ## we fix nodes for which the penalty for moving to non (locally) optimal copy state
      ## is greater than k / slack.prior penalty (where k is some copy difference
      ## that we would never imagine a "reasonable" slack to have to over-rule      
      fix = as.integer(which(residual.diff>(8/slack.prior))) ## 8 is a constant that is conservative, but basically assumes that no node will have more than 4 neighbors (todo: make adjustable per node)
      
      if (verbose)
        cat('Fixing', length(fix), 'nodes that are unmovable by slack\n')

      ##
      ## now we will create a graph of unfixed nodes and fixed node "halves"
      ## i.e. we split each fixed node to a node that is receiving edges
      ## and a node that is sending edges
      ## 
      unfix = as.numeric(setdiff(1:length(segstats), fix))
      G = graph(as.numeric(t(which(adj!=0, arr.ind = T))), n = length(segstats), directed = T)
      V(G)$name = as.numeric(V(G)) ##  1:length(V(G)) ## igraph vertex naming is a mystery

      if (length(fix)>0)
          G.unfix = induced.subgraph(G, unfix) + vertices(c(paste('from', fix), paste('to', fix)))
      else
          G.unfix = induced.subgraph(G, unfix)

      if (length(fix)>0 & length(unfix)>0)
        node.map = structure(c(unfix, fix, fix), names = c(as.character(unfix), paste('from', fix), paste('to', fix)))
      else if (length(fix)>0)
        node.map = structure(c(fix, fix), names = c(paste('from', fix), paste('to', fix)))
      else
        node.map = structure(c(unfix), names = c(as.character(unfix)))
      
      ## add nodes representing the "receiving" and "sending" side of fixed nodes
      if (length(fix)>0 & length(unfix)>0)
        {
          tofix = which(adj[unfix, fix]!=0, arr.ind = T)
          fromfix = which(adj[fix, unfix]!=0, arr.ind = T)
        }
      else
        {
          tofix = c()
          fromfix = c()            
        }

      if (length(fix)>0)
        fixtofix = which(adj[fix, fix]!=0, arr.ind = T)
      else
        fixtofix = c()
      
      if (length(tofix)>0)
        e.tofix = edges(as.vector(rbind(unfix[tofix[,1]], paste('to', fix[tofix[,2]]))))
      else
        e.tofix = edges()

      if (length(fromfix)>0)        
        e.fromfix = edges(rbind(paste('from', fix[fromfix[,1]]), unfix[fromfix[,2]]))
      else
        e.fromfix = edges()

      if (length(fixtofix)>0)
        e.fixtofix = edges(rbind(paste('from', fix[fixtofix[,1]]), paste('to', fix[fixtofix[,2]])))
      else
        e.fixtofix = edges()
      
      ## add edges to graph from fixed to unfixed, unfixed to fix, and fixed to fixed node sides
      G.unfix = G.unfix + e.tofix + e.fromfix + e.fixtofix      

      ## find connected components in these graphs
      cl = clusters(G.unfix, 'weak')  
      cll = split(V(G.unfix)$name, cl$membership) ## keep augmented graph names, use node.map later
      
      ## combine components with their reverse complement components
      ## (only intervals that have a (fixed node free) path from their positive to their negative strand
      ## will be part of the same component .. all other intervals will be separated from their
      ## reverse complement.  However, in the MIP we always optimize
      ## over both strands, and thus must merge components with their reverse complement
      pos.ix = which( as.logical( strand(segstats)=='+') )
      neg.ix = which( as.logical( strand(segstats)=='-') )

      ## maps segments and reverse complements
      seg.map = c(1:length(pos.ix), suppressWarnings(pos.ix[match(segstats[neg.ix], gr.flipstrand(segstats[pos.ix]))]))
      
      cll.m = sapply(cll, function(x) paste(sort(seg.map[node.map[x]]), collapse = ' '))
      dup.ix = match(cll.m, unique(cll.m))
#      cll = lapply(split(1:length(dup.ix), dup.ix), function(x) sort(unique(do.call('c', cll[x]))))
      cll = lapply(split(1:length(dup.ix), dup.ix), function(x) c(cll[[x[1]]], cll[[x[2]]]))
      
      ord.ix = order(-sapply(cll, length))
      cll = cll[ord.ix]
      
      if (verbose)
        cat('Partitioned graph into ', length(cll), ' connected components with the size of the highest 10 components being:\n',
            paste(sapply(cll[1:min(10, length(cll))], length), collapse = ','), '\n')
      
      cn.fix = ifelse(1:length(segstats) %in% fix, cnmle, NA)

      ## force "non lazy" evaluation of args in order to avoid weird R ghosts (WTF) downstream in do.call
      args = as.list(match.call())[-1]
      args = structure(lapply(names(args), function(x) eval(parse(text = x))), names = names(args))      

      sols = mclapply(1:length(cll), function(k, args)
        {
          if (debug)
            browser()
          
          ix = node.map[cll[[k]]] ## indices in the original graph
          uix = unique(ix)
          fr.ix = grepl('from', cll[[k]])
          to.ix = grepl('to', cll[[k]])

          ## we want to make sure that fixed nodes that straddle
          ## two clusters will only have the "correct"
          ## half included in this run          
          fronly.ix = setdiff(ix[fr.ix], ix[to.ix])
          toonly.ix = setdiff(ix[to.ix], ix[fr.ix])

          ## now we want to make sure that fronly.ix don't have incoming edges
          tmp.adj = adj[uix, uix, drop = F]

          if (length(fronly.ix)>0)
            tmp.adj[, as.character(as.integer(fronly.ix))] = 0

          ## and toonly.ix don't have outgoing edges
          if (length(toonly.ix)>0)
            tmp.adj[as.character(as.integer(toonly.ix)), ] = 0
          
          args$adj = tmp.adj
          args$adj.nudge = adj.nudge[uix, uix, drop = F]
          args$na.node.nudge = na.node.nudge
          args$adj.lb = adj.lb[uix, uix, drop = F]
          args$segstats = segstats[uix]
          args$cn.fix = cn.fix[uix]
          args$cn.lb = cn.lb[uix]
          args$cn.ub = cn.ub[uix]
          args$partition = F
          args$nsolutions = 1 
          args$ploidy.min = 0 ## no ploidy constraints          
          args$ploidy.max = max(c(100, cnmle[ix]), na.rm = T)*1.5          
          
          if (verbose)
            cat('Starting cluster ', k, 'of', length(cll), 'which has', length(uix), 'nodes comprising',
                sum(as.numeric(width(segstats[uix])))/2/1e6, 'MB and', length(unique(seqnames((segstats[uix])))),
                'chromosomes, including', paste(names(sort(-table(as.character(seqnames((segstats[uix])))))[1:min(4,
                     length(unique(seqnames((segstats[uix])))))]), collapse = ', '), '\n')

          
          out = do.call('jbaMIP', args)

          gc() ## garbage collect .. not sure why this needs to be done
          
          return(out)
        }, args, mc.cores = mc.cores)

      out = list()
      for (f in c('residual', 'nll.cn', 'nll.opt', 'gap.cn', 'cn.prior', 'slack.prior')) ## scalar fields --> length(cluster) vector
        out[[paste('component', f, sep = '')]] = sapply(sols, function(x) x[[f]])
      
      for (f in c('ploidy.constraints', 'beta.constraints')) ## length 2 fields --> length(cluster) x 2 matrix
        out[[paste('component', f, sep = '')]] = do.call('rbind', lapply(sols, function(x) x[[f]]))
            
      ## adjacency matrix
      out$adj = 0 * adj      
      for (i in 1:length(sols))
        {
          ix1 = as.numeric(rownames(sols[[i]]$adj))
          out$adj[ix1, ix1] = out$adj[ix1, ix1] + sols[[i]]$adj
        }

      ## segstats
      sol.ix = lapply(sols, function(x) as.numeric(rownames(x$adj)))
      out$segstats = do.call('grbind', lapply(sols, function(x) x$segstats))[match(1:length(segstats), unlist(sol.ix))]

      ## annotate segstats keep to keep track and "fixed nodes"
      out$segstats$fixed = 1:length(out$segstats) %in% fix
      out$segstats$cn.fix = cn.fix
      out$segstats$cl  = NA
      out$segstats$id = 1:length(out$segstats)
      
      ## keep track of which clusters segments originated
      sol.ixul = munlist(sol.ix)
      tmp = vaggregate(sol.ixul[,1], by = list(sol.ixul[,3]), FUN = paste, collapse = ',')
      out$segstats$cl = NA
      out$segstats$cl[as.numeric(names(tmp))] = tmp

      out$purity = 2/(2+gamma.guess)
      v = out$segstats$cn; w = as.numeric(width(out$segstats))
      out$ploidy = sum((v*w)[!is.na(v)]) / sum(w[!is.na(v)])
      out$beta = beta.guess;
      out$gamma = gamma.guess;
      
      target.less = Matrix::rowSums(adj, na.rm = T)==0
      source.less = Matrix::colSums(adj, na.rm = T)==0
      out$segstats$eslack.out[!target.less] = out$segstats$cn[!target.less] - Matrix::rowSums(out$adj)[!target.less]
      out$segstats$eslack.in[!source.less] =  out$segstats$cn[!source.less] - Matrix::colSums(out$adj)[!source.less]

      out$segstats$ecn.out =  Matrix::rowSums(out$adj)
      out$segstats$ecn.in =  Matrix::colSums(out$adj)

      out$segstats$edges.in = sapply(1:length(out$segstats),
        function(x) {ix = which(adj[,x]!=0); paste(ix, '(', out$adj[ix,x], ')', '->', sep = '', collapse = ',')})
      out$segstats$edges.out = sapply(1:length(out$segstats),
        function(x) {ix = which(adj[x, ]!=0); paste('->', ix, '(', out$adj[x,ix], ')', sep = '', collapse = ',')})

      ###
      ###
      ncn = rep(2, length(segstats))
      if (!is.null(field.ncn))
        if (field.ncn %in% names(values(segstats)))
          ncn = values(segstats)[, field.ncn]
            
      
      nnix = !is.na(out$segstats$mean) & !is.na(out$segstats$sd) & !is.na(out$segstats$cn)

      ##       out$obj = 1/4*sum(((out$segstats$cn[nnix] + out$gamma - out$beta*out$segstats$mean[nnix])/out$segstats$sd[nnix])^2) + 
#        1/slack.prior * (sum(out$segstats$eslack.in + out$segstats$eslack.out, na.rm = T)) ## 1/4 because our original objective is 1/2 for pos strand intervals only      
      ##

      ### new obj allowing variable normal copy number
      out$obj = 1/4*sum(((out$segstats$cn[nnix] + ncn[nnix]/2*out$gamma - out$beta*out$segstats$mean[nnix])/out$segstats$sd[nnix])^2) + 
        1/slack.prior * (sum(out$segstats$eslack.in + out$segstats$eslack.out, na.rm = T)) ## 1/4 because our original objective is 1/2 for pos strand intervals only      
      out$nll.cn = (1/2*sum(((out$segstats$cn[nnix] + out$gamma - out$beta*out$segstats$mean[nnix])/out$segstats$sd[nnix])^2))
      out$nll.opt = (1/2*sum(((cnmle[nnix] + out$gamma - out$beta*out$segstats$mean[nnix])/out$segstats$sd[nnix])^2))
      out$gap.cn = as.numeric(1 - out$nll.opt / out$nll.cn)
      out$sols = sols
      
      return(out)
    }
        
  # map intervals to their reverse complement to couple their copy number (and edge variables)  
  pos.ix = which( as.logical( strand(segstats)=='+') )
  neg.ix = which( as.logical( strand(segstats)=='-') )

  ## "original vertices"
  og.ix = pos.ix 

  ## "duplicates" of og.ix i.e. revcomp vertices
  dup.ix = suppressWarnings(neg.ix[match(segstats[og.ix], gr.flipstrand(segstats[neg.ix]))])

  if (!identical(segstats$mean[og.ix] , segstats$mean[dup.ix]) & !identical(segstats$sd[og.ix] , segstats$sd[dup.ix]))
    stop('Segstats mean or sd not identical for all pos / neg strand interval pairs: check segstats computation')
  
  edges = which(adj!=0, arr.ind = T)

  if (verbose)
    cat('Setting up matrices .. \n')

  varmeta = data.frame() ## store meta data about variables to keep track
  consmeta = data.frame() ## store meta data about constraints to keep track
  n = 2*nrow(adj) + nrow(edges) + 2;  # number of vertices + slack variables, number of edges, and beta, gamma parameters.
  v.ix = 1:nrow(adj)
  s.ix = length(v.ix) + v.ix

  varmeta = data.frame(id = v.ix, subid = 1:length(v.ix), label = paste('interval', 1:length(v.ix), sep = ''), type = 'interval',
    stringsAsFactors = F)
  varmeta = rbind(varmeta, data.frame(id = s.ix, subid = 1:length(s.ix), label = paste('residual', 1:length(s.ix), sep = ''),
    type = 'residual', stringsAsFactors = F))
  
#  s.ix = length(v.ix) + v.ix + 1
  if (nrow(edges)>0)
    {
      e.ix = max(s.ix) + (1:nrow(edges))
      varmeta = rbind(varmeta, data.frame(id = e.ix, subid = 1:length(e.ix), label = paste('edge', 1:length(e.ix), sep = ''),
        type = 'edge', stringsAsFactors = F))
    }
  else
    e.ix = integer();
  
  gamma.ix = max(c(s.ix, e.ix))+1;
  beta.ix = max(c(s.ix, e.ix))+2;
  varmeta = rbind(varmeta, data.frame(id = c(gamma.ix, beta.ix), subid = rep(1, 2), label = c('gamma', 'beta'), type = 'global', stringsAsFactors = F))
  
  ## add cn prior variables
  if (!is.na(cn.prior))
    {
      d.ix = v.ix + n;
      varmeta = rbind(varmeta, data.frame(id = d.ix, subid = 1:length(d.ix), label = paste('cn.prior', 1:length(d.ix), sep = ''), type = 'cn.prior', stringsAsFactors = F))      
      ploidy.ix = max(d.ix)+1;
      n = length(v.ix)+n+1;
      varmeta = rbind(varmeta, data.frame(id = ploidy.ix, subid = 1, label = c('ploidy.prior'), type = 'global', stringsAsFactors = F))     
    }
  
  if (!any(is.na(purity.prior)))
    {
      pd.ix = n + 1; ## measure our deviation from the purity prior target
      n = n + 1
      varmeta = rbind(varmeta, data.frame(id = pd.ix, subid = 1, label = c('purity.prior'), type = 'global', stringsAsFactors = F))
    }

  if (edge.slack) # slack on edge consistency constraints
    {
      es.s.ix = n+(1:length(v.ix)) ## "source slack" variable
      varmeta = rbind(varmeta, data.frame(id = es.s.ix, subid = 1:length(es.s.ix), label = paste('source.slack', 1:length(es.s.ix), sep = ''), type = 'source.slack', stringsAsFactors = F))      
      es.t.ix = n+(length(v.ix) + 1:length(v.ix))
      varmeta = rbind(varmeta, data.frame(id = es.t.ix, subid = 1:length(es.t.ix), label = paste('target.slack', 1:length(es.t.ix), sep = ''), type = 'target.slack', stringsAsFactors = F))
      n = n+2*length(v.ix);
    }
  
  vtype = rep('C', n); vtype[c(v.ix, e.ix)] = 'I'
  lb = rep(0, n); lb[s.ix] = -Inf;

  if (nrow(edges)>0)
    lb[e.ix] = adj.lb[edges]

  if (any(ix <<- !is.na(cn.fix)))
    cat('Fixing copy states on', sum(ix), 'vertices\n')
  
### implement lower bounds and fixes
  if (any(!is.na(cn.lb)))
    lb[v.ix[!is.na(cn.lb)]] = cn.lb[!is.na(cn.lb)]

  ub = rep(Inf, n);

  if (any(!is.na(cn.ub)))
    ub[v.ix[!is.na(cn.ub)]] = cn.ub[!is.na(cn.ub)]

  ## add these vars to varmeta TODO: convert everything to data frame
  varmeta$vtype = vtype
  varmeta$lb = vtype
  varmeta$ub = ub
  
##   if (!is.na(purity.guess))
##     {
##       cat('applying purity guess .. ', purity.guess, ' \n')
##       gamma.guess = 2 * (1-purity.guess) / purity.guess
##       ub[gamma.ix] = lb[gamma.ix] = gamma.guess
##     }

##   if (!is.na(ploidy.guess))
##     {
##       cat('applying ploidy guess .. ', ploidy.guess, ' \n')

##       ## if we have ploidy and purity, then let's solve the problem, i.e. compute beta
##       ## (unless we have already, i.e. doing local analysis)
##       if (!is.na(purity.guess) & is.na(beta.guess)) 
##         {
##           mu = segstats$mean
##           w = as.numeric(width(segstats))
##           nna = !is.na(mu)
##           beta.guess = (2*(1-purity.guess) + purity.guess * ploidy.guess)/((purity.guess * sum(w[nna] * mu[nna]))/sum(w[nna]))
##         }
##       else ## otherwise we can only explicitly constrain ploidy, and implicitly constrain beta (beware on subgraphs)
##         ploidy.min = ploidy.max = ploidy.guess

##       ploidy.min = ploidy.max = ploidy.guess

##       ## override ignore.cons
## #      ignore.cons = T

## #      print(ignore.cons)
##     }
  
  if (!is.na(gamma.guess) & !is.na(beta.guess))
    cn.prior = cn.sd = NA    
  
  if (edge.slack)
    vtype[c(es.s.ix, es.t.ix)] = 'I'
  
  Zero = sparseMatrix(1, 1, x = 0, dims = c(n, n))
  
  # vertices that will actually have constraints (i.e. those that have non NA segstats )
  v.ix.c = setdiff(v.ix[!is.na(segstats$mean) & !is.na(segstats$sd)], dup.ix)

  v.ix.na = which(is.na(segstats$mean) | is.na(segstats$sd))
  
  # weighted mean across vertices contributing to mean
  if (length(v.ix.c))
    mu.all = (width(segstats)[v.ix.c] %*% segstats$mean[v.ix.c]) / sum(as.numeric(width(segstats)[v.ix.c]))
  else
    mu.all = NA

  if (verbose)
    cat('cn constraints .. \n')

  if (length(v.ix.c)>0)
    {
      ## take into account (variable) normal cn 
      ncn = rep(2, length(segstats)) 
      if (!is.null(field.ncn))
        if (field.ncn %in% names(values(segstats)))
          ncn = values(segstats)[, field.ncn]

      normal_ploidy = sum(width(segstats)[v.ix.c]*ncn[v.ix.c]) / sum(as.numeric(width(segstats))[v.ix.c])
      
      ## copy number constraints  
      Acn = Zero[rep(1, length(v.ix.c)+1), ]
      Acn[cbind(1:length(v.ix.c), v.ix.c)] = 1;
      Acn[cbind(1:length(v.ix.c), s.ix[v.ix.c])] = 1
##      Acn[cbind(1:length(v.ix.c), gamma.ix)] = 1  ## replacing with below
      Acn[cbind(1:length(v.ix.c), gamma.ix)] = ncn[v.ix.c]/2 ## taking into account (normal) variable cn
      Acn[cbind(1:length(v.ix.c), beta.ix)] = -segstats$mean[v.ix.c]

      ## final "conservation" constraint
      Acn[length(v.ix.c)+1, v.ix] = width(segstats)/sum(as.numeric(width(segstats)));
                                        #  Acn[length(v.ix.c)+1, s.ix[length(s.ix)]] = 1
##      Acn[length(v.ix.c)+1, gamma.ix] = 1; ## replacing with below
      Acn[length(v.ix.c)+1, gamma.ix] = normal_ploidy/2; ## taking into account (normal) variable cn
      Acn[length(v.ix.c)+1, beta.ix] = -mu.all;
      bcn = rep(0, nrow(Acn))
      sensecn = rep("E", length(bcn))

      consmeta = rbind(consmeta, data.frame(type = 'Copy', label = paste('Copy', 1:nrow(Acn)), sense = 'E', b = bcn, stringsAsFactors = F)) 
      
      if (ignore.cons | T)
        {
          Acn[nrow(Acn), ] = 0
          bcn[length(bcn)] = 0
        }      
    }
  else
    { ## abort abort!
      sol = list()
      sol$residual = NA
      sol$beta = beta.guess
      sol$gamma = gamma.guess
      sol$purity = NA
      sol$ploidy = NA
      sol$adj = adj*NA
      sol$nll.cn = NA           
      sol$nll.opt = NA           
      sol$gap.cn = NA
      sol$segstats = segstats[, c('mean', 'sd')]
      sol$segstats$cn = NA
      sol$segstats$ecn.in = NA
      sol$segstats$ecn.out = NA
      segstats$ncn = NA
      sol$segstats$edges.out = sol$segstats$edges.in = rep('', length(segstats))           
      if (edge.slack)
        {
          sol$segstats$eslack.in = NA
          sol$segstats$eslack.out = NA          
        }           
      sol$ploidy.constraints = c(ploidy.min, ploidy.max)
      sol$beta.constraints = c(beta.min, beta.max)
      sol$cn.prior = cn.prior
      sol$slack.prior = slack.prior
      return(sol)
    }
    
  # dup constraints on vertices
  # constrain every vertex to get the same copy number as its reverse complement
  Dcn = Zero[rep(1, length(dup.ix)),, drop = F];
  Dcn[cbind(1:nrow(Dcn), dup.ix)] = 1
  Dcn[cbind(1:nrow(Dcn), og.ix)] = -1
  dcn = rep(0, nrow(Dcn))
  sensedcn = rep("E", nrow(Dcn))
  
  consmeta = rbind(consmeta, data.frame(type = 'Dup', label = paste('Dup', 1:nrow(Dcn)), sense = 'E', b = dcn, stringsAsFactors = F))  
  
  if (edge.slack)
    {
      if (verbose)
        cat('edge slack .. \n')
      
      # dup constraints on (reverse complement) edge.slack
      # (these make sure that reverse complement edge.slacks are
      # given the same solution as their reverse complement)
      Ecn = Zero[rep(1, length(dup.ix)*2),, drop = F];
      Ecn[cbind(1:nrow(Ecn), c(es.s.ix[dup.ix], es.t.ix[dup.ix]))] = 1
      Ecn[cbind(1:nrow(Ecn), c(es.t.ix[og.ix], es.s.ix[og.ix]))] = -1
      ecn = rep(0, nrow(Ecn))
      Dcn = rBind(Dcn, Ecn)
      dcn = c(dcn, ecn);
      sensedcn = c(sensedcn, rep("E", nrow(Ecn)))

      consmeta = rbind(consmeta, data.frame(type = 'EdgeSlack', label = paste('EdgeSlack', 1:nrow(Ecn)), sense = 'E', b = ecn, stringsAsFactors = F)) 
    }
  
  Acn = rBind(Acn, Dcn)
  bcn = c(bcn, dcn)
  sensecn = c(sensecn, sensedcn)
  
  if (!is.na(cn.prior))
    {
      if (verbose)
        cat('cn prior .. \n')            
      Pcn = Zero[rep(1, length(v.ix)+1), ]
      Pcn[cbind(v.ix, v.ix)] = 1
      Pcn[v.ix, ploidy.ix] = -1
      Pcn[cbind(v.ix, d.ix)] = -1
      Pcn[length(v.ix)+1, v.ix] = width(segstats)/sum(as.numeric(width(segstats)))
      Pcn[length(v.ix)+1, ploidy.ix] = -1;
      bpcn = rep(0, nrow(Pcn))
      Acn = rBind(Acn, Pcn)
      bcn = c(bcn, bpcn)
      sensecn = c(sensecn, rep("E", length(bpcn)))
      lb[d.ix] = -Inf;

      consmeta = rbind(consmeta, data.frame(type = 'CNPrior', label = paste('CNPrior', 1:nrow(Pcn)), sense = 'E', b = bpcn, stringsAsFactors = F)) 
    }
  
  if (!any(is.na(purity.prior)))
    {
      if (verbose)
        cat('purity prior .. \n')
      
      Ppd = Zero[1, , drop = FALSE]
      Ppd[1, gamma.ix] = 1
      Ppd[1, pd.ix] = -1
      bpd = 2/purity.prior[1] - 2
      Acn = rBind(Acn, Ppd)
      bcn = c(bcn, bpd)
      sensecn = c(sensecn, 'E')
      lb[pd.ix] = -Inf;
      ub[pd.ix] = -Inf;

      consmeta = rbind(consmeta, data.frame(type = 'PurityPrior', label = paste('PurityPrior', 1:nrow(Ppd)), sense = 'E', b = bpd, stringsAsFactors = F)) 
    }
  
  # ploidy constraints
  Aineq = Zero[rep(1, 2), , drop = F];
  Aineq[1:2 , v.ix] = rbind(width(segstats)/sum(as.numeric(width(segstats))), width(segstats)/sum(as.numeric(width(segstats))))
  bineq = c(pmax(0, ploidy.min), pmin(Inf, ploidy.max))
  senseineq = c("G", "L")

      
  # override gamma 
  if (!is.na(gamma))
    lb[gamma.ix] = ub[gamma.ix] = gamma;

  # override beta
  if (!is.na(beta))
    lb[beta.ix] = ub[beta.ix] = beta.guess;

  # beta.max
  if (!is.na(beta.max))
    ub[beta.ix] = beta.max 

  # beta.max
  if (!is.na(beta.min))
    lb[beta.ix] = beta.min

  # gamma.max
  if (!is.na(gamma.min))
    lb[gamma.ix] = gamma.min

  # gamma.max
  if (!is.na(gamma.max))
    ub[gamma.ix] = gamma.max 
  
  if (!ignore.edge)
    {
      if (verbose)
        cat('edge consistency matrix .. \n')
      
      ## add edge consistency criteria
      ## for every node that is source of an edge
      ## ensure that sum of weights on outgoing edges
      ## = node weight
      ## do the same for nodes that are targets of edges

      v.ix.s = unique(edges[,1])
      v.ix.t = unique(edges[,2])
      Bs = Zero[v.ix.s, , drop = F]
      Bt = Zero[v.ix.t, , drop = F]

      if (nrow(edges)>0)
        {
          Bs[cbind(1:nrow(Bs), v.ix.s)] = 1
          Bt[cbind(1:nrow(Bt), v.ix.t)] = 1
          Bs[cbind(match(edges[,1], v.ix.s), e.ix)] = -1
          Bt[cbind(match(edges[,2], v.ix.t), e.ix)] = -1
          
          if (edge.slack)
            {
              Bs[cbind(1:nrow(Bs), es.s.ix[v.ix.s])] = -1  # provide "fake" edges bringing flux into and out of vertex
              Bt[cbind(1:nrow(Bt), es.t.ix[v.ix.t])] = -1
            }
        }
      # B = rbind(as.matrix(Bs), as.matrix(Bt))
      B = rBind(Bs, Bt)
      
      if (verbose)
        cat('populating linear constraints .. \n')

      consmeta = rbind(consmeta,
        data.frame(type = 'EdgeSource', label = paste('EdgeSource', 1:nrow(Bs)), sense = 'E', b = 0, stringsAsFactors = F),
        data.frame(type = 'EdgeTarget', label = paste('EdgeSource', 1:nrow(Bt)), sense = 'E', b = 0, stringsAsFactors = F)
        )
      
      ## populate linear constraints
      Aed = B;
      bed = rep(0, nrow(B))
      senseed = rep("E", length(bed))      
      
#      Amat = rbind(as.matrix(Acn), as.matrix(Aed), as.matrix(Aineq));
      Amat = rBind(Acn, Aed, Aineq);
      b = c(bcn, bed, bineq);
      sense = c(sensecn, senseed, senseineq);
    }
  else
    {
      Amat = rBind(Acn, Aineq);
      b = c(bcn, bineq);
      sense = c(sensecn, senseineq);
    }

  consmeta = rbind(consmeta, data.frame(type = 'PloidyConstraint', label = paste('PloidyConstraint', 1:nrow(Aineq)), sense = senseineq, b = bineq, stringsAsFactors = F)) 

  ## ecn.out.ub constraints (if any)
  if (any(!is.na(ecn.out.ub)))
    {
      ix = which(!is.na(ecn.out.ub))
      
      Aineq_eout = Zero[ix, , drop = F]
      bineq_eout = ecn.out.ub[ix]
      senseineq_eout = rep('L', length(ix))
      
      for (i in 1:length(ix))
        if (any(iy <<- edges[,1] %in% ix[i]))
          Aineq_eout[i, e.ix[iy]] = 1 ## constrain the sum of edges exiting this vertex
      
      consmeta = rbind(consmeta, data.frame(type = 'EdgeOutUB', label = paste('EdgeOutUB', ix), sense = 'L', b = bineq_eout, stringsAsFactors = F))

      Amat = rBind(Amat, Aineq_eout)
      b = c(b, bineq_eout)
      sense = c(sense, senseineq_eout)
    }

  ## ecn.in.ub constraints (if any)
  if (any(!is.na(ecn.in.ub)))
    {
      ix = which(!is.na(ecn.in.ub))
      
      Aineq_ein = Zero[ix, , drop = F]
      bineq_ein = ecn.in.ub[ix]
      senseineq_ein = rep('L', length(ix))
      
      for (i in 1:length(ix))
        if (any(iy <<- edges[,2] %in% ix[i]))
          Aineq_ein[i, e.ix[iy]] = 1  ## constrain the sum of edges entering this vertex

      consmeta = rbind(consmeta, data.frame(type = 'EdgeInUB', label = paste('EdgeInUB', ix), sense = 'L', b = bineq_ein, stringsAsFactors = F))
      
      Amat = rBind(Amat, Aineq_ein)
      b = c(b, bineq_ein)
      sense = c(sense, senseineq_ein)
    }
  
  # quadratic portion of objective function
  Qobj = Zero;

  if (length(v.ix.c>0))
      Qobj[cbind(s.ix[v.ix.c], s.ix[v.ix.c])] = 1/segstats$sd[v.ix.c]^2;

  EPS = 0.000000000001
  if (na.node.nudge) ## added Monday, Oct 23, 2017 05:17:25 PM to prevent unconstrained nodes from blowing up
      if (length(v.ix.na)>0)
      {
          if (verbose)
              message('NA nudging node!!')
          Qobj[cbind(v.ix.na, v.ix.na)] = EPS
      }

#  if (!ignore.cons)
#    Qobj[s.ix[length(s.ix)], s.ix[length(s.ix)]] = 1

  # linear portion of objective function
  cvec = Zero[,1]

  if (nrow(edges)>0)
    {
      if (verbose)
        cat('Adding ', sum(adj.nudge[edges]), "of edge nudge across", sum(adj.nudge[edges]!=0), "edges\n")

      cvec[e.ix] = -adj.nudge[edges] ### reward each edge use in proportion to position in edge nudge
    }

  if (!is.na(cn.prior))
    Qobj[cbind(d.ix, d.ix)] = 1/cn.prior^2 

  ## the slack prior will determine the degree of "coupling" enforced between neighboring copy states
  ## this should be high if we think that our rearrangement annotation is quite complete
  ## in the end, there will be tension between enforcing edge consistency and consistency with means / sd
  ## abundances at intervals
  if (edge.slack)
    {      
      cvec[c(es.s.ix, es.t.ix)] = 1/slack.prior

      ## let any specified "loose ends" have unpenalized slack
      if (length(loose.ends)>0)
        {
          cat('Relaxing slack penalty on', length(loose.ends), 'loose ends\n')
          cvec[c(es.s.ix[loose.ends], es.t.ix[loose.ends])] = 0
        }
      
      cat(sprintf('Total mass on cn portion of objective function: %s. Total mass on edge slack: %s\n', sum(Qobj[cbind(s.ix, s.ix)]), sum(cvec[cbind(es.s.ix, es.t.ix)])))

      
    }
  
  if (!any(is.na(purity.prior)))
    Qobj[cbind(pd.ix, pd.ix)] = 1/purity.prior[2]^2
  
  if (verbose)
    cat('Beginning optimization .. \n')
  
  # setup MIP
  if (gurobi) # translate into gurobi
    {
      print('Running gurobi!')

      model = list()
      model$A = Amat
      model$rhs = b;
      model$sense = c('E'='=', 'G'='>=', 'L'='<=')[sense]
      model$Q = Qobj;
      model$obj = cvec;
      model$lb = lb;
      model$ub = ub;
      model$vtype = vtype;      
      model$modelsense = 'min';

      if (!is.na(gamma))
        {
          mu_hat = as.vector(round(((segstats$mean-gamma)/beta))); 
          model$start = rep(NA, n);
          model$start[v.ix] = mu_hat;
          model$start[s.ix] = segstats$mean-(mu_hat*beta+gamma);
          model$start[gamma.ix] = gamma;
          model$start[is.infinite(model$start)] = NA;
        }
      
      sol = gurobi(model, params = c(list(TimeLimit=tilim), list(...)));
      sol$xopt = sol$x;
    }
  else
    sol = Rcplex(cvec = cvec, Amat = Amat, bvec = b, sense = sense, Qmat = Qobj, lb = lb, ub = ub, n = nsolutions, objsense = "min", vtype = vtype, control = c(list(...), list(tilim = tilim, epgap = epgap ,mipemphasis = mipemphasis)))
  
  if (is.null(sol$xopt))
    sol.l = sol
  else
    sol.l = list(sol);
  
  adj = as(adj, 'sparseMatrix');
  mu = segstats$mean
  sd = segstats$sd
  segstats = segstats[, c()]
  segstats$mean = mu
  segstats$sd = sd
  segstats$ncn = ncn

  if (verbose)
    cat('Post processing .. \n')
  
  sol.l = lapply(sol.l, function(sol)
  {
           vcn = round(sol$xopt[v.ix])
           ecn = round(sol$xopt[e.ix])
           sol$residual = round(sol$xopt[s.ix])
           sol$beta = sol$xopt[beta.ix]
           sol$gamma = sol$xopt[gamma.ix]
           sol$purity = 2/(2+sol$gamma)
           sol$ploidy = (vcn%*%width(segstats))/sum(as.numeric(width(segstats)))
#           sol$mu.all = mu.all
           sol$adj = adj*0;
           if (length(v.ix.c)>0)
             sol$nll.cn = (sol$xopt[s.ix[v.ix.c]]%*%Qobj[s.ix[v.ix.c], s.ix[v.ix.c]])%*%sol$xopt[s.ix[v.ix.c]]
           else
             sol$nll.cn = NA
           
#           sol$nll.opt = pp.nll(segstats[v.ix.c], sol$purity, sol$ploidy, field = 'mean')$NLL
           if (length(v.ix.c)>0)
             sol$nll.opt = pp.nll(segstats[v.ix.c], gamma = sol$gamma, beta = sol$beta, field = 'mean', field.ncn = field.ncn)$NLL
           else
             sol$nll.opt = NA
           
           sol$gap.cn = as.numeric(1 - sol$nll.opt / sol$nll.cn)
           sol$adj[edges] = ecn;
           sol$segstats = segstats
           sol$segstats$cn = round(vcn)
           sol$segstats$ecn.in = round(Matrix::colSums(sol$adj))
           sol$segstats$ecn.out = round(Matrix::rowSums(sol$adj))
           sol$segstats$edges.in = sapply(1:length(sol$segstats),
             function(x) {ix = which(adj[,x]!=0); paste(ix, '(', sol$adj[ix,x], ')', '->', sep = '', collapse = ',')})
           sol$segstats$edges.out = sapply(1:length(sol$segstats),
             function(x) {ix = which(adj[x, ]!=0); paste('->', ix, '(', sol$adj[x,ix], ')', sep = '', collapse = ',')})
           
           if (edge.slack)
             {
##                sol$segstats$eslack.s = round(sol$xopt[es.s.ix])
##                sol$segstats$eslack.t = round(sol$xopt[es.t.ix])
##                sol$eslack.s = round(sol$xopt[es.s.ix])
##                sol$eslack.t = round(sol$xopt[es.t.ix])
               sol$segstats$eslack.in = round(sol$xopt[es.t.ix])
               sol$segstats$eslack.out = round(sol$xopt[es.s.ix])
               sol$eslack.in = round(sol$xopt[es.t.ix])
               sol$eslack.out = round(sol$xopt[es.s.ix])
             }
           
           sol$ploidy.constraints = c(ploidy.min, ploidy.max)
           sol$beta.constraints = c(beta.min, beta.max)
           sol$cn.prior = cn.prior
           sol$slack.prior = slack.prior

#           if (any(!is.na(sol$eslack.in)))
 #            sol = .correct.slack(sol)

#           print('did not correct slack')

           if (debug)
             browser()
           
           return(sol)
         });
  
  sol.l = sol.l[order(sapply(sol.l, function(x) x$obj))]

  if (length(sol.l)==1)
    sol.l = sol.l[[1]]

  return(sol.l)
}


##############################################################
#' jbaMIP.summarize
#' 
#' summarizes miqp result (i.e. multiple solutions) outputting data frame
#' with summary info
#' 
##############################################################
jbaMIP.summarize = function(sol)
{
  if (!is.null(sol$obj))
    sol = list(sol)
  
  df = data.frame(purity = sapply(sol, function(x) x$purity),
    ploidy = sapply(sol, function(x) x$ploidy),
    obj = sapply(sol, function(x) x$obj),
    max.cn = sapply(sol, function(x) max(round(x$vcn), na.rm = T)),
    max.ecn = sapply(sol, function(x) max(round(x$ecn), na.rm = T)),
    tot.eslack = sapply(sol, function(x) c(sum(x$segstats$eslack.out) + sum(x$segstats$eslack.in), NA)[1]),
    num.eslack = sapply(sol, function(x) c(sum(x$segstats$eslack.out!=0) + sum(x$segstats$eslack.in!=0), NA)[1]),
    max.eslack = sapply(sol, function(x) pmax(max(c(x$segstats$eslack.out, 0)), max(c(x$segstats$eslack.in, 0)))))

  return(df)    
}


####################################################################
#' jabba.melt
#' 
#' @details
#' melt JaBbA graph into "events" that decompose the total ploidy into amplicons (or deleticons, if anti = TRUE)
#' Each amplicons / deleticon is flanked by either (1) junctions (2) loose ends or (3) chromosome ends / telomeres
#'
#' 
#' @param jab JaBbA object "undigested"
#' @param kag karyograph (original karyograph input to JaBbA), if NULL then will "redigest" JaBbA object
#' @param verbose logical flag
#' @param keep.all keep.all (default TRUE) whether to keep 0 copy junctions or collapse segments across these as well
#' @export
####################################################################
jabba.melt = function(jab, anti = FALSE, verbose = FALSE, mc.cores = 1, max.del = 10)
    {
        abs = rbind(jab$ab.edges[,,1], jab$ab.edges[,,2])[, 1:2]
        abs = abs[!is.na(abs[,1]), ]
        adj = data.table(i = abs[,1], j = abs[,2])
        adj[, cn := jab$adj[cbind(i, j)]]
        setkeyv(adj, c('i', 'j'))
        junc.right = adj[, sum(cn), keyby = i]
        junc.left = adj[, sum(cn), keyby = j]

        gr = gr2dt(jab$segstats)[, seg.id := 1:length(seqnames)][loose == FALSE & strand == '+' & !is.na(cn), ]
        gr[, cluster := {
            tmp = rle(cn)
            rep(paste(seqnames[1], 1:length(tmp$values), sep = '.'), tmp$length)
        }, by = seqnames]
        
        lix = jab$segstats$loose        
        gr[, loose.right := Matrix::rowSums(as.matrix(jab$adj[seg.id, lix, drop = FALSE]))]
        gr[, loose.left := Matrix::colSums(as.matrix(jab$adj[lix, seg.id, drop = FALSE]))]

        gr[, loose.right := Matrix::rowSums(as.matrix(jab$adj[seg.id, lix, drop = FALSE]))]
        gr[, loose.left := Matrix::colSums(as.matrix(jab$adj[lix, seg.id, drop = FALSE]))]
        
        if (nrow(junc.left)>0)
            {
                gr[, ab.left := ifelse(is.na(junc.left[list(seg.id), V1]), 0, junc.left[list(seg.id), V1])]
                gr[, ab.left.ix := gr[ , adj[, ][cn>0, i[1], keyby = j][list(seg.id), V1]]]
            }
        else
            gr[, ab.left := 0]

        if (nrow(junc.right)>0)
            {
                gr[, ab.right := ifelse(is.na(junc.right[list(seg.id), V1]), 0, junc.right[list(seg.id), V1])]
                gr[, ab.right.ix := gr[ , adj[, ][cn>0, j[1], keyby = i][list(seg.id), V1]]]
            }
        
        else
            gr[, ab.right := 0]

        gr[, id := 1:nrow(gr)]
        setkey(gr, id)

        .flip = function(gr)
            {
                gr = as.data.table(as.data.frame(gr))
                if (nrow(gr) <= 1)
                    return(gr)
                max.cn = gr[, pmin(max(cn), max.del)] ## cap max cn for "anti analysis" since we may not care about dels in high copy contexts, and want to avoid edge jabba cases where a region got infinite copy number 

                ## flip copy number upside down
                gr[, cn := max.cn-cn]

                ## assign right.ab[i] and right.loose[i] to left.ab[i+1]

                for (i in 1:(nrow(gr)-1))
                    {
                        tmp.loose = gr[i+1, loose.left]
                        tmp.ab = gr[i+1, ab.left]
                        gr[i+1, loose.left := gr$loose.right[i]]
                        gr[i+1, ab.left := gr$ab.right[i]]
                        gr[i, loose.right := tmp.loose]
                        gr[i, ab.right := tmp.ab]
                    }
                return(gr)
            }
        
        if (anti)
            gr = rbindlist(lapply(split(gr, gr$seqnames), .flip))
            
        .fun = function(gr)
            {                
                ## this assumes the input is from one seqneme
                gr = as.data.table(as.data.frame(gr))
                cl = split(gr$id, gr$cluster)
                ix = order(-gr[, cn])
                
                if (nrow(gr)==0)
                    return(NULL)

                gr[, done := FALSE]
                gr[1, loose.left := cn] ## give telomeres temporary loose ends                
                gr[nrow(gr), loose.right := cn]
                
                out = data.table(seqnames = rep(as.character(NA), nrow(gr)), start = as.numeric(NA), end = as.numeric(NA), cn = as.numeric(NA),
                    cn.og = as.numeric(NA),
                    left.ab = as.numeric(NA),
                    right.ab = as.numeric(NA), 
                    left.loose = as.numeric(NA),
                    right.loose = as.numeric(NA),
                    left.ix = as.numeric(NA),
                    right.ix = as.numeric(NA))
                out.i = 1
                
                for (i in gr$id[ix])
                    {
                        if (verbose)
                            cat('.')
                        setkey(gr, id)
                        if (!gr[list(i), done])
                            {
#                                if (i==21)
 #                                   browser()
                                this.cl = cl[[gr[list(i), cluster]]]
                                out[out.i, ":="(left.ix = this.cl[1],
                                                right.ix = this.cl[length(this.cl)])]
                                left.cl = out[out.i, gr[list(c(left.ix-1, left.ix)), setdiff(cluster, NA)[1]]]
                                right.cl = out[out.i, gr[list(c(right.ix+1, right.ix)), setdiff(cluster, NA)[1]]]
                                ## neighbor CN is max of left and right
                                left.cn = out[out.i, gr[list(left.ix-1), cn]]
                                right.cn = out[out.i, gr[list(right.ix+1), cn]]
                                neighbor.cn = max(c(0, left.cn,  right.cn), na.rm = TRUE)

                                out.cn = gr[list(i), cn - neighbor.cn]
                                
                                out[out.i, ":="(cn = out.cn, seqnames = gr[list(i), seqnames], start = gr[list(left.ix), start], end = gr[list(right.ix), end], cn.og = gr[list(i), cn])]
                                tmp = out[out.i, pmin(gr[list(left.ix), ab.left], out.cn)]
                                out[out.i, ':='(left.ab = tmp, left.loose = out.cn - tmp)]
                                
                                tmp = out[out.i, pmin(gr[list(right.ix), ab.right], out.cn)]
                                out[out.i, ':='(right.ab = tmp, right.loose = out.cn - tmp)]
                                
                                        # update gr
                                        # subtracting the copies we have assigned to this outgoing event

                                gr[list(out[out.i, left.ix]), ab.left := ab.left - out[out.i, left.ab]]
                                gr[list(out[out.i, left.ix]), loose.left := loose.left - out[out.i, left.loose]]
                                gr[list(out[out.i, right.ix]), ab.right := ab.right - out[out.i, right.ab]]
                                gr[list(out[out.i, right.ix]), loose.right := loose.right - out[out.i, right.loose]]
                                gr[list(this.cl), ":="(cn = cn - out.cn)]

                                ## update clusters - merge left and right depending on which cn are equals
                                ## the new copy number for this interval should be equal to the left or right cn
                                this.cli = gr[list(i), cluster]
                                left.cn = ifelse(is.na(left.cn), -Inf, left.cn)
                                right.cn = ifelse(is.na(right.cn), +Inf, right.cn)
                                if (left.cn == right.cn)
                                    {
                                        gr[list(this.cl), ":="(cluster = left.cl)]
                                        gr[list(cl[[right.cl]]), ":="(cluster = left.cl)]
                                        cl[[left.cl]] = do.call('c', cl[unique(c(left.cl, this.cli, right.cl))])
                                    }
                                else if (left.cn == gr[list(i), cn])
                                    {
                                        gr[list(this.cl), ":="(cluster = left.cl)]
                                        cl[[left.cl]] = do.call('c', cl[unique(c(left.cl, this.cli))])
                                    }
                                else                                
                                    {
                                        gr[list(this.cl), ":=" (cluster = right.cl)]
                                        cl[[right.cl]] = do.call('c', cl[unique(c(this.cli, right.cl))])
                                    }                                   
                                gr[list(this.cl), done := TRUE] ## don't go back to any interval in this cluster
                                ## advance output
                                
                                out.i = out.i +1
                            }
                    }
                
                out[, left.tel := 0]
                out[, right.tel := 0]
                out[left.ix == gr[, id[1]], ":="(left.tel = left.loose, left.loose = 0)]
                out[right.ix == gr[, id[nrow(gr)]], ":="(right.tel = right.loose, right.loose = 0)]
                return(out)
            }

        tmp = mclapply(split(gr, gr$seqnames), .fun, mc.cores = mc.cores)
        out = rbindlist(tmp)[!is.na(cn), ][cn>0, ]
        out = seg2gr(out, seqlengths = seqlengths(jab$segstats))
        out$cn.max = gr.val(out, jab$segstats, 'cn', FUN = max, weighted = FALSE)$cn
        out$cn.min = gr.val(out, jab$segstats, 'cn', FUN = min, weighted = FALSE)$cn
        return(out)
    }


####################################################################
#' JaBbA.digest
#' 
#' @details
#' processes JaBbA object 
#' (1) collapsing segments with same copy number that lack loose ends 
#' (2) (optional) +/- adds segments correponding to loose ends
#' (3) outputting edges data frame with colors, and other formatting information
#' (4) outputting junctions GRangesList with copy number, color, lty and other plotting components
#'
#' TODO: replace with proper object instantiator
#'
#' @param jab JaBbA object "undigested"
#' @param kag karyograph (original karyograph input to JaBbA), if NULL then will "redigest" JaBbA object
#' @param verbose logical flag
#' @param keep.all keep.all (default TRUE) whether to keep 0 copy junctions or collapse segments across these as well
#' @export
############################################
JaBbA.digest = function(jab, kag = NULL, verbose = T, keep.all = T)
  {
      if (is.null(kag))
          {
              if (verbose)
                  cat('Redigesting..\n')
              ## redigesting
              ix = rep(TRUE, length(jab$segstats))
              if (!is.null(jab$segstats$loose))
                  ix = !jab$segstats$loose
              
              jab = list(segstats = jab$segstats[ix], adj = jab$adj[ix, ix], ab.edges = jab$ab.edges, purity = jab$purity, ploidy = jab$ploidy)

              tmp.adj = jab$adj*0
              jab$segstats$ref.child = gr.match(gr.flipstrand(flank(gr.flipstrand(jab$segstats),1)), jab$segstats, ignore.strand = F)
              jab$segstats$ref.parent = gr.match(flank(jab$segstats,1), jab$segstats, ignore.strand = F)
              ix = suppressWarnings(cbind(1:length(jab$segstats), jab$segstats$ref.child))
              nnaix = Matrix::rowSums(is.na(ix))==0
              if (any(nnaix))
                  tmp.adj[ix[nnaix,, drop = F]] = 1

              nnab = !ifelse(is.na(jab$ab.edges[,3,1]), TRUE, jab$ab.edges[,3,1]==0)
              
              if (any(nnab))
                  tmp.adj[rbind(jab$ab.edges[nnab, 1:2, 1], jab$ab.edges[nnab, 1:2, 2])] = 1
                                
              kag = list(adj = tmp.adj, ab.edges = jab$ab.edges, segstats = jab$segstats)
          }
      
    if (any(dim(jab$adj) != dim(kag$adj)))
      stop('JaBbA and karyograph object mismatch')
    
      bk.adj = kag$adj ## the background graph will make sure our collapsed paths are reference adjacent
      
      nnab = !ifelse(is.na(kag$ab.edges[,3,1]), TRUE, kag$ab.edges[,3,1]==0)
      
    if (!keep.all & nrow(kag$ab.edges)>0) ## we can throw out unused aberrant edges by throwing them out of the background graph (used ones will still be in jab$adj)
      {
        bk.adj = kag$adj
        bk.adj[rbind(kag$ab.edges[nnab,1:2, 1])] = 0 
        bk.adj[rbind(kag$ab.edges[nnab,1:2, 2])] = 0
      }
    
      adj = sign(bk.adj)*0.01 + jab$adj ## keep a hint of 0 copy edges
      segstats = jab$segstats
      
      lends = loose.ends(jab, kag)
      
      if (!is.null(lends))
          lends = lends[lends$type != 'type4'] ## don't include type 4 loose ends (i.e. unused rearrangements)
      
    if (length(lends)>0)
      {
        strand(lends) = '+'        
        lends = c(lends, gr.flipstrand(lends))
        lends$partner.id = gr.match((lends), jab$segstats, ignore.strand = F)
        lends$id = nrow(adj) + c(1:length(lends))
        lends$right = end(lends) == end(jab$segstats)[lends$partner.id]            
        adj.plus = rBind(cBind(adj, sparseMatrix(1,1,x = 0, dims = c(nrow(adj), length(lends)))),
          cBind(sparseMatrix(1,1,x = 0, dims = c(length(lends), ncol(adj))), sparseMatrix(1,1,x = 0, dims = c(length(lends), length(lends)))))
        
        ## right side ends of '+' and left side ends of '-' are sinks
        sink.ix = as.logical((lends$right & as.logical( strand(lends)=='+') )| (!lends$right & as.logical( strand(lends)=='-')) )
        adj.plus[cbind(lends$partner.id, lends$id)[sink.ix, , drop = F]] = lends$cn[sink.ix]+0.01
        adj.plus[cbind(lends$id, lends$partner.id)[!sink.ix, , drop = F]] = lends$cn[!sink.ix]+0.01
        adj = adj.plus
        lends$loose = T
        segstats$loose = F    
        segstats = grbind(jab$segstats, lends)
        values(segstats) = rrbind(values(jab$segstats), values(lends))
      }
    else
      segstats$loose = F    
    
    out = list()
    
    ## now we have augmented adjacency matrix with loose ends, let's simplify the structure
    ## by collapsing all simple paths

    collapsed = collapse.paths(adj, verbose = verbose)

    ## new segstats formed by reducing "collapsed' sets
    segstats$set.id = collapsed$map
    out$segstats = gr.fix(gr.simplify(segstats[unlist(lapply(collapsed$sets, sort))], val = rep(1:length(collapsed$sets), sapply(collapsed$sets, length))), segstats)

    tmp.ss = gr.string(gr.stripstrand(out$segstats), other.col = 'loose')
    check1 = all(table(match(tmp.ss, tmp.ss)))
    check2 = identical(1:length(collapsed$sets), sort(out$segstats$set.id))
      
    if (!check1 | !check2) ## quick sanity check to make sure we didn't screw up collapsing 
      stop('collapse yielded funny / missing segments')          
    else
      out$segstats = out$segstats[match(1:length(collapsed$sets), out$segstats$set.id), c('cn')]
    
#    out$segstats$og.ix = sapply(collapsed$sets, function(x) paste(sort(x), collapse = ','))

    tmp.start.ix = sapply(collapsed$sets, function(x) sort(x)[1])
    tmp.end.ix = sapply(collapsed$sets, function(x) -sort(-x)[1])    
    out$segstats$start.ix = ifelse(as.logical(strand(out$segstats)=='+'), tmp.start.ix, tmp.end.ix)
    out$segstats$end.ix = ifelse(as.logical(strand(out$segstats)=='+'), tmp.end.ix, tmp.start.ix)
    out$segstats$eslack.in = segstats$eslack.in[out$segstats$start.ix]
    out$segstats$eslack.out = segstats$eslack.out[out$segstats$end.ix]
    out$segstats$loose = F
    if (any(loose.ix <- which(segstats$loose)))
      out$segstats$loose = out$segstats$start.ix %in% loose.ix
            
    adj.new.ix = out$adj = collapsed$adj*0 ## rewire copy numbers and edge indices according to collapsed
    edge.ix = which(collapsed$adj!=0, arr.ind = T)

    if (nrow(edge.ix)>0)
      {
        out$adj[edge.ix] = round(adj[cbind(out$segstats$end.ix[edge.ix[,1]], out$segstats$start.ix[edge.ix[,2]])])
        adj.new.ix[edge.ix] = 1:nrow(edge.ix)
      }
    
    out$ab.edges = array(NA, dim = c(nrow(kag$ab.edges), 3, 2), dimnames = list(NULL, c('from', 'to', 'edge.ix'), c('+', '-')))

    ## match ab edges to new graph, excluding any edges that aren't included in graph (i.e. not given >0 copy number)
    ## (tmp.ix may map some ab.edges to "internal" vertices, so need to weed these out via keep TODO:cleanup)
      if (any(nnab))
          {
              ## og junctions
              tmp.ix = cbind(rep(NA, nrow(kag$ab.edges)), rep(NA, nrow(kag$ab.edges)))
              tmp.ix[nnab,] = cbind(collapsed$map[kag$ab.edges[nnab,1,1]], collapsed$map[kag$ab.edges[nnab,2,1]])
              keep = rep(FALSE, length(nnab))
              keep[nnab] = (out$segstats$end.ix[collapsed$map[kag$ab.edges[nnab,1,1]]] == kag$ab.edges[nnab,1,1]) & (out$segstats$start.ix[collapsed$map[kag$ab.edges[nnab,2,1]]] == kag$ab.edges[nnab,2,1])
              tmp.ix[!keep, ] = NA ## not really needed but let's keep it
              if (any(keep))
                  out$ab.edges[keep,,1] = cbind(tmp.ix[keep, , drop = F], adj.new.ix[tmp.ix[keep, , drop = F]])

              ## rev comp junctions
              tmp.ix = cbind(rep(NA, nrow(kag$ab.edges)), rep(NA, nrow(kag$ab.edges)))              
              tmp.ix[nnab] = cbind(collapsed$map[kag$ab.edges[nnab,1,2]], collapsed$map[kag$ab.edges[nnab,2,2]])
              keep = rep(FALSE, length(nnab))
              keep[nnab] = (out$segstats$end.ix[collapsed$map[kag$ab.edges[nnab,1,2]]] == kag$ab.edges[nnab,1,2]) & (out$segstats$start.ix[collapsed$map[kag$ab.edges[nnab,2,2]]] == kag$ab.edges[nnab,2,2])
              tmp.ix[!keep, ] = NA ## not really needed, but let's keep it    
              if (any(keep))
                  out$ab.edges[keep,,2] = cbind(tmp.ix[keep, , drop = F], adj.new.ix[tmp.ix[keep, , drop = F]])
          }
      
    ## convert to "simplified form"
    out$edges = data.frame(from = edge.ix[,1], to = edge.ix[,2], cn = out$adj[edge.ix])

    estr = paste(edge.ix[,1], edge.ix[,2])
    abestr = paste(out$ab.edges[,1,1:2], out$ab.edges[,2,1:2])

    if (nrow(out$edges)>0)
        {
            out$edges$type = 'reference'
            if (any(ix <- estr %in% abestr))
                out$edges$type[ix] = 'aberrant'
            if (any(ix <- out$segstats$loose[out$edges[,'from']] | out$segstats$loose[out$edges[,'to']]))
                out$edges$type[ix] = 'loose'
    
            out$edges$col = ifelse(out$edges$type == 'aberrant', ifelse(out$edges$cn>0, alpha('red', 0.4), alpha('purple', 0.3)), alpha('gray', 0.2))

            loose.ix = which(out$edges$type == 'loose')
            out$edges$h = 1
        }
      
    if (length(loose.ix)>0)
      {
        seg.map = match(out$segstats, gr.flipstrand(out$segstats)) ## maps segs to rev comp
        ed.map = match(paste(out$edges[, 'from'], out$edges[, 'to']), paste(seg.map[out$edges[, 'to']], seg.map[out$edges[, 'from']])) ## maps edges to rev comp
        temp.ix = which(ed.map>(1:length(ed.map)));
        ed.id = ed.map
        ed.id[temp.ix] = temp.ix  ## edges whose rev comp has higher id we name with their index, and we name their rev comp with their index        
        out$edges$col[loose.ix] = alpha('blue', 0.6)
        rh = 0.5 + runif(length(loose.ix)/2)
        out$edges$h[loose.ix] = rh[match(ed.id[loose.ix], unique(ed.id[loose.ix]))]
        # out$edges$h = ifelse(out$edges$type == 'loose', rand(nrow(out$edges)), 1)
      }

    if (nrow(out$edges)>0)
        {
            out$edges$lwd = ifelse(out$edges$type == 'aberrant', 1 + log2(0.2*out$edges$cn+2), 1)
            out$edges$lwd[out$edges$cn==0] = 0.2    
            out$edges$lty = ifelse(out$edges$type == 'loose', 3, ifelse(out$edges$cn==0, 2, 1))
            out$edges$lty[out$edges$cn==0] = 3
            out$edges$col[out$edges$cn==0 & out$edges$type == 'loose'] = alpha('purple', 0.6)
            out$edges$cex.arrow = 0
            out$edges$v = 1
            out$edges$not.flat = out$edges$type == 'aberrant'
            out$edges$v[out$edges$type == 'aberrant'] = 2
            out$edges$h[out$edges$type == 'aberrant'] = 2
            out$edges$dangle.w = 0.5
        }
        
    out$G = graph.adjacency(adj.new.ix, weighted = 'edge.ix')
    
    out$segstats$edges.in = sapply(1:length(out$segstats),
      function(x) {ix = which(out$adj[,x]!=0); paste(ix, '(', out$adj[ix,x], ')', '->', sep = '', collapse = ',')})
    out$segstats$edges.out = sapply(1:length(out$segstats),
      function(x) {ix = which(out$adj[x, ]!=0); paste('->', ix, '(', out$adj[x,ix], ')', sep = '', collapse = ',')})

    pos.ix = which( as.logical( strand(out$segstats)=='+') )
    out$segstats$tile.id = match(gr.stripstrand(out$segstats), gr.stripstrand(out$segstats[pos.ix]))
    
    ss = out$segstats
    ss$right = segstats$right[ss$start.ix]
    ss$partner.id = segstats$partner.id[ss$start.ix]
    ss$col = ifelse(ss$loose, alpha('blue', 0.3), alpha('gray', 0.5))
    ss$border = ifelse(ss$loose, alpha('blue', 0.3), alpha('black', 0.5))
    ss$ywid = 0.8

    if (any(ss$loose))
      {
#        ss$cn[ss$loose] = ss$cn[ss$loose]+0.5
        ss$cn[ss$loose] = ifelse(ss$right[ss$loose], segstats$cn[ss$partner.id[ss$loose]]*1.2, segstats$cn[ss$partner.id[ss$loose]]*1.2)
        ss[ss$loose] = GenomicRanges::shift(ss[ss$loose], ifelse(ss[ss$loose]$right, -500, 500)) + 500
#        ss[ss$loose] = gr.flipstrand(ss[ss$loose]) + 100
        ss[ss$loose]$ywid = 0.001
        ss[ss$loose]$col = alpha('white', 0)
        ss[ss$loose]$border = alpha('white', 0)
      }
    
    out$td = gTrack(ss, y.field = 'cn', edges = out$edges[order(out$edges$cn), ], name ='JaBbA', angle = 0)

    out$purity = jab$purity
    out$ploidy = jab$ploidy

    return(out)
  }

####################################################################
#' jbaMIP.process
#' 
#' process jbaMIP solution "sol" given original graph "g" (karyograph() list output)
#' into JaBbA object
#' 
#' output is 
#' 
#' @param sol JaBbA object
#' @param allelic logical flag specifying whether object is allelic
#' @return
#' list with items:
#' $B incidence matrix of augmented graph (including slack vertices) (vertices x edges)
#' rownames of $B are vertex names of $G and colnames of B are named with character version of their $G indices
#' (i.e. column order of B  respects the original edge order in the solution)
#'
#' $e edge constraints for downstream karyoMIP, i.e the copy numbers at the edges
#' $e.ij numedges x 2 vertex pair matrix denoting what are the vertex pairs corresponding to the cols of $B and entries of $e, $eclass, $etype etc
#' $eclass id for each unique edge / anti-edge equivalence class
#' $etype specifies whether edge is slack or nonslack
#' @export
###################################################################
jbaMIP.process = function(
  ## output of jbaMIP, sol$segstats needs to have field $tile.id whose unique values appear exactly twice in the object, 
  ## corresponding to + and - strands of the same interval
  sol,
  allelic = F
  )
  {
    if (allelic)
      sol = list(segstats = sol$asegstats, adj = sol$aadj)
    
    if (!all(c('segstats', 'adj') %in% names(sol)))
      stop('sol must be output of jbaMIP()')

    if (is.null(sol$segstats$tile.id))
      stop('sol$segstats must be populated with tile.id')
    else
      {
        if (!all(table(sol$segstats$tile.id)==2))
          stop('sol$segstats$tile.id are malformed, there should be exactly two instances of each tile.id in sol$segstats, one for the positive and one for the negative strand of the same interval')

        tmp = lapply(split(1:length(sol$segstats$tile.id), sol$segstats$tile.id), rev)

        recip.ix = rep(NA, length(sol$segstats))
        recip.ix[order(sol$segstats$tile.id)] = unlist(tmp)
      }
    
    if (is.null(sol$segstats$eslack.in))
      sol$segstats$eslack.in = sol$segstats$slack.in

    if (is.null(sol$segstats$eslack.out))
      sol$segstats$eslack.out = sol$segstats$slack.out

    ed.ij = which(sol$adj!=0, arr.ind = T)

    ## B is vertices x edges (i.e. signed incidence matrix)
    B = sparseMatrix(c(ed.ij[,1], ed.ij[,2]), rep(1:nrow(ed.ij), 2), x = rep(c(-1.00001, 1), each = nrow(ed.ij)), dims = c(nrow(sol$adj), nrow(ed.ij)))
    
    rownames(B) = 1:nrow(B)
    
    tmp.ix = which(abs(B)>=1)
    B[tmp.ix] = round(B[tmp.ix]) ## "0.00001" hack to take care of eclass matching below, these are length 1 self loop edge cases
    
    ix.tel.5 = which(Matrix::colSums(sol$adj!=0)==0)  ## make fake slacks for telomeres    
    sol$segstats$eslack.in[ix.tel.5] = sol$segstats$cn[ix.tel.5]  

    ix.tel.3 = which(Matrix::rowSums(sol$adj!=0)==0)
    sol$segstats$eslack.out[ix.tel.3] = sol$segstats$cn[ix.tel.3]  ## make fake slacks for telomeres
    
    ix.eslack.out = which(sol$segstats$eslack.out!=0);
    names(ix.eslack.out) = paste('out slack', ix.eslack.out)
    ix.eslack.in = which(sol$segstats$eslack.in!=0);
    names(ix.eslack.in) = paste('in slack', ix.eslack.in)

    names(ix.eslack.in)[ix.eslack.in %in% ix.tel.3] = paste(names(ix.eslack.in)[ix.eslack.in %in% ix.tel.3], 'tel')
    names(ix.eslack.out)[ix.eslack.out %in% ix.tel.5] = paste(names(ix.eslack.out)[ix.eslack.out %in% ix.tel.5], 'tel')

    ## we add "slack edges" and "slack nodes" to incidence matrix
    Zero = sparseMatrix(1, 1, x = 0, dims = c(length(ix.eslack.in) + length(ix.eslack.out), ncol(B)))

    if (nrow(Zero)>0)
      rownames(Zero) = c(paste('slack in', 1:length(ix.eslack.in)), paste('slack out', 1:length(ix.eslack.out)))

    Bs = rBind(B, Zero)
    ed.ij = rbind(ed.ij, cbind(ix.eslack.out, NA), cbind(NA, ix.eslack.in))
    
    Is = Diagonal(n = nrow(Bs), rep(1, nrow(Bs)))
    
    Bs = cBind(Bs, -Is[, ix.eslack.out], Is[, ix.eslack.in])
    colnames(Bs) = c(as.character(1:ncol(B)), names(ix.eslack.out), names(ix.eslack.in))
    
    ## map new "slack nodes" to their reciprocals
    recip.ix = c(recip.ix,
      nrow(B) + length(ix.eslack.out) +  match(recip.ix[ix.eslack.out], ix.eslack.in),
      nrow(B) + match(recip.ix[ix.eslack.in], ix.eslack.out)
      )

    ## match matrix against its reverse complement (i.e. rotation) to find reciprocal edges
    erecip.ix = mmatch(t(Bs), t(-Bs[recip.ix, ])) ## maps edges to their reciprocals

    tmp.na = which(is.na(erecip.ix))
    if (length(tmp.na)>0) ## fix the self loops so that they match
      erecip.ix[tmp.na] = tmp.na[mmatch(t(Bs[1:nrow(Bs), tmp.na]), t(Bs[recip.ix,tmp.na, drop = F]))]
    
    ## now use this mapping to define edge equivalence classes
##    rmat = t(apply(cbind(erecip.ix, erecip.ix[erecip.ix]), 1, sort)) ## length(erecip.ix) x 2 matrix of edge ids and their reciprocal, sorted

    rmat = cbind(pmin(erecip.ix, erecip.ix[erecip.ix]), pmax(erecip.ix, erecip.ix[erecip.ix]))
    
    ## eclass will map length(erecip.ix) edges to length(erecip.ix)/2 edge equivalence class ids
    eclass = mmatch(rmat, rmat[!duplicated(rmat), ])

    Bs = round(Bs) ## remove the 0.0001 dummy coefficients (i.e. for self loops)
    
    ## e will store observed copy states corresponding to edges (i.e. columns of Bs)
    e = c(sol$adj[which(sol$adj!=0)], sol$segstats$eslack.out[ix.eslack.out],  sol$segstats$eslack.in[ix.eslack.in])
    
    return(list(e = e, e.ij = ed.ij, B = Bs, eclass = eclass, etype = c(ifelse(grepl('slack', colnames(Bs)), 'slack', 'nonslack'))))
  }

#################################################
#' jbaMIP.allelic
#' 
#' Takes adj and segstats from output from jbaMIP and
#' a granges of het.sites with $ref.count and $alt.count
#' 
#' assumes segstats has fields $cn populated
#' and adj has copy states
#'
#' @param adj adjacency matrix populated with total copy counts on junctions
#' @param segstats granges tiling genome populated with total copy counts on interval, mu_high, sd_high, mu_low, sd_low variables on alleles
#' @param purity purity from solution
#' @param gamma gamma param from jbaMIP
#' @param slack.prior 1/slack.prior = penalty for each slack i.e. loose end copy in solution
#' @return
#'
#' list with fields
#' $segstats = input segstats annotated with fitted cn.high, cn.low columns
#' $asegstats = output "allelic" segstats, with $cn, $parent.node, $eslack.in, $eslack.out, $phased fields filled in
#' $adj = output length(segstats) x length(segstats) x 2  "allelic" adjacency matrix with inferred allelic copy numbers on edges
#' $aadj = flattened output 2*length(segstats) x 2* length(segstats)  "allelic" adjacency matrix with inferred allelic copy numbers on edges
#' 
#' @export
#################################################
jbaMIP.allelic = function(
  adj, ## adjacency matrix populated with total copy counts on junctions
  segstats,  ## granges tiling genome populated with total copy counts on interval, mu_high, sd_high, mu_low, sd_low variables on alleles
  purity, ## purity from solution
  gamma, ## gamma param from jbaMIP
  partition = T, 
  slack.prior = 0.001  
  )
  {    
    ploidy = sum(segstats$cn)/sum(as.numeric(width(segstats)))
    
    mu = c(segstats$mu_high, segstats$mu_low)
    ix = !is.na(mu)
    total = sum((mu * 2 * width(segstats))[ix])
    sw = sum(as.numeric(2*width(segstats))[ix])
    
    gamma = 2*(1-purity)/purity  
    beta = (2*(1-purity)*sw + purity*ploidy*sw) / (purity * total)
    
    mu_high = segstats$mu_high*beta + gamma
    sd_high = segstats$sd_high*beta
    mu_low = segstats$mu_low*beta + gamma
    sd_low = segstats$sd_low*beta
    
    ## find the reference junctions
    ord.ix = order(segstats)
    rev.ix = as.logical(strand(segstats[ord.ix]) == '-')
    ord.ix = c(ord.ix[!rev.ix], rev(ord.ix[rev.ix]))
    
    ref.jun = cbind(ord.ix[-length(ord.ix)], ord.ix[-1])
    ref.jun = ref.jun[adj[ref.jun]>0, ]
    ab.adj = adj
    ab.adj[ref.jun] = 0
    ab.jun = which(ab.adj!=0, arr.ind = T) ## aberrant junctions

    ## this will map vertices to their (positive) duplicate
    ## doing this will help contain some of the dimensionality
    dup.vmap = 1:length(segstats)   
    pos.ix = which(as.logical( strand(segstats)=='+') )# "neg vertices" duplicates of pos vertices
    neg.ix = which( as.logical( strand(segstats)=='-') )
    dup.ix = suppressWarnings(neg.ix[gr.match(segstats[pos.ix], segstats[neg.ix])])
    dup.vmap[dup.ix] = pos.ix ## map neg vertices to their positive parent

    ## find dup ref and ab junctions
    ## these are neg-neg junctions (dup of pos-pos)
    ## (all neg-pos and pos-neg junctions are unique)    
    dup.ref.emap = 1:nrow(ref.jun)
    ref.pos.ix = ref.jun[,1] %in% pos.ix & ref.jun[,2] %in% pos.ix
    ref.neg.ix = ref.jun[,1] %in% neg.ix & ref.jun[,2] %in% neg.ix
    ref.dup.ix = mmatch(ref.jun[ref.pos.ix, ], cbind(dup.vmap[ref.jun[ref.neg.ix,2]], dup.vmap[ref.jun[ref.neg.ix,1]])) ## ij ~ n(j)n(i)
    dup.ref.emap[ref.dup.ix] = ref.pos.ix
    
    dup.ab.emap = 1:nrow(ab.jun)
    ab.pos.ix = ab.jun[,1] %in% pos.ix & ab.jun[,2] %in% pos.ix
    ab.neg.ix = ab.jun[,1] %in% neg.ix & ab.jun[,2] %in% neg.ix
    ab.dup.ix = mmatch(ab.jun[ab.pos.ix, ], cbind(dup.vmap[ab.jun[ab.neg.ix,2]], dup.vmap[ab.jun[ab.neg.ix,1]])) ## ij ~ n(j)n(i)
    dup.ab.emap[ab.dup.ix] = ab.pos.ix

    ## we will need the following variables
    ## a1, a2 = length(segstats) allelic copy states for low (1) and high (2) state
    ## is1, is2 = incoming allele slack for low vs high state
    ## os1, os2 = outgoing allele slack for low vs high state    
    ## r11, r12, r21, r22 = allelic copies on reference junctions for {low, high} x {low, high} combos
    ## n11, n12, n21, n22 = allelic copies on non-reference junctions for {low, high} x {low, high} combos
    ## i_r11, i_r12, i_r21, i_r22 = binary indicator variables representing positivity of the allelic copy state for ref junctions
    ## i_n11, i_n12, i_n21, i_n22 = binary indicator variables representing positivity of the allelic copy state for nonref junctions
    ## ns = (linearly penalized) non reference allelic junction slack (length non reference edges)
    ## eps1, eps2 = (quadratic penalized) epsilon residual between observed allelic segment means and the integer fit 
    varmeta = data.frame( ## meta data of all variables we will be using
      var = c(
        rep(c('a1', 'a2'), each = length(segstats)),
        rep(c('is1', 'is2'), each = length(segstats)),
        rep(c('os1', 'os2'), each = length(segstats)),
        rep(c('r11', 'r12', 'r21', 'r22'), each = nrow(ref.jun)),
        rep(c('n11', 'n12', 'n21', 'n22'), each = nrow(ab.jun)),
        rep(c('i_r11', 'i_r12', 'i_r21', 'i_r22'), each = nrow(ref.jun)),
        rep(c('i_n11', 'i_n12', 'i_n21', 'i_n22'), each = nrow(ab.jun)),
        rep(c('i_is1', 'i_is2', 'i_os1', 'i_os2'), each = length(pos.ix)),
        rep('ns', nrow(ab.jun)),
        rep(c('eps1', 'eps2'), each = length(pos.ix))
        ),
      parent = c(
        rep(1:length(segstats), 2),
        rep(1:length(segstats), 2),
        rep(1:length(segstats), 2),
        rep(1:nrow(ref.jun), 4),
        rep(1:nrow(ab.jun), 4),
        rep(1:nrow(ref.jun), 4),
        rep(1:nrow(ab.jun), 4),
        rep(pos.ix, 4),
        1:nrow(ab.jun), 
        rep(pos.ix, 2)
        ),
      stringsAsFactors = F)
    varmeta$vtype = 'I'
    varmeta$vtype[grepl('i_', varmeta$var)] = 'B'
    varmeta$vtype[grepl('eps', varmeta$var)] = 'C'    
    varmeta$lb = 0
    varmeta$lb[varmeta$vtype == 'C'] = -Inf
    varmeta$ub = Inf
    varmeta$id = 1:nrow(varmeta)
    rownames(varmeta) = paste(varmeta$var, varmeta$parent)
    
    var = split(1:nrow(varmeta), varmeta$var) ## handy structure to keep track of variables
    Zero = sparseMatrix(1, 1, x = 0, dims = c(nrow(varmeta), nrow(varmeta)))
    
    ## we will have the following sets of constraints:
    ##
    ## copy1, copy2 = constraints linking counts to copy numbers via alpha, beta, and residual
    ## acopysum = constraints constraining allelic reference copy sums to equal total sums
    ## rcopysum, ncopysum = constraints constraining allelic junction sums to equal total sums
    ## iscopysum = constraints constraining incoming slack sums to equal total sums
    ## oscopysum = constraints constraining outoing slack sums to equal total sums
    ## rphase** = reference edge phase constraints that limit only one pair of r11, r12, r21, r22 to be positive
    ## nphase = aberrant edge phase constraints that limit only one allele pair to be non-negative
    ## oedge1, oedge2 = outgoing allelic edge conservation on each allele
    ## iedge1, iedge2 = incoming allelic edge conservation on each allele
    ## adup1, adup2 = duplicate constraints coupling positive interval copy to reverse complement interval copy
    ## rdup11, rdup12, rdup21, rdup22 = duplicate constraints coupling positive reference edge to reverse complement reference edge copy
    ## adup11, adup12, adup21, adup22 = duplicate constraints coupling positive aberrant edge copy to reverse complement aberrant edge copy
    ## isphase, osphase = phasing for incoming and outgoing slack
    ## i_r**, i_n**, i_is*, i_os* = "big M" constraints instantiating indicator variables 
    consmeta = data.frame(
      cons = c(
        rep(c('copy1', 'copy2'), each = length(pos.ix)),
        rep(c('adup1', 'adup2'), each = length(pos.ix)),
        rep(c('acopysum'), each = length(pos.ix)),
        rep('rcopysum', each = length(ref.pos.ix)),
        rep('ncopysum', each = length(ab.pos.ix)),
        rep(c('rdup11', 'rdup12', 'rdup21', 'rdup22'), each = length(ref.pos.ix)),
        rep(c('ndup11', 'ndup12', 'ndup21', 'ndup22'), each = length(ab.pos.ix)),
        rep(c('oedge1', 'oedge2', 'iedge1', 'iedge2'), each = length(pos.ix)),       
        rep(c('rphase*1', 'rphase*2', 'rphase1*', 'rphase2*'), each = length(ref.pos.ix)),
        rep(c('nphase'), each = length(ab.pos.ix)),
        rep(c('isphase', 'osphase'), each = length(pos.ix)),
        rep(c('Mlb_r11', 'Mlb_r12', 'Mlb_r21', 'Mlb_r22'), each = length(ref.pos.ix)),
        rep(c('Mub_r11', 'Mub_r12', 'Mub_r21', 'Mub_r22'), each = length(ref.pos.ix)),
        rep(c('Mlb_n11', 'Mlb_n12', 'Mlb_n21', 'Mlb_n22'), each = length(ab.pos.ix)),
        rep(c('Mub_n11', 'Mub_n12', 'Mub_n21', 'Mub_n22'), each = length(ab.pos.ix)),
        rep(c('Mlb_is1', 'Mlb_is2', 'Mlb_os1', 'Mlb_os2'), each = length(pos.ix)),
        rep(c('Mub_is1', 'Mub_is2', 'Mub_os1', 'Mub_os2'), each = length(pos.ix))
        ),
      num = c(
        rep(1:length(pos.ix), 2),
        rep(1:length(pos.ix), 2),
        rep(1:length(pos.ix), 1),
        rep(1:length(ref.pos.ix), 1),
        rep(1:length(ab.pos.ix), 1),
        rep(1:length(ref.pos.ix), 4),
        rep(1:length(ab.pos.ix), 4),
        rep(1:length(pos.ix), 4),
        rep(1:length(ref.pos.ix), 4),
        rep(1:length(ab.pos.ix), 1),
        rep(1:length(pos.ix), 2),
        rep(1:length(ref.pos.ix), 4),
        rep(1:length(ref.pos.ix), 4),
        rep(1:length(ab.pos.ix), 4),
        rep(1:length(ab.pos.ix), 4),
        rep(1:length(pos.ix), 4),
        rep(1:length(pos.ix), 4)
        ),
      stringsAsFactors = F
      )
    consmeta$sense = 'E'    
    consmeta$sense[grepl('(phase)|(Mlb)|(Mub)', consmeta$cons)]= 'L'
    
    ## populate constraints
    ##
    n = nrow(varmeta)
    m = nrow(consmeta)
    A = sparseMatrix(1, 1, x = 0, dims = c(m, n)) #
    consmeta$b = rep(NA, length(n))
        
    M = 1e7;

    browser()
    
    ## copy state + eps constraints
    ## a1 = mu1 + eps
    cix = which(consmeta$cons == 'copy1')    
    A[cbind(cix, varmeta[paste('a1', pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('eps1', pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = mu_high[varmeta[paste('eps1', pos.ix[consmeta$num[cix]]), 'parent']]
    
    ## a2 = mu2 + eps
    cix = which(consmeta$cons == 'copy2')    
    A[cbind(cix, varmeta[paste('a2', pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('eps2', pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = mu_low[varmeta[paste('eps2', pos.ix[consmeta$num[cix]]), 'parent']]
    
    ## dup vertex constraints
    ## a1[neg.ix] = a1[dup.vmap[neg.ix]]
    cix = which(consmeta$cons == 'adup1')    
    A[cbind(cix, varmeta[paste('a1', neg.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('a1', dup.vmap[neg.ix[consmeta$num[cix]]]), 'id'])] = -1
    consmeta$b[cix] = 0

    ## a2[neg.ix] = a2[dup.vmap[neg.ix]]
    cix = which(consmeta$cons == 'adup2')    
    A[cbind(cix, varmeta[paste('a2', neg.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('a2', dup.vmap[neg.ix[consmeta$num[cix]]]), 'id'])] = -1
    consmeta$b[cix] = 0

    ## vertex allelic copy sum constraints
    ## a = a1 + a2
    cix = which(consmeta$cons == 'acopysum')    
    A[cbind(cix, varmeta[paste('a1', pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('a2', pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = segstats$cn[varmeta[paste('a1', pos.ix[consmeta$num[cix]]), 'parent']]
    
    ## reference junction copy sum constraints
    ## r = r11 + r12 + r21 + r22
    cix = which(consmeta$cons == 'rcopysum')    
    A[cbind(cix, varmeta[paste('r11', consmeta$num[cix]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r12', consmeta$num[cix]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r21', consmeta$num[cix]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r22', consmeta$num[cix]), 'id'])] = 1
    consmeta$b[cix] = adj[ref.jun[ref.pos.ix[varmeta[paste('r22', consmeta$num[cix]), 'parent']], ]]
    
    ## aberrant junction copy sum constraints
    ## n = n11 + n12 + n21 + n22
    cix = which(consmeta$cons == 'ncopysum')    
    A[cbind(cix, varmeta[paste('n11', consmeta$num[cix]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('n12', consmeta$num[cix]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('n21', consmeta$num[cix]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('n22', consmeta$num[cix]), 'id'])] = 1
    consmeta$b[cix] = adj[ab.jun[ab.pos.ix[varmeta[paste('n22', consmeta$num[cix]), 'parent']], ]]

    ## reference junction dup constraints
    ##
    cix = which(consmeta$cons == 'rdup11')    
    A[cbind(cix, varmeta[paste('r11', ref.neg.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r11', dup.ref.emap[ref.neg.ix[consmeta$num[cix]]]), 'id'])] = 1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'rdup12')    
    A[cbind(cix, varmeta[paste('r12', ref.neg.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r12', dup.ref.emap[ref.neg.ix[consmeta$num[cix]]]), 'id'])] = 1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'rdup21')    
    A[cbind(cix, varmeta[paste('r21', ref.neg.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r21', dup.ref.emap[ref.neg.ix[consmeta$num[cix]]]), 'id'])] = 1
    consmeta$b[cix] = 0
    
    cix = which(consmeta$cons == 'rdup22')    
    A[cbind(cix, varmeta[paste('r22', ref.neg.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r22', dup.ref.emap[ref.neg.ix[consmeta$num[cix]]]), 'id'])] = 1
    consmeta$b[cix] = 0

    ## aberrant junction dup constraints
    ##
    cix = which(consmeta$cons == 'rdup11')    
    A[cbind(cix, varmeta[paste('r11', ab.neg.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r11', dup.ab.emap[ab.neg.ix[consmeta$num[cix]]]), 'id'])] = 1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'rdup12')    
    A[cbind(cix, varmeta[paste('r12', ab.neg.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r12', dup.ab.emap[ab.neg.ix[consmeta$num[cix]]]), 'id'])] = 1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'rdup21')    
    A[cbind(cix, varmeta[paste('r21', ab.neg.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r21', dup.ab.emap[ab.neg.ix[consmeta$num[cix]]]), 'id'])] = 1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'rdup22')    
    A[cbind(cix, varmeta[paste('r22', ab.neg.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r22', dup.ab.emap[ab.neg.ix[consmeta$num[cix]]]), 'id'])] = 1
    consmeta$b[cix] = 0
    
    ## outgoing allelic edge constraints
    ##
    ## a1 = r11 + r12 + sum_k {n11}_k + sum_k {n12}_k + os1
    ##
    cix = which(consmeta$cons == 'oedge1')
    rj.ix = lapply(pos.ix[consmeta$num[cix]], function(x) which(ref.jun[,1] %in% x))
    aj.ix = lapply(pos.ix[consmeta$num[cix]], function(x) which(ab.jun[,1] %in% x))
    A[cbind(cix, varmeta[paste('a1', pos.ix[consmeta$num[cix]]), 'id'])] = -1
    for (i in 1:length(rj.ix))
      if (length(rj.ix[[i]])>0)
        {
          A[cbind(cix[i], varmeta[paste('r11', rj.ix[[i]]), 'id'])] = 1
          A[cbind(cix[i], varmeta[paste('r12', rj.ix[[i]]), 'id'])] = 1
        }
    
    for (i in length(aj.ix))
      if (length(aj.ix[[i]])>0)
        {
          A[cbind(cix[i], varmeta[paste('n11', aj.ix[[i]]), 'id'])] = 1
          A[cbind(cix[i], varmeta[paste('n12', aj.ix[[i]]), 'id'])] = 1
        }    
    
    ## a2 = r21 + r22 + sum_k {n21}_k + sum_k {n22}_k + os2
    cix = which(consmeta$cons == 'oedge2')
    rj.ix = lapply(pos.ix[consmeta$num[cix]], function(x) which(ref.jun[,1] %in% x))
    aj.ix = lapply(pos.ix[consmeta$num[cix]], function(x) which(ab.jun[,1] %in% x))
    A[cbind(cix, varmeta[paste('a2', pos.ix[consmeta$num[cix]]), 'id'])] = -1
    for (i in 1:length(rj.ix))
      if (length(rj.ix[[i]])>0)
        {
          A[cbind(cix[i], varmeta[paste('r21', rj.ix[[i]]), 'id'])] = 1
          A[cbind(cix[i], varmeta[paste('r22', rj.ix[[i]]), 'id'])] = 1
        }    
    for (i in length(aj.ix))
      if (length(aj.ix[[i]])>0)
        {
          A[cbind(cix[i], varmeta[paste('n21', aj.ix[[i]]), 'id'])] = 1
          A[cbind(cix[i], varmeta[paste('n22', aj.ix[[i]]), 'id'])] = 1
        }
        
    ## incoming allelic edge constraints
    ##
    ## a1 = r11 + r21 + sum_k {n11}_k + sum_k {n21}_k + is1
    ##
    cix = which(consmeta$cons == 'iedge1')
    rj.ix = lapply(pos.ix[consmeta$num[cix]], function(x) which(ref.jun[,2] %in% x))
    aj.ix = lapply(pos.ix[consmeta$num[cix]], function(x) which(ab.jun[,2] %in% x))
    A[cbind(cix, varmeta[paste('a1', pos.ix[consmeta$num[cix]]), 'id'])] = -1
    for (i in 1:length(rj.ix))
      if (length(rj.ix[[i]])>0)
        {
          A[cbind(cix[i], varmeta[paste('r11', rj.ix[[i]]), 'id'])] = 1
          A[cbind(cix[i], varmeta[paste('r21', rj.ix[[i]]), 'id'])] = 1
        }
    
    for (i in length(aj.ix))
      if (length(aj.ix[[i]])>0)
        {
          A[cbind(cix[i], varmeta[paste('n11', aj.ix[[i]]), 'id'])] = 1
          A[cbind(cix[i], varmeta[paste('n21', aj.ix[[i]]), 'id'])] = 1
        }    
    
    ## a2 = r12 + r22 + sum_k {n12}_k + sum_k {n22}_k + is2
    cix = which(consmeta$cons == 'oedge2')
    rj.ix = lapply(pos.ix[consmeta$num[cix]], function(x) which(ref.jun[,2] %in% x))
    aj.ix = lapply(pos.ix[consmeta$num[cix]], function(x) which(ab.jun[,2] %in% x))
    A[cbind(cix, varmeta[paste('a2', pos.ix[consmeta$num[cix]]), 'id'])] = -1
    for (i in 1:length(rj.ix))
      if (length(rj.ix[[i]])>0)
        {
          A[cbind(cix[i], varmeta[paste('r12', rj.ix[[i]]), 'id'])] = 1
          A[cbind(cix[i], varmeta[paste('r22', rj.ix[[i]]), 'id'])] = 1
        }    
    for (i in length(aj.ix))
      if (length(aj.ix[[i]])>0)
        {
          A[cbind(cix[i], varmeta[paste('n12', aj.ix[[i]]), 'id'])] = 1
          A[cbind(cix[i], varmeta[paste('n22', aj.ix[[i]]), 'id'])] = 1
        }

    ## reference phase constraints
    ##

    ## i_r11 + i_r12 <=1
    cix = which(consmeta$cons == 'rphase*1')
    A[cbind(cix, varmeta[paste('i_r11', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('i_r21', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 1

    ## i_r12 + i_r22 <=1
    cix = which(consmeta$cons == 'rphase*2')
    A[cbind(cix, varmeta[paste('i_r12', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('i_r22', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 1

    ## i_r11 + i_r12 <=1
    cix = which(consmeta$cons == 'rphase1*')
    A[cbind(cix, varmeta[paste('i_r11', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('i_r12', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 1

    ## i_r21 + i_r22 <=1
    cix = which(consmeta$cons == 'rphase2*')
    A[cbind(cix, varmeta[paste('i_r21', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('i_r22', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 1
    
    ## aberrant phase constraints
    ## i_n11 + i_n12 + i_n21 + i_n22 <= 1
    cix = which(consmeta$cons == 'nphase')
    A[cbind(cix, varmeta[paste('i_n11', ab.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('i_n12', ab.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('i_n21', ab.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('i_n22', ab.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 1

    ## incoming slack constraints
    cix = which(consmeta$cons == 'isphase')    
    A[cbind(cix, varmeta[paste('i_is1', pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('i_is2', pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 1

    ## outgoing slack constraints
    cix = which(consmeta$cons == 'osphase')    
    A[cbind(cix, varmeta[paste('i_os1', pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('i_os2', pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 1

    ## "Big M" lower and upper bound constraints
    ##
    ## -M*i_r11 <= r11 <= M*i_r11
    cix = which(consmeta$cons == 'Mlb_r11')
    A[cbind(cix, varmeta[paste('i_r11', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('r11', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'Mub_r11')
    A[cbind(cix, varmeta[paste('i_r11', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('r11', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 0

    ## -M*i_r12 <= r12 <= M*i_r12
    cix = which(consmeta$cons == 'Mlb_r12')
    A[cbind(cix, varmeta[paste('i_r12', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('r12', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'Mub_r12')
    A[cbind(cix, varmeta[paste('i_r12', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('r12', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1

    ## -M*i_r21 <= r21 <= M*i_r21
    cix = which(consmeta$cons == 'Mlb_r21')
    A[cbind(cix, varmeta[paste('i_r21', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('r21', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'Mub_r21')
    A[cbind(cix, varmeta[paste('i_r21', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('r21', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1

    ## -M*i_r22 <= r22 <= M*i_r22
    cix = which(consmeta$cons == 'Mlb_r22')
    A[cbind(cix, varmeta[paste('i_r22', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('r22', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    sconsmeta$b[cix] = 0
    
    cix = which(consmeta$cons == 'Mub_r22')
    A[cbind(cix, varmeta[paste('i_r22', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('r22', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1

    ## -M*i_n11 <= n11 <= M*i_n11
    cix = which(consmeta$cons == 'Mlb_n11')
    A[cbind(cix, varmeta[paste('i_n11', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('n11', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'Mub_n11')
    A[cbind(cix, varmeta[paste('i_n11', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('n11', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 0

    ## -M*i_n12 <= n12 <= M*i_n12
    cix = which(consmeta$cons == 'Mlb_n12')
    A[cbind(cix, varmeta[paste('i_n12', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('n12', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'Mub_n12')
    A[cbind(cix, varmeta[paste('i_n12', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('n12', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1

    ## -M*i_n21 <= n21 <= M*i_n21
    cix = which(consmeta$cons == 'Mlb_n21')
    A[cbind(cix, varmeta[paste('i_n21', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('n21', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'Mub_n21')
    A[cbind(cix, varmeta[paste('i_n21', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('n21', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1

    ## -M*i_n22 <= n22 <= M*i_n22
    cix = which(consmeta$cons == 'Mlb_n22')
    A[cbind(cix, varmeta[paste('i_n22', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('n22', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    sconsmeta$b[cix] = 0
    
    cix = which(consmeta$cons == 'Mub_n22')
    A[cbind(cix, varmeta[paste('i_n22', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('n22', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    
    ## -M*i_is1 <= is1 <= M*i_is1
    cix = which(consmeta$cons == 'Mlb_is1')
    A[cbind(cix, varmeta[paste('i_is1', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('is1', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'Mub_is1')
    A[cbind(cix, varmeta[paste('i_is1', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('is1', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 0

    ## -M*i_is2 <= is2 <= M*i_is2
    cix = which(consmeta$cons == 'Mlb_is2')
    A[cbind(cix, varmeta[paste('i_is2', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('is2', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'Mub_is2')
    A[cbind(cix, varmeta[paste('i_is2', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('is2', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 0

    ## -M*i_os1 <= os1 <= M*i_os1
    cix = which(consmeta$cons == 'Mlb_os1')
    A[cbind(cix, varmeta[paste('i_os1', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('os1', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'Mub_os1')
    A[cbind(cix, varmeta[paste('i_os1', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('os1', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 0

    ## -M*i_os2 <= os2 <= M*i_os2
    cix = which(consmeta$cons == 'Mlb_os2')
    A[cbind(cix, varmeta[paste('i_os2', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('os2', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    consmeta$b[cix] = 0
    
    cix = which(consmeta$cons == 'Mub_os2')
    A[cbind(cix, varmeta[paste('i_os2', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('os2', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 0
    
  }

####################
#' jabba.alleles
#'
#' Populates allelic value s for JaBbA object.  This does not explicitly impose junction balance constraints on alleles, but rather just computes
#' the maximum likelihood estimate given allelic counts and the inferred total copy number on a given segment according to JaBbA
#'
#' @param jab JaBbA object 
#' @param het.sites GRanges with meta data fields (see below) for alt and rref count
#' @param alt.count.field character specifying alt.count meta data field in input het.sites (default alt.count.t)
#' @param ref.count.field character specifying ref.count meta data field in input het.sites (default ref.count.t)
#' @param split.ab logical flag whether to split aberrant segmetns (segmentss with ab edge entering or leaving prior to computing allelic states (default FALSE) 
#' @param uncoupled logical flag whether to not collapse segments after inferring MLE estimate (default FALSE), if FALSE will try to merge adjacent segments and populate allele-specific junctions with copy numbers on the basis of the MLE fit on individual allelic segments
#' @param conservative if TRUE then will leave certain allelic segments "unphased" if one cannot sync the high / low interval state with the incoming and / or outgoing junction state
#' @param verbose logical flag
#' @return
#' list with following fields:
#' $segstats = GRanges of input segments with $cn.high and $cn.low segments populated
#' $asegstats = GRanges of allelic segments (length is 2*length(segstats)) with high and low segments each having $cn, this is a "melted" segstats GRAnges
#' $atd = gTrack of allelic segments and supporting input het.sites
#' $aadj = allelic adjacency matrix of allele specific junctions
#' $ab.ix = indices of aberrant edges in $aadj
#' $ref.ix = indices of reference edges in $aadj
#' @export
############################################
jabba.alleles = function(
    jab,
    het.sites, ## granges with meta data fields for alt.count and 
    alt.count.field = 'alt.count.t',
    ref.count.field = 'ref.count.t',
    baf.field = 'baf.t',
    split.ab = F, ## if split.ab == T, then will split across any "aberrant" segment (i.e. segment with ab edge entering or leaving prior to computing allelic states (note: this might create gaps)
    uncoupled = FALSE, ## if uncoupled, we just assign each high low allele the MLE conditioning on the total copy number
    conservative = FALSE, ## if TRUE then will leave certain allelic segments "unphased" if one cannot sync the high / low interval state with the incoming and / or outgoing junction state
    verbose = F  
  )
  {
      if (!all(c(alt.count.field, ref.count.field) %in% names(values(het.sites)))){
          cat('count fields not found in meta data of het.sites input, trying BAF...')
          if (!(baf.field %in% names(values(het.sites))))
              stop('BAF field not found in meta data of het.sites input either!')
          else{
              ## TODO: change the stats model to beta
              ## outputs are re.seg$low and re.seg$high
              ## test deviations of observed BAF from expected by beta distribution
              if (verbose)
                  cat('Processing', length(het.sites),
                      'het sites using fields', baf.field, '\n')

############ below are copied from old poisson model
              
    ## ## now test deviation from each absolute copy combo using poisson model
    ## ## i.e. counts ~ poisson(expected mean)
    ## ##
    ## re.seg$low = sapply(1:length(re.seg), function(i)
    ##   {
    ##     if (verbose)
    ##       cat('.')    
    ##     x = lows[[i]]
    ##     if (length(x)==0)
    ##       return(NA)
    ##     y = highs[[i]]
    ##     tot.cn = cn[i]
    ##     ll = sapply(0:(floor(tot.cn/2)), function(j) sum(ppois(x,centers[j+1], log.p = T) + ppois(y,centers[tot.cn-j+1],log.p = T)))
    ##     ll = ll - min(ll)
    ##     return(which.max(ll)-1)
    ##   })

    ## if (verbose)
    ##   cat('\n')

    ## re.seg$high = re.seg$cn-re.seg$low

              
          }
      } else {
          ## stop('count fields not found in meta data of het.sites input')

    if (verbose)
      cat('Processing', length(het.sites), 'het sites using fields', alt.count.field, 'and', ref.count.field, '\n')
    
    het.sites$low.count = pmin(values(het.sites)[, alt.count.field], values(het.sites)[, ref.count.field])
    het.sites$high.count = pmax(values(het.sites)[, alt.count.field], values(het.sites)[, ref.count.field])

    het.sites = het.sites[!is.na(het.sites$low.count) & !is.na(het.sites$high.count)]
    
    ss.p = jab$segstats[ as.logical( strand(jab$segstats)=='+' ) ] 

    ## find the reference junctions    
    ord.ix = order(jab$segstats)
    rev.ix = as.logical(strand(jab$segstats[ord.ix]) == '-')
    ord.ix = c(ord.ix[!rev.ix], rev(ord.ix[rev.ix]))

    ref.jun = cbind(ord.ix[-length(ord.ix)], ord.ix[-1])
    ref.jun = ref.jun[which(jab$adj[ref.jun]>0), ]

    has.ab.rand = 0
    if (split.ab)
      {
        ab.adj = jab$adj
        ab.adj[ref.jun] = 0        
        has.ab = as.numeric(Matrix::rowSums(ab.adj!=0)!=0 | Matrix::colSums(ab.adj!=0)!=0)[which( as.logical( strand(jab$segstats)=='+')) ]
        has.ab.rand = runif(length(ss.p)) * 1e-6 * has.ab
      }

    ss.p = ss.p[!is.na(ss.p$cn)]
    re.seg = as(coverage(ss.p, weight = ss.p$cn + has.ab.rand), 'GRanges')
    re.seg$cn = round(re.seg$score)
      
    het.sites$ix = gr.match(het.sites, re.seg)

    if (verbose)
      cat('Computed high / low counts and matched to segs\n')
   
    highs = split(het.sites$high.count, het.sites$ix)[as.character(1:length(re.seg))]
    lows = split(het.sites$low.count, het.sites$ix)[as.character(1:length(re.seg))]
    
    het.sites$cn = re.seg$cn[het.sites$ix]
    purity = jab$purity
    ploidy = mean(het.sites$cn, na.rm = T) ## ploidy may be slightly different from "global ploidy" depending on the distribution of sites
    
    sw = length(het.sites)
    total = sum(as.numeric(c(het.sites$high.count, het.sites$low.count)))

    cn = re.seg$cn
    ## gamma = 2*(1-purity)/purity  ## gammas and betas need to be recomputed for 
    ## beta = (2*(1-purity)*sw + purity*ploidy*sw) / (purity * total)
    gamma = 1*(1-purity)/purity  ## gammas and betas need to be recomputed for  (1 since we are looking at het alleles)
    beta = (1*(1-purity)*sw + purity*ploidy*sw) / (purity * total)
    centers = (0:(max(cn)) + gamma)/beta

    if (verbose)
      cat('Computed SNP ploidy and allelic copy centers\n')
    
    ## now test deviation from each absolute copy combo using poisson model
    ## i.e. counts ~ poisson(expected mean)
    ##
    re.seg$low = sapply(1:length(re.seg), function(i)
      {
        if (verbose)
          cat('.')    
        x = lows[[i]]
        if (length(x)==0)
          return(NA)
        y = highs[[i]]
        tot.cn = cn[i]
        ll = sapply(0:(floor(tot.cn/2)), function(j) sum(ppois(x,centers[j+1], log.p = T) + ppois(y,centers[tot.cn-j+1],log.p = T)))
        ll = ll - min(ll)
        return(which.max(ll)-1)
      })

    if (verbose)
      cat('\n')

    re.seg$high = re.seg$cn-re.seg$low
          }
###########################################################################
    ## borderline, below are common to both methods    
    jab$segstats$cn.low = round(gr.val(jab$segstats, re.seg, 'low', na.rm = TRUE)$low)
    jab$segstats$cn.high = round(gr.val(jab$segstats, re.seg, 'high', na.rm = TRUE)$high)
    na.ix = !gr.val(jab$segstats, re.seg, 'low', FUN = function(x,w,na.rm) any(!is.na(x)))$low | !gr.val(jab$segstats, re.seg, 'high', FUN = function(x,w,na.rm) any(!is.na(x)))$high
    jab$segstats$cn.low[na.ix] = jab$segstats$cn.high[na.ix] = NA
    
    #############
    # phasing
    #
    #############
    
    ## iterate through all reference junctions and apply (wishful thinking) heuristic
    ##
    ## populate n x n x 2 adjacency matrix, which we will later expand to a bigger matrix
    
    adj.ab = jab$adj    
    adj.ab[ref.jun] = 0
    adj.ref = jab$adj*0
    adj.ref[ref.jun] = jab$adj[ref.jun]
    high = low = jab$segstats[, c()]
    high$cn = jab$segstats$cn.high
    low$cn = jab$segstats$cn.low
    high$parent = low$parent = 1:length(jab$segstats)
    high$type = 'high'
    low$type = 'low'
    high$id = 1:length(jab$segstats)
    low$id = length(jab$segstats) + 1:length(jab$segstats)
    asegstats = c(high, low)
    amap = cbind(high$id, low$id) ## maps segstats id x allele combos to asegstats id
    
    aadj = sparseMatrix(1, 1, x = 0, dims = c(length(asegstats), length(asegstats)))
    
    .flip = function(x) x %% 2+1

    asegstats = c(high, low)
    acn = cbind(high$cn, low$cn)

    phased.out = phased.in = rep(TRUE, length(asegstats))

    str = strand(asegstats)
    
    if (verbose)
      cat('Starting phasing \n')
    
    for (k in 1:nrow(ref.jun))
      {
        if (verbose)
          cat('.')
        i = ref.jun[k, 1]
        j = ref.jun[k, 2]        
        a = acn[ref.jun[k,1],]
        b = acn[ref.jun[k,2],]

        phased.out[amap[i, ]] = FALSE
        phased.in[amap[j, ]] = FALSE
        
        pairs.ij = cbind(rep(c(1:2), 2), rep(c(1:2), each = 2))
        m = setdiff(which(a[pairs.ij[,1]] == b[pairs.ij[,2]]), NA)
        
##         cat(any(c(amap[i,], amap[j, ]) %in% wtf), '\n')
        
##         if (any(c(amap[i,], amap[j, ]) %in% wtf))
##              browser()
        
        if (!(length(m) %in% c(0, 4))) ## 1,2, and 3 matches are fine (3 matches occur if one interval is in allelic balance, and the other not
          {
            if (length(m)==2) ## pick the phase that the alleles can handle
              m = rev(m[order(as.numeric(sum(adj.ab[i, ])<=a[pairs.ij[m,1]]) + as.numeric(sum(adj.ab[, j])<=b[pairs.ij[m,2]]))])

            m.ij = pairs.ij[m[1], ]
            fm.ij = .flip(m.ij)
            aadj[amap[i, m.ij[1]], amap[j, m.ij[2]]] = min(a[m.ij[1]], jab$adj[i, j])
            aadj[amap[i, fm.ij[1]], amap[j, fm.ij[2]]] = jab$adj[i, j] - aadj[amap[i, m.ij[1]], amap[j, m.ij[2]]]

            phased.out[amap[i, ]] = TRUE
            phased.in[amap[j, ]] = TRUE
            
            if (length(a.ab <- which(adj.ab[i,]!=0))>0)
              {
                ## if a.ab (partner) is already phased then unpopulate the non-ab allelic junction, otherwise populate both alleles of partner
                if (any(ph <- aadj[amap[i, fm.ij[1]], amap[a.ab, ]] !=0))
                  {
                    aadj[amap[i, fm.ij[1]], amap[a.ab, ph]] = adj.ab[i, a.ab]
                    aadj[amap[i, m.ij[1]], amap[a.ab, ph]] = 0
                  }
                else
                  ## otherwise diffuse copy into both alleles of the partner (will be resolved when we resolve phase for the partner interval)
                  ## or collapse unphased nodes back
                  aadj[amap[i, fm.ij[1]], amap[a.ab, ]] = adj.ab[i, a.ab]/2

                if (!conservative)
                  if (a[fm.ij[1]] < adj.ab[i, a.ab]) # if the allelic node can't handle the outgoing allelic edge flux, so unphase
                    phased.out[amap[i, ]] = FALSE
              }
            
            if (length(b.ab <- which(adj.ab[,j]!=0))>0)
              {
                ## if b.ab (partner) is already phased then concentrate all of the junction copy into the aberrant allele of this interval
                if (any(ph <- aadj[amap[b.ab, ], amap[j, fm.ij[2]]] !=0))
                  {
                    aadj[amap[b.ab, ph], amap[j, fm.ij[2]]] = adj.ab[b.ab, j]
                    aadj[amap[b.ab, ph], amap[j, m.ij[2]]] = 0
                  }
                else
                  ## otherwise diffuse copy into both alleles of the partner (will be resolved when we resolve phase for the partner interval)
                  ## or collapse unphased nodes back
                  aadj[amap[b.ab,], amap[j, fm.ij[2]]] = adj.ab[b.ab, j]/2
                
                if (!conservative)
                  if (b[fm.ij[2]] < adj.ab[b.ab, j]) # the allelic node cn can't handle the incoming allelic edge flux, so unphase
                    phased.in[amap[j, ]] = FALSE
              }
          }
      }

    if (verbose)
      cat('\nFinished phasing, finalizing \n')
    
    asegstats$phased.in = phased.in    
    asegstats$phased.out = phased.out

    if (uncoupled)
        unphased = rep(FALSE, length(asegstats))
    else
        unphased = !asegstats$phased.out | !asegstats$phased.in
            
    unphased.parents = unique(asegstats$parent[unphased])
    aadj.unphunph = jab$adj[unphased.parents, unphased.parents]
    aadj.phph = aadj[!unphased, !unphased]

    asegstats$new.ind = NA
    asegstats$new.ind[!unphased] = 1:sum(!unphased)
    asegstats$new.ind[unphased] = as.integer(factor(asegstats$parent[unphased], unphased.parents))
    mat.collapse = sparseMatrix(which(unphased), asegstats$new.ind[unphased], x = 1, dims = c(nrow(aadj), length(unphased.parents)))
    
    aadj.phunph = aadj[!unphased, ] %*% mat.collapse    
    aadj.unphph = t(mat.collapse) %*% aadj[, which(!unphased)]

    aadj.final = rBind(
      cBind(aadj.phph, aadj.phunph),
      cBind(aadj.unphph, aadj.unphunph)
      )
    
    asegstats.unphased = asegstats[match(unphased.parents, asegstats$parent)]
    asegstats.unphased$cn = jab$segstats$cn[asegstats.unphased$parent]
    asegstats.final = c(asegstats[!unphased], asegstats.unphased)
    asegstats.final$phased = c(rep(T, sum(!unphased)), rep(F, length(asegstats.unphased)))
    asegstats.final$type[!asegstats.final$phased] = 'total'

    tmp.str = gr.string(gr.stripstrand(asegstats), mb = F, other.cols = 'type');
    asegstats$tile.id = as.integer(factor(tmp.str, unique(tmp.str)))
    
    ix = order(asegstats.final)
    asegstats.final = asegstats.final[ix]
    aadj.final = aadj.final[ix, ix]

    if (verbose)
      cat('Annotating allelic vertices\n')

    tmp.string = gr.string(asegstats, mb = F, other.col = 'type'); tmp.string2 = gr.string(gr.flipstrand(asegstats), mb = F, other.col = 'type')
    asegstats$flip.ix = match(tmp.string, tmp.string2)
    asegstats$phased = !unphased
        
    asegstats.final$edges.in = sapply(1:length(asegstats.final),
      function(x) {ix = which(aadj.final[,x]!=0); paste(ix, '(', aadj.final[ix,x], ')', '->', sep = '', collapse = ',')})
    asegstats.final$edges.out = sapply(1:length(asegstats.final),
      function(x) {ix = which(aadj.final[x, ]!=0); paste('->', ix, '(', aadj.final[x,ix], ')', sep = '', collapse = ',')})

    asegstats.final$slack.in = asegstats.final$cn - Matrix::colSums(aadj.final)
    asegstats.final$slack.out = asegstats.final$cn - Matrix::rowSums(aadj.final)

    asegstats.final$new.ind = asegstats.final$phased.out = asegstats.final$phased.in = asegstats.final$id = NULL
    asegstats.final$tile.id = as.integer(factor(gr.string(gr.stripstrand(asegstats.final), mb = F, other.cols = 'type')))
    
    m = sparseMatrix(1:length(asegstats.final), asegstats.final$parent, x = 1);

    hh = rep(het.sites[, c()], 2)
    hh$count = c(het.sites$high.count, het.sites$low.count)
    hh$type = rep(c('high', 'low'), each = length(het.sites))

    hh$ywid = 0.5  
    atd = c(
        gTrack(hh, angle = 0, y.field = 'count', y0 = 0, 
                  colormap = list(type = c('high' = alpha('red', 0.3), 'low' = alpha('blue', 0.3))), name = 'hets', y.quantile = 0.001, lwd.border = 2),
        gTrack(asegstats.final, angle = 0, y.field = 'cn', y0 = 0,
              colormap = list(type = c('high' = alpha('red', 0.3), 'low' = alpha('blue', 0.3), 'total' = alpha('purple', 0.3))), name = 'alleles')
        )
      
    out = list(
        segstats = jab$segstats,
        asegstats = asegstats.final,
        atd = atd, 
        aadj = aadj.final,
        ab.ix = which((m %*% adj.ab %*% t(m))!=0, arr.ind = T),
        ref.ix = which((m %*% adj.ref %*% t(m))!=0, arr.ind = T)
      )

    return(out)
}



#################################################
#' junction.paths
#' 
#' Applies "pigeonhole principle" to enumerate all junction paths
#' in a karyograph that can be proven to have copy number greater than 0
#' 
#' Takes as input adjacency matrix specifying junction copy numbers
#' and numeric vector specifying node copy numbers. 
#'
#' @param cn length n vector of integer copy numbers
#' @param adj nxn matrix of junction copy numbers
#' @return  list
#' with fields
#' $paths  list of n paths, each path i consisting of an n_i x 2 matrix specifying sequences of n_i junctions (each an ij node pair)
#' $mcn minimum copy number associated with path i
#' @export
#################################################
junction.paths = function(cn, adj)
  {
    ## preallocate, preallocate, preallocate
    ed = which(adj!=0)
    NMAX = length(cn)*3 ## should be larger than the number of anticipated paths
    EMAX = 1000
    BOOSTER.ROW = 1e4
    BOOSTER.COL = 500    
    paths = array(NA, dim = c(NMAX, EMAX))
    firstnode = lastnode = lastix = mcn = rep(NA, NMAX)
    numpaths = 0
    numrows = nrow(paths) ## yes this is a very slow query for giant arrays
    numcols = ncol(paths)
    torem = rep(F, NMAX)
        
    .sub2ind = function(dim, r, c) (c-1)*dim[1] + r
          
    if (nrow(adj) != ncol(adj) | nrow(adj)!=length(cn))
      stop('Adjacency matrix must be n x n and cn must be a length n vector')
    
    for (i in which(!is.na(cn)))
      {
        outgoing.nodes = which(adj[i, ]>0)        
        incoming.nodes = which(adj[, i]>0)        
        outgoing.cn = adj[i, outgoing.nodes]
        incoming.cn = adj[incoming.nodes, i]
        outgoing.edges = .sub2ind(dim(adj), i, outgoing.nodes)
        incoming.edges = .sub2ind(dim(adj), incoming.nodes, i)
        
        ## augment existing paths and adjust their minimum cn
        if (length(outgoing.edges)>0)
          {
            if (numpaths>0)
              {
                ix =  which(lastnode[1:numpaths] == i)
                if( length(ix)>0 )
                  for (j in ix)
                    {
                      new.lastix = rep(lastix[j]+1, length(outgoing.edges))
                      new.paths = paths[rep(j, length(outgoing.edges)), 1:new.lastix[1], drop = F]
                      new.paths[, new.lastix] = outgoing.edges                  
                      new.mcn =  mcn[j]-(cn[i]-outgoing.cn)
                      new.paths = new.paths[new.mcn>0,, drop = F]                  
                      new.lastix = new.lastix[new.mcn>0]
                      new.firstnode = rep(firstnode[j], sum(new.mcn>0))
                      new.lastnode = outgoing.nodes[new.mcn>0]
                      new.mcn = new.mcn[new.mcn>0]

                      if (any(is.na(new.lastix)))                        
                        browser()
                      
                      if (length(new.mcn)>0)
                        {
                          if (any(mcn[j]==new.mcn)) ## only throw away an old path if the mcn of the new path did not decay
                            {
                              torem[j] = T
                              firstnode[j] = lastnode[j] = lastix[j] = mcn[j] = NA
                            }

                          new.ix = numpaths + (1:nrow(new.paths))
                          paths[new.ix, 1:new.lastix[1]] = new.paths[, 1:new.lastix[1]]
                          mcn[new.ix] = new.mcn
                          lastix[new.ix] = new.lastix
                          lastnode[new.ix] = new.lastnode
                          firstnode[new.ix] = new.firstnode
                          numpaths = numpaths + nrow(new.paths)
                        }
                    }
              }
            
            ## add brand new paths
            new.ix = numpaths + (1:length(outgoing.edges))
            paths[new.ix, 1] = outgoing.edges
            mcn[new.ix] = outgoing.cn
            firstnode[new.ix] = i
            lastnode[new.ix] = outgoing.nodes
            lastix[new.ix] = 1
            numpaths = numpaths + length(outgoing.edges)
          }

        ## augment existing paths and adjust their minimum cn
        if (length(incoming.edges)>0)
          {
            if (numpaths>0)
              {
                ix =  which(firstnode[1:numpaths] == i)
                if( length(ix)>0 )
                  for (j in ix)
                    {
                      now = Sys.time()
                      new.paths = cbind(0, paths[rep(j, length(incoming.edges)), 1:lastix[j], drop = F])
                      new.lastix = rep(lastix[j]+1, length(incoming.edges))
                      new.paths[, 1] = incoming.edges
                      new.mcn =  mcn[j]-(cn[i]-incoming.cn)
                      new.paths = new.paths[new.mcn>0,, drop = F]                  
                      new.lastix = new.lastix[new.mcn>0]
                      new.firstnode = incoming.nodes[new.mcn>0]
                      new.lastnode = rep(lastnode[j], sum(new.mcn>0))
                      new.mcn = new.mcn[new.mcn>0]

                      if (any(is.na(new.lastix)))                        
                        browser()
                      
                      if (length(new.mcn)>0)
                        {                      
                          if (any(mcn[j]==new.mcn)) ## only throw away an old path if the mcn of a new path did not decay
                            {
                              torem[j] = T
                              firstnode[j] = lastnode[j] = lastix[j] = mcn[j] = NA
                            }

                          new.ix = numpaths + (1:nrow(new.paths))
                          paths[new.ix, 1:new.lastix[1]] = new.paths[, 1:new.lastix[1]]
                          mcn[new.ix] = new.mcn
                          lastix[new.ix] = new.lastix
                          lastnode[new.ix] = new.lastnode
                          firstnode[new.ix] = new.firstnode
                          numpaths = numpaths + nrow(new.paths)
                        }
                    }
              }
            
            ## add brand new paths
            new.ix = numpaths + (1:length(incoming.edges))
            paths[new.ix, 1] = incoming.edges
            mcn[new.ix] = incoming.cn
            firstnode[new.ix] = incoming.nodes
            lastnode[new.ix] = i
            lastix[new.ix] = 1
            numpaths = numpaths + length(incoming.edges)            
          }
        
        ## augment if necessary
        if (numpaths>(numrows-BOOSTER.ROW))
          {
            cat('Allocating more row space\n')
            paths = rbind(paths, array(NA, dim = c(BOOSTER.ROW, ncol(paths))))
            numrows = numrows + BOOSTER.ROW
            torem = c(torem, rep(NA, BOOSTER.ROW))
            firstnode = c(firstnode, rep(NA, BOOSTER.ROW))
            lastnode = c(lastnode, rep(NA, BOOSTER.ROW))
            lastix = c(lastix, rep(NA, BOOSTER.ROW))
            mcn = c(mcn, rep(NA, BOOSTER.ROW))             
          }
        
        if (max(c(0, lastix[1:numpaths]), na.rm = T)>(numcols-BOOSTER.COL))
          {
            cat('Allocating more column space\n')
            paths = cbind(paths, array(NA, dim = c(nrow(paths), BOOSTER.COL)))
            numcols = numcols + BOOSTER.COL
          }
        
        if ((i %% 500)==0)
          {
            cat(i, numpaths, '\n')            
#            print(table(rowSums(!is.na(paths[1:numpaths,]))))
            cat('all last nodes\n')
            print(sort(table(seqnames(this.asol$asegstats[setdiff(lastnode, NA)]))))

            cat('all traversed nodes\n')
            print(sort(table(seqnames(this.asol$asegstats[1:i]))))
            rc = ind2sub(dim(adj), setdiff(as.vector(paths[1:numpaths, 1:max(lastix, na.rm = T)]), NA))

            cat('all path nodes\n')
            print(sort(table(seqnames(this.asol$asegstats[unique(as.numeric(rc))]))))

            cat('adj\n')
            print(sort(table(adj[rc])))

            keep.ix = which(!torem[1:numpaths])
            tmp.lastix = lastix[keep.ix]
            tmp.paths = paths[keep.ix, 1:max(tmp.lastix), drop = F]
            tmp.firstnode = firstnode[keep.ix]
            tmp.lastnode = lastnode[keep.ix]
            tmp.mcn = mcn[keep.ix]
            tmp.numpaths = length(keep.ix)
              
            paths[1:numpaths, 1:max(tmp.lastix)] = NA
            lastix[1:numpaths] = mcn[1:numpaths] = firstnode[1:numpaths] = lastnode[1:numpaths] = NA
            torem[1:numpaths] = F

            paths[1:tmp.numpaths, 1:max(tmp.lastix)] = tmp.paths[, 1:max(tmp.lastix), drop = F]
            firstnode[1:tmp.numpaths] = tmp.firstnode
            lastnode[1:tmp.numpaths] = tmp.lastnode
            mcn[1:tmp.numpaths] = tmp.mcn
            lastix[1:tmp.numpaths] = tmp.lastix
            numpaths = tmp.numpaths
 #           print(lapply(order(-rowSums(!is.na(paths[1:numpaths, ])))[1:2], function(x) paths[x, !is.na(paths[x, ])]))

            saveRDS(list(paths = paths[1:numpaths,], mcn = mcn[1:numpaths]), 'paths.rds')
            cat(i, numpaths, '\n')            
          }
      }

    keep = which(rowSums(!is.na(paths)) !=0 & !torem)

    paths = paths[keep, 1:max(lastix, na.rm = T)]
    paths = lapply(1:nrow(paths), function(x) paths[x, !is.na(paths[x, ])])
    mcn = mcn[keep]
    
    return(list(paths = paths, mcn = mcn))
  }


#################################################
#' loose.ends
#' 
#' takes jbaMIP output and outputs a vector of ranges
#' on the right or left end of the intervals
#' that have type 1-4 labels where
#' 
#' type1 = cn drop, no junction on this side, slack
#' type2 = cn drop, no junction on this side slack
#' type3 = no cn diff, used junction on other side, slack
#' type4 = no cn diff, unused junction on other side, no slack
#'
#'
#' @param sol JaBbA object
#' @param kag karyograph object
#' @return vector of ranges
#' on the right or left end of the intervals
#' that have type 1-4 labels where
#' @export
#################################################
loose.ends = function(sol, kag)
  {
    if(!any(sol$segstats$eslack.in>0 | sol$segstats$eslack.out>0, na.rm = T))
      return(GRanges(seqinfo = seqinfo(kag$segstats)))

    nnab = !ifelse(is.na(sol$ab.edges[,3,1]), TRUE, sol$ab.edges[,3,1]==0)

    if (any(nnab))
        adj.ab = sparseMatrix(as.numeric(kag$ab.edges[nnab,1,]), as.numeric(kag$ab.edges[nnab,2,]),
            x = sol$adj[cbind(as.numeric(kag$ab.edges[nnab,1,]), as.numeric(kag$ab.edges[nnab,2,]))], dims = dim(sol$adj))
    else
        adj.ab = sol$adj*0
    
    ss = sol$segstats
    ss$num = 1:length(ss)
    n = length(ss)
    ss$left.ab = ss$cn.diff = ss$right.ab = -1;
    neg.ix = which( as.logical( strand(ss)=='-') )
    adj.ab[neg.ix, ] = adj.ab[rev(neg.ix), ]
    adj.ab[ ,neg.ix] = adj.ab[, rev(neg.ix)]
    ix.right = 1:n %in% as.numeric(kag$ab.edges[,1,])
    ix.left = 1:n %in% as.numeric(kag$ab.edges[,2,])
    ix.right[neg.ix] = ix.right[rev(neg.ix)]
    ix.left[neg.ix] = ix.left[rev(neg.ix)]
    
    ss[ as.logical( strand(ss)=='-' )] = rev(ss[ as.logical( strand(ss)=='-' ) ])
    tmp.right = Matrix::rowSums(adj.ab)
    tmp.left = Matrix::colSums(adj.ab)
    mask = c(as.numeric(diff(as.numeric(as.factor(seqnames(ss))))==0 & diff(as.numeric(as.factor(strand(ss))))==0), 0)
    ss$left.ab[ix.left] = tmp.left[ix.left]
    ss$right.ab[ix.right] = tmp.right[ix.right]
    ss$cn.diff = c(diff(ss$cn), 0)* mask
    
    ## now classify loose ends
    
    ## (1)
    ix.next = 1+(1:n)
    
    type1 = (rowSums(cbind(ss$eslack.out>0, ss$eslack.in[ix.next]>0))>0 &
             rowSums(cbind(ss$right.ab>0, ss$left.ab[ix.next]>0))==0 &
             ss$cn.diff != 0) * mask
        
    type2 = (rowSums(cbind(ss$eslack.out>0, ss$eslack.in[ix.next]>0))>0 &
             rowSums(cbind(ss$right.ab>0, ss$left.ab[ix.next]>0))>0 &
             ss$cn.diff != 0) * mask

    type3 = (rowSums(cbind(ss$eslack.out>0, ss$eslack.in[ix.next]>0))>0 &
             rowSums(cbind(ss$right.ab>0, ss$left.ab[ix.next]>0))==1 &
             ss$cn.diff == 0) * mask
            
    type4 = (rowSums(cbind(ss$right.ab==0, ss$left.ab[ix.next]==0))>0 &
             rowSums(cbind(ss$right.ab>0, ss$left.ab[ix.next]>0))==0 &
             ss$cn.diff == 0) * mask
        
    ss.p = ss[ as.logical( strand(ss)=='+' ) ]
    win.size = 1
    
    ss.p$num = 1:length(ss.p)
    
    slacks.tmp1 = gr.end(ss.p[which(ss.p$eslack.out>0)], win.size, force = T)
    if (length(slacks.tmp1)>0)
      {
        slacks.tmp1$type = '?'
        slacks.tmp1$type[which(ss.p$eslack.out>0) %in% which(type1!=0)] = 'type1'
        slacks.tmp1$type[which(ss.p$eslack.out>0) %in% which(type2!=0)] = 'type2'
        slacks.tmp1$type[which(ss.p$eslack.out>0) %in% which(type3!=0)] = 'type3'
        slacks.tmp1$cn = slacks.tmp1$eslack.out
        slacks.tmp1$sink = TRUE
      }

    slacks.tmp2 = gr.flipstrand(gr.start(ss.p[which(ss.p$eslack.in>0)], win.size, force = T))
    if (length(slacks.tmp2)>0)
      {        
        slacks.tmp2$type = '?'
        slacks.tmp2$type[which(ss.p$eslack.in>0) %in% (which(type1!=0)+1)] = 'type1'
        slacks.tmp2$type[which(ss.p$eslack.in>0) %in% (which(type2!=0)+1)] = 'type2'
        slacks.tmp2$type[which(ss.p$eslack.in>0) %in% (which(type3!=0)+1)] = 'type3'
        slacks.tmp2$cn = slacks.tmp2$eslack.in
        slacks.tmp2$sink = FALSE ## i.e. source
      }
    
    ## with type 4 there is either a right facing ab w copy 0 or a left facing ab w copy 0 (from the "next" interval)
    t4.ix = intersect(which(type4!=0), 1:length(ss.p))
    t4.ix = ifelse(ss.p[t4.ix]$right.ab==0, -(t4.ix+1), t4.ix) ##
          
    slacks.t4 = c(gr.end(ss.p[t4.ix[t4.ix>0]], win.size, force = T), gr.flipstrand(gr.start(ss.p[-t4.ix[t4.ix<0]], win.size, force = T)))

    if (length(slacks.t4)>0)
      {
        slacks.t4$type = 'type4'
        slacks.t4$cn = 0
        slacks.t4$sink = t4.ix>0 ## t4.ix>0 means an incoming (i.e. left facing) ab edge w copy 0, hence the slack is a "sink" node leaving the other side of the breakpoint
      }
    
    loose.ends = grbind(slacks.tmp1, slacks.tmp2, slacks.t4)[, c('cn', 'num', 'type', 'sink')]
    
    return(loose.ends)    
  }

#################################################
#' slack.breaks
#' 
#' makes "breaks" GRL from slack edges in jbaMIP segstats output (GRanges object)
#' 
#' makes fake slack chromosomes corresponding to each nonzero slack edge
#' and connects points on the genome in the following orientation:
#' 
#' + segment source slack --> - granges at endpoint of segment
#' + segment target slack --> + granges at start point of segment
#' - segment source slack --> + granges at start point of segment
#' - segment target slack --> - granges at endpoint of segment
#' 
#################################################
slack.breaks = function(segstats)
  {
    if (any(!(c('eslack.out', 'eslack.in') %in% names(values(segstats)))))
      stop('Input segstats should be $segstats field output of jbaMIP, i.e. be a GRAnges with fields eslack.out and eslack.in')

    ix.pos.source = which(as.logical( strand(segstats)== '+' ) & segstats$eslack.out>0)
    ix.pos.target = which( as.logical( strand(segstats)== '+' ) & segstats$eslack.in>0)
    ix.neg.source = which( as.logical( strand(segstats)== '-' ) & segstats$eslack.out>0)
    ix.neg.target = which( as.logical( strand(segstats)== '-' ) & segstats$eslack.in>0)

    if (length(ix.pos.source)>0)
      names(ix.pos.source) = paste('source.slack.', ix.pos.source, sep = '')

    if (length(ix.neg.source)>0)
      names(ix.neg.source) = paste('source.slack.', ix.neg.source, sep = '')

    if (length(ix.pos.target)>0)
      names(ix.pos.target) = paste('target.slack.', ix.pos.target, sep = '')

    if (length(ix.neg.target)>0)
      names(ix.neg.target) = paste('target.slack.', ix.neg.target, sep = '')
       
    all.names = c(names(ix.pos.source), names(ix.neg.source), names(ix.pos.target), names(ix.neg.target))
    sl = c(seqlengths(segstats), structure(rep(1, length(all.names)), names = all.names))
    
    out1 = GRanges(seqnames(segstats)[c(ix.pos.source, ix.pos.target, ix.neg.source, ix.neg.target)],
      IRanges(c(end(segstats)[ix.pos.source], start(segstats)[ix.pos.target], start(segstats)[ix.neg.source], end(segstats)[ix.neg.target]),
              width = 1), str = c(rep('-', length(ix.pos.source)), rep('+', length(ix.pos.target)), rep('+', length(ix.neg.source)), rep('-', length(ix.neg.target))), seqlengths = sl)

    out2 = GRanges(names(c(ix.pos.source, ix.pos.target, ix.neg.source, ix.neg.target)), IRanges(1, width = 1), str = '+', seqlengths = sl)
    
    out = split(c(out1, out2), names(c(ix.pos.source, ix.pos.target, ix.neg.source, ix.neg.target)))
    values(out)$cn = c(segstats$eslack.out[ix.pos.source], segstats$eslack.in[ix.pos.target],
                 segstats$eslack.out[ix.neg.source], segstats$eslack.in[ix.neg.target])
    
    return(out)     
  }


###############################################################
#' karyoMIP
#'
#' MIP to locally compute walks in an existing JaBbA reconstruction, note: usually many optimal solutions to a given run.
#' Used by jabba.walk.
#'
#' TODO: Make user friendly, still pretty raw
#' 
#' takes |E| x k matrix of k extreme paths (i.e. contigs) across e edges of the karyograph
#' and length |E| vector of edge copy numbers (eclass), length |E| vector of edge equivalence classes (both outputs of jbaMIP.process)
#' and computes most likely karyotypes that fit the edge copy number profile subject to some prior likelihood
#' over the k extreme paths
#'
#' @param K  |E| x k binary matrix of k "extreme" contigs across |E| edges
#' @param e  edge copy numbers across |E| edges
#' @param eclass  edge equivalence classes, used to constrain strand flipped contigs to appear in solutions together, each class can have at most 2 members
#' @param prior  prior log likelihood of a given contig being in the karyotype
#' @param cpenalty karyotype complexity penalty - log likelihood penalty given to having a novel contig in the karyotype, should be calibrated to prior, i.e. higher than the contig-contig variance in the prior, otherwise complex karyotypes may be favored
#' @param tilim time limit to optimizatoin
#' @param nsolutions how many equivalent solutions to report
#' @return
#' Rcplex solution list object with additional field $kcn for path copy number, $kclass for k class id, $mval for mval
#' @export 
###############################################################
karyoMIP = function(K, # |E| x k binary matrix of k "extreme" contigs across |E| edges
  e, # edge copy numbers across |E| edges
  eclass = 1:length(e), # edge equivalence classes, used to constrain strand flipped contigs to appear in solutions together,
                        # each class can have at most 2 members
  kclass = NULL,
  prior = rep(0, ncol(K)), # prior log likelihood of a given contig being in the karyotype
  mprior = NULL, # matrix prior which should be a binary matrix of m x k, eg mapping contigs to their read / barcode support
                 # will result in the addition of a quadratic objective in addition to the complexity penalty 
  cpenalty = 1, # karyotype complexity penalty - log likelihood penalty given to having a novel contig in the karyotype,
                # should be calibrated to prior, i.e. higher than the contig-contig variance in the prior,
                # otherwise complex karyotypes may be favored
  tilim = 100, epgap = 1, nsolutions = 50, objsense = 'max', ...)
  {
    require(Rcplex)
    
    M = 1e7;
    K = as(K, 'sparseMatrix')

    if (length(prior)!=ncol(K))
      stop('prior must be of the same length as number of columns in K')
        
    # variable indices
    v.ix = 1:ncol(K)
    M.ix = max(v.ix) + (1:ncol(K))
    n = max(M.ix);
    
    # add big M constraints
    Zero = sparseMatrix(1, 1, x = 0, dims = c(n, n)) # upper bound is "infinity" if indicator is positive 
    Amub = Zero[1:length(M.ix), ]
    Amub[cbind(1:length(M.ix), v.ix)] = 1
    Amub[cbind(1:length(M.ix), M.ix)] = -M
    
    Amlb = Zero[1:length(M.ix), ] # lower bound a little > 0 if indicator is positive 
    Amlb[cbind(1:length(M.ix), v.ix)] = 1
    Amlb[cbind(1:length(M.ix), M.ix)] = -0.1

    if (is.null(kclass))
      kclass = .e2class(K, eclass)
    
    kclass.counts = table(kclass)
    if (any(kclass.counts>1)) ## any equiv i.e. strand flipped contig pairs? then make sure they appear in solutions togethrer
      {
        bikclass = which(kclass.counts>1)
        Ac = Zero[1:length(bikclass), ]
        pairs = matrix(unlist(split(1:length(kclass), kclass)[as.character(bikclass)]), ncol = 2, byrow = T)
        Ac[cbind(1:nrow(pairs), pairs[,1])] = 1
        Ac[cbind(1:nrow(pairs), pairs[,2])] = -1        
      }  
    else
      Ac = Zero[1,,drop = FALSE]
                  
    # combine constraints
    A = rBind(cBind(K, Zero[rep(1, nrow(K)), M.ix]), Amub, Amlb, Ac);
    b = c(e, rep(0, nrow(Amlb)*2), rep(0, nrow(Ac)));
    sense = c(rep('E', nrow(K)), rep('L', nrow(Amlb)), rep('G', nrow(Amlb)), rep('E', nrow(Ac)))
    vtype = c(rep('I', length(v.ix)), rep('B', length(M.ix)))
    cvec = c(rep(0, length(v.ix)), prior-cpenalty*rep(1, length(M.ix)))

    if (is.null(mprior))
       sol = Rcplex(cvec = cvec, Amat = A, bvec = b, sense = sense, Qmat = NULL, lb = 0, ub = Inf, n = nsolutions, objsense = objsense, vtype = vtype, control = c(list(...), list(tilim = tilim, epgap = epgap)))
    else
    {
      if (!is.matrix(mprior))
        stop('mprior must be matrix')

      if (ncol(mprior) != ncol(K))
        stop('mprior must be matrix with as many columns as there are walks')

      m = nrow(mprior)
      message('Adding mprior to karyoMIP')

      ## column cat matrices with blank Zero matrix on left, with
      ## constraints acting on binary variables, then identity matrix on right to capture
      ## the new "barcode residual" variables and their associated indicator variables

      ## barcode constraints + 2*m additional variabes
      Ap = cbind(Zero[rep(1, nrow(mprior)), rep(1, length(M.ix))], sign(mprior), -diag(rep(1, nrow(mprior))), 0*diag(rep(1, nrow(mprior))))
      prix = n + 1:m ## indices of prior residuals
      iprix = n + m + 1:m ## indices of indicators of prior residuals
      pb = rep(0, nrow(mprior))
      psense = rep('E', nrow(mprior))

      Mplb = Mpub = 0*Ap ## upper and lower bounds for indictors
      Mpub[cbind(1:length(prix), prix)] = 1
      Mpub[cbind(1:length(prix), iprix)] = -M
      Mplb[cbind(1:length(prix), prix)] = 1
      Mplb[cbind(1:length(prix), iprix)] = -0.1
      pmb = rep(0, 2*nrow(Mpub))
      pmsense = c(rep('L', nrow(Mpub)), rep('G', nrow(Mplb)))
      
      ## define additional variables
      pvtype = c(rep('C', nrow(mprior)), rep('B', nrow(mprior)))

      ## objective function weighs the rows of mprior (barcodes) according to their max weight
      ## so user can weigh importance of individual barcodes
      ## or tune the overall importance of barcodes vs parsimony
      pcvec = c(rep(0, m), apply(mprior, 1, max)) 

      A = rBind(cBind(A, sparseMatrix(1, 1, x = 0, dims = c(nrow(A), 2*m))), Ap, Mpub, Mplb)
      b = c(b, pb, pmb)
      sense = c(sense, psense, pmsense)
      vtype = c(vtype, pvtype)
      cvec = c(cvec, pcvec)      
      
      message('Solving optimization with additional ', m, ' matrix prior terms')      

      sol = Rcplex(cvec = cvec, Amat = A, bvec = b, sense = sense, Qmat = NULL, lb = 0, ub = Inf, n = nsolutions, objsense = objsense, vtype = vtype, control = c(list(...), list(tilim = tilim, epgap = epgap)))
    }
    
    if (!is.null(sol$xopt))
      sol = list(sol)
    
    sol = lapply(sol, function(x)
      {
        x$kcn = round(x$xopt[v.ix])
        x$kclass = kclass
        x$mval= round(x$xopt[M.ix])
        return(x)
      })
        
    return(sol)
  }

##############################################################
#' karyoMIP.to.path 
#' 
#' for a karyoMIP solution and associated K matrix of n x e elementary paths  (input to karyoMIP), and v x e edge signed incidence matrix
#' 
#'
#' @param sol solution to karyoMIP
#' @param K matrix of elementary paths (input to karyoMIP)
#' @param e nrow(K) x 2 edge matrix representing vertex pairs (i.e. edges to which K is referring to)
#' @param gr optional GRanges whose names are indexed by rownames of B
#' @param mc.cores integer number of cores 
#' @param verbose flag
#' @return
#' A list with following items: 
#' $path length k list of paths, cycles (each item i is vector denoting sequence of vertices in G )
#' $is.cycle length k logical vector whose component i denotes whether path i is cyclic
#' $cn  length k integer vector whose component i denotes copy number of contig i
#' $path.grl if path.grl == T
#' @export
##############################################################
karyoMIP.to.path = function(sol, ## karyoMIP solutions, i.e. list with $kcn, $kclass (edges vectors)
  K, ## K matrix input to karyomip (edges x paths)
  e, ## nrow(K) x 2 edge matrix representing vertex pairs (i.e. edges to which K is referring to)
  gr = NULL, ## optional GRanges who names are indexed by <<rownames>> of B
  mc.cores = 1,
  verbose = T
  )
{
  contigs = which(sol$kcn!=0)
  c1 =  contigs[!duplicated(sol$kclass[contigs])]
  c2 = setdiff(contigs, c1)
  c2 = c2[match(sol$kclass[c1], sol$kclass[c2])]
  contigs = c1
  contigs2 = c2
  
  nm.gr = names(gr)
  names(gr) = NULL
  
  if (is.null(nm.gr))
    nm.gr  = 1:length(gr)
  
  if (any(duplicated(nm.gr)))
    nm.gr = 1:length(gr)
  
  if (!is.character(e))
    e = matrix(as.character(e), ncol = 2)

  out = list();

  i1 = which(!is.na(e[,1]))
  i2 = which(!is.na(e[,2]))
  B = sparseMatrix(as.numeric(c(e[i1,1], e[i2,2])),  c(i1,i2), x = c(rep(-1, length(i1)), rep(1, length(i1))))
  rownames(B) = 1:nrow(B)

  ## tells us whether the given contig is a cycle .. cycles represent any path lacking net flow in a 
  ## non-slack vertex

  is.slack = rowSums(is.na(e))!=0
  
  out$is.cyc = Matrix::colSums(K[is.slack, contigs, drop = F])==0 & Matrix::colSums((B %*% K[, contigs, drop = F])!=0)==0
  out$cn = sol$kcn[contigs]
  out$kix = contigs;
  out$kix2 = contigs2;

  K = K[, contigs, drop = F]
  out$paths = mclapply(1:length(contigs),
    function(i)
    {
      if (verbose)
        cat('contig', i, 'of', length(contigs), '\n')

      k = K[, i]
      v.all = setdiff(as.vector(e[k!=0,]), NA)
##      v.all = rownames(B)[which(rowSums(abs(B) %*% k)>0)]  ## vertices associated with edges in path / cycle  k
      
      if (length(v.all)==1) ## this is a slack to slack path involving 1 node
        return(v.all)
      
      ## make subgraph corresponding to edges in this path / cycle
##       B.tmp = B[, which(!is.slack)[k[!is.slack]!=0], drop = F] ##      
##       so = rownames(B.tmp)[apply(B.tmp, 2, function(x) which(x<0))]
##       si = rownames(B.tmp)[apply(B.tmp, 2, function(x) which(x>0))]
##       sG = graph(rbind(so, si))
##       sG = graph(rbind(so, si))

      tmp.e = e[k!=0, ,drop = F]
      tmp.e = tmp.e[rowSums(is.na(tmp.e))==0,,drop = F]      
      sG = graph(t(tmp.e))
            
      if (out$is.cyc[i])
        {
          p.fwd = names(get.shortest.paths(sG, v.all[1], v.all[pmin(length(v.all), 2)])$vpath[[1]])
          p.bwd = names(get.shortest.paths(sG, v.all[pmin(length(v.all), 2)], v.all[1])$vpath[[1]])
          return(unique(unlist(c(p.fwd, p.bwd))))
        }
      else
        {
          io = as.numeric(B[, !is.slack, drop = F] %*% k[!is.slack])
          v.in = rownames(B)[io<0][1]
          v.out = rownames(B)[io>0][1]
          return(names(get.shortest.paths(sG, v.in, v.out)$vpath[[1]]))
        }
    }, mc.cores = mc.cores)
  
  if (!is.null(gr))
      {
      if (is.null(nm.gr))
        nm.gr = names(B)
      names(gr) = NULL
      out$grl = do.call('GRangesList', lapply(out$paths, function(x) gr[match(x, nm.gr), c()]))  ## match non-slack vertices
      names(out$grl) = paste('Contig ', out$kix, ' (CN = ', out$cn, ')', sep = '')
      values(out$grl)$is.cycle = out$is.cyc      
    }
   
  return(out)
}


####################################################
#' jabba.hood
#' 
#' Given JaBbA  object
#' and seed window "win", outputs a reduced set of windows within neighborhoof of n coordinate (ork nodes)
#' within the seed region(s) on the graph (only includes edges with weight !=0)
#'
#' @param jab JaBbA object
#' @param win GRanges of window of interest
#' @param d = distance in coordinates on graph
#' @param k Neighborhood on graph around window of interest to query
#' @param pad pad level at which to collapse nearly reference adjacent intervals
#' @return a reduced set of windows within neighborhood k
#' of seed on the graph (only includes edges with weight !=0)
#' @export
#########x############################################
jabba.hood = function(jab, win, d = 0, k = NULL, pad = 0, ignore.strand = TRUE, bagel = FALSE, verbose = FALSE)
{

    if (!is(win, 'GRanges'))
        win = dt2gr(gr2dt(win))
        
    if (ignore.strand)
        win = gr.stripstrand(win)
    
    if (is.null(k)) ## use distance
        {            
            ss = tryCatch(c(jab$segstats[jab$segstats$loose == FALSE, c()], win[, c()]), error = function(e) NULL)
            
            if (is.null(ss))
                ss = grbind(c(jab$segstats[jab$segstats$loose == FALSE, c()], win[, c()]))
            
            if (ignore.strand)
                ss = gr.stripstrand(ss)

            ss = disjoin(ss)            
            win = gr.findoverlaps(ss, win, ignore.strand = ignore.strand)
                         
            seg.s = suppressWarnings(gr.start(ss, ignore.strand = TRUE))
            seg.e = suppressWarnings(gr.end(ss, ignore.strand = TRUE))
            D.s = suppressWarnings(jabba.dist(jab, win, seg.s, verbose = verbose))
            D.e = suppressWarnings(jabba.dist(jab, win, seg.e, verbose = verbose))
            
            min.s = apply(D.s, 2, min, na.rm = TRUE)
            min.e = apply(D.e, 2, min, na.rm = TRUE)
            s.close = min.s<=d
            e.close = min.e<=d

            ## now for all "left close" starts we add whatever distance to that point + pad
            gr.start(ss)[s.close]
            

            out = GRanges()
            if (any(s.close))
                out = c(out, GenomicRanges::flank(seg.s[s.close], -(d-min.s[s.close])))
    
            if (any(e.close))
                out = c(out, GenomicRanges::shift(flank(seg.e[e.close], d-min.e[e.close]),1))

            if (!bagel)
                out = streduce(c(win[, c()], out[, c()]))

            return(streduce(out, pad))
        }
    else ## use graph connections
        {      
            G = tryCatch(graph.adjacency(jab$adj!=0), error = function(e) NULL)

            ix = which(jab$segstats %^% win)
            if (is.null(G)) ## sometimes igraph doesn't like Matrix
                G = graph.edgelist(which(jab$adj!=0, arr.ind = TRUE))
            vix = unique(unlist(neighborhood(G, ix, order = k)))              
            return(streduce(jab$segstats[vix], pad))
        }
}


######################################################
#' jabba.dist
#'
#' Computes distance between pairs of intervals on JaBbA graph
#' 
#' Given "jabba" object and input granges gr1 and gr2 of (signed) intervals
#'
#'
#' @param jab JaBbA object
#' @param gr1 interval set 1 GRanges
#' @param gr2 interval set 2 GRanges
#' @param matrix flag whteher to output a matrix
#' @param max.dist numeric (default = Inf), if non-infinity then output will be a sparse matrix with all entries that are greater than max.dist set to zero
#' @return a length(gr1) x length(gr2) matrix whose entries ij store the distance between
#' the 3' end of gr1[i] and 5' end of gr2[j]
#' @export
#######################################################
jabba.dist = function(jab, gr1, gr2,
                      matrix = T, ## if false then will return a data frame with fields $i $j $dist specifying distance between ij pairs
                      directed= FALSE, ## flag specifying whether we are computing a "directed distance" across only paths FROM gr1 TO gr2 on graph (ie gr2-->gr1 paths do not count
                      max.dist = Inf, ## if max.dist is not Inf then a sparse matrix will be returned that has 0 at all locations greater than max.dist
                      include.internal = TRUE, ## includes internal connections eg if a junction lies inside a feature then that feature is "close" to another feature
                      verbose = FALSE,
                      EPS = 1e-9  ## the value used for "real 0" if a sparse matrix is returned  
  )
{
    if (verbose)
        now = Sys.time()
    
    intersect.ix = gr.findoverlaps(gr1, gr2, ignore.strand = FALSE)
    
    ngr1 = length(gr1)
    ngr2 = length(gr2)
    
    if (is.null(jab$segstats))
      tiles = jab$tile
    else
      tiles = jab$segstats;
    
    if (is.null(jab$G))
      G = graph.adjacency(jab$adj!=0)
    else
      G = jab$G

    ## keep track of original ids when we collapse 
    gr1$id = 1:length(gr1)
    gr2$id = 1:length(gr2)
    
    ## check for double stranded intervals
    ## add corresponding nodes if present
    if (any(ix <- as.logical( strand(gr1)=='*')) )
    {
        strand(gr1)[ix] = '+'
        gr1 = c(gr1, gr.flipstrand(gr1[ix]))
    }
    
    if (any(ix <- as.logical( strand(gr2)=='*')))
    {
        strand(gr2)[ix] = '+'
        gr2 = c(gr2, gr.flipstrand(gr2[ix]))
    }

    ## expand nodes by jabba model to get internal connectivity
    if (include.internal)
    {
        gr1 = gr1[, 'id'] %**% jab$segstats
        gr2 = gr2[, 'id'] %**% jab$segstats	
    }

    if (verbose)
        {
            message('Finished making gr objects')
            print(Sys.time() -now)
        }
    
    tmp = get.edges(G, E(G))
    E(G)$from = tmp[,1]
    E(G)$to = tmp[,2]
    E(G)$weight = width(tiles)[E(G)$to]
               
    gr1.e = gr.end(gr1, ignore.strand = FALSE)
    gr2.s = gr.start(gr2, ignore.strand = FALSE)


    if (!directed)
        {
            gr1.s = gr.start(gr1, ignore.strand = FALSE)
            gr2.e = gr.end(gr2, ignore.strand = FALSE)
        }

    gr1.e$ix = gr.match(gr1.e, tiles, ignore.strand = F) ## graph node corresponding to end of gr1.ew
    gr2.s$ix= gr.match(gr2.s, tiles, ignore.strand = F) ## graph node corresponding to beginning of gr2

    if (!directed)
        {
            gr1.s$ix = gr.match(gr1.s, tiles, ignore.strand = F) ## graph node corresponding to end of gr1.ew
            gr2.e$ix= gr.match(gr2.e, tiles, ignore.strand = F) ## graph node corresponding to beginning of gr2
        }
            
    ## 3' offset from 3' end of query intervals to ends of jabba segs  to add / subtract to distance when query is in middle of a node
    off1 = ifelse(as.logical(strand(gr1.e)=='+'), end(tiles)[gr1.e$ix]-end(gr1.e), start(gr1.e) - start(tiles)[gr1.e$ix])
    off2 = ifelse(as.logical(strand(gr2.s)=='+'), end(tiles)[gr2.s$ix]-end(gr2.s), start(gr2.s) - start(tiles)[gr2.s$ix])

    ## reverse offset now calculate 3' offset from 5' of intervals
    if (!directed)
        {
            off1r = ifelse(as.logical(strand(gr1.s)=='+'), end(tiles)[gr1.s$ix]-start(gr1.s), end(gr1.s) - start(tiles)[gr1.s$ix])
            off2r = ifelse(as.logical(strand(gr2.e)=='+'), end(tiles)[gr2.e$ix]-start(gr2.e), end(gr2.e) - start(tiles)[gr2.e$ix])
        }
    
    ## compute unique indices for forward and reverse analyses
    uix1 = unique(gr1.e$ix)
    uix2 = unique(gr2.s$ix)

    if (!directed)
        {
            uix1r = unique(gr1.s$ix)
            uix2r = unique(gr2.e$ix)
        }

    ## and map back to original indices
    uix1map = match(gr1.e$ix, uix1)
    uix2map = match(gr2.s$ix, uix2)

    if (!directed)
        {
            uix1mapr = match(gr1.s$ix, uix1r)
            uix2mapr = match(gr2.e$ix, uix2r)
        }
    
    self.l = which(diag(jab$adj)>0)

    if (verbose)
        {
            message('Finished mapping gr1 and gr2 objects to jabba graph')
            print(Sys.time() -now)
        }
    
    if (is.infinite(max.dist)) ## in this case we do not bother making sparse matrix and can compute distances very quickly with one call to shortest.paths
    {        
        ## need to take into account forward and reverse scenarios of "distance" here
        ## ie upstream and downstream connections between query and target
        ## edges are annotated with width of target

        ## so for "downstream distance"  we are getting matrix of shortest paths between from uix1 and uix2 node pair 
        ## and then correcting those distances by (1) adding the 3' offset of uix1 from its node
        ## and (2) subtracting the 3' offset of uix2        
        Df = sweep(
            sweep(
                shortest.paths(G, uix1, uix2, weights = E(G)$weight, mode = 'out')[uix1map, uix2map, drop = F],
                1, off1, '+'), ## add uix1 3' offset to all distances 
            2, off2, '-') ## subtract uix2 3' offset to all distances 


        if (!directed)
            {
                ## now looking upstream - ie essentially flipping edges on our graph - the edge weights
                ## now represent "source" node widths (ie of the flipped edges)
                                        # need to correct these distances by (1) subtracting 3' offset of uix1 from its node
                ## and (2) adding the 3' offset of uix2
                ## and using the reverse indices
                Dr = sweep(
                    sweep(
                        t(shortest.paths(G, uix2r, uix1r, weights = E(G)$weight, mode = 'out'))[uix1mapr, uix2mapr, drop = F],
                        1, off1r, '-'), ## substract  uix1 offset to all distances but subtract weight of <first> node
                    2, off2r , '+') ## add uix2 offset to all distances

                Df2 = sweep(
                    sweep(
                        shortest.paths(G, uix1r, uix2, weights = E(G)$weight, mode = 'out')[uix1mapr, uix2map, drop = F],
                        1, off1r, '+'), ## add uix1 3' offset to all distances 
                    2, off2, '-') ## subtract uix2 3' offset to all distances 

                Dr2 = sweep(
                    sweep(
                        t(shortest.paths(G, uix2r, uix1, weights = E(G)$weight, mode = 'out'))[uix1map, uix2mapr, drop = F],
                        1, off1, '-'), ## substract  uix1 offset to all distances but subtract weight of <first> node
                    2, off2r , '+') ## add uix2 offset to all distances
                D = pmin(abs(Df), abs(Dr), abs(Df2), abs(Dr2))
            }
        else
            D = Df
        
        # then we do the same thing but flipping uix1r vs uix        


        if (verbose)
            {
                message('Finished computing distances')
                print(Sys.time() -now)
            }
    
        
        ## take care of edge cases where ranges land on the same node, since igraph will just give them "0" distance
        ## ij contains pairs of gr1 and gr2 indices that map to the same node
        ij = as.matrix(merge(cbind(i = 1:length(gr1.e), nid = gr1.e$ix), cbind(j = 1:length(gr2.s), nid = gr2.s$ix)))
        
        ## among ij pairs that land on the same (strand of the same) node
        ##
        ## several possibilities:
        ## (1) if gr1.e[i] < gr2.s[j] then keep original distance (i.e. was correctly calculated)
        ## (2) if gr1.e[i] > gr2.s[j] then either
        ##   (a) check if there is a self loop and adjust accordingly (i.e. add back width of current tile)
        ##   (b) PITA case, compute shortest distance from i's child(ren) to j

        if (nrow(ij)>0)
          {
            ## rix are present 
              rix = as.logical((
                  (as.logical( strand(gr1)[ij[,'i']] == '+' ) &
                   as.logical( strand(gr2)[ij[,'j']] == '+' ) &
                   end(gr1)[ij[,'i']] <= start(gr2[ij[,'j']])) |
                  ( as.logical( strand(gr1)[ij[,'i']] == '-' ) &
                    as.logical( strand(gr2)[ij[,'j']] == '-' ) &
                    start(gr1)[ij[,'i']] >= end(gr2)[ij[,'j']])))

              ij = ij[!rix, , drop = F] ## NTD with rix == TRUE these since they are calculated correctly
            
            if (nrow(ij)>0) ## any remaining will either be self loops or complicated loops              
              {
                selfix = (ij[, 'nid'] %in% self.l)
                
                if (any(selfix)) ## correct distance for direct self loops (add back width of current node)
                  D[ij[selfix, c('i', 'j'), drop = F]]  = D[ij[selfix, c('i', 'j'), drop = F]] + width(tiles)[ij[selfix, 'nid']] 
                  
                ij = ij[!selfix, , drop = F]

                if (nrow(ij)>0) ## remaining are pain in the ass indirect self loops
                  {                  
                    ch = G[[ij[, 'nid']]] ## list of i nodes children for all remaining ij pairs
                    chu = munlist(ch) ## unlisted children, third column are the child id's, first column is the position of nrix

                    if (ncol(chu)>1)
                        {
                    
                    ## now find paths from children to corresponding j
                            epaths = suppressWarnings(get.shortest.paths(G, chu[, 3], ij[chu[,'ix'], 'nid'], weights = E(G)$weight, mode = 'out', output = 'epath')$epath)
                            epathw = sapply(epaths, function(x,w) if (length(x)==0) Inf else sum(w[x]), E(G)$weight) ## calculate the path weights                    
                            epathw = epathw + width(tiles)[chu[, 3]] + off1[ij[chu[, 'ix'], 'i']] + off2[ij[chu[,'ix'], 'j']] - width(tiles)[ij[chu[, 'ix'], 'nid']]
                            
                            ## aggregate (i.e. in case there are multiple children per node) by taking min width
                            D[ij[, c('i', 'j'), drop = F]] = vaggregate(epathw, by = list(chu[, 'ix']), min)[as.character(1:nrow(ij))]
                        }
                    }
              }
          }
        
        if (verbose)
            {
                message('Finished correcting distances')
                print(Sys.time() -now)
            }
      }

    ## need to collapse matrix ie if there were "*" strand inputs and if we are counting internal
    ## connections inside our queries ..
    ## collapsing +/- rows and columns by max value based on their id mapping to their original "*" interval


    ## melt distance matrix into ij
    ij = as.matrix(expand.grid(1:nrow(D), 1:ncol(D)))
    dt = data.table(i = ij[,1], j = ij[,2], value = D[ij])[, id1 := gr1$id[i]][, id2 := gr2$id[j]]
    
    tmp = dcast.data.table(dt, id1 ~ id2, fun.aggregate = function(x) min(as.numeric(x)))
    setkey(tmp, id1)    
    Dtmp = as.matrix(tmp[list(1:ngr1), -1, with = FALSE])
    D = matrix(NA, nrow = ngr1, ncol = ngr2, dimnames = list(NULL,
                                                             1:as.character(ngr2)))
    D[1:nrow(Dtmp), colnames(Dtmp)] = Dtmp
   

    ## finally zero out any intervals that actually intersect
    ## (edge case not captured when we just examine ends)
    if (length(intersect.ix)>0)
        D[cbind(intersect.ix$query.id, intersect.ix$subject.id)] = 0    

    if (verbose)
        {
            message('Finished aggregating distances to original object')
            print(Sys.time() -now)
        }
    return(D)    
  }


########################################
#' jabba2vcf
#' 
#' Converts jabba output to vcf file according to 4.2 "BND" syntax
#' 
#' 
#' @param jab JaBbA object
#' @param fn output file name
#' @param sampleid sample id
#' @param hg human genome as BSgenome or ffTrack
#' @param cnv flag whether to dump in CNV format
#' @return returns character string or writes to file if specified
#' @export
#########################################
jabba2vcf = function(jab, fn = NULL, sampleid = 'sample', hg = skidb::read_hg(fft = T), cnv = FALSE)
    {
        require(data.table)
        vcffields = c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'GENO')
        
        ## convert all aberrant connections into pairs of VCF rows
        if (!cnv)
            {
                jix = which(!is.na(jab$ab.edges[,3,1])) ## these are the only junctions with breaks in the reconstruction
                abs = rbind(jab$ab.edges[jix,1:2,1])        
                rabs = rbind(jab$ab.edges[jix,1:2,2])
                rcix = match(jab$segstats, gr.flipstrand(jab$segstats)) ## map of seg to its reverse complement

                adj.ref = jab$adj ## reference graph has reference copy numbers, we obtain by zeroing out all ab.edges and loose end edges
                adj.ref[rbind(jab$ab.edges[jix,1:2,1])] = 0
                adj.ref[rbind(jab$ab.edges[jix,1:2,2])] = 0

                if (any(jab$segstatsloose))
                    {
                        adj.ref[jab$segstats$loose, ] = 0
                        adj.ref[,jab$segstats$loose] = 0
                    }

                if (length(jix)>0)
                    {
                        jcn = jab$adj[abs]
                        gr1 = gr.end(jab$segstats[abs[,1]], ignore.strand = F)[, 'cn']
                        gr1$jid = jix
                        gr1$nid = abs[,1]
                        gr1$acn = jcn
                        gr1$rcn = Matrix::rowSums(adj.ref[gr1$nid, , drop = FALSE])        
                        gr1$ID = paste(sampleid, '_seg', abs[,1], ifelse(as.logical(strand(gr1)=='+'), '_R', '_L'), sep = '')
                        
                        gr2 = gr.start(jab$segstats[abs[,2]], ignore.strand = F)[, 'cn']
                        gr2$jid = jix
                        gr2$nid = abs[,2]
                        gr2$acn = jcn
                        gr2$rcn = Matrix::colSums(adj.ref[,gr2$nid, drop = FALSE])
                        gr2$ID = paste(sampleid, '_seg', abs[,2], ifelse(as.logical(strand(gr2)=='+'), '_L', '_R'), sep = '')
                        
                        gr1$mid = gr2$ID
                        gr2$mid = gr1$ID
                        
                        gr1$REF = tryCatch(as.character(ffTrack::get_seq(hg, gr.stripstrand(gr1))), error = function(e) 'N')
                        gr1$ALT = ifelse(as.logical(strand(gr1)=='+') & as.logical(strand(gr2)=='+'), paste(gr1$REF, '[', seqnames(gr2), ':', start(gr2), '[', sep = ''),
                            ifelse(as.logical(strand(gr1)=='+') & as.logical(strand(gr2)=='-'), paste(gr1$REF, ']', seqnames(gr2), ':', start(gr2), ']', sep = ''),
                                   ifelse(as.logical(strand(gr1)=='-') & as.logical(strand(gr2)=='+'), paste('[', seqnames(gr2), ':', start(gr2), '[', gr1$REF, sep = ''),
                                          paste(']', seqnames(gr2), ':', start(gr2), ']', gr1$REF, sep = '')))) ## last one is A- --> B-
                        
                        gr2$REF = tryCatch(as.character(ffTrack::get_seq(hg, gr.stripstrand(gr2))), error = function(e) 'N')
                        gr2$ALT = ifelse(as.logical(strand(gr1)=='+') & as.logical(strand(gr2)=='+'), paste(']', seqnames(gr1), ':', start(gr1), ']', gr2$REF, sep = ''),
                            ifelse(as.logical(strand(gr1)=='-') & as.logical(strand(gr2)=='+'), paste('[', seqnames(gr1), ':', start(gr1), '[', gr2$REF, sep = ''),
                                   ifelse(as.logical(strand(gr1)=='+') & as.logical(strand(gr2)=='-'), paste(gr2$REF, ']', seqnames(gr1), ':', start(gr1), ']', sep = ''),
                                          paste(gr2$REF, '[', seqnames(gr1), ':', start(gr1), '[', sep = '')))) ## last one is B- --> A-
                        
                        
                        gr1$FILTER = gr2$FILTER = ifelse(gr2$acn != 0, "PASS", "NINC")
                        gr2$FORMAT = gr1$FORMAT = "GT:CN:RCN:SCN"
                        gr1$GENO = paste(ifelse(gr1$rcn>0, '0/1', '1'), gr1$acn, gr1$rcn, gr1$cn, sep = ":")
                        gr2$GENO = paste(ifelse(gr2$rcn>0, '0/1', '1'), gr2$acn, gr2$rcn, gr2$cn, sep = ":")
                        
                        gr2$CHROM = as.character(seqnames(gr2))
                        gr1$CHROM = as.character(seqnames(gr1))

                        gr2$POS = as.character(start(gr2))
                        gr1$POS = as.character(start(gr1))
                        
                        gr1$QUAL = gr2$QUAL = '.'
                        gr1$INFO = paste("SVTYPE=BND", ";MATEID=", gr1$mid, ";CNADJ=", gr1$acn, ";CNRADJ=", gr1$rcn, ";CN=", gr1$cn,
                            ";JABID=", abs[,1], ";RJABID=", rcix[abs[,1]], ";JUNCID=", 1:length(gr1), sep = '')
                        gr2$INFO = paste("SVTYPE=BND", ";MATEID=", gr2$mid, ";CNADJ=", gr2$acn, ";CNRADJ=", gr2$rcn, ";CN=", gr2$cn,
                            ";JABID=", abs[,2], ";RJABID=", rcix[abs[,2]], ";JUNCID=", 1:length(gr2), sep = '')

                        gr1 = gr1[, vcffields]
                        gr2 = gr2[, vcffields]
                    }
                else
                    {
                        gr1 = GRanges()
                        gr2 = GRanges()
                    }
                
                ## now loose ends 
                lix = which(jab$segstats$loose & as.logical(strand(jab$segstats)=="+"))
                if (length(lix)>0)
                    { 
                        gr.loose = gr.start(jab$segstats[lix, 'cn']) ## loose ends should be width 1, but just in case
                        gr.loose$jid = NA
                        gr.loose$nid = lix
                        gr.loose$acn = gr.loose$cn                    

                        ## query parents / children of loose ends
                        tmp1 = apply(jab$adj[gr.loose$nid, , drop = FALSE],1, function(x) which(x!=0)[1]) ## only non NA if loose end has children (i.e. is a parent)
                        tmp2 = apply(jab$adj[, gr.loose$nid, drop = FALSE],2, function(x) which(x!=0)[1]) ## only non NA if loose end has parents (i.e. is child)
                        
                        isp = apply(cbind(tmp1, tmp2), 1, function(x) which(is.na(x))[1])==2 ## is the loose end a parent or child of a seg? 
                        pcid = pmax(tmp1, tmp2, na.rm = T) ## which seg is the parent / child of the loose end

##                        gr.loose$rcn = ifelse(isp, colSums(adj.ref[,pcid, drop = FALSE]), Matrix::rowSums(adj.ref[pcid,, drop = FALSE])) ## if loose end is parent of seg, we want num copies of that segs reference parent, 

                        gr.loose$rcn = ifelse(isp, Matrix::colSums(adj.ref[,pcid, drop = FALSE]), Matrix::rowSums(adj.ref[pcid,, drop = FALSE])) ## if loose end is parent of seg, we want num copies of that segs reference parent, 

                        
                        gr.loose$cn = jab$segstats$cn[pcid]
                        ## if loose end is the parent of a seg, then it is a "left" loose end (since + strand) otherwise "right"
                        gr.loose$ID = paste(sampleid, '_looseend', pcid, ifelse(isp, '_L', '_R'), sep = '')
                        gr.loose$mid = NA
                        gr.loose$REF = tryCatch(as.character(ffTrack::get_seq(hg, gr.stripstrand(gr.loose))), error = function(e) 'N')

                        ## again, same rationale as above, parent = left loose end, child = right loose end
                        gr.loose$ALT = ifelse(isp, paste('.', gr.loose$REF, sep = ''), paste(gr.loose$REF, '.', sep = ''))
                        gr.loose$FORMAT = "GT:CN:RCN:SCN"
                        gr.loose$GENO = paste(ifelse(gr.loose$rcn>0, '0/1', '1'), gr.loose$acn, gr.loose$rcn, gr.loose$cn, sep = ":")
                        gr.loose$CHROM = as.character(seqnames(gr.loose))
                        gr.loose$POS = start(gr.loose)
                        gr.loose$FILTER = "LOOSEEND"
                        gr.loose$QUAL = '.'
                        gr.loose$INFO = paste("SVTYPE=BND", ";CNADJ=", gr.loose$acn, ";CNRADJ=", gr.loose$rcn, ";CN=", gr.loose$cn, ";JABID=", lix, ";RJABID=", rcix[lix], sep = '')
                        gr.loose = gr.loose[, vcffields]
                    }
                else
                    gr.loose = GRanges()
                
                                        # make header
                sl = seqlengths(jab$segstats)
                header = '##fileformat=VCFv4.2'
                header = c(header, sprintf('##fileDate=%s', format(Sys.Date(), '%Y%m%d')))
                header = c(header, '##source=JaBbAV0.1')        
                if (inherits(hg, "BSgenome"))
                    {
                        header = c(header, sprintf("##reference=%s", sourceUrl(hg)))
                        header = c(header, unlist(mapply(function(x, y) sprintf('##contig=<ID=%s,length=%s,assembly=%s,species="%s">', x, y, providerVersion(hg), organism(hg)), names(sl), sl)))
                    }
                else
                    {
                        header = c(header, sprintf("##reference=%s", filename(hg)['rds']))
                        header = c(header, unlist(mapply(function(x, y) sprintf('##contig=<ID=%s,length=%s>', x, y), names(sl), sl)))
                    }
                
                header = c(header,
                    '##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">', 
                    '##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">',
                    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
                    '##INFO=<ID=CNADJ,Number=.,Type=Integer,Description="Copy number of variant adjacency at breakend">',
                    '##INFO=<ID=CNRADJ,Number=1,Type=Integer,Description="Copy number of reference adjacency at breakend">',
                    '##INFO=<ID=CN,Number=1,Type=Integer,Description="Copy number of segment containing breakend">',
                    '##INFO=<ID=JUNCID,Number=.,Type=Integer,Description="Index of allele(s) in JaBbA junction input pile">',
                    '##INFO=<ID=JABID,Number=.,Type=Integer,Description="Index of the interval containing allele(s) in JaBbA reconstruction">',
                    '##INFO=<ID=RJABID,Number=.,Type=Integer,Description="Index of the reverse complement containing interval allele(s) in JaBbA reconstruction">',
                    '##FILTER=<ID=PASS,Description="Junction incorporated into JaBbA MIP reconstruction at nonzero copy number">',
                    '##FILTER=<ID=LOOSEEND,Description="Loose end incorporated into JaBbA MIP reconstruction at unexplained copy break">',
                    '##FILTER=<ID=NINC,Description="Not Incorporated into JaBbA MIP reconstruction: either subclonal or false positive">',
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                    '##FORMAT=<ID=CN,Number=.,Type=String,Description="Copy number of variant adjacencies">',
                    '##FORMAT=<ID=RCN,Number=1,Type=String,Description="Copy number of reference adjacency at breakend">',
                    '##FORMAT=<ID=SCN,Number=1,Type=String,Description="Copy number of segment containing breakend">'
                           )
               
                if ((length(gr1) + length(gr.loose))>0)
                    {
                        body = as.data.frame(values(c(gr1, gr2, gr.loose)))
                        body$ord = NA
                        body$ord[order(gr.stripstrand(c(gr1, gr2, gr.loose)))] = 1:nrow(body)
                        body = as.data.table(body)
                        setkeyv(body, c('CHROM', 'POS'))
                        
                        ## useful for merging info records in painful coordinate deduping process below
                        .infomelt = function(str)
                            {
                                z = lapply(strsplit(str, ";"), function(x) {y = matrix(unlist(strsplit(x, "=")), ncol = 2, byrow = T); return(structure(y[,2], names = y[,1]))})
                                unames = unique(unlist(sapply(z, names)))
                                mfields = c("MATEID", "CNADJ", "JUNCID", "JABID", "RJABID")
                                mergeix = unames %in% mfields
                                out = sapply(1:length(unames), function(i) {
                                    x = unames[i]
                                    tmp = sapply(z, function(y) y[x])
                                    if (mergeix[i])
                                        {                        
                                            if (any(is.na(tmp))) tmp[is.na(tmp)] = '.';
                                            return(paste(tmp, collapse = ','))
                                        }
                                    else
                                        return(tmp[!is.na(tmp)][1])})
                                names(out) = unames
                                return(paste(names(out), "=", out, collapse = ";", sep = ''))
                            }

                        .genomelt = function(geno, format)
                            {
                                if (length(geno)==1)
                                    return(genos)
                                genos = strsplit(geno, ":")
                                formats = strsplit(format[1], ":")[[1]]
                                
                                mergeix = formats == "CN"
                                gtix = formats == "GT" ## should be only one such entry
                                out = genos[[1]] ## pick the first one as the default - this should be a vector length(formats)
                                if (any(gtix))
                                    out[gtix] = paste(genos[[1]][gtix][1], '/', paste(2:length(genos), collapse = '/'), sep = '') ## add "fake allele names"

                                
                                if (any(mergeix))
                                    out[mergeix] = sapply(which(mergeix), function(x)
                                        paste(sapply(genos, function(y) y[[x]]), collapse = ','))
                                return(paste(out, collapse = ":"))
                            }
                        
                        ## dedup and collapse breakends that share coordinates into single variant sites with several alleles
                        body = body[, list(
                            ID = ID[1],
                            REF = REF[1],
                            ALT = if (length(ALT)>1) paste(ALT, collapse = ',') else ALT,            
                            QUAL = QUAL[1],
                            FILTER = if (length(FILTER)>1) paste(FILTER, collapse = ";") else FILTER,
                            INFO = if (length(INFO)>1) .infomelt(INFO) else INFO,
                            FORMAT = FORMAT[1], ## assume everything is same format
                            GENO = if (length(GENO)>1) .genomelt(GENO, FORMAT) else GENO,
                            ord = ord[1]
                        ), by = c("CHROM", "POS")]
                        
                        setkey(body, "ord")

                        body = as.data.frame(body)[, vcffields]
                        names(body)[ncol(body)] = sampleid
                        names(body)[1] = '#CHROM'
                    }
                else
                    {
                        body = as.data.frame(matrix(NA, ncol = length(vcffields), dimnames = list(c(), vcffields), ))[c(), ]
                        names(body)[1] = '#CHROM'
                        names(body)[ncol(body)] = sampleid
                    }
            }
        else ## CNV mode
            {
                sl = seqlengths(jab$segstats)
                header = '##fileformat=VCFv4.2'
                header = c(header, sprintf('##fileDate=%s', format(Sys.Date(), '%Y%m%d')))
                header = c(header, '##source=JaBbAV0.1 CNV')

                if (inherits(hg, "BSgenome"))
                    {
                        header = c(header, sprintf("##reference=%s", sourceUrl(hg)))
                        header = c(header, unlist(mapply(function(x, y) sprintf('##contig=<ID=%s,length=%s,assembly=%s,species="%s">', x, y, providerVersion(hg), organism(hg)), names(sl), sl)))
                    }
                else if (inherits(hg, "ffTrack"))
                {
                    header = c(header, sprintf("##reference=%s", filename(hg)['rds']))
                    header = c(header, unlist(mapply(function(x, y) sprintf('##contig=<ID=%s,length=%s>', x, y), names(sl), sl)))
                }                
                else
                {
                    header = c(header, sprintf("##reference=NA"))
                    header = c(header, unlist(mapply(function(x, y) sprintf('##contig=<ID=%s,length=%s>', x, y), names(sl), sl)))
                }
                
                header = c(header,
                    '##ALT=<ID=DEL,Description="Decreased copy number relative to reference">',
                    '##ALT=<ID=DUP,Description="Increased copy number relative to reference">',
                    '##ALT=<ID=DIP,Description="Normal diploid copy number">',
                    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
                    '##INFO=<ID=START,Number=1,Type=Integer,Description="Start of copy number annotated interval">',
                    '##INFO=<ID=END,Number=1,Type=Integer,Description="End of copy number annotated interval">',
                    '##INFO=<ID=JABID,Number=1,Type=Integer,Description="Index of the interval containing allele(s) in JaBbA reconstruction">',
                    '##INFO=<ID=RJABID,Number=1,Type=Integer,Description="Index of the reverse complement  interval allele(s) in JaBbA reconstruction">',                    
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                    '##FORMAT=<ID=CN,Number=1,Type=String,Description="Copy number">'
                           )

                rcix = match(jab$segstats, gr.flipstrand(jab$segstats)) ## map of seg to its reverse complement
                six = which(!is.na(jab$segstats$cn) & !jab$segstats$loose & as.character(strand(jab$segstats))=='+')
                ss = jab$segstats[six]
                
                if (length(ss)>0)
                    {
                        body = data.frame("CHROM" = seqnames(ss), POS = start(ss),
                            ID = paste(sampleid, '_seg', six, sep = ''),
                            REF = as.character(ffTrack::get_seq(hg, gr.stripstrand(gr.start(ss,1)))),
                            ALT = ifelse(ss$cn==2, "<DIP>", ifelse(ss$cn>2, "<DUP>", "<DEL>")),
                            QUAL = ".",
                            FILTER = "PASS",
                            INFO = paste("SVTYPE=", ifelse(ss$cn>2, "DUP", "DEL"),
                                ";START=", start(ss), ";END=", end(ss), ";JABID=", six, ";RJABID=", rcix[six], sep = ''),
                            FORMAT = "GT:CN",
                            GENO = paste("./.", ss$cn, sep = ':'))
                        names(body)[1] = '#CHROM'
                        names(body)[ncol(body)] = sampleid
                    }
                else
                    {
                        body = as.data.frame(matrix(NA, ncol = length(vcffields), dimnames = list(c(), vcffields), ))[c(), ]
                        names(body)[1] = '#CHROM'
                        names(body)[ncol(body)] = sampleid
                    }
            }
        
        if (!is.null(fn))
            {
                writeLines(header, fn)
                suppressWarnings(write.table(body, fn, sep = '\t', quote = F, row.names = F, append = T))
            }
        else
            {
                t = textConnection("out", 'w')
                writeLines(header, t)
                suppressWarnings(write.table(body, t, sep = '\t', quote = F, row.names = F, append = T))
                close(t)
                return(out)
            }        
    }


get.constrained.shortest.path = function(cn.adj, ## copy number matrix
                                         G, ## graph with distances as weights
                                         allD=NULL, ## shortest path between all nodes in graph
                                         v,
                                         to,
                                         weight,
                                         edges,
                                         verbose = TRUE,
                                         mip = TRUE
                                         )
{

    if (is.null(allD)) allD = shortest.paths(G, mode="out", weights = weight)

    v = as.numeric(v)
    to = as.numeric(to)
    
    if (is.infinite(allD[v, to]) | allD[v, to]==0) return(NULL)

    edges$cn = cn.adj[cbind(edges$from, edges$to)]
           
    ## ASSUME: from, to are scalars, within node range, to is reachable from from
    ## ASSUME edges contains eid key and eclass mapping
    tmp.p = as.numeric(get.shortest.paths(G, from=v, to=to, "out", weights=weight)$vpath[[1]])
    tmp.e = cbind(tmp.p[-length(tmp.p)], tmp.p[-1])
    tmp.eid = paste(tmp.e[, 1], tmp.e[, 2])
    tmp.eclass = edges[.(tmp.eid), eclass]


    ## the cn of this path is the max number of copies that the network will allow
    ## here we have to group by eclass, i.e. so if there are two edges from an eclass
    ## in a given path then we need to halve the "remaining copies" constraint
    tmp.pcn = edges[.(tmp.eid), if (length(cn)>1) cn/2 else cn, by = eclass][, floor(min(V1))]

    
    edges[, rationed := cn<(tmp.pcn*2)]

    D.totarget = allD[, as.numeric(to)]
    edges[, distance_to_target :=  D.totarget[to]]
    edges = edges[!is.infinite(distance_to_target) & cn>0, ]
    
    rationed.edges = edges[rationed == TRUE, ]

    ## find overdrafted eclasses - meaning two instances in this path but only one remaining copy
    overdrafts.eclass = intersect(names(which(table(tmp.eclass)==2)), rationed.edges$eclass)
     
    first.overdraft = which(tmp.eclass %in% overdrafts.eclass & duplicated(tmp.eclass))[1]

    ## no overdrafts?, then return
    if (is.na(first.overdraft) & tmp.pcn>0)
    {
        if (verbose)
            message('Shortest path is good enough!')
        return(tmp.p)
    }

    if (!mip)
        return(NULL)


    ## use MIP to find constrained path
    edges[, enum := 1:length(eid)]

    ## incidence matrix constraints + 1 for tmp.pcn
    A = sparseMatrix(edges$to, edges$enum, x = 1, dims = c(nrow(cn.adj), nrow(edges))) -
        sparseMatrix(edges$from, edges$enum, x = 1, dims = c(nrow(cn.adj), nrow(edges))) 
    b = rep(0, nrow(A))
    b[v] = -1
    b[to] = 1

    ix = which(rowSums(A!=0)!=0) ## remove zero constraints

    ## "ration" or reverse complementarity constraints
    tmp.constraints = edges[, list(e1 = enum[1], e2 = enum[2], ub = cn[1]), by = eclass]
    tmp.constraints = tmp.constraints[!is.na(e1) & !is.na(e2), ]

    R = sparseMatrix(rep(1:nrow(tmp.constraints), 2),
                     c(tmp.constraints$e1, tmp.constraints$e2),
                     x = 1, dims = c(nrow(tmp.constraints), nrow(edges)))
    Rb = tmp.constraints$ub

    ## minimize weight of path
    c = edges$weight

    res = Rcplex(c, rbind(A[ix,], R), c(b[ix], Rb), sense = c(rep('E', length(ix)), rep('L', length(Rb))),
                 lb = 0, vtype = "B",
                 objsense = 'min')


    if (verbose)
        message('YES WE ARE DOING PROPER MIP!!!!')

    if (!(res$status %in% c(101, 102)))
    {
        if (verbose)
            message('No solution to MIP!')

        return(NULL)
    }
    
    ## use igraph to sort these edges into a path, i.e. make simple graph with one path and extract it using igraph (lazy :)
    tmp.p = as.numeric(get.shortest.paths(graph_from_edgelist(edges[res$xopt!=0, cbind(from, to)]), v, to)$vpath[[1]])    
    
    ## check if overdrafted
    if (verbose)
    {
        tmp.e = cbind(tmp.p[-length(tmp.p)], tmp.p[-1])
        tmp.eid = paste(tmp.e[, 1], tmp.e[, 2])
        tmp.eclass = edges[.(tmp.eid), eclass]
        tmp.pcn = edges[.(tmp.eid), if (length(cn)>1) cn/2 else cn, by = eclass][, min(V1)]
        overdrafts.eclass = intersect(names(which(table(tmp.eclass)==2)), rationed.edges$eclass)
        if (length(overdrafts.eclass)==0)
            message('No overdrafts after MIP')
        else
        {
            stop('Still overdraft!')
        }
    }

#    browser()
    return(tmp.p)
}



#' @name jabba.gwalk
#' @title jabba.gwalk
#' @description
#'
#' Computes greedy collection (i.e. assembly) of genome-wide walks (graphs and cycles) by finding shortest paths in JaBbA graph.
#' 
#' @param jab JaBbA object
#' #
#' @return GRangesList of walks with copy number as field $cn, cyclic walks denoted as field $is.cycle == TRUE, and $wid (width) and $len (segment length) of walks as additional metadata
#' @export
#' @import igraph
#' @export
#' @author Marcin Imielinski
#' @author Xiaotong Yao
jabba.gwalk = function(jab, verbose = FALSE)
{
    cn.adj = jab$adj
    adj = as.matrix(cn.adj)
    adj.new = adj*0
    ## ALERT!!! see below
    adj[which(adj!=0, arr.ind = TRUE)] = width(jab$segstats)[which(adj!=0, arr.ind = TRUE)[,2]] ## make all edges a large number by default
    ## adj[which(adj!=0, arr.ind = TRUE)] = width(jab$segstats)[which(adj!=0, arr.ind = TRUE)[,1]] ## make all edges a large number by default
    if (verbose){
        ## ALERT!!! I'm gonna switch to source node width for default weight of edges
        message('Setting edge weights to destination widths for reference edges and 1 for aberrant edges')
        ## message('Setting default edge weights to SOURCE widths for edges and 1% less for aberrant edges')
    }
    
    ab.edges = rbind(jab$ab.edges[,1:2, 1], jab$ab.edges[,1:2, 2])
    ab.edges = ab.edges[rowSums(is.na(ab.edges))==0, ]
    ## ALERT!!!
    adj[ab.edges] = sign(cn.adj[ab.edges]) ## make ab.edges = 1
    ## adj[ab.edges] = adj[ab.edges] * 0.99 ## make ab.edges 1 bp shorter than ref!
    adj[is.na(adj)] = 0
    cn.adj[which(is.na(cn.adj))] = 0

    ## ALERT!!! major change
    ## adjj = adj/as.matrix(cn.adj)
    ## adjj[which(is.nan(adjj))] = 0
    ## adjj[which(adjj<0)] = 0
    ## G = graph.adjacency(adjj, weighted = 'weight')
    ## esl = which(adj != 0, arr.ind=T)
    ## eids = paste(esl[,1], esl[,2])
    ## weights = adj[esl]
    ## eclasses = ed[.(eids), eclass]
    G = graph.adjacency(adj, weighted = 'weight')
    ## G = make_graph(t(esl), )

    ## DD = shortest.paths(G, mode="out")
    ## IJ = which(!is.infinite(DD), arr.ind=T)
    
    ## define ends not using degree (old method) but using either telomeres or loose ends
    ## (otherwise lots of fake ends at homozygous deleted segments)
    ss = gr2dt(jab$segstats)[ , vid:= 1:length(seqnames)]
    ss[loose == TRUE, is.end := TRUE]
    ss[loose == FALSE, is.end := 1:length(loose) %in% c(which.min(start), which.max(end)), by = list(seqnames, strand)]
    ends = which(ss$is.end)
    src = (colSums(adj)[ends]==0) ## indicate which are sources

    ## sanity check
    unb = which(!ss$is.end & rowSums(jab$adj, na.rm = TRUE) != colSums(jab$adj, na.rm = TRUE))

    if (length(unb)>0)
    {
        message(sprintf('JaBbA model not junction balanced at %s non-ends! Adding these to "ends"', length(unb)))
        ends = c(ends, unb)         ## shameless HACK ... TOFIX
    }
    
    ## ends = which(degree(G, mode = 'out')==0 | degree(G, mode = 'in')==0)    
    i = 0
    ## adjust weight just before creating D
    ## assign lighter weight to higher copy 
    ## D = shortest.paths(G, v = ends, mode = 'out', weight = E(G)$weight)[, ends]

    ## D records distance from ends to every node
    D = shortest.paths(G, v = ends, mode = 'out', weight = E(G)$weight)[, ends]
    
    ## sort shortest paths
    ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row != col, ][order(dist), ][, row := ends[row]][, col := ends[col]]

    ## ij only record end to end distance
    ## ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[col %in% ends, ][, dist := D[cbind(row, col)]][, row := ends[row]][row != col, ][order(dist), ]
    
    maxrow = length(ends)*max(cn.adj[ends, ends], na.rm = TRUE)
    vpaths = rep(list(NA), maxrow)
    epaths = rep(list(NA), maxrow)
    cns = rep(NA, maxrow)
    palindromic.path = rep(FALSE, maxrow)
    palindromic.cycle = rep(FALSE, maxrow)

    nb.all = which(rowSums(cn.adj) != colSums(cn.adj))
    cn.adj0 = cn.adj
    G0 = G
    D0 = D

    #' first peel off "simple" paths i.e. zero degree
    #' ends with >0 copy number
    psimp =  which(degree(G, mode = 'out')==0 & degree(G, mode = 'in')==0 & jab$segstats$cn>0)
    i = 0
    if (length(psimp)>0)
    {
        vpaths[1:length(psimp)] = split(psimp, 1:length(psimp))
        epaths[1:length(psimp)] = lapply(psimp, function(x) cbind(NA, NA)) ## there is no "edge" associated with a zero total degree node 
        cns[1:length(psimp)] = jab$segstats$cn[psimp]
        i = length(psimp)
    }
    
    ## now iterate from shortest to longest path
    ## peel that path off and see if it is still there ..
    ## and see if it is still there 
    ## peel off top path and add to stack, then update cn.adj

    jab$segstats$tile.id = jab$segstats$tile.id + as.numeric(jab$segstats$loose)*0.5
    
    tile.map =
        gr2dt(jab$segstats)[, .(id = 1:length(tile.id),
                                tile.id = ifelse(strand == '+', 1, -1)*tile.id)]    
    rtile.map =
        gr2dt(jab$segstats)[, .(id = 1:length(tile.id),
                                tile.id = ifelse(strand == '+', 1, -1)*tile.id)]
    setkey(tile.map, id)
    setkey(rtile.map, tile.id)

    ## unique pair of edge ids: rev comp of a foldback edge will be identical to itself!!!
    ed = data.table(jab$edges)[cn>0, .(from, to , cn)]

    if (nrow(ed)==0)
        return(GRangesList())
    
    ed[, ":="(fromss = tile.map[ .(from), tile.id],
              toss = tile.map[ .(to), tile.id]),
       by = 1:nrow(ed)]
    ed[, weight :=  adj[cbind(from, to)]]
    print(ed)
    ed[fromss*toss > 0, eclass := ifelse(fromss>0, paste(fromss, toss), paste(-toss, -fromss))]
    ed[fromss*toss < 0, eclass := ifelse(abs(fromss)<=abs(toss),
                                         paste(fromss, toss), paste(-toss, -fromss))]
    ed[, eclass := as.numeric(as.factor(eclass))]
    ed[, eid := paste(from, to)]
    setkey(ed, "eid")
    eclass.cn = ed[!duplicated(eclass), setNames(cn, eclass)]

    cleanup_mode = FALSE

    
    while (nrow(ij)>0)
    {
        if (verbose)
            message('Path peeling iteration ', i, ' with ', sum(adj!=0, na.rm = TRUE), ' edges left and ', nrow(ij), ' end-pairs to resolve' )
        i = i+1
        ## swap this
        ##        vpaths[[i]] = p = as.numeric(get.shortest.paths(G, ij[1, 1], ij[1, 2], mode = 'out', weight = E(G)$weight)$vpath[[1]])


        p = get.constrained.shortest.path(cn.adj, G, v = ij[1, 1], to = ij[1, 2], weight = E(G)$weight, edges = ed, verbose = TRUE, mip = cleanup_mode)
         
        if (is.null(p)){
            message('Came up empty!')
            i = i -1
            ij = ij[-1, , drop = FALSE]
        }
        else
        {
            ## Don't forget to update ed here
            ed$cn = cn.adj[cbind(ed$from, ed$to)]
            
            vpaths[[i]] = p
            epaths[[i]] = cbind(p[-length(p)], p[-1])                          
            eids = paste(epaths[[i]][,1], epaths[[i]][,2])            
            cns[i] = ed[.(eids), if (length(cn)>1) cn/2 else cn, by = eclass][, floor(min(V1))] ## update cn correctly, adjusting constraints for palinrdomic edges by 1/2
                        
            rvpath = rtile.map[list(tile.map[list(vpaths[[i]]), -rev(tile.id)]), id]
            repath = cbind(rvpath[-length(rvpath)], rvpath[-1])
            plen = length(rvpath)
            hplen = floor(length(rvpath)/2)

            ## (awkward) check for palindromicity for odd and even length palindromes
            ## if (all((vpaths[[i]]==rvpath)[c(1:hplen,(plen-hplen+1):plen)]))
            if (ed[eids, any(table(eclass)>1)])
                palindromic.path[i] = TRUE         
            ## else
            ## {
            vpaths[[i+1]] = rvpath
            epaths[[i+1]] = repath
            cns[i+1] = cns[i]
            palindromic.path[i+1] = TRUE
            ## }
            ##        palindromic = TRUE ## set to true while we "figure things out"


            #' so now we want to subtract that cn units of that path from the graph
            #' so we want to update the current adjacency matrix to remove that path
            #' while keeping track of of the paths on the stack
            cn.adj[epaths[[i]]] = cn.adj[epaths[[i]]]-cns[i]
            
            ## if (!palindromic) ## update reverse complement unless palindromic
            cn.adj[epaths[[i+1]]] = cn.adj[epaths[[i+1]]]-cns[i+1]
            
            if (!all(cn.adj[epaths[[i]]]>=0)) ## something wrong, backtrack
            {
                message('backtracking ...') ## maybe we got stuck in a quasi-palindrome and need to backtrack
                                        #            browser()
                cn.adj[epaths[[i]]] = cn.adj[epaths[[i]]]+cns[i]           
                ## if (!palindromic) ## update reverse complement unless palindromic
                cn.adj[epaths[[i+1]]] = cn.adj[epaths[[i+1]]]+cns[i+1]
                i = i-1
                ij = ij[-1, , drop = FALSE]
            }
            else ## continue, reduce
            {
                adj.new[epaths[[i]]] = adj.new[epaths[[i]]] + cns[i]
                ## if (!palindromic)
                adj.new[epaths[[i+1]]] = adj.new[epaths[[i+1]]] + cns[i]

                ## ## make sure I didn't overuse any edge
                ## if (nrow(overdue <- which((as.matrix(jab$adj)-adj.new)<0, arr.ind=T))>0) {
                ##     print("Edge copy deficit!")
                ##     browser()
                ## }

                ## intermediate check
                ## if (length(which(((adj.new + cn.adj) - jab$adj)!=0, arr.ind = TRUE)))
                ##     browser()
                
                to.rm = epaths[[i]][which(cn.adj[epaths[[i]]]==0), ,drop = FALSE]
                ## if (!palindromic) ## update reverse complement
                to.rm = rbind(to.rm, epaths[[i+1]][which(cn.adj[epaths[[i+1]]]==0), ,drop = FALSE])
                
                if (nrow(to.rm)>0)
                {
                    adj[to.rm] = 0
                    ## ALERT!!! major change
                    ## adjj = adj/as.matrix(cn.adj)
                    ## adjj[which(is.nan(adjj))] = 0
                    ## adjj[which(adjj<0)] = 0
                    G = graph.adjacency(adj, weighted = 'weight')
                    ## G = graph.adjacency(adjj, weighted = 'weight')
                    new.ends = setdiff(which(
                    (degree(G, mode = 'out')==0 | degree(G, mode = 'in')==0)
                    & degree(G)>0), ends)

                    ## ## check if cn.adj out of balance
                    ## if (any((colSums(cn.adj)*rowSums(cn.adj) != 0) & (colSums(cn.adj) != rowSums(cn.adj)))){
                    ##     print("Junction OUT OF BALANCE!")
                    ##     browser()
                    ## }

                    ## ## should be no new ends
                    ## if (length(new.ends)>0){
                    ##     print("Please, no new ends!")
                    ##     browser()
                    ## }

                    ## remain = as.matrix(jab$adj) - adj.new
                    ## nb <- which(colSums(remain) != rowSums(remain))
                    ## if (any(!is.element(nb, nb.all)))
                    ##     browser()
                    
                    D = shortest.paths(G, v = ends, mode = 'out', weight = E(G)$weight)[, ends]
                    ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row != col, ][order(dist), ][, row := ends[row]][, col := ends[col]]                    
                }
                else
                    ij = ij[-1, , drop = FALSE]

                ## if (!palindromic) ## increase extra counter to account for reverse complement
                ## TOFIX: just update counter by 2 above, since we are just doing every path and its rc 
                i = i+1            
            }
        }


        ## DEBUG DEBUG DEBUG
        seg.ix = which(strand(jab$segstats)=='+'); seg.rix = which(strand(jab$segstats)=='-');
        
        
        if (nrow(ij)==0 & cleanup_mode == FALSE)
        {
            message('!!!!!!!!!!!!!!!!!!!!!!!!!!STARTING CLEANUP MODE FOR PATHS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row != col, ][order(dist), ][, row := ends[row]][, col := ends[col]]
            cleanup_mode = TRUE
        }
    }
    if (verbose)
        message('Path peeling iteration ', i, ' with ', sum(adj!=0, na.rm = TRUE), ' edges left ', nrow(ij) )

    ## ## record G, D, remaining edges at the end of path peeling
    ## G1 = G
    ## D1 = D
    ## remain1 = remain
    
    vpaths = vpaths[1:i]
    epaths = epaths[1:i]
    cns = cns[1:i]
    palindromic.path = palindromic.path[1:i]

    vcycles = rep(list(NA), maxrow) 
    ecycles = rep(list(NA), maxrow)
    ccns = rep(NA, maxrow)

    csimp = which(diag(cn.adj)!=0)
    ipath = i
    i = 0
    if (length(csimp)>0)
    {
        vcycles[1:length(csimp)] = split(csimp, 1:length(csimp))
        ecycles[1:length(csimp)] = lapply(csimp, function(x) cbind(x, x))
        ccns[1:length(csimp)] = diag(cn.adj)[csimp]
        cn.adj[cbind(csimp, csimp)] = 0
        adj[cbind(csimp, csimp)] = 0
        i = length(csimp)
        
        for (j in 1:length(csimp))
            adj.new[ecycles[[j]]] = adj.new[ecycles[[j]]] + ccns[j]
    }
    
    ## sort shortest paths and find which connect a node to its ancestor (i.e. is a cycle)
    .parents = function(adj)
    {
        tmp = apply(adj, 2, function(x) which(x!=0))
        ix = which(sapply(tmp, length)>0)
        if (length(ix)>0)
        {
            parents = rbindlist(lapply(ix, function(x) data.table(x, tmp[[x]])))
            setnames(parents, c('node', 'parent'))
            setkey(parents, node)
        } else {
            parents = data.table(node = 0, parent = NA)
            setkey(parents, node)
        }
    }

    parents = .parents(adj)
    
    #' then find paths that begin at a node and end at (one of its) immediate upstream neighbors
    #' this will be a path for whom col index is = parent(row) for one of the rows
    ## ALERT!!! major change
    ## adjj = adj/as.matrix(cn.adj)
    ## adjj[which(is.nan(adjj))] = 0
    ## adjj[which(adjj<0)] = 0
    G = graph.adjacency(adj, weighted = 'weight')
    ## G = graph.adjacency(adjj, weighted = 'weight')
    D = shortest.paths(G, mode = 'out', weight = E(G)$weight)

    ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row %in% parents$parent & row != col, ][order(dist), ][, is.cycle := parents[list(row), col %in% parent], by = row][is.cycle == TRUE, ]
    

    ## now iterate from shortest to longest path
    ## peel that path off and see if it is still there ..
    ## and see if it is still there 

    ## peel off top path and add to stack, then update cn.adj

    cleanup_mode = FALSE
    while (nrow(ij)>0)
    {
        if (verbose)
            message('Cycle-peeling iteration ', i, ' with ', sum(adj!=0, na.rm = TRUE), ' edges left ', nrow(ij) )
        i = i+1
                                        #        p = as.numeric(get.shortest.paths(G, ij[1, 1], ij[1, 2], mode = 'out', weight = E(G)$weight)$vpath[[1]])
        
        p = get.constrained.shortest.path(cn.adj, G, allD = D, v = ij[1, 1], to = ij[1, 2], weight = E(G)$weight, edges = ed, verbose = TRUE, mip = cleanup_mode)

        if (is.null(p)){
            message('Came up empty!')
            i = i -1
            ij = ij[-1, , drop = FALSE]
        } else
        {        

            ed$cn = cn.adj[cbind(ed$from, ed$to)]
            vcycles[[i]] = p
            ecycles[[i]] = cbind(p, c(p[-1], p[1]))
            eids = paste(ecycles[[i]][,1], ecycles[[i]][,2])
            ccns[i] = ed[.(eids), if (length(cn)>1) cn/2 else cn, by = eclass][, floor(min(V1))] ## update cn correctly, adjusting constraints for palindromic edges by 1/2
            
            rvcycle = rtile.map[list(tile.map[list(vcycles[[i]]), -rev(tile.id)]), id]
            recycle = cbind(rvcycle, c(rvcycle[-1], rvcycle[1]))
            clen = length(rvcycle)
            hclen = floor(length(rvcycle)/2)
            ## (awkward) check for palindromicity for odd and even length palindromes

            ## if (all((vcycles[[i]]==rvcycle)[c(1:hclen,(clen-hclen+1):clen)]))
            if (ed[eids, any(table(eclass)>1)])
                palindromic.cycle[i] = TRUE
            ## else
            ## {
            vcycles[[i+1]] = rvcycle
            ecycles[[i+1]] = recycle
            ccns[i+1] = ccns[i]
            palindromic.cycle[i+1] = TRUE
            ##     palindromic = FALSE
            ## }
            ##        palindromic = TRUE ## set to true while we "figure things out"
            
            #' so now we want to subtract that cn units of that path from the graph
            #' so we want to update the current adjacency matrix to remove that path
            #' while keeping track of of the cycles on the stack        
            cn.adj[ecycles[[i]]] = cn.adj[ecycles[[i]]]-ccns[i]
            ## if (!palindromic) ## update reverse complement unless palindromic
            cn.adj[ecycles[[i+1]]] = cn.adj[ecycles[[i+1]]]-ccns[i+1]
            
            if (!all(cn.adj[ecycles[[i]]]>=0))
            {
                message('backtracking')
                ## browser()
                cn.adj[ecycles[[i]]] = cn.adj[ecycles[[i]]]+ccns[i]           
                ## if (!palindromic) ## update reverse complement unless palindromic
                cn.adj[ecycles[[i+1]]] = cn.adj[ecycles[[i+1]]]+ccns[i+1]
                i = i-1
                ij = ij[-1, , drop = FALSE]                
            }
            else
            {
                adj.new[ecycles[[i]]] = adj.new[ecycles[[i]]] + ccns[i]
                
                ## ## if (!palindromic)
                ##     adj.new[ecycles[[i+1]]] = adj.new[ecycles[[i+1]]] + ccns[i]

                ## ## ## make sure I didn't overuse any edge
                ## ## if (length(overdue <- which((as.matrix(jab$adj)-adj.new)<0))) {
                ## ##     print("Edge copy deficit!")
                ## ##     browser()
                ## ## }

                ## ## ## intermediate cross check
                ## ## if (length(which(((adj.new + cn.adj) - jab$adj)!=0, arr.ind = TRUE)))
                ## ##     browser()
                
                to.rm = ecycles[[i]][which(cn.adj[ecycles[[i]]]==0), ,drop = FALSE]
                
                ## if (!palindromic) ## update reverse complement
                to.rm = rbind(to.rm, ecycles[[i+1]][which(cn.adj[ecycles[[i+1]]]==0), ,drop = FALSE])
                
                if (nrow(to.rm)>0)
                {
                    adj[to.rm] = 0
                    parents = .parents(adj)
                    ## G = graph.adjacency(adj, weighted = 'weight')

                    ## ALERT!!! major change
                    ## adjj = adj/as.matrix(cn.adj)
                    ## adjj[which(is.nan(adjj))] = 0
                    ## adjj[which(adjj<0)] = 0
                    G = graph.adjacency(adj, weighted = 'weight')
                    ## G = graph.adjacency(adjj, weighted = 'weight')
                    
                    ## if (any((colSums(cn.adj)*rowSums(cn.adj) != 0) & (colSums(cn.adj) != rowSums(cn.adj)))){
                    ##     print("Junction OUT OF BALANCE!")
                    ##     browser()
                    ## }

                    ## remain = as.matrix(jab$adj) - adj.new
                    ## nb <- which(colSums(remain) != rowSums(remain))
                    ## if (any(!is.element(nb, nb.all)))
                    ##     browser()
                    
                    D = shortest.paths(G, mode = 'out', weight = E(G)$weight)
                    ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row %in% parents$parent & row != col, ][order(dist), ][, is.cycle := parents[list(row), col %in% parent], by = row][is.cycle == TRUE, ]
                }
                else
                    ij = ij[-1, ,drop = FALSE]
                
                ## if (!palindromic) ## increase extra counter to account for reverse complement
                i = i+1
            }
        }
        
        if (nrow(ij)==0 & cleanup_mode == FALSE)
        {
            message('!!!!!!!!!!!!!!!!!!!!!!!!!!STARTING CLEANUP MODE FOR CYCLES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row %in% parents$parent & row != col, ][order(dist), ][, is.cycle := parents[list(row), col %in% parent], by = row][is.cycle == TRUE, ]
            
            cleanup_mode = TRUE
        }                        
    }
    
    if (verbose)
        message('Cycle peeling iteration ', i, ' with ', sum(adj!=0, na.rm = TRUE), ' edges left ', nrow(ij) )
    
    
    if (i>0)
    {
        vcycles = vcycles[1:i]
        ecycles = ecycles[1:i]
        ccns = ccns[1:i]
    }
    else
    {
        vcycles = NULL
        ecycles = NULL
        ccns = NULL
    }
    
    vall = c(vpaths, vcycles)
    eall = c(epaths, ecycles)
    ecn = c(cns, ccns)

    ## ## record G, D, remaining edges at the end of cycle peeling
    ## G2 = G
    ## D2 = D
    ## remain2 = remain
    remain = as.matrix(jab$adj) - adj.new
    remain.ends = which(colSums(remain)*rowSums(remain)==0 & colSums(remain)-rowSums(remain)!=0)
    if (length(remain.ends)>0){
        if (verbose)
            message(length(remain.ends), "ends were not properly assigned a path. Do them.")
    }
    
    tmp = cbind(do.call(rbind, eall), rep(ecn, sapply(eall, nrow)), munlist(eall))
    ix = which(rowSums(is.na(tmp[, 1:2]))==0)

    if (length(ix)>0)
        adj.new = sparseMatrix(tmp[ix,1], tmp[ix,2], x = tmp[ix,3], dims = dim(adj))
    else
        adj.new = sparseMatrix(1, 1, x = 0, dims = dim(adj))
    vix = munlist(vall)

    jab$segstats$node.id = 1:length(jab$segstats)
    pathsegs = jab$segstats[vix[,3]]
    pathsegs$grl.ix = vix[,1]    
    abjuncs =  as.data.table(rbind(jab$ab.edges[, 1:2, '+'], jab$ab.edges[, 1:2, '-']))[, 
                                   id := rep(1:nrow(jab$ab.edges),2)*
                                       rep(c(1, -1), each = nrow(jab$ab.edges))][!is.na(from), ]    
    abjuncs = abjuncs[, tag := structure(paste(from, to), names = id)]
    setkey(abjuncs, tag)           

    pathsegs$ab.id = gr2dt(pathsegs)[ , .(ab.id = c(abjuncs[paste(node.id[-length(node.id)], node.id[-1]), id], NA)), by = grl.ix][, ab.id]

    paths = split(pathsegs, vix[,1])        
    values(paths)$ogid = 1:length(paths)
    values(paths)$cn = ecn[as.numeric(names(paths))]
    values(paths)$label = paste('CN=', ecn[as.numeric(names(paths))], sep = '')
    values(paths)$is.cycle = !(as.numeric(names(paths)) %in% 1:length(vpaths))
    values(paths)$numsegs = elementNROWS(paths)
    values(paths)$num.ab = sapply(paths, function(x) sum(!is.na(x$ab.id)))
    values(paths)$wid = sapply(lapply(paths, width), sum)

    check = which((adj.new - jab$adj) !=0, arr.ind = TRUE)
    
    if (length(check)>0)
        stop('Alleles do not add up to marginal copy number profile!')
    else if (verbose)
        message('Cross check successful: sum of walk copy numbers = marginal JaBbA edge set!')

    ## match up paths and their reverse complements
    psig = lapply(paths, function(x) ifelse(as.logical(strand(x)=='+'), 1, -1)*x$tile.id)
    psig.flip = sapply(psig, function(x) -rev(x))

    unmix = data.table(
        ix = 1:length(paths),
        mix = match(sapply(psig, paste, collapse = ','), sapply(psig.flip, paste, collapse = ',')))[, pos := 1:length(mix)<mix][order(!pos), ]
    setkey(unmix, ix)
    unmix[is.na(mix), pos := TRUE] ## if we have paths with no reverse complement i.e. NA mix, then call "+" for now

    remix = rbind(
        unmix[pos == TRUE, ][, id := 1:length(ix)],
        unmix[list(unmix[pos == TRUE, mix]), ][, id := 1:length(ix)][!is.na(ix), ]
    )
    
    paths = paths[remix$ix]
    names(paths) = paste(remix$id, ifelse(remix$pos, '+', '-'), sep = '')
    values(paths)$id = remix$id
    values(paths)$str = ifelse(remix$pos, '+', '-')

    if (length(setdiff(values(paths)$ogid, 1:length(paths))))
        message('Warning!!! Some paths missing!')
    
    return(paths)
}



#' @name jabba.kid  
#' @title jabba.kid
#' @description
#'
#' Given a set of gwalks (grangeslist outputs of jabba.gwalk) identifies strings "kidnapped" fragments i.e.
#' strings of tempaled insertions, which are outputted as a grangeslist
#' 
#' @export
#' @param gwalks grangeslist of walks (e.g.outputted from jabba.gwalks)
#' @param pad how much bp neighboring material around each kidnapped string to include in outputs
#' @param min.ab minimal bp distance a junction needs to be considered aberrant (e.g. to exclude very local deletions)
#' @param min.run how many aberrant junctions to  require in the outputted kidnapped fragments
#' @export
#' @author Marcin Imielinski
jabba.kid = function(gwalks, pad = 5e5, min.ab = 5e5, min.run = 2)
{
    gw = gr2dt(grl.unlist(gwalks))    
    gw[, ab.chunk := cumsum(!is.na(ab.id)),  by = grl.ix]
    gw[, dist := c(ifelse((seqnames[-1] != seqnames[-length(seqnames)]) |
                          (strand[-1] != strand[-length(strand)]), Inf,
                   ifelse(strand[-1]=='+',
                          start[-1]-end[-length(end)],
                          start[-length(end)]-end[-1]))
                 , Inf), by = grl.ix]
    gw[, dist.nostrand := c(ifelse((seqnames[-1] != seqnames[-length(seqnames)]), Inf,
                   ifelse(strand[-1]=='+',
                          start[-1]-end[-length(end)],
                          start[-length(end)]-end[-1]))
                 , Inf), by = grl.ix]
    
    
    gw[, dist := ifelse((1:length(grl.iix) %in% length(grl.iix)), as.numeric(NA), dist), by = grl.ix]
    gw[, dist.nostrand := ifelse((1:length(grl.iix) %in% length(grl.iix)), as.numeric(NA), dist),
       by = grl.ix]

    ## get rid of little dels
    gw[which(dist<min.ab & dist>0), ab.id := NA]

    gwu = gw[, .(wid = sum(width), ab.id = ab.id[!is.na(ab.id)][1],
                           start = grl.iix[1], end = grl.iix[length(grl.iix)]), by = .(ab.chunk, grl.ix)]

    .labrun = function(x) ifelse(x, cumsum(diff(as.numeric(c(FALSE, x)))>0), NA)    
    
    ## every run of "trues" i.e. wid<something is labeled
    gwu[, runtag := as.numeric(.labrun(wid<pad)), by = grl.ix]
    gwu = gwu[!is.na(runtag), ]
    gwu[, runlab := paste(grl.ix, runtag), by = grl.ix]
    gwu[, runlen := sum(!is.na(ab.id)), by = .(grl.ix, runlab)]
    

    ## choose only chunk s with min run of abs
    kidnapped = gwu[!is.na(runlab) & runlen>=min.run, ]

    if (nrow(kidnapped) == 0)
        return(GRangesList())


    ## expand chunks back out
    kidnapped[, first := (1:length(start))==1, by = .(runlab)]
    kidnapped[, last := (1:length(start))==length(start), by = .(runlab)]
    ix = kidnapped[, .(runlab, frag.id = paste(grl.ix, ab.chunk), 
                       grl.iix = unique(c(start-first, start:end, end+last))), by = .(ab.chunk, grl.ix)]
    setkeyv(ix, c('grl.ix', 'grl.iix'))
    setkeyv(gw, c('grl.ix', 'grl.iix'))
    
    kidnapped = gw[ix, ][!is.na(start),]
    kidnapped[, frag.iid := 1:length(grl.iix), by = frag.id]
    kidnapped[, frag.pos := cumsum(width), by = frag.id]

    kidnapped[, first := (1:length(grl.iix)) == 1, by = runlab]
    kidnapped[, last := (1:length(grl.iix)) == length(grl.iix), by = runlab]
    
    gr.kn = dt2gr(kidnapped)

    gr.kn$width = NULL
    ## trim ends leading to and out of segment
    gr.kn[gr.kn$first] = gr.end(gr.kn[gr.kn$first], pad, ignore.strand = FALSE)
    gr.kn[gr.kn$last] = gr.start(gr.kn[gr.kn$last], pad, ignore.strand = FALSE)

    setkey(kidnapped, frag.id)

    walks.kn = split(gr.kn[, c('ab.id','grl.iix', 'cn', 'cn.1', 'frag.id', 'frag.iid', 'frag.pos', 'dist', 'dist.nostrand')], gr.kn$runlab)

    kidnapped$runlab = as.character(kidnapped$runlab)
    setkey(kidnapped, 'runlab')    
    values(walks.kn) = as.data.frame(kidnapped[ ,.(pair = pair[1], grl.ix = grl.ix[1], 
                                                   len = length(setdiff(unique(ab.id), NA))), keyby = runlab][names(walks.kn), ])
    return(walks.kn)
}


#' @name anno.hops
#' @title anno.hops
#' @description
#'
#' Adds simple annotations to GRangesList of walks including
#' distance along each reference fragment and distance
#' between "hops"
#' 
#' @export
#' @param walks walks to annotate
#' 
#' @author Marcin Imielinski
anno.hop = function(walks)
{
  gw = gr2dt(grl.unlist(walks))
  gw[, ab.chunk := cumsum(!is.na(ab.id)),  by = grl.ix]
  gw[, dist := c(ifelse((seqnames[-1] != seqnames[-length(seqnames)]) |
                        (strand[-1] != strand[-length(strand)]), Inf,
                 ifelse(strand[-1]=='+',
                        start[-1]-end[-length(end)],
                        start[-length(end)]-end[-1]))
               , Inf), by = grl.ix]

  gw[, dist.nostrand := c(ifelse((seqnames[-1] != seqnames[-length(seqnames)]), Inf,
                          ifelse(strand[-1]=='+',
                                 start[-1]-end[-length(end)],
                                 start[-length(end)]-end[-1]))
                        , Inf), by = grl.ix]

  gw[, dist := ifelse((1:length(grl.iix) %in% length(grl.iix)), as.numeric(NA), dist), by = grl.ix]
  gw[, dist.nostrand := ifelse((1:length(grl.iix) %in% length(grl.iix)), as.numeric(NA), dist),
     by = grl.ix]

  gw[, ":="(frag.id = paste(grl.ix, ab.chunk), 
            frag.iid = 1:length(grl.ix),
            frag.pos = cumsum(width)
            ), by = .(ab.chunk, grl.ix)]

  gr.out = dt2gr(gw)
  gr.out$width = NULL

  setkey(gw, frag.id)

  grl.out = split(gr.out[, c('ab.id','grl.iix', 'cn', 'cn.1', 'frag.id', 'frag.iid', 'frag.pos', 'dist', 'dist.nostrand')], gr.out$grl.ix)[as.character(1:length(walks))]

  names(grl.out) = names(walks)
  values(grl.out) = values(walks)

    return(grl.out)
}

####################################################
#' jabba.walk
#'
#' Computes walks around all aberrant edges in JABbA object
#' 
#' Takes in JaBbA solution and computes local
#' reconstructions around all aberrant edges (default).  Reconstructions (i.e. Huts) consists
#' of collections of walks, each walk associated with a copy number, and a given
#' region (collection of genomic windows).  The interval sum of walks in a given region, weighted
#' by copy numbers will recapitulate the marginal copy profile (as estimated by JaBbA).
#' The reconstruction is chosen to maximize parsimony.
#' 
#' Optional flags allow making huts around specific junctions or specified loci (GRangesList)
#' 
#' Walks are reconstructed locally within "clustersize" nodes of each aberrant edge, where
#' clustersize is measured by the number of total edges.  Larger cluster sizes may fail to be
#' computationally tractable, i.e. with a highly rearranged genome in an area of dense interconnectivity. 
#'
#' @param sol JaBbA object
#' @param outdir output directory
#' @param junction.ix junction indices around which to build walks (default is all junctions)
#' @param loci  loci around which to build walks (over-rides junction.ix), alternatively can be a list of  "all.paths" objects (i.e. each a list utput of initial all.paths = TRUE run  +/- field $prior for walk to re-eval a given all.paths combo
#' @param clustersize size of the cluster to output around the locus or junction of interest
#' @param trim logical flag whether trim in neighborhood of junction (only applicable if loci = NULL, default = TRUE)
#' @param trim.w integer width to which to trim to
#' @param prune flag whether to prune trivial walks for whom a path can be drawn from first to last interval in a graph linking intervals with pairwise distance < d1 on the walk or distance < d2 on the reference
#' @param prune.d1 local distance threshold for walk pruning 
#' @param prune.d2 referenc distance threshold for walk pruning 
#' @param mc.cores number of cores to use, default 1
#' @param genes character vector of gene symbols with which to annotate walk (eg cancer genes)
#' @param verbose logical flag
#' @return list of walk set around each locus or junction that is inputted to analysis, each list item is a list with the following fields
#' $win = input locus of interest, $grl = GRangesList of walks, $grs is a collapsed footprint of all walks in the walk list for this locu
#' $td gTrack of of the output, additional outputs for debugging: $sol, $K, $Bc, $eix, $vix, $h
#' @export
####################################################
jabba.walk = function(sol, kag = NULL, digested = T, outdir = 'temp.walk', junction.ix = NULL, loci = NULL, clustersize = 100,
  trim = FALSE, ## whether to trim around junction (only applicable when loci = NULL)
  trim.w = 1e6, ## how far to trim in neighborhood of junction (only applicable when loci = NULL
  prune = FALSE, ## whether to prune trivial walks i.e. those for whom a path can be drawn from first to last interval in a graph linking intervals with pairwise distance < d1 on the walk or distance < d2 on the reference
  prune.d1 = 1e5, ## local distance threshold for walk pruning 
  prune.d2 = 1e5, ## reference distance threshold for walk pruning
  maxiterations = Inf, mc.cores = 1, genes = read.delim('~/DB/COSMIC/cancer_gene_census.tsv', strings = F)$Symbol, verbose = T, max.threads = 4, customparams = T, mem = 6, all.paths = FALSE, nomip = F, tilim = 100, nsolutions = 100, cb.interval = 1e4, cb.chunksize = 1e4, cb.maxchunks = 1e10)
{
  system(paste('mkdir -p', outdir))
  ## awkward workaround to limit the number of processors Cplex will gobble up
  ##
  
  if (customparams)
    {
      out.file = paste(outdir, 'tmp.prm', sep = '/')
      max.threads = Sys.getenv("LSB_DJOB_NUMPROC")
      if (nchar(max.threads) == 0)
        max.threads = Inf
      else
        max.threads = as.numeric(max.threads)
      max.threads = min(max.threads, mc.cores)
      if (is.infinite(max.threads))
        max.threads = 0 
      
      param.file = paste(out.file, '.prm', sep = '')
      .cplex_customparams(param.file, max.threads, treememlim = mem * 1e3)
      
      Sys.setenv(ILOG_CPLEX_PARAMETER_FILE = normalizePath(param.file))
      print(Sys.getenv('ILOG_CPLEX_PARAMETER_FILE'))
    }    


   if (is.null(sol))
      sol = kag

  if (is.null(sol$segstats))
      {
          sol$segstats = sol$tile
          sol$segstats$cn = 2
          sol$segstats$eslack.out = 0
          sol$segstats$eslack.in = 0
      }
  
  if (is.null(kag))
      kag = sol
  
  
  out = list()
  tmp.adj = sol$adj
  if (digested)  ## if input is already "digested", then don't need to bother with slacks
      {
      sol$segstats$eslack.in = 0
      sol$segstats$eslack.out = 0
      G = sol$G
    }
  else ## soon to be deprecated
    {
      ix = which(sol$segstats$eslack.in!=0 | sol$segstats$eslack.out!=0)
      tmp.adj[ix, ix] = 0
      pos.ix = which(strand(sol$segstats)=='+')
      sol$segstats$tile.id = match(gr.stripstrand(sol$segstats), gr.stripstrand(sol$segstats[pos.ix]))
      G = graph.adjacency(tmp.adj!=0)      
    }

  if (verbose)
    message(paste('Processing JaBbA'))

  h = jbaMIP.process(sol)
  
  if (!is.null(genes))
    td.rg = track.gencode(genes = genes, height = 3)

  if (is.null(junction.ix) & is.null(loci))
    junction.ix = 1:nrow(kag$ab.edges)

  if (!is.null(junction.ix))
    if (is.null(names(junction.ix)))
      names(junction.ix) = 1:length(junction.ix)

  if (is.null(loci)) ## junction.ix should be not null here 
    {
      loci = do.call('GRangesList', mclapply(junction.ix, function(i)
        {
             if (verbose)
               cat(paste('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nDefining subgraph around junction', i, '\n'))
             vix = vix.i = setdiff(kag$ab.edges[i, 1:2, ], NA)
             if (length(vix)==0)
                 return(GRanges())
             k = 0
             last.clustersize = 0
             while (length(vix)<clustersize & k < maxiterations & length(vix)>last.clustersize)
               {
                 k = k + 1
                 last.clustersize = length(vix)
                 vix = unique(unlist(neighborhood(G, vix.i, order = k)))
               }
             if (verbose)
               cat(paste('Outputting', length(vix), 'vertices around junction', i, '\n'))
             
             return(kag$segstats[vix])
           }        
        , mc.cores = mc.cores))
      
      names(loci) = names(junction.ix)
      loci = loci[sapply(loci, length)>0]
    }
  else ## if loci are provided (i.e. not junction centric) then we will not trim or prune
    {
      trim = F 
      prune = F
    }

  if (verbose)
    cat(paste('Finished defining subgraphs\n'))
  
  starts = gr.start(sol$segstats, ignore.strand = F)
  ends = gr.end(sol$segstats, ignore.strand = F)

  names(sol$segstats) = 1:length(sol$segstats)

  if (is.null(names(loci)))
    lnames =  paste('locus', 1:length(loci), sep = '')
  else
    lnames = names(loci)

  all.junc.pair = c(paste(sol$ab.edges[, 1, 1], sol$ab.edges[, 2, 1], sep = ','), paste(sol$ab.edges[, 1, 2], sol$ab.edges[, 2, 2], sep = ','))
  names(all.junc.pair) = c(1:nrow(sol$ab.edges), -c(1:nrow(sol$ab.edges)))
  
  if (length(loci)>0)
    {
      out = mclapply(1:length(loci), function(i)
          {
            label = lnames[i]
            mprior = NULL
              outfile.rds = sprintf('%s/%s.rds', outdir, label)
              outfile.pdf = sprintf('%s/%s.pdf', outdir, label)
              outfile.txt = sprintf('%s/%s.txt', outdir, label)
              outfile.allpaths.txt = sprintf('%s/%s.allpaths.txt', outdir, label)
              if (is(loci[[i]], 'GRanges'))
                  {
                      vix = which(gr.in(kag$segstats, loci[[i]]))
                      cat('Number of vertices:', length(vix), '\n')
                      eix = which((h$e.ij[,1] %in% vix | h$e.ij[,2] %in% vix) & h$e>0)
                      Bc = as.matrix(h$B)[vix, eix]
                      K = tryCatch(convex.basis(Bc, interval = cb.interval, chunksize = cb.chunksize, verbose = T, maxchunks = cb.maxchunks), error = function(e) as.character(e))
                      if (is.character(K))
                          return(list(README = K))
                      prior = rep(1, ncol(K))
                  }
              else ## assume we are re-heating a previous all.paths = TRUE output (and presumably adding a prior)
                  {
                      K = loci[[i]]$K
                      h = loci[[i]]$h
                      eix = loci[[i]]$eix
                      Bc = loci[[i]]$Bc
                      vix = loci[[i]]$vix
                      prior = rep(1, ncol(K))

                      if (is.matrix(loci[[i]]$mprior))
                      {
                        if (verbose)
                          cat(paste('Adding a matrix prior!!!!!!\n'))

                        ## initialize matrix 
                        mprior = array(0, dim = c(nrow(loci[[i]]$mprior), ncol(K)))
                        loci
                        mprior[, values(loci[[i]]$allpaths.og)$kix] = loci[[i]]$mprior
                        mprior[, values(loci[[i]]$allpaths.og)$kix2] = loci[[i]]$mprior
                        colnames(mprior) = 1:ncol(mprior)
                        colnames(mprior)[values(loci[[i]]$allpaths.og)$kix] = as.character(1:ncol(mprior))
                        colnames(mprior)[values(loci[[i]]$allpaths.og)$kix2] = as.character(-(1:ncol(mprior)))
                      }

                      if (!is.null(loci[[i]]$prior))
                      {
                        if (verbose)
                          cat(paste('Adding a prior!!!!!!\n'))
                        prior[c(values(loci[[i]]$allpaths.og)$kix,values(loci[[i]]$allpaths.og)$kix2)]  = loci[[i]]$prior
                      }
                      loci[[i]] = loci[[i]]$win
                  }
              
              is.cyc = Matrix::colSums(K[h$etype[eix] == 'slack', ])==0 & Matrix::colSums((h$B[, eix, drop = F] %*% K)!=0)==0
              
              karyo.sol = karyoMIP(K, h$e[eix], h$eclass[eix], nsolutions = nsolutions, tilim = tilim, cpenalty = 1/prior, mprior = mprior)
              kag.sol = karyo.sol[[1]]
          p = karyoMIP.to.path(kag.sol, K, h$e.ij[eix, ], sol$segstats, mc.cores = pmin(4, mc.cores))
          values(p$grl)$cn = p$cn
          values(p$grl)$is.cyc = p$is.cyc
          td.rg$stack.gap = 5e6
          
          if (!is.null(kag$junctions))
            {
              values(kag$junctions)$lwd = sol$adj[kag$ab.edges[,1:2, 1]]
              values(kag$junctions)$lty = 1
              values(kag$junctions)$label = ifelse(sol$adj[kag$ab.edges[,1:2, 1]]>0, sol$adj[kag$ab.edges[,1:2, 1]], '')
              values(kag$junctions)$col = ifelse(sol$adj[kag$ab.edges[,1:2, 1]]>0, alpha('red', 0.3), alpha('white', 0))
            }
          win = streduce(sol$segstats[vix], 1e4)
          
          y1 = max(sol$segstats$cn[gr.in(sol$segstats, win)], na.rm = T)*1.1
          pdf(outfile.pdf, height = 30, width = 24)
          grs = gr.simplify(grl.unlist(p$grl), 'grl.ix', split = T)
          values(grs) = values(p$grl)
          names(grs) = names(p$grl)

          if (!is.null(sol$td))
            {
              td.seg = sol$td
              td.seg$y1 = y1
              td = c(td.seg, td.rg)
            }
          else
              {
                  td.seg = gTrack(sol$segstats, y.field = 'cn', angle = 0, col ='black', height = 6, labels.suppress = T, y1 = y1)
                  
                                        #          td = c(gTrack(grs, draw.paths = T, path.cex.arrow = 0, border = NA, angle = 0, ywid = 0.5, path.stack.x.gap = 1e6, height = 20, labels.suppress.gr = T),
                                    
                  gt.walk = gTrack(grs, draw.paths = T, border = NA, angle = 0, ywid = 0.5, height = 20, labels.suppress.gr = T)
                  gt.walk$path.cex.arrow = 0
                  gt.walk$path.stack.x.gap = 1e6
                  td = c(
                      gt.walk,
                      td.seg,
                      td.rg)
                  plot(td,
                       windows = win, links = kag$junctions)
                  dev.off()
              }
          
          df = data.frame(label = label, cn = p$cn, walk = sapply(grs, function(x) paste(gr.string(x, mb = F), collapse = ',')), widths = sapply(grs, function(x) paste(width(x), collapse = ',')), width = sapply(grs, function(x) sum(width(x))), numpieces = sapply(grs, length), type = 'walk') 
          df = rbind(data.frame(label = label, cn = NA, walk = paste(gr.string(win, mb = F), collapse = ','), widths = paste(width(win), collapse = ','), width = sum(width(win)), type = 'window', numpieces = length(win)), df)
          write.tab(df, outfile.txt)
          out = list(
            win = win, grl = p$grl, grls = grs, td = td, sol = karyo.sol,
            K = K, Bc = Bc, eix = eix, vix = vix, h = h, 
            README = 'win=windows, grl = raw granges list corresponding to paths, grls = simplified granges list corresponding to paths, td = gTrack object plotting walks, sol = solution object from karyoMIP of local walks, K = incidence matrix input to karyomip, Bc = input to convex.basis, eix = eix input to karyomip, vix = vix input corresponding to rows of Bc, h = h input to karyomip')

          if (all.paths)
            {
              outfile.allpaths.pdf = sprintf('%s/%s.allpaths.pdf', outdir, label)
              
              if (verbose)
                cat('Generating all walks\n')

              ## repurpose karyoMIP.to.path to generate all paths using "fake solution" i.e. all 1 weights,  to karyoMIP as input
              pallp = karyoMIP.to.path(list(kcn = kag.sol$kcn*0 + 1, kclass = kag.sol$kclass), K, h$e.ij[eix, ], sol$segstats, mc.cores = pmin(4, mc.cores), verbose = verbose)
              allp = pallp$grl

              allps = gr.simplify(grl.unlist(allp), 'grl.ix', split = T)
              allps[values(allp)$is.cycle] = do.call('GRangesList', lapply(which(values(allp)$is.cycle), function(x) c(allps[[x]], allps[[x]])))
              allps.og = allps; ## save for later
              values(allps.og)$kix = pallp$kix
              values(allps.og)$kix2 = pallp$kix2

              ## text encoding of junctions
              if (!is.null(junction.ix))
                junc.pair = paste(sol$ab.edges[junction.ix[i], 1, ], sol$ab.edges[junction.ix[i], 2, ], sep = ',')
              
              if (trim | prune) ## junction.ix should be not null here (i.e. they were provided as input or loci = NULL)
                {
                  allps.u = grl.unlist(allps)
                  allps.u$ix.s = gr.match(gr.start(allps.u, ignore.strand = F), starts, ignore.strand = F)
                  allps.u$ix.e = gr.match(gr.end(allps.u, ignore.strand = F), ends, ignore.strand = F)
                  allps = split(allps.u, allps.u$grl.ix)
                  allps.ixs = split(allps.u$ix.s, allps.u$grl.ix) ## start indices of walk intervals in sol$segstats
                  allps.ixe = split(allps.u$ix.e, allps.u$grl.ix) ## end indices of walks intervals in sol$segstats
                  allps.w = split(width(allps.u), allps.u$grl.ix)
                  allps.endc = split(levapply(width(allps.u), by = list(allps.u$grl.ix), FUN = cumsum), allps.u$grl.ix)
                
                  if (trim) ## only include windows around the junction of interest
                    {
                      ## allps.ix.pairs tells us what junction indices are present in a walk collection                        
                      allps.ix.pairs = mapply(function(x,y) if (length(x)<=1) NULL else which(paste(x[-length(x)], y[-1], sep = ',') %in% junc.pair), allps.ixe, allps.ixs, SIMPLIFY = F)
                      ## first, which windows contain the junction

                      wix = which(sapply(allps.ix.pairs, length)>0)
                      allps = allps[wix]

                      if (length(allps)>0)
                        {
                          allps.ixs = allps.ixs[wix] ## start interval id of kth interval in ith walk
                          allps.ixe = allps.ixe[wix] ## end interval id of kth interval in ith walk
                          allps.endc = allps.endc[wix] ## end walk coordinate of kth interval in ith walk
                          allps.w = allps.w[wix]
                          allps.ix.pairs = allps.ix.pairs[wix]
                          
                          ## start window for trimming
                          values(allps)$allps.junc.first =
                            pmax(0, mapply(function(x, y) y[x[1]], allps.ix.pairs, allps.endc)) ## walk position of first junction
                          values(allps)$allps.junc.last =
                            pmax(0, mapply(function(x, y) y[x[length(x)]], allps.ix.pairs, allps.endc)) ## walk position of last junction
                          
                          ## check for any quasi-palindromic walks that contain both orientations of a junction
                          ## split each of these into two so we can maintain the width limit
                          pal.wix = which(values(allps)$allps.win.firstix != values(allps)$allps.win.lastix)
                          if (length(pal.wix)>0)
                            {                              
                              allps.dup = allps[pal.wix]
                              values(allps.dup)$allps.junc.first = values(allps)$allps.junc.last
                              allps = c(allps, allps.dup)
                              allps.endc = c(allps.endc, allps.endc[pal.wix])
                              allps.w = c(allps.w, allps.w[pal.wix])
                            }
                          
                          values(allps)$allps.win.first =
                            pmax(0, values(allps)$allps.junc.first - trim.w) ## walk coordinate of new window start
                          values(allps)$allps.win.last =
                            pmin(sapply(allps.endc, function(x) x[length(x)]), values(allps)$allps.junc.first + trim.w) ## walk coordinate of new window end
                          values(allps)$allps.win.firstix = ## first walk interval to trim to
                            mapply(function(x, y) setdiff(c(which(x>y)[1], 1), NA)[1], allps.endc, values(allps)$allps.win.first)
                          values(allps)$allps.win.lastix = ## last walk interval to trim to 
                            mapply(function(x, y) setdiff(c(which(x>y)[1], length(x)), NA)[1], allps.endc, values(allps)$allps.win.last)
                          values(allps)$allps.win.first.keep =
                            mapply(function(p,e,i) e[i] - p, values(allps)$allps.win.first, allps.endc, values(allps)$allps.win.firstix)
                          values(allps)$allps.win.last.keep =
                            mapply(function(p,e,i,w) w[i] - (e[i] - p), values(allps)$allps.win.last, allps.endc, values(allps)$allps.win.lastix, allps.w)                          
                          ## apply trimming
                          ## we are trimming walks so that they are within trim.w bases of junction                          
                          allps.u = grl.unlist(allps)
                          iix = mapply(function(x,y) y %in% values(allps)$allps.win.firstix[x]:values(allps)$allps.win.lastix[x], allps.u$grl.ix, allps.u$grl.iix)
                          allps.u = allps.u[iix]
                          allps.u$keep.end = mapply(function(x, y)
                            ifelse(y == values(allps)$allps.win.firstix[x], values(allps)$allps.win.first.keep[x], NA), allps.u$grl.ix, allps.u$grl.iix)
                          allps.u$keep.start = mapply(function(x, y)
                            ifelse(y == values(allps)$allps.win.lastix[x], values(allps)$allps.win.last.keep[x], NA), allps.u$grl.ix, allps.u$grl.iix)
                                                    
                          if (any(tmp.ix <- !is.na(allps.u$keep.start))) ## we keep the end of the first segment
                            allps.u[tmp.ix] = gr.start(allps.u[tmp.ix], allps.u$keep.start[tmp.ix], ignore.strand = F)

                          if (any(tmp.ix <- !is.na(allps.u$keep.end))) ## we keep the beginning of the last segment
                            allps.u[tmp.ix] = gr.end(allps.u[tmp.ix], allps.u$keep.end[tmp.ix], ignore.strand = F)
                          
                          ## if there are multiple walks with the same aberrant junction set, then pick the longest of these
                          
                          ## first need to find the aberrant walks in each set
                          ij = paste(allps.u$ix.e[-length(allps.u)], allps.u$ix.s[-1], sep = ',') ## indices of all walk adjacent interval pairs
                          names(ij) = 1:length(ij)
                          ij = ij[diff(allps.u$grl.ix)==0] ## only pick intra-walk interval pairs
                          ij.ix = names(all.junc.pair)[match(ij, all.junc.pair)]
                          ## then compute the width of each walk

                          allps = split(allps.u, allps.u$grl.ix)
                          ij.ix.l = split(ij.ix, allps.u$grl.ix[as.numeric(names(ij))])[names(allps)]
                          values(allps)$ab.junc = lapply(ij.ix.l, paste, collapse = ',')
                          values(allps)$wid = vaggregate(width(allps.u), by = list(allps.u$grl.ix), FUN = sum)[names(allps)]
                          ix.w = order(-values(allps)$wid) 
                          allps = allps[ix.w[which(!duplicated(values(allps)$ab.junc[ix.w]))]] ## only keep the longest non-duplicate walks
                        }
                    }
                  
                  ## now dedup and trim contigs to locus (mainly useful if loci was provided as argument)
                  if (length(allps)>0)
                    {
                      win = reduce(gr.stripstrand(loci[[i]]))
                      allps.u = grl.unlist(allps)
                      
                      ## trim to locus
                      ix = gr.match(allps.u, win)
                      allps.u = allps.u[!is.na(ix)]
                      ix = ix[!is.na(ix)]
                      start(allps.u) = pmax(start(allps.u), start(win)[ix])
                      end(allps.u) = pmin(end(allps.u), end(win)[ix])
                      
                      allps.u$ix.s = gr.match(gr.start(allps.u, ignore.strand = F), starts, ignore.strand = F)
                      allps.u$ix.e = gr.match(gr.end(allps.u, ignore.strand = F), ends, ignore.strand = F)
                      
                      ## remove dups
                      allps.ixs = split(allps.u$ix.s, allps.u$grl.ix) ## start indices of intervals
                      allps.ixe = split(allps.u$ix.e, allps.u$grl.ix) ## end indices of intervals
                      
                      allps.u = allps.u[allps.u$grl.ix %in% which(!duplicated(paste(sapply(allps.ixs, paste, collapse = ','), sapply(allps.ixe, paste, collapse = ','))))]
                      allps = split(allps.u, allps.u$grl.ix)
                    }
                  
                  
                  if (prune & length(allps)>0)
                    ## this is to prune pseudo-aberrant walks that basically consist of short insertions of non-reference
                    ## sequences in a big reference chunk
                    {
                      ## for each walk create graph of intervals by determining whether pair ij is BOTH near on the walk (<= d1)
                      ## and near on the refernce (<= d2)
                      allps.u = grl.unlist(allps)

                      ## what are the ij pairs we want to test from this collapsed list
                      ij = merge(cbind(i = 1:length(allps.u), ix = allps.u$grl.ix), cbind(j = 1:length(allps.u), ix = allps.u$grl.ix))[, c('i', 'j')]

                      tmp = levapply(width(allps.u), by = list(allps.u$grl.ix), FUN = cumsum)
                      allps.u.ir = IRanges(tmp - width(allps.u) + 1, tmp)

                      ## distance on the walk
                      D1 = sparseMatrix(ij[, 'i'],ij[, 'j'],
                        x = suppressWarnings(
                          distance(IRanges(start = end(allps.u.ir[ij[,'i']]), width = 1),
                                   IRanges(start(allps.u.ir[ij[,'j']]), width = 1))) + 1e-5, dims = rep(length(allps.u.ir), 2))

                      ## distance on the reference
                      D2 = sparseMatrix(ij[, 'i'],ij[, 'j'],
                        x = suppressWarnings(
                          distance(gr.end(allps.u[ij[,'i']], ignore.strand = F),
                                   gr.start(allps.u[ij[,'j']], ignore.strand = F))) + 1e-5, dims = rep(length(allps.u.ir), 2))

                      D1 = pmin(as.matrix(D1), as.matrix(t(D1)))
                      D2 = pmin(as.matrix(D2), as.matrix(t(D2)))

                      tmp = D1>0 & D1<prune.d1 & D2>0 & D2<prune.d2
                      tmp[which(is.na(tmp))] = FALSE
                      G = graph.adjacency(tmp)
                      cl = clusters(G, 'weak')$membership ## clusters based on this adjacency relationship
                      cls = split(1:length(cl), cl)
                      lens = sapply(allps, length)

                      ## check if there any clusters that contain both the first and last member  of a walk
                      cls.fl = cls[mapply(function(x) all(c(1,lens[allps.u$grl.ix[x[1]]]) %in% allps.u$grl.iix[x]), cls)]

                      if (length(cls.fl)>0)
                        {
                          toprune = allps.u$grl.ix[sapply(cls.fl, function(x) x[1])]
                          if (length(toprune)>0)
                            cat('Pruning', length(toprune), 'walks\n')
                          allps = allps[-toprune]
                        }                          
                    }
                }
              
              if (length(allps)>0)
                win = streduce(unlist(allps), 0)
#                win = streduce(unlist(allps), sum(width(unlist(allps)))*0)

              values(allps) = NULL
              out$allpaths = allps
              out$allpaths.og = allps.og ## untouched all.paths if we want to reheat eg after computing 10X support
              gt.walk = gTrack(out$allpaths, draw.paths = T,border = NA, angle = 0, ywid = 0.5, height = 20, labels.suppress.gr = T)
              gt.walk$path.cex.arrow = 0
              gt.walk$path.stack.x.gap = 1e6
              out$td.allpaths = c(
                  gt.walk,
                  td.seg,
                  td.rg)
              pdf(outfile.allpaths.pdf, height = 30, width = 24)
              plot(out$td.allpaths,
                      windows = win, links = kag$junctions)
              dev.off()                                  
              out$README = paste(out$README, 'allpaths= all paths through windows (not just optimal ones), td.allpaths = gTrack object of plot of all paths')
            }
          
          ## if junction.ix was specified then label which positions in the walks represent the rearrangement junction
          if (!is.null(junction.ix) & length(out$allpaths)>0)
            {
              allps = out$allpaths
              allps.u = grl.unlist(allps)
              allps.u$ix.s = gr.match(gr.start(allps.u, ignore.strand = F), starts, ignore.strand = F)
              allps.u$ix.e = gr.match(gr.end(allps.u, ignore.strand = F), ends, ignore.strand = F)
              allps.ixs = split(allps.u$ix.s, allps.u$grl.ix) ## start indices of walk intervals in sol$segstats
              allps.ixe = split(allps.u$ix.e, allps.u$grl.ix) ## end indices of walks intervals in sol$segstats
              allps.ix.pairs = sapply(mapply(function(x,y) if (length(x)<=1) NULL else which(paste(x[-length(x)], y[-1], sep = ',') %in% junc.pair), allps.ixe, allps.ixs, SIMPLIFY = F), paste, collapse = ',')
              values(allps)$junction.id = names(junction.ix)[i]
              values(allps)$junction.ix = allps.ix.pairs
              out$allpaths = allps
            }

          if (length(out$allpaths)>0)
            {
              values(out$allpaths)$string = grl.string(out$allpaths)
              values(out$allpaths)$wid = sapply(out$allpaths, function(x) sum(width(x)))
              values(out$allpaths)$wids = sapply(out$allpaths, function(x) paste(width(x), collapse = ','))              
              write.tab(as.data.frame(values(out$allpaths)), outfile.allpaths.txt)
            }

          saveRDS(out, outfile.rds)             
          return(out)
        }, mc.cores = mc.cores)
    }
  
  ## awkward workaround to limit the number of processors Cplex will gobble up
  if (customparams)
    {
      system(paste('rm', param.file))
      Sys.setenv(ILOG_CPLEX_PARAMETER_FILE='')    
      cat('Finished\n')
    }
  
  return(out)
}


##############################################
#' reads.to.walks
#'
#' Utility function to realign reads to walks.
#' 
#' Takes bam and collection of walks (GRanges list of signed intervals on hg19 or other BSgenome)
#' pulls reads in regions of walks, then uses bwa mem to realign reads to walks, returns paths to new bams
#' These bams are in "walk coordinates"
#'
#' Assumes bwa mem installed an runnable from command line.
#' 
#' @param bam path to bam file
#' @param walks GRangesList of walks
#' @param outdir outdir to compute into 
#' @param hg human genome sequence BSGenome or ffTrack
#' @param mc.cores number of cores
#' @param insert.size >= max insert size of library so that all relevant read pairs are extracted from the original bam
#' @return paths to new bams
#' These bams are in "walk coordinates"
#' @export
##############################################
reads.to.walks = function(bam, walks, outdir = './test', hg = skidb::read_hg(fft = T), mc.cores = 1, insert.size = 1e3, verbose = T)
  {
    system(paste('mkdir -p', outdir))
    
    .pairs2fq = function(x)
      {
        x.gr = grl.unlist(x)
        x.gr = x.gr[!is.na(x.gr$seq)]
        x.gr$first = bamflag(x.gr$flag)[,'isFirstMateRead']!=0
        x.gr$unmapped = bamflag(x.gr$flag)[,'isUnmappedQuery']!=0
        x.gr1 = x.gr[x.gr$first]
        x.gr2 = x.gr[!x.gr$first]
        x.gr1 = x.gr1[x.gr1$qname %in% x.gr2$qname]
        x.gr2 = x.gr2[match(x.gr1$qname, x.gr2$qname), ]
        if (any(ix <- as.logical(strand(x.gr1)=='-') & !x.gr1$unmapped)) ## rev comp sequence and reverse quality scores
          {
            x.gr1$seq[ix] = as.character(reverseComplement(DNAStringSet(x.gr1$seq[ix])))
            x.gr1$qual[ix] = sapply(x.gr1$qual[ix], function(x) rawToChar(rev(charToRaw(x))))
          }
        
        if (any(ix <- as.logical(strand(x.gr2)=='-') & !x.gr1$unmapped)) ## rev comp sequence and reverse quality scores
          {
            x.gr2$seq[ix] = as.character(reverseComplement(DNAStringSet(x.gr2$seq[ix])))
            x.gr2$qual[ix] = sapply(x.gr2$qual[ix], function(x) rawToChar(rev(charToRaw(x))))
          }                   
        fq1 = as.vector(t(cbind(paste('@', x.gr1$qname, sep = ''), x.gr1$seq, '+', x.gr1$qual)))
        fq2 = as.vector(t(cbind(paste('@', x.gr2$qname, sep = ''), x.gr2$seq, '+', x.gr2$qual)))
        return(list(fq1 = fq1, fq2 = fq2))
      }

    if (is.null(names(walks)))
      names(walks) = paste('walk', 1:length(walks))

    outdir = normalizePath(outdir)
    
    reads.fq = paste(outdir, '/reads.', 1:2, '.fq', sep = '')
    walk.fa = paste(outdir, '/', names(walks), '.fa', sep = '')
    walks.gff = paste(outdir, '/walks.gff', sep = '')
      
    if (!all(file.exists(walk.fa)))
      {
        if (verbose)
          cat('Fetching walk sequences\n')
        
        walk.seq = ffTrack::get_seq(hg, walks, mc.cores = mc.cores)
        names(walk.seq) = names(walks)
                
        sapply(1:length(walks), function(x) writeXStringSet(walk.seq[x], walk.fa[x]))
                
        if (is.list(walk.seq))
          {
            c(walk.seq[[1]], walk.seq[[2]]) ## weird R ghost
            walk.seq = do.call('c', walk.seq)
          }

        ## write compiled fasta
        writeXStringSet(walk.seq, paste(outdir, '/walks.fa', sep = ''))
      }

    if (!all(file.exists(walks.gff)))
      {
        tmp = spChain(walks)$y;
        export.gff(split(tmp, seqnames(tmp)), walks.gff)
      }
       
    ## grab reads from region and output to fq
    if (!all(file.exists(reads.fq)))
      {
        if (verbose)
          cat(sprintf('Fetching mapped and unmapped reads associated with region from %s\n', bam))
        
        reads = read.bam(bam, intervals = streduce(unlist(walks)+insert.size))
        reads.seq = .pairs2fq(reads)
        
        ## write read fastqs
        if (verbose)
          cat(sprintf('Writing fastqs for %s read pairs\n', length(reads[[1]])))
        
        bla = mapply(function(x, y) writeLines(x, y), reads.seq, reads.fq)
      }
    
    if (verbose)
      cat('Indexing walk fastas')

    walks.faidx = paste(walk.fa, '.bwt', sep = '')

    if (any(ix <- !file.exists(walks.faidx)))      
      mclapply(walk.fa[ix], function(x) system(paste('bwa index', x)), mc.cores = mc.cores)
    
    if (verbose)
      cat('Running bwa mem\n')
      
    mclapply(1:length(walk.bam), function(x)
             {
               cmd = sprintf('bwa mem %s %s %s | samtools view -bSh -F4 - > %s; samtools sort %s %s; samtools index %s',
                     walk.fa[x], reads.fq[1], reads.fq[2], walk.bam[x], walk.bam[x], gsub('.bam$', '', walk.bam[x]), walk.bam[x])
               if (verbose)
                 cat(cmd, '\n')
               system(cmd)
               }, mc.cores = mc.cores)

    if (verbose)
      cat('Done\n')
    
    return(walk.bam)    
  }




###############################################
#' annotate.walks
#'
#' Low level function to annotate walks (GRanges list) with cds / promoter annotations
#' 
#' given:
#' walks: input grl of walks on the genome
#' tx:  gr annotating transcript boundaries
#' or grl annotating exon level transcripts e.g. refgene or grl with 
#' grl-level meta data fields $s1, $s2, $e1, $e2 annotating start and end positions of transcipt and cds
#' respectively, $gene_sym representing gene label, $chr - chromosome, $str strand
#' and gr level features (exon_frame) annotating the frame (0,1,2) of the first exon position in a + transcript
#' and last exon position in a - transcript.
#' Assumes that exons are ordered in each grl item in the order of transcription (i.e. right most exon for negative strand transcripts)
#' (e.g. output of read_refGene(grl = T))
#' 
#' @param walks GRangesList of walks to query (eg from traversal of JaBbA object)
#' @param cds  GRangesList of CDS annotation, one per transcript (each a coding region of an exon)
#' @param promoters GRanges of promoters (same length as transcript)
#' @param filter.splice flag whether to filter out splice variants of a given gene
#' @param verbose  flag 
#' @param prom.window window to use around each transcript to identify putative promoter if promoter is NULL
#' @param mc.cores number of cores to use
#' @return
#' a grl of putative fusions with annotations in values field:
#' $label
#' $type  e.g. promoter-fusion, 5-UTR fusion, In-frame fusion, 3' truncated fusion, 5' truncated fusion, in-frame poly-fusion, 3' truncated poly-fusion,
#'             5' truncated poly-fusion
#' $genes genes involved in fusion
#' $transcripts transcripts involved in fusion
#' 
#' annotates every possible altered transcript in region including
#'   - transcripts with truncated 
#' @export
############################################
annotate.walks = function(walks, cds, promoters = NULL, filter.splice = T, verbose = F, prom.window = 1e3, max.chunk = 1e9, mc.cores = 1, exhaustive = FALSE)
    {
        require(igraph)

        if (inherits(walks, 'GRanges'))
            walks = GRangesList(walks)
        
        if (is(walks, 'list'))
            walks = do.call(GRangesList, walks)

        if (!is(cds, 'GRangesList'))
            {
                if (verbose)
                    cat('splitting cds\n')
                cds = .gencode_split(cds, by = 'transcript_id')
            }
        
        tx.span = seg2gr(values(cds)) ## assumed that transcript span is encoded in the cds metadata (i.e. beginning end including UTR)

        cdsu = gr2dt(grl.unlist(cds)[, c('grl.ix')])
        setkey(cdsu, grl.ix)

        cds.span = cdsu[, list(start = min(start), end = max(end)), keyby = grl.ix][list(1:length(tx.span)), ]

        
        utr.left = seg2gr(gr2dt(tx.span)[, list(seqnames = seqnames, start = start, strand = strand, end = cds.span[list(1:length(start)), start], transcript_id = Transcript_id, transcript_name = Transcript_name, gene_name = Gene_name)])
        utr.right = seg2gr(gr2dt(tx.span)[, list(seqnames = seqnames, start = cds.span[list(1:length(start)), end], strand = strand, end = end, transcript_id = Transcript_id, transcript_name = Transcript_name, gene_name = Gene_name)])
        utr = c(utr.left, utr.right)
        
        names(values(tx.span))[match(c('Transcript_id', 'Transcript_name', 'Gene_name'), names(values(tx.span)))] = c('transcript_id', 'transcript_name', 'gene_name')
        
        if (is.null(promoters))
            {
                promoters = flank(tx.span, prom.window)                
                values(tx.span) = values(promoters)
            }
        
        tx.span$type = 'cds'
        promoters$type = 'gene'
        utr.left$type = 'gene'
        utr.right$type = 'gene'
        tx.span$cds.id = 1:length(tx.span)
        tx.span = seg2gr(as.data.table(rrbind(as.data.frame(tx.span), as.data.frame(promoters)))[, list(seqnames = seqnames[1], start = min(start), end = max(end), strand = strand[1], gene_name = gene_name[1], transcript_id = transcript_id[1], transcript_name = transcript_name[1], cds.id = cds.id), keyby = cds.id][!is.na(cds.id), ], seqlengths = seqlengths(tx.span))
                        
        # match up tx.span to walks
        walks.u = grl.unlist(walks)

        ## these are fragments of transcripts that overlap walks
        this.tx.span = gr.findoverlaps(tx.span, walks.u, qcol = c('transcript_id', 'transcript_name', 'gene_name'), verbose = verbose, max.chunk = max.chunk)
        this.tx.span$tx.id = this.tx.span$query.id

        strand(this.tx.span) = strand(tx.span)[this.tx.span$query.id]
        this.tx.span$left.broken = start(this.tx.span) != start(tx.span)[this.tx.span$query.id]
        this.tx.span$right.broken = end(this.tx.span) != end(tx.span)[this.tx.span$query.id]
        
        # remove elements that are "unbroken" by the window boundaries
        this.tx.span = this.tx.span[this.tx.span$left.broken | this.tx.span$right.broken]
        this.tx.span$cds.sign = c('-'= -1, '+' = 1)[as.character(strand(this.tx.span))]
        this.tx.span$window.sign = c('-'= -1, '+' = 1)[as.character(strand(walks.u)[this.tx.span$subject.id])]

        # annotate left and right ends (if information available)
        # i.e. UTR, CDS, Promoter

        ## we trim cds by 1 nucleotide at both ends so that 'cds' annotation only refers to a cds fragment (not including start and stop)
                ## annotate ends with feature types so that positions internal to cds bases will be 'cds',
        ## external to cds but in cds will 'utr'
        ## and 5' flank (with respect to cds orientation) will be promoter
        this.tx.span$left.feat = NA
        this.tx.span$left.feat[gr.findoverlaps(gr.start(this.tx.span), tx.span, by = 'transcript_id', max.chunk = max.chunk)$query.id] = 'cds'
        this.tx.span$left.feat[gr.findoverlaps(gr.start(this.tx.span), utr, by = 'transcript_id', max.chunk = max.chunk)$query.id] = 'utr'
        this.tx.span$left.feat[gr.findoverlaps(gr.start(this.tx.span), promoters, by = 'transcript_id', max.chunk = max.chunk)$query.id] = 'promoter'
                   
        this.tx.span$right.feat = NA
        this.tx.span$right.feat[gr.findoverlaps(gr.end(this.tx.span), tx.span, by = 'transcript_id', max.chunk = max.chunk)$query.id] = 'cds'
        this.tx.span$right.feat[gr.findoverlaps(gr.end(this.tx.span), utr, by = 'transcript_id', max.chunk = max.chunk)$query.id] = 'utr'
        this.tx.span$right.feat[gr.findoverlaps(gr.end(this.tx.span), promoters, by = 'transcript_id', max.chunk = max.chunk)$query.id] = 'promoter'
                
        ## if lands in CDS annotate this.tx.span ends with right and/or left exon frame
        ## (if lands in intron then annotate left end with frame of next exon on right and right end
        ## with frame of next exon on left)
        ## otherwise annotate as NA
        ## we will eventually integrate frames across walks and call a transition "in frame" if the
        ## frame of the right (left) end of the previous + (-) interval
        
        ## now we want to find the first exon to the right of the left boundary and
        ## the first exon to the left of the right boundary for each fragment
        ##tix = match(this.tx.span$transcript_id, tx.span$transcript_id)

        cds.u = grl.unlist(cds[this.tx.span$tx.id])
                                        #        ranges(cds.u) =  ranges(pintersect(cds.u, tx.span[this.tx.span$tx.id[cds.u$grl.ix]], resolve.empty = 'start.x'))
        ranges(cds.u) =  ranges(pintersect(cds.u, tx.span[this.tx.span$tx.id[cds.u$grl.ix]]))
        
        tmp = gr.findoverlaps(this.tx.span, cds.u, scol = c('start.local', 'end.local', 'exon_number'), by = 'transcript_id', verbose = verbose, max.chunk = max.chunk)
        
        leftmost.cds.exon = data.table(id = 1:length(tmp), qid = tmp$query.id, start = start(tmp))[, id[which.min(start)], by = qid][, V1]
        rightmost.cds.exon = data.table(id = 1:length(tmp), qid = tmp$query.id, end = end(tmp))[, id[which.max(end)], by = qid][, V1]

        ## now we want to get the frame of the left and right base
        ## of each leftmost and rightmost exon (etc.
        ## this will depend on orientation of the exon and side that we are querying
        ## for left side of - exon, (exonFrame + width) %% 3
        ## for right side of - exon, (exonFrame + width(og.exon) - width %% 3
        ## for left side of + exon, (exonFrame + width(og.exon) - width) %%3
        ## for right side of - exon, (exonFrame + int.exon) %% 3

        leftmost.coord = ifelse(as.logical(strand(cds.u[tmp$subject.id[leftmost.cds.exon]])=='+'),
            (start(tmp)[leftmost.cds.exon] - start(cds.u)[tmp$subject.id[leftmost.cds.exon]] + cds.u$start.local[tmp$subject.id[leftmost.cds.exon]]),
            (end(cds.u)[tmp$subject.id[leftmost.cds.exon]] - start(tmp)[leftmost.cds.exon] + cds.u$start.local[tmp$subject.id[leftmost.cds.exon]]))
        
        rightmost.coord = ifelse(as.logical(strand(cds.u[tmp$subject.id[rightmost.cds.exon]])=='+'),
            (end(tmp)[rightmost.cds.exon] - start(cds.u)[tmp$subject.id[rightmost.cds.exon]] + cds.u$start.local[tmp$subject.id[rightmost.cds.exon]]),
            (end(cds.u)[tmp$subject.id[rightmost.cds.exon]] - end(tmp)[rightmost.cds.exon] + cds.u$start.local[tmp$subject.id[rightmost.cds.exon]]))

        leftmost.frame = ifelse(as.logical(strand(cds.u[tmp$subject.id[leftmost.cds.exon]])=='+'),
            (start(tmp)[leftmost.cds.exon] - start(cds.u)[tmp$subject.id[leftmost.cds.exon]] + cds.u$phase[tmp$subject.id[leftmost.cds.exon]]) %% 3,
            (end(cds.u)[tmp$subject.id[leftmost.cds.exon]] - start(tmp)[leftmost.cds.exon] + cds.u$phase[tmp$subject.id[leftmost.cds.exon]]) %% 3)
        
        rightmost.frame = ifelse(as.logical(strand(cds.u[tmp$subject.id[rightmost.cds.exon]])=='+'),
            (end(tmp)[rightmost.cds.exon] - start(cds.u)[tmp$subject.id[rightmost.cds.exon]] + cds.u$phase[tmp$subject.id[rightmost.cds.exon]]) %% 3,
          (end(cds.u)[tmp$subject.id[rightmost.cds.exon]] - end(tmp)[rightmost.cds.exon] + cds.u$phase[tmp$subject.id[rightmost.cds.exon]]) %% 3)
                
        leftmost.frame = leftmost.coord %% 3
        rightmost.frame = rightmost.coord %% 3

        this.tx.span$left.coord = this.tx.span$left.boundary = this.tx.span$right.coord = this.tx.span$right.boundary = NA
        
        this.tx.span$left.coord[tmp$query.id[leftmost.cds.exon]] = leftmost.coord
        this.tx.span$right.coord[tmp$query.id[rightmost.cds.exon]] = rightmost.coord
        
        this.tx.span$left.boundary[tmp$query.id[leftmost.cds.exon]] =
            ifelse(strand(this.tx.span)[tmp$query.id[leftmost.cds.exon]]=="+",
                   tmp$start.local[leftmost.cds.exon], tmp$end.local[leftmost.cds.exon])
        
        this.tx.span$right.boundary[tmp$query.id[rightmost.cds.exon]] =
            ifelse(strand(this.tx.span)[tmp$query.id[rightmost.cds.exon]]=="+",
                   tmp$end.local[rightmost.cds.exon], tmp$start.local[rightmost.cds.exon])
                
        this.tx.span$right.exon_del = this.tx.span$left.exon_del = NA;
        this.tx.span$right.frame = this.tx.span$left.frame = NA;
        this.tx.span$right.exon_id= this.tx.span$left.exon_id = NA;
        
        ## keep track of exon frames to left and right
        this.tx.span$left.frame[tmp$query.id[leftmost.cds.exon]] = leftmost.frame
        this.tx.span$right.frame[tmp$query.id[rightmost.cds.exon]] = rightmost.frame
               
        this.tx.span$left.exon_id[tmp$query.id[leftmost.cds.exon]] = cds.u[tmp$subject.id[leftmost.cds.exon]]$exon_number
        this.tx.span$right.exon_id[tmp$query.id[rightmost.cds.exon]] = cds.u[tmp$subject.id[rightmost.cds.exon]]$exon_number

        ## keep track of which exons have any sort of sequence deletion
        this.tx.span$left.exon_del[tmp$query.id[leftmost.cds.exon]] =
            start(this.tx.span)[tmp$query.id[leftmost.cds.exon]] > start(cds.u)[tmp$subject.id[leftmost.cds.exon]]
        this.tx.span$right.exon_del[tmp$query.id[rightmost.cds.exon]] =
            end(this.tx.span)[tmp$query.id[rightmost.cds.exon]] <  end(cds.u)[tmp$subject.id[rightmost.cds.exon]]
        
        
        ## now traverse each walk and connect "broken transcripts"
        ## each walk is a list of windows, connections will be made with respect to the window walk
        ## paying attention to (1) side of breakage, (2) orientation of adjacent windows in walk (3) orientation of transcript
        ##                               
        ## right broken + cds upstream of ++ junction attaches to left  broken + cds in next window        
        ##              - cds upstream of ++ junction attached to left  broken - cds in next window
        ##              + cds upstream of +- junction attaches to right broken - cds in next window
        ##              - cds upstream of +- junction attaches to right broken + cds in next window
        ##
        ## left  broken + cds upstream of -- junction attaches to right broken + cds in next window        
        ##              - cds upstream of -- junction attached to right broken - cds in next window
        ##              + cds upstream of -+ junction attaches to left  broken - cds in next window
        ##              - cds upstream of -+ junction attaches to left  broken + cds in next window
        ##
        ##
        ## 
        
        # to achieve this create a graphs on elements of this.tx.span
        # connecting elements i and j if a window pair k k+1 in (the corresponding) input grl
        # produces a connection with the correct orientation

        # match up on the basis of the above factors to determine edges in graph

        ##

        edges = merge(
            data.table(
                i = 1:length(this.tx.span),
                key1 = walks.u$grl.ix[this.tx.span$subject.id],
                key2 = walks.u$grl.iix[this.tx.span$subject.id],
                key3 = this.tx.span$cds.sign*this.tx.span$window.sign,
                key4 = ifelse(this.tx.span$window.sign>0, this.tx.span$right.broken, this.tx.span$left.broken)),
            data.table(
                j = 1:length(this.tx.span),
                key1 = walks.u$grl.ix[this.tx.span$subject.id],
                key2 = walks.u$grl.iix[this.tx.span$subject.id]-1,
                key3 = this.tx.span$cds.sign*this.tx.span$window.sign,
                key4 = ifelse(this.tx.span$window.sign>0, this.tx.span$left.broken, this.tx.span$right.broken)), by = c('key1', 'key2', 'key3', 'key4'), allow.cartesian = TRUE
            )[key4 == TRUE, ]
        
        ## remove edges that link different transcripts of same gene

        if (filter.splice)
            edges = edges[!(this.tx.span$gene_name[edges$i] == this.tx.span$gene_name[edges$j] & this.tx.span$transcript_id[edges$i] != this.tx.span$transcript_id[edges$j]), ]

        if (nrow(edges)==0)
            return(GRangesList())

        require(Matrix)
        A = sparseMatrix(edges$i, edges$j, x = 1, dims = rep(length(this.tx.span),2))
        sources = which(Matrix::colSums(A!=0)==0)
        sinks = which(Matrix::rowSums(A!=0)==0)

        G = graph.adjacency(A)
        C = clusters(G, 'weak')
        vL = split(1:nrow(A), C$membership)

        ## collate all paths through this graph
        paths = do.call('c', mclapply(1:length(vL), function(i) {
            if (verbose & (i %% 10)==0)
                cat(i, ' of ', length(vL), '\n')
            x = vL[[i]]
            tmp.source = setdiff(match(sources, x), NA)
            tmp.sink = setdiff(match(sinks, x), NA)
            tmp.mat = A[x, x, drop = FALSE]!=0
            if (length(x)<=1)
                return(NULL)
            if (length(x)==2)
                list(x[c(tmp.source, tmp.sink)])            
            else if (all(Matrix::rowSums(tmp.mat)<=1) & all(Matrix::colSums(tmp.mat)<=1))
                get.shortest.paths(G, from = intersect(x, sources), intersect(x, sinks))$vpath
            else
                {
                    if (exhaustive)                    
                        lapply(all.paths(A[x,x, drop = FALSE], source.vertices = tmp.source, sink.vertices = tmp.sink, verbose = FALSE)$paths, function(y) x[y])
                    else
                        {
                            out = do.call('c', lapply(intersect(x, sources),
                                function(x, sinks) suppressWarnings(get.shortest.paths(G, from = x, to = sinks)$vpath), sinks = intersect(x, sinks)))
                            out = out[sapply(out, length)!=0]
                            if (length(out)>0)
                                out = out[!duplicated(sapply(out, paste, collapse = ','))]
                            return(out)                       
                        }
                }
        }, mc.cores = mc.cores))
                    
        fus.sign = this.tx.span$cds.sign * this.tx.span$window.sign
        paths = lapply(paths, function(x) if (fus.sign[x][1]<0) rev(x) else x) ## reverse "backward paths" (i.e. those producing fusions in backward order)
        paths.first = sapply(paths, function(x) x[1])
        paths.last = sapply(paths, function(x) x[length(x)])
        paths.broken.start = ifelse(as.logical(strand(this.tx.span)[paths.first] == '+'), this.tx.span$left.broken[paths.first], this.tx.span$right.broken[paths.first])
        paths.broken.end = ifelse(as.logical(strand(this.tx.span)[paths.last] == '+'), this.tx.span$right.broken[paths.last], this.tx.span$left.broken[paths.last])
        paths.u = unlist(paths)
        paths.i = unlist(lapply(1:length(paths), function(x) rep(x, length(paths[[x]]))))
        tmp.gr = this.tx.span[paths.u]
        
        ## annotate steps of walk with out of frame vs in frame (if cds is a grl)
                
        ## left and right exon frame
        paths.u.lec = tmp.gr$left.coord
        paths.u.rec = tmp.gr$right.coord
        paths.u.lef = tmp.gr$left.frame
        paths.u.ref = tmp.gr$right.frame
        paths.u.str = as.character(strand(tmp.gr))
        paths.u.lcds = tmp.gr$left.feat == 'cds'
        paths.u.rcds = tmp.gr$right.feat == 'cds'
        paths.u.lout = tmp.gr$left.feat %in% c('promoter', 'utr')
        paths.u.rout = tmp.gr$right.feat %in% c('promoter', 'utr')
        
                                        # a fragment is in frame either if (1) it begins a walk at frame 0 outside of the cds
                                        # or if (2) its frame is concordant with previous cds
        paths.u.inframe = rep(NA, length(paths.u))
        paths.u.cdsend = paths.u.cdsstart = rep(FALSE, length(paths.u)) # this keeps track of cds starts and cds ends in the fusion
        
        outside = TRUE        
        for (i in 1:length(paths.u.inframe))
            {
                if (i == 1)
                    outside = TRUE
                else if (paths.i[i] != paths.i[i-1] | (paths.u.str[i] == '+' & paths.u.lout[i]) | (paths.u.str[i] == '-' & paths.u.rout[i]))
                    outside = TRUE
                
                if (outside)
                    {                  
                        if (paths.u.str[i] == '+')
                            paths.u.inframe[i] = paths.u.lout[i] & paths.u.lef[i] == 0
                        else
                            paths.u.inframe[i] = paths.u.rout[i] & paths.u.ref[i] == 0
                        
                        paths.u.cdsstart[i] = paths.u.inframe[i]                  
                        outside = F
                    }
                else
                    {
                        if (paths.u.str[i] == '+' & paths.u.str[i-1] == '+')
                            paths.u.inframe[i] = paths.u.lec[i] != 1 & paths.u.lef[i] == ((paths.u.ref[i-1]+1) %% 3) & paths.u.lcds[i] & paths.u.rcds[i-1]
                        else if (paths.u.str[i] == '+' & paths.u.str[i-1] == '-')
                            paths.u.inframe[i] = paths.u.lec[i] != 1 & paths.u.ref[i]  == ((paths.u.ref[i-1]+1) %% 3) & paths.u.rcds[i] & paths.u.rcds[i-1]
                        else if (paths.u.str[i] == '-' & paths.u.str[i-1] == '-') 
                            paths.u.inframe[i] = paths.u.rec[i] != 1 & paths.u.ref[i] == ((paths.u.lef[i-1]+1) %% 3) & paths.u.rcds[i] & paths.u.lcds[i-1]
                        else if (paths.u.str[i] == '-' & paths.u.str[i-1] == '+')
                            paths.u.inframe[i] = paths.u.rec[i] != 1 & paths.u.lef[i] == ((paths.u.lef[i-1]+1) %% 3) & paths.u.lcds[i] & paths.u.lcds[i-1]
                    }

                if ((paths.u.str[i] == '+' & paths.u.rout[i]) | (paths.u.str[i] == '-' & paths.u.lout[i]))
                    {
                        paths.u.cdsend[i] = paths.u.inframe[i]
                        outside = T
                    }
            }

        tmp.gr$in.frame = paths.u.inframe;
        tmp.gr$cds.start = paths.u.cdsstart
        tmp.gr$cds.end = paths.u.cdsend
        tmp.gr$del5 = ifelse(paths.u.str == '+', tmp.gr$left.exon_del, tmp.gr$right.exon_del)
        tmp.gr$del3 = ifelse(paths.u.str == '+', tmp.gr$right.exon_del, tmp.gr$left.exon_del)
        tmp.gr$first.coord = ifelse(paths.u.str == '+', tmp.gr$left.coord, tmp.gr$right.coord)
        tmp.gr$last.coord  = ifelse(paths.u.str == '+', tmp.gr$right.coord, tmp.gr$left.coord)
        tmp.gr$first.boundary = ifelse(paths.u.str == '+', tmp.gr$left.boundary, tmp.gr$right.boundary)
        tmp.gr$last.boundary  = ifelse(paths.u.str == '+', tmp.gr$right.boundary, tmp.gr$left.boundary)                
        tmp.gr$first.exon = ifelse(paths.u.str == '+', tmp.gr$left.exon_id, tmp.gr$right.exon_id)        
        tmp.gr$last.exon = ifelse(paths.u.str == '+', tmp.gr$right.exon_id, tmp.gr$left.exon_id)                
        fusions = split(tmp.gr, paths.i)

                                        # now annotate fusions

        values(fusions)$walk.id = data.table(wid = walks.u$grl.ix[tmp.gr$subject.id], fid = paths.i)[, wid[1], keyby = fid][, V1]
 #    values(fusions)$walk.id = vaggregate(walks.u$grl.ix[tmp.gr$subject.id], by = list(paths.i), FUN = function(x) x[1])
                
        tmp.g = tmp.gr$gene_name
        tmp.cds = tmp.gr$transcript_id
        tmp.fe = as.numeric(tmp.gr$first.exon)
        tmp.le = as.numeric(tmp.gr$last.exon)
        tmp.fb = tmp.gr$first.boundary
        tmp.lb = tmp.gr$last.boundary
        tmp.fc = tmp.gr$first.coord
        tmp.lc = tmp.gr$last.coord
        tmp.5d = tmp.gr$del5
        tmp.3d = tmp.gr$del3
        tmp.5d[is.na(tmp.5d)] = FALSE
        tmp.3d[is.na(tmp.3d)] = FALSE
        paths.u.cdsend[is.na(paths.u.cdsend)] = FALSE
        paths.u.cdsstart[is.na(paths.u.cdsstart)] = FALSE
        paths.u.inframe[is.na(paths.u.inframe)] = FALSE

        totpaths = max(paths.i)

        if (verbose)
            cat('Populating coordinates\n')
        
        values(fusions)[, 'coords'] = mcmapply(function(x) paste(unique(x), collapse = '; '),
                           split(gr.string(this.tx.span[paths.u], mb = TRUE, round = 1), paths.i), mc.cores = mc.cores)

        if (verbose)
            cat('Populating transcript names\n')
        values(fusions)[, 'transcript_names'] = mcmapply(function(x, y) paste(x, ' (', y, ')', sep = '', collapse = '; '),
                           split(values(tx.span)[, 'gene_name'][this.tx.span$query.id[paths.u]], paths.i), 
                           split(values(tx.span)[, 'transcript_name'][this.tx.span$query.id[paths.u]], paths.i), mc.cores = mc.cores)

        if (verbose)
            cat('Populating transcript ids\n')
        values(fusions)[, 'transcript_ids'] = mcmapply(function(x, y) paste(x, ' (', y, ')', sep = '', collapse = '; '),
                           split(values(tx.span)[, 'gene_name'][this.tx.span$query.id[paths.u]], paths.i), 
                           split(values(tx.span)[, 'transcript_id'][this.tx.span$query.id[paths.u]], paths.i), mc.cores = mc.cores)

        if (verbose)
            cat('Populating gene names\n')
        values(fusions)[, 'genes'] = mcmapply(function(x) paste(unique(x), collapse = '; '),
                           split(values(tx.span)[, 'gene_name'][this.tx.span$query.id[paths.u]], paths.i), mc.cores = mc.cores)

        if (verbose)
            cat('Populating alteration\n')
        values(fusions)$alteration =  vaggregate(1:length(paths.i), by = list(paths.i),
                           FUN = function(x)
                               {
                                   if (verbose & (x[1] %% 10)==0)
                                       cat('Path', unique(paths.i[x]), 'of', totpaths, '\n')
                                   if (length(unique((tmp.cds[x])))==1) ## single transcript event
                                       {
                                           out = NULL
                                           x = x[!is.na(tmp.fe[x]) & !is.na(tmp.le[x])]                                           
                                           if (length(x)>0)
                                               {
#                                                   browser()
                                                   ir = IRanges(pmin(tmp.le[x], tmp.fe[x]), pmax(tmp.fe[x], tmp.le[x]))
                                                   if (length(del <- setdiff(IRanges(min(tmp.fe[x]), max(tmp.le[x])), ir))>0)
                                                       {
                                                           del.fc = pmax(tmp.lc[x[match(start(del)-1, tmp.le[x])]]+1, 1, na.rm = TRUE)
                                                           del.lc = pmin(tmp.fc[x[match(end(del)+1, tmp.fe[x])]]-1, max(tmp.lc[x]), na.rm = TRUE)
                                                           out = c(out,## some portion deleted 
                                                               ifelse(start(del)==end(del),
                                                                      paste('deletion of exon ', start(del),
                                                                            ' [', del.fc, '-', del.lc, 'bp]', 
                                                                            sep = '', collapse = ', '),
                                                                      paste('deletion of exons ', start(del), '-', end(del),
                                                                            ' [', del.fc, '-', del.lc, 'bp]', 
                                                                            sep = '', collapse = ', ')))
                                                       }
                                                   
                                                   if (length(amp <- IRanges(coverage(ir)>1))>0)
                                                       {
                                                           amp.fc = tmp.lc[x[match(start(amp), tmp.le[x])]]
                                                           amp.lc = tmp.fc[x[match(end(amp), tmp.fe[x])]]
                                                           out = c(out,   ## some portion duplicated
                                                               ifelse(start(amp)==end(amp),
                                                                      paste('duplication of exon ', end(amp), 
                                                                            '[', amp.fc, '-', amp.lc, 'bp]', 
                                                                            sep = '', collapse = ', '),
                                                           paste('duplication of exons ', start(amp), '-', end(amp),
                                                                 ' [', amp.fc, '-', amp.lc, 'bp]',
                                                                 sep = '', collapse = ', ')))
                                                       }
                                                   
                                                   if (any(ix <- tmp.5d[x]))
                                                       {
                                                           out = c(out, paste("partial 5' deletion of exon ", tmp.fe[x[ix]],
                                                               ' [', tmp.fb[x[ix]], '-', tmp.fc[x[ix]], 'bp]', 
                                                               sep = '', collapse = ', '))  ## some portion duplicated
                                                       }
                                                   if (any(ix <- tmp.3d[x]))
                                                       {
                                                           del.fc = pmax(tmp.lc[x[ix-1]] + 1, 1, na.rm = TRUE)
                                                           del.lc = pmin(tmp.fc[x[ix+1]]-1, max(tmp.lc[x]), na.rm = TRUE)
                                                           out = c(out, paste("partial 3' deletion of exon ", tmp.fe[x[ix]],
                                                                ' [', tmp.lc[x[ix]], '-', tmp.lb[x[ix]], 'bp]',
                                                               sep = '', collapse = ', '))  ## some portion duplicated
                                                       }
                                               }
                                           
                                           if (length(out)>0)
                                               paste(out, collapse = '; ')
                                           else
                                               ''                   
                                       }
                                   else
                                       {
                                           return(paste(tmp.g[x], ' ', ifelse(paths.u.cdsstart[x], 'S', ''), ifelse(tmp.5d[x],  'tr', ''),
                                                        ifelse(is.na(tmp.fe[x]), 'UTR', 
                                                               ifelse(tmp.le[x]==tmp.fe[x],
                                                                      paste('exon ', tmp.le[x], sep = ''),
                                                                      paste('exons ', tmp.fe[x], '-', tmp.le[x], sep = ''))),
                                                               ' [', tmp.fc[x], '-', tmp.lc[x], 'bp]', 
                                                               ifelse(tmp.3d[x],  'tr', ''), ifelse(paths.u.cdsend[x], 'E', ''), ' ',
                                                               ifelse(c(paths.u.inframe[x[-1]], FALSE), '-',
                                                                      ifelse((1:length(x))!=length(x), '-X', '')), sep = '', collapse = '-> '))
                                       }
                               })

        values(fusions)$max.inframe = vaggregate(paths.u.inframe, by = list(paths.i),
                           FUN = function(x) return(max(c(0, rle(x)$lengths[which(rle(x)$values == T)]))))
        values(fusions)$num.win = vaggregate(paths.u.inframe, by = list(paths.i), length)
#        values(fusions)$broken.start = paths.broken.start[as.numeric(names(fusions))]
#        values(fusions)$broken.end = paths.broken.end[as.numeric(names(fusions))]

        values(fusions) = cbind(values(walks)[values(fusions)$walk.id, ], values(fusions))

        fusions = fusions[nchar(values(fusions)$alteration)>0, ]
        return(fusions)        
    }

#####################################
#' fusions
#' 
#' annotates all gene fusions given an n x n adjacency matrix A of n genomic segments seg and grl of transcripts (eg output of read_RefGene)
#' seg must be (1) a tiling of the genome and (2) have copies of both + and - intervals for each genomic range (eg output of karyograph)
#' 
#' alternate input is a pile of junctions of ranges with strand orientation pointing AWAY from breakend
#' 
#' cds = gencode cds GRanges gff3 / gtf input
#' 
#' "gene_name" GRangesList meta data field is used in annotation and in not creating "splice" fusions that arise from different transcripts of the same gene.
#' 
#' @param junctions GRangesList of junctions (each a length 2 GRanges)
#' @param jab  JaBbA object (overrides junctions input)
#' @param cds  CDS annotations (GrangesList of transcript composed of coordinates coding regions of exons)
#' @param promoters GRanges of promoters (same length as transcript)
#' @param query optional query limiting walks to specific regions of interest
#' @param prom.window window to use around each transcript to identify putative promoter if promoter is NULL
#' @return GRangesList of walks corresponding to transcript boundaires
#' @export
################################
fusions = function(junctions = NULL, jab = NULL, cds = NULL, promoters = NULL, query = NULL, prom.window = 1e3, verbose = T, max.chunk = 1e10, cb.interval = 1e4, cb.chunksize = 1e4, cb.maxchunks = 1e10, exhaustive = FALSE, debug = NULL, mc.cores = 1)
    {
        if (!is.null(junctions))
            {
                jab = karyograph(junctions)                
                A = jab$adj
                seg = jab$tile
                if (verbose)
                    cat('made karyograph\n')
            }
        else if (!is.null(jab))
            {
                if (!is.null(jab$segstats))
                    seg = jab$segstats
                else
                    seg = jab$tile
                A = jab$adj
            }
        else
            stop('either jab or junctions input must be non NULL')                    
        
        if (is.null(A) | is.null(seg))
            stop('Some essential args are NULL')

        if (!is(cds, 'GRangesList'))
            cds = NULL
        else if (!all(c('Transcript_id', 'Gene_name') %in% names(values(cds))))
            cds = NULL                               
        
        if (is.null(cds))
            {
                if (verbose)
                    message('CDS object missing or malformed (e.g. does not contain Transcript_id and Gene_name GRangesList metadata fields\nReading in from gencode CDS via skidb::read_gencode("cds"))')                
                cds = read_gencode('cds')
            }
        
        tx.span = seg2gr(values(cds))

        names(values(tx.span))[match(c('transcript_id', 'gene_name'), tolower(names(values(tx.span))))] = c('transcript_id', 'gene_name')
        
        if (verbose)
            cat('got transcript boundaries\n')
        
        ## determine set of transcript fragments
        ## these correspond to transcripts that intersect a segment boundary
        
        cds.frag.left = gr.findoverlaps(tx.span, gr.start(seg), qcol = c('gene_name', 'transcript_id'), ignore.strand = F, max.chunk = max.chunk)
        strand(cds.frag.left) = strand(tx.span)[cds.frag.left$query.id]
        cds.frag.right = gr.findoverlaps(tx.span, gr.end(seg), qcol = c('gene_name', 'transcript_id'),  ignore.strand = F, max.chunk = max.chunk)
        strand(cds.frag.right) = strand(tx.span)[cds.frag.right$query.id]

        ## I want to find all unique walks that involve tx fragments
        
        if (length(cds.frag.left)>0 & length(cds.frag.right) > 0 )
            tmp = merge(data.frame(i = 1:length(cds.frag.left), key1 = cds.frag.left$query.id, key2 = cds.frag.left$subject.id),
                data.frame(j = 1:length(cds.frag.right), key1 = cds.frag.right$query.id, key2 = cds.frag.right$subject.id), all = T)
        else
            return(GRangesList())

        pos.right = which( as.logical( strand(cds.frag.right)=='+'))
        pos.left = which( as.logical( strand(cds.frag.left)=='+') )
        neg.right = which(as.logical( strand(cds.frag.right)=='-') )
        neg.left = which(as.logical( strand(cds.frag.left)=='-'))
        
        ## positive start fragments will be "right" fragments
        cds.start.frag.pos = cds.frag.right[tmp[is.na(tmp$i) & tmp$j %in% pos.right, ]$j]
        start(cds.start.frag.pos) = start(tx.span)[cds.start.frag.pos$query.id]
        if (length(cds.start.frag.pos)>0)
            cds.start.frag.pos$type = 'start'

        ## positive end fragments will be "left" fragments
        cds.end.frag.pos = cds.frag.left[tmp[is.na(tmp$j) & tmp$i %in% pos.left, ]$i]
        end(cds.end.frag.pos) = end(tx.span)[cds.end.frag.pos$query.id]
        if (length(cds.end.frag.pos)>0)
            cds.end.frag.pos$type = 'end'

        ## negative start fragments will be "right" fragments
        cds.start.frag.neg = cds.frag.left[tmp[is.na(tmp$j) & tmp$i %in% neg.left, ]$i]
        end(cds.start.frag.neg) = end(tx.span)[cds.start.frag.neg$query.id]
        if (length(cds.start.frag.neg)>0)
            cds.start.frag.neg$type = 'start'

        ## negative end fragments will be "left" fragments
        cds.end.frag.neg = cds.frag.right[tmp[is.na(tmp$i) & tmp$j %in% neg.right, ]$j]
        start(cds.end.frag.neg) = start(tx.span)[cds.end.frag.neg$query.id]
        if (length(cds.end.frag.neg)>0)
            cds.end.frag.neg$type = 'end'
        
        ## remaining will be "middle" fragments
        middle.frag = cds.frag.left[tmp[!is.na(tmp$i) & !is.na(tmp$j),]$i]
        end(middle.frag) = end(cds.frag.right[tmp[!is.na(tmp$i) & !is.na(tmp$j),]$j])
        if (length(middle.frag)>0)
            middle.frag$type = 'middle'

        ## concatenate fragments
        ## subject.id of frags is the id of the node on the graph

        all.frags = c(cds.start.frag.pos, cds.end.frag.pos, cds.start.frag.neg, cds.end.frag.neg, middle.frag)

        ##
        ## now connect all.frags according to A
        ## i.e. apply A connections to our fragments, so draw an edge between fragments
        ## if
        ## (1) there exists an edge connecting segment and
        ## (2) only allowable connections are 'start' --> 'middle' --> 'middle' --> 'end'
        ##
        
        seg.edges = as.data.frame(which(A!=0, arr.ind = T))
        colnames(seg.edges) = c('from.seg', 'to.seg')
        edges = merge(merge(data.frame(i = 1:length(all.frags), from.seg = all.frags$subject.id),
            seg.edges), data.frame(j = 1:length(all.frags), to.seg = all.frags$subject.id))

        edges = edges[all.frags$type[edges$i] == 'start' & all.frags$type[edges$j] == 'middle' |
                          all.frags$type[edges$i] == 'start' & all.frags$type[edges$j] == 'end' |
                              all.frags$type[edges$i] == 'middle' & all.frags$type[edges$j] == 'middle' |
                                  all.frags$type[edges$i] == 'middle' & all.frags$type[edges$j] == 'end', ]

        ## this removes splice variants .. keeping only links that fuse different genes or same transcripts
        edges = edges[tx.span$gene_name[all.frags$query.id[edges$i]] != tx.span$gene_name[all.frags$query.id[edges$j]] |
                          all.frags$query.id[edges$i] == all.frags$query.id[edges$j],]

        if (nrow(edges)==0)
            return(GRangesList())

        if (verbose)
            cat('computed subgraph\n')
        
        A.frag = sparseMatrix(edges$i, edges$j, x = 1, dims = rep(length(all.frags),2))
        keep.nodes = which(Matrix::rowSums(A.frag)>0 | Matrix::colSums(A.frag)>0)
        A.frag = A.frag[keep.nodes, keep.nodes]
        all.frags = all.frags[keep.nodes]
        
        sources = which(all.frags$type == 'start')
        sinks = which(all.frags$type == 'end')

        G = graph.adjacency(A.frag)
        C = clusters(G, 'weak')      
        vL = split(1:nrow(A.frag), C$membership)
       
        paths = do.call('c', mclapply(1:length(vL), function(i) {
            if (verbose & (i %% 10)==0)
                cat(i, ' of ', length(vL), '\n')
            x = vL[[i]]
            if (!is.null(debug))
                if (i %in% debug) 
                    browser()
            tmp.source = setdiff(match(sources, x), NA)
            tmp.sink = setdiff(match(sinks, x), NA)
            tmp.mat = A.frag[x, x, drop = FALSE]!=0
            if (length(x)<=1)
                return(NULL)
            if (length(x)==2)
                list(x[c(tmp.source, tmp.sink)])            
            else if (all(Matrix::rowSums(tmp.mat)<=1) & all(Matrix::colSums(tmp.mat)<=1))
                get.shortest.paths(G, from = intersect(x, sources), intersect(x, sinks))$vpath
            else
                {
                    if (exhaustive)                    
                        lapply(all.paths(A.frag[x,x, drop = FALSE], source.vertices = tmp.source, sink.vertices = tmp.sink, verbose = verbose)$paths, function(y) x[y])
                    else
                        {
                            out = do.call('c', lapply(intersect(x, sources),
                                function(x, sinks) suppressWarnings(get.shortest.paths(G, from = x, to = sinks)$vpath), sinks = intersect(x, sinks)))
                            out = out[sapply(out, length)!=0]
                            if (length(out)>0)
                                out = out[!duplicated(sapply(out, paste, collapse = ','))]
                            return(out)                       
                        }
                }
        }, mc.cores = mc.cores))        

        if (verbose)
            cat('computed paths\n')

        paths.u = unlist(paths)
        paths.i = unlist(lapply(1:length(paths), function(x) rep(x, length(paths[[x]]))))
        walks = split(seg[all.frags$subject.id[paths.u]], paths.i)
        values(walks)$seg.id = split(all.frags$subject.id[paths.u], paths.i)

        ## for now just want to pick the non duplicated paths on the original graph and send these to the walk annotation module
        walks = walks[!duplicated(sapply(values(walks)$seg.id, function(x) paste(x, collapse = ',')))]

        ## note aberrant junction ids and filter out trivial walks that don't employ any ab junctions
        A.ab = sparseMatrix(1, 1, x = as.numeric(NA), dims = dim(jab$adj))
        ab.ix = !is.na(rowSums(rbind(jab$ab.edges[, 1:2,1])))

        A.ref = sign(jab$adj)-sign(A.ab)
        
        
        if (any(ab.ix))
            A.ab[rbind(jab$ab.edges[ab.ix, 1:2, 1])] = A.ab[rbind(jab$ab.edges[ab.ix, 1:2, 2])] = which(ab.ix)

        values(walks)$ab.id = lapply(values(walks)$seg.id, function(x)
            if (length(x)==1) c() else setdiff(A.ab[cbind(x[-length(x)], x[-1])], NA))

        values(walks)$junc.type = lapply(values(walks)$seg.id, function(x)
            if (length(x)==1) c() else sign(A.ref[cbind(x[-length(x)], x[-1])]) + 2*sign(A.ab[cbind(x[-length(x)], x[-1])]))

        walks = walks[!sapply(values(walks)$junc.type, function(x) all(x == 1))]
        
        values(walks)$seg.id = sapply(values(walks)$seg.id, paste, collapse = ',')
        values(walks)$ab.id = sapply(values(walks)$ab.id, paste, collapse = ',')
        values(walks)$junc.type = NULL
        
        if (verbose)
            cat(sprintf('Annotating %s walks\n', length(walks)))        
        
        if (length(walks)==0)
            return(walks)
        else
            {
                names(walks) = 1:length(walks)
                if (!is.null(query))
                    walks = walks[grl.in(walks, query, some = TRUE)]
                return(annotate.walks(walks, cds, promoters, verbose = verbose, exhaustive = FALSE, mc.cores = mc.cores))
            }
    }


####################################
#' .e2class
#' 
#' edge to contig class conversion
#' 
#' given matrix K of k contigs over e edges, each belonging to cardinality 1 or cardinality 2 equivalence classes,
#' assigns id's to equivalent contigs
#' 
####################################
.e2class = function(K, eclass)
  {
    eclass = factor(as.character(eclass))

    if (length(eclass)!=nrow(K))
      stop('eclass must be of the same length as number of rows in K')
    
    eclass = factor(as.character(eclass))
    class.count = table(eclass);
    
    if (any(class.count)>2)
      stop('Edge equivalence classes can have at most 2 members')
    
    biclasses = names(class.count)[class.count==2];  # classes with two edges        
    
    if (length(biclasses)>0)
      {
        # edge class rotation matrix
        R = diag(!(eclass %in% biclasses));  ## edges belonging to classes of cardinality 1 are on the diagonal

        ix = matrix(unlist(split(1:length(eclass), eclass)[biclasses]), ncol = 2, byrow = T); # index pairs corresponding to edges in biclasses
        R[ix[, 1:2]] = 1
        R[ix[, 2:1]] = 1
        
        Kr = R %*% K        
        eix = mmatch(t(Kr), t(K))
        eix[is.na(eix)] = 0
        pairs = t(apply(cbind(1:length(eix), eix), 1, sort))
        pairs = pairs[!duplicated(pairs) & rowSums(pairs==0)==0, , drop = FALSE]

        kclass = rep(NA, ncol(K))
        kclass[pairs[,1]] = 1:nrow(pairs);
        kclass[pairs[,2]] = 1:nrow(pairs);
        kclass[is.na(kclass)] = nrow(pairs) + (1:sum(is.na(kclass)))
      }
    else
      kclass = 1:ncol(K)

    return(kclass)
  }


################################################################################
#' karyoMIP utility functions
#' 
#' 
#' 
################################################################################

##################################
#' convex.basis
#' 
#' Outputs a matrix K of the convex basis of matrix A
#' 
#' i.e. each column x = K[,i] is a minimal solution (with respect to sparsity) to 
#' Ax = 0, x>=0
#' 
#' exclude.basis =  0, 1 matrix of dimension k x ncol(A) specifying k sparsity patterns that we would
#' like to exclude from the convex.basis.  This can speed up computation since any non-negative 
#' combinations of vectors that satisfy an exclusion property will also be excludable, and thus
#' we can remove such vectors as soon as we detect them..
#' 
#' exclude.range = 9, 1 matrix of dimension k x nrow(A) specifying k sparsity patterns that we would like
#' exclude, but these are specified in terms of the range of abs(A) .. i.e. we want to exclude all
#' basis vectors v such that nz(exclude.ranges[i, ]) C  nz(abs(A)*v)) for some pattern i.  Again
#' any non-neg linear comb of any intermediate-basis vector that satisfies this property will satisfy it,
#' as a result we can exclude these vectors when we see them. 
#' 
#' 
#' 
##################################
convex.basis = function(A, interval = 80, chunksize = 100, exclude.basis = NULL, exclude.range = NULL, maxchunks = Inf, 
  verbose = F)
  {
    ZERO = 1e-8;
    remaining = 1:nrow(A);
    iter = 0;
    i = 0;
#    order = c()
    numelmos = c()
    K_i = I = as(diag(rep(1, ncol(A))), 'sparseMatrix');
#    A_i = as(A %*% K_i, 'sparseMatrix');
    K_i = I = diag(rep(1, ncol(A)))
    A_i = A %*% K_i

    if (!is.null(exclude.basis))
      {
        exclude.basis = sign(exclude.basis)
        exclude.basis = exclude.basis[rowSums(exclude.basis)>0, ]
        if (nrow(exclude.basis) == 0)
          exclude.basis = NULL
      }
    
    if (!is.null(exclude.range))
      {
        exclude.range = sign(exclude.range)
        exclude.range = exclude.range[rowSums(exclude.range)>0, ]
        if (nrow(exclude.range) == 0)
          exclude.range = NULL
      }
    
    # vector to help rescale matrix (avoid numerical issues)
    mp  = apply(abs(A), 1, min); # minimum value of each column
    mp[mp[ZERO]] = 1; # columns with zero minimum get scale "1"

    st = Sys.time()
    # iterate through rows of A, "canceling" them out
    while (length(remaining)>0)
      {
        if (nrow(K_i)==0 | ncol(K_i)==0) ## TODO figure out why we have to check this so many times          
          return(matrix())
        
        iter = iter+1;
        K_last = K_i;

        if (verbose)
          print(Sys.time() - st)

        if (verbose)
          cat('Iter ', iter, '(of',  nrow(A_i),  ') Num basis vectors: ', nrow(K_i), " Num active components: ", sum(Matrix::rowSums(K_i!=0)), "\n")

        i = remaining[which.min(Matrix::rowSums(A_i[remaining,, drop = FALSE]>=ZERO)*Matrix::rowSums(A_i[remaining,, drop = FALSE]<=(-ZERO)))]  # chose "cheapest" rows

        remaining = setdiff(remaining, i);
#        order = c(order, i);

        zero_elements = which(abs(A_i[i, ]) <= ZERO);
        K_i1 = K_last[zero_elements, , drop = FALSE];  ## K_i1 = rows of K_last that are already orthogonal to row i of A
        K_i2 = NULL; ## K_i1 = will store positive combs of K_last rows that are orthogonal to row i of A (will compute these below)

        pos_elements = which(A_i[i, ]>ZERO)
        neg_elements = which(A_i[i, ]<(-ZERO))

        if (verbose)
          cat('Iter ', iter, " Row ", i, ":", length(zero_elements), " zero elements ", length(pos_elements), " pos elements ", length(neg_elements), " neg elements \n")
        
        if (length(pos_elements)>0 & length(neg_elements)>0)
          for (m in seq(1, length(pos_elements), interval))
            for (l in seq(1, length(neg_elements), interval))
              {
                ind_pos = c(m:min(c(m+interval, length(pos_elements))))
                ind_neg = c(l:min(c(l+interval, length(neg_elements))))

                indpairs = cbind(rep(pos_elements[ind_pos], length(ind_neg)),
                  rep(neg_elements[ind_neg], each = length(ind_pos))); # cartesian product of ind_pos and ind_neg
                pix = rep(1:nrow(indpairs), 2)
                ix = c(indpairs[,1], indpairs[,2])
#                coeff = c(-A_i[i, indpairs[,2]], A_i[i, indpairs[,1]])  ## dealing with Matrix ghost
                coeff = c(-A_i[i, ][indpairs[,2]], A_i[i, ][indpairs[,1]])  ## 
                combs = sparseMatrix(pix, ix, x = coeff, dims = c(nrow(indpairs), nrow(K_last)))
                combs[cbind(pix, ix)] = coeff;

                H = combs %*% K_last;

                # remove duplicated rows in H (with respect to sparsity)
                H = H[!duplicated(as.matrix(H)>ZERO), ];

                # remove rows in H that have subsets in H (with respect to sparsity) ..
                if ((as.numeric(nrow(H))*as.numeric(nrow(H)))>maxchunks)
                  {
                    print('Exceeding maximum number of chunks in convex.basis computation')
                    stop('Exceeding maximum number of chunks in convex.basis computation')
                  }
                keep = which(Matrix::colSums(sparse_subset(abs(H)>ZERO, abs(H)>ZERO, chunksize = chunksize, quiet = !verbose))<=1) # <=1 since every H is its own subset
                H = H[keep, , drop = FALSE]
                
                # remove rows in H that have subsets in K_i2
                if (!is.null(K_i2))
                  if (nrow(K_i2)>0)
                    {
                      if ((as.numeric(nrow(K_i2))*as.numeric(nrow(H)))>maxchunks)
                        {
                          print('Exceeding maximum number of chunks in convex.basis computation')
                          stop('Exceeding maximum number of chunks in convex.basis computation')
                        }
                      keep = which(Matrix::colSums(sparse_subset(abs(K_i2)>ZERO, abs(H)>ZERO, chunksize = chunksize, quiet = !verbose))==0) 
                      H = H[keep, , drop = FALSE]
                    }
                
                # remove rows in H that have subsets in K_i1
                if (!is.null(K_i1))
                  if (nrow(K_i1)>0)
                    {
                      if ((as.numeric(nrow(K_i1))*as.numeric(nrow(H)))>maxchunks)
                        {
                          print('Exceeding maximum number of chunks in convex.basis computation')
                          stop('Exceeding maximum number of chunks in convex.basis computation')
                        }
                      keep = which(Matrix::colSums(sparse_subset(abs(K_i1)>ZERO, abs(H)>ZERO, chunksize = chunksize, quiet = !verbose))==0) 
                      H = H[keep, , drop = FALSE]
                    }
                
                # maintain numerical stability
                if ((iter %% 10)==0)
                  H = diag(1/apply(abs(H), 1, max)) %*% H
                
#                K_i2 = rBind(K_i2, H)
                K_i2 = rbind(K_i2, as.matrix(H))                
              }

#        K_i = rBind(K_i1, K_i2)
        K_i = rbind(K_i1, K_i2) ## new basis set

        if (nrow(K_i)==0)
          return(matrix())
        
        if (!is.null(exclude.basis)) ## only keep vectors that fail to intersect all vectors "exclude" in matrix
          {
            if ((as.numeric(nrow(exclude.basis))*as.numeric(nrow(K_i)))>maxchunks)
              {
                print('Exceeding maximum number of chunks in convex.basis computation')
                stop('Exceeding maximum number of chunks in convex.basis computation')
              } 
            keep = Matrix::colSums(sparse_subset(exclude.basis>0, K_i>ZERO))==0
            if (verbose)
              cat('Applying basis exclusion and removing', sum(keep==0), 'basis vectors\n')
            K_i = K_i[keep, , drop = F]            
          }
          
        if (!is.null(exclude.range)) ## only keep vectors that fail to intersect all vectors "exclude" in matrix
          {
            A_i_abs = abs(A) %*% t(K_i)
            if ((as.numeric(nrow(exclude.range))*as.numeric*ncol(A_i_abs))>maxchunks)
              {
                print('Exceeding maximum number of chunks in convex.basis computation')
                stop('Exceeding maximum number of chunks in convex.basis computation')
              } 
            keep = Matrix::colSums(sparse_subset(exclude.range>0, t(A_i_abs), quiet = !verbose))==0
            if (verbose)
              cat('Applying range exclusion and removing', sum(keep==0), 'basis vectors\n')
            K_i = K_i[keep, , drop = F]            
          }
        
        A_i = A %*% t(K_i)       
      }
    
    return(t(K_i))       
  }


###############################################
#' collapse.paths 
#' 
#' collapse simple paths in a graph G (adjacency matrix or igraph object)
#' returns m x m new adjacency matrix and map of old vertex id's to new ones
#' $adj = m x m matrix
#' #map = length n with indices 1 .. m
#' 
###############################################
collapse.paths = function(G, verbose = T)
{
  if (inherits(G, 'igraph'))
      G = G[,]

  out = G!=0

  if (verbose)
      cat('graph size:', nrow(out), 'nodes\n')
  
  ## first identify all nodes with exactly one parent and child to do initial collapsing of graph
  singletons = which(Matrix::rowSums(out)==1 & Matrix::colSums(out)==1)

  if (verbose)
      cat('Collapsing simple paths..\n')

  sets = split(1:nrow(G), 1:nrow(G))
  if (length(singletons)>0)
      {
          tmp = out[singletons, singletons]
          cl = clusters(graph(as.numeric(t(which(tmp, arr.ind = TRUE))), n = nrow(tmp)), 'weak')$membership
          dix = unique(cl)
          if (length(dix)>0)
              {
                  for (j in dix)
                      {
                          if (verbose)
                              cat('.')

                          ## grab nodes in this cluster
                          setj = singletons[which(cl == j)]

                          ## move all members into a single set
                          sets[setj[1]] = list(setj)
                          sets[setj[-1]] = list(NULL)

                          ## connect this node to the parent and child of the set                         
                          parent = setdiff(which(Matrix::rowSums(out[, setj, drop = FALSE])>0), setj)
                          child = setdiff(which(Matrix::colSums(out[setj, , drop = FALSE])>0), setj)
                          out[setj, c(setj, child)] = FALSE
                          out[c(setj, parent), setj] = FALSE
                          out[parent, setj[1]] = TRUE
                          out[setj[1], child] = TRUE
                      }
              }
      }

  if (verbose)
      cat('done\nnow fixing branches\n')
  
  todo = rep(FALSE, nrow(G))
  todo[Matrix::rowSums(out)==1 | Matrix::colSums(out)==1] = TRUE

  while (sum(todo)>0)
      {
        sets.last = sets
        out.last = out
          
        if (verbose)
            if ((sum(todo) %% 200)==0)
                cat('todo:', sum(todo), 'num sets:', sum(!sapply(sets, is.null)), '\n')
        
        i = which(todo)[1]
        
        todo[i] = F
        
        child = which(out[i, ])
        parent = which(out[,i])
 
        if (length(child)==1 & length(parent)==1) ## if there is exactly one child and one parent then we want to merge with one or both
            {
                ## if i-child has no other parents and i-parent has no other child
                ## then merge i, i-parent and i-child              
                if (sum(out[,  child])==1 & sum(out[parent, ])==1)  
                    {                      
                        grandch = which(out[child, ])
                        if (length(grandch)>0)
                            {
                                out[parent, grandch] = TRUE  ## parent inherits grandchildren of i
                                out[child, grandch] = FALSE
                            }
                        out[parent, i] = FALSE ## remove node i's edges
                        out[i, child] = FALSE
                        sets[[parent]] = c(sets[[parent]], sets[[child]], sets[[i]])
                        sets[c(i, child)] = list(NULL)
                        todo[child] = F ## no longer have to do i-child, since they have already been merged with parent
                    }
                ## otherwise if either i-child has no other parent or i-parent has no other children (but not both)
                ## then connect i-parent to i-child, but do not merge them (but merge ONE of them with i)
                else if (sum(out[,  child])==1 | sum(out[parent, ])==1)
                    {
                        ## if parent has no other children then merge with him
                        if (sum(out[parent, ])==1)                            
                            sets[[parent]] = c(sets[[parent]], sets[[i]])
                        else
                            sets[[child]] = c(sets[[child]], sets[[i]])

                        out[parent, child] = TRUE
                        out[parent, i] = FALSE ## remove node i's edges
                        out[i, child] = FALSE
                        sets[i] = list(NULL)                        
                    }
            }
        else if (length(child)==1 & length(parent)>1) ## if i has more than one parent but one child, we merge with child if child has no other parents
            {
                if (sum(out[, child])==1)
                    {
                        sets[[child]] = c(sets[[child]], sets[[i]])
                        out[parent, child] = TRUE
                        out[parent, i] = FALSE ## remove node i's edges
                        out[i, child] = FALSE ## remove node i's edges
                        sets[i] = list(NULL)
                    }

        
            }        
        else if (length(child)>1 & length(parent)==1) ## if i has more than one child but one parent, then merge with parent if parent has no other children
            {
                if (sum(out[parent, ])==1)
                    {
                        sets[[parent]] = c(sets[[parent]], sets[[i]])
                        out[parent, child] = TRUE
                        out[parent, i] = FALSE ## remove node i's edges
                        out[i, child] = FALSE ## remove node i's edges
                        sets[i] = list(NULL)
                    }
            }
    }          

  slen = sapply(sets, length)
  ix = which(slen>0)
  map = rep(NA, nrow(G))
  map[unlist(sets)] = match(rep(1:length(sets), slen), ix)
  out = out[ix, ix]
  colnames(out) = rownames(out) = NULL  

  return(list(adj = out, map = map, sets = split(1:length(map), map)))
}


###############################################
#' sparse_subset
#' 
#' given k1 x n matrix A and k2 x n matrix B
#' returns k1 x k2 matrix C whose entries ij = 1 if the set of nonzero components of row i of A is
#' a (+/- strict) subset of the nonzero components of row j of B
#' 
###############################################
sparse_subset = function(A, B, strict = FALSE, chunksize = 100, quiet = FALSE)
  {
    nz = Matrix::colSums(as.matrix(A)!=0, 1)>0
    
    if (is.null(dim(A)) | is.null(dim(B)))
      return(NULL)

    C = sparseMatrix(i = c(), j = c(), dims = c(nrow(A), nrow(B)))
    
    for (i in seq(1, nrow(A), chunksize))
      {
        ixA = i:min(nrow(A), i+chunksize-1)
        for (j in seq(1, nrow(B), chunksize))
          {
            ixB = j:min(nrow(B), j+chunksize-1)

            if (length(ixA)>0 & length(ixB)>0 & !quiet)
              cat(sprintf('\t interval A %s to %s (%d) \t interval B %d to %d (%d)\n', ixA[1], ixA[length(ixA)], nrow(A), ixB[1], ixB[length(ixB)], nrow(B)))
            if (strict)
              C[ixA, ixB] = (sign((A[ixA, , drop = FALSE]!=0)) %*% sign(t(B[ixB, , drop = FALSE]!=0))) * (sign((A[ixA, , drop = FALSE]==0)) %*% sign(t(B[ixB, , drop = FALSE]!=0))>0)
            else
              C[ixA, ixB] = (sign(A[ixA, nz, drop = FALSE]!=0) %*% sign(t(B[ixB, nz, drop = FALSE]==0)))==0
          }
      }
    
    return(C)
  }


######################################################
#' mmatch
#' 
#' Low level utility function to match rows of matrix A to matrix B
#'
######################################################
mmatch = function(A, B, dir = 1, default.value = 0)
{
  nzix = which(A!=default.value, arr.ind = TRUE)
  Adt = as.data.table(nzix)[, v := A[nzix]]
  if (dir == 2)
    setnames(Adt, c('row', 'col'), c('col', 'row'))
  sA = Adt[, paste(col, v, collapse = ' '), by = row]
  setkey(sA, row)
  
  nzix = which(B!=default.value, arr.ind = TRUE)
  Bdt = as.data.table(nzix)[, v := B[nzix]]
  if (dir == 2)
    setnames(Bdt, c('row', 'col'), c('col', 'row'))
  sB = Bdt[, paste(col, v, collapse = ' '), by = row]      
  setkey(sB, V1)

  ix = sB[.(sA[.(1:nrow(A)), ]$V1), unname(row)]

  return(ix)
}

mmatch.og = function(A, B, dir = 1)
  {
    SEP = '|';

    if (is.null(dim(A)))
        A = rbind(A)

    if (is.null(dim(B)))
        B = rbind(B)
    
    if (dim(A)[(dir %% 2)+1] != dim(B)[(dir %% 2)+1])
      stop('Dimensions of A and B matrices mismatch')
      

    if (inherits(A, 'sparseMatrix') | inherits(B, 'sparseMatrix'))
      {
          if (dir == 2)
              {
                  A = t(A)
                  B = t(B)             
              }        
          ixA = which(A!=0, arr.ind = T)
          ixB = which(B!=0, arr.ind = T)

          if (nrow(ixA)>0)
              ixAl = split(1:nrow(ixA), ixA[,1])
          else
              ixAl = c()
          
          if (nrow(ixB)>0)
              ixBl = split(1:nrow(ixB), ixB[,1])
          else
              ixBl = c()
                            
          Atxt = rep('', nrow(A))
          Btxt = rep('', nrow(B))

          if (length(ixAl))
              {
                  tmp.ix = as.numeric(names(ixAl))
                  Atxt[tmp.ix] = sapply(1:length(ixAl), function(x) paste(ixA[ixAl[[x]], 2], ':',
                          as.character(A[tmp.ix[x], ixA[ixAl[[x]], 2], drop = FALSE]), collapse = SEP))
              }
          
          if (length(ixBl)>0)
              {
                  tmp.ix = as.numeric(names(ixBl))
                  Btxt[tmp.ix] = sapply(1:length(ixBl), function(x) paste(ixB[ixBl[[x]], 2], ':',
                          as.character(B[tmp.ix[x], ixB[ixBl[[x]], 2], drop = FALSE]), collapse = SEP))
              }          
      }
    else
      {
        Atxt = apply(A, dir, function(x) paste(x, collapse = SEP))
        Btxt = apply(B, dir, function(x) paste(x, collapse = SEP))
      }
        
    return(match(Atxt, Btxt))
  }


#################################################
#' adj2inc
#' 
#' converts adjacency matrix (of directed graph) into incidence matrix - ie
#' an nodes x edges matrix, for each edge i connecting node j to k, column i will have -1 at position
#' j and +1 at position k
#' 
#################################################
adj2inc = function(A)
  {
    ij = which(A!=0, arr.ind = T)   
    return(sparseMatrix(c(ij[,1], ij[,2]), rep(1:nrow(ij), 2), x = rep(c(-1, 1), each = nrow(ij)), dims = c(nrow(A), nrow(ij))))
  }

#####################################################
#' all.paths
#'
#' Low level function to enumerate all elementary paths and cycles through graph
#' 
#' takes directed graph represented by n x n binary adjacency matrix  A and outputs all cycles and paths between source.vertices, sink.vertices
#' 
#' 
#' @param A nxn adjacency matrix
#' @param all logical flag, if all = T, will include all sources (parentless vertices) and sinks (childless vertices) in path computati
#' @param ALL logical flag, if ALL = T, will also include vertices without outgoing and incoming edges in paths
#' @param sources graph indices to treat as sources (by default is empty)
#' @param sinks graph indices to treat as sinks (by default is empty)
#' @param verbose logical flag
#' @return list of integer vectors corresponding to indices in A (i.e. vertices)
#' $paths = paths indices
#' $cycles = cycle indices
#' @export
#####################################################
all.paths = function(A, all = F, ALL = F, sources = c(), sinks = c(), source.vertices = sources, sink.vertices = sinks,
  exclude = NULL, ## specifies illegal subpaths, all such paths / cycles and
                  ## their supersets will be excluded, specified as k x nrow(A) matrix of vertex sets
  verbose = FALSE,...)
  {
    require(igraph)

    blank.vertices = which(Matrix::rowSums(A)==0 & Matrix::colSums(A)==0)
    
    if (ALL)
      all = T
    
    if (all)
      {
        source.vertices = which(Matrix::rowSums(A)>0 & Matrix::colSums(A)==0)
        sink.vertices = which(Matrix::colSums(A)>0 & Matrix::rowSums(A)==0)
      }

    out = list(cycles = NULL, paths = NULL)
    
    node.ix = which(Matrix::rowSums(A!=0)>0 | Matrix::colSums(A!=0)>0)
    if (length(node.ix)==0)
      return(out)

    A = A[node.ix, node.ix]
    
    if (!is.null(exclude))
      exclude = sign(abs(exclude[, node.ix]))
    
    ij = which(A!=0, arr.ind = T)   
    B = sparseMatrix(c(ij[,1], ij[,2]), rep(1:nrow(ij), 2), x = rep(c(-1, 1), each = nrow(ij)), dims = c(nrow(A), nrow(ij)))    
    I = diag(rep(1, nrow(A)))

    source.vertices = setdiff(match(source.vertices, node.ix), NA)
    sink.vertices = setdiff(match(sink.vertices, node.ix), NA)

    B2 = cBind(B, I[, source.vertices, drop = FALSE], -I[, sink.vertices, drop = FALSE])
    
    if (verbose)
      cat(sprintf('Computing paths for %s vertices and %s edges\n', nrow(B2), ncol(B2)))
    
    K = convex.basis(B2, verbose = verbose, exclude.range = exclude, ...)

    if (all(is.na(K)))
      return(out)
    
    K = K[, Matrix::colSums(K[1:ncol(B), ,drop = FALSE])!=0, drop = FALSE] ## remove any pure source to sink paths
    
    is.cyc = Matrix::colSums(B %*% K[1:ncol(B), ,drop = FALSE]!=0)==0
    
    out$cycles = lapply(which(is.cyc),
      function(i)
      {
        k = which(K[1:ncol(B), i]!=0)
        v.all = unique(as.vector(ij[k, , drop = FALSE]))
        sG = graph.edgelist(ij[k, , drop = FALSE])
        tmp.v = v.all[c(1,length(v.all))]
        p.fwd = get.shortest.paths(sG, tmp.v[1], tmp.v[2])
        p.bwd = get.shortest.paths(sG, tmp.v[2], tmp.v[1])
        return(node.ix[unique(unlist(c(p.fwd, p.bwd)))])
      })

    out$paths = lapply(which(!is.cyc),
      function(i)
      {        
        k = K[1:ncol(B), i]
        eix = which(k!=0)
        v.all = unique(as.vector(ij[eix, , drop = FALSE]))
        sG = graph.edgelist(ij[eix, , drop = FALSE])
        io = B %*% k
        v.in = which(io<0)[1]
        v.out = which(io>0)[1]
        return(node.ix[unlist(get.shortest.paths(sG, v.in, v.out))])
      })

    if (length(out$cycles)>0)
      {
        tmp.cix = cbind(unlist(lapply(1:length(out$cycles), function(x) rep(x, length(out$cycles[[x]])))), unlist(out$cycles))
        out$cycles = out$cycles[!duplicated(as.matrix(sparseMatrix(tmp.cix[,1], tmp.cix[,2], x = 1)))]
      }

    if (length(out$paths)>0)
      {
        tmp.pix = cbind(unlist(lapply(1:length(out$paths), function(x) rep(x, length(out$paths[[x]])))), unlist(out$paths))
        out$paths = out$paths[!duplicated(as.matrix(sparseMatrix(tmp.pix[,1], tmp.pix[,2], x = 1)))]
      }
    
    if (ALL & length(blank.vertices)>0)
      out$paths = c(out$paths, lapply(blank.vertices, identity))
    
    return(out)
  }



#################################################
#' chromoplexy
#'
#' Determines chromoplexy paths from standard JaBbA output
#' 
#' Outputs all chromoplexy paths and cycles
#' (i.e. paths and cycles in breakpoint graph) allowing quasi-reciprocal
#' rearrangements with amplification / deletion bridge distance threshold "dist" 
#' 
#' @param jab JaBbA object
#' @param all logical flag whteher to enumerate all possible cycles, otherwise will return (an arbitrary) minimal decompositoin into the shortest "chains" of balanced rearrangements
#' @param dist maximum distance at which to cluster junctions, i.e. in which to consider a deletion or amplification bridge
#' @param cn.dist minimum distance at which to enforce copy concordance for deletion and amplification bridges
#' @param paths logical flag if paths = T, will also try to compute paths (in addition to cycles), default = FALSE
#' @return paths and cycles of as list of vectors of aberrant edge index sequences,
#' aberrant edges refer to edges described in kag$ab.edges matrix
#' @export
#####################################################
chromoplexy = function(kag = NULL, # output of karyograph
  jab = NULL, ## optional alternate input, if NOT null then this will be used in place of kag
  sol = NULL, ## if sol is null, then copy state is ignored when determining amp or del bridges
  all = F, ## if TRUE, will try to enumerate all possible cycles, otherwise will return (an arbitrary) minimal decomposition into the shortest "chains" of balanced rearrangements
  ref.only = F, ## if T will only compute distance criteria on reference (i.e. won't use any subsequent rearrangements)
  filt.jab = T, ## filter out 0 copy edges if input is a jabba object
  reciprocal = TRUE, ## aka deletion bridge
  hijacked = TRUE,  ## aka amplification bridge
  paths = F, dist = 1e3, cn.dist = dist, verbose = F, interval = 400, chunksize = 5000)
  {
    if (!is.null(jab))
      {
        if (filt.jab)
          {
              nnab = which(rowSums(is.na(rbind(jab$ab.edges[, 1:2, 1])))==0)
              edge.ix = which(jab$adj[rbind(jab$ab.edges[nnab, 1:2, 1])]>0)
              jab$ab.edges = jab$ab.edges[edge.ix, ,,drop = F]
          }
        else
          edge.ix = 1:nrow(kag$ab.edges)
        kag = jab
        sol = jab
      }
    else
      edge.ix = 1:nrow(kag$ab.edges)
        
    G = kag$G
    
    if (is.null(kag$tile))
      kag$tile = kag$segstats

    nnab = which(rowSums(is.na(rbind(kag$ab.edges[, 1:2, 1])))==0)
    if (ref.only)
      {
          adj2 = kag$adj
          adj2[kag$ab.edges[nnab, 1:2, 1]] = 0
          adj2[kag$ab.edges[nnab, 1:2, 2]] = 0
          G = graph(as.numeric(t(which(adj2!=0, arr.ind = T))), n = length(kag$segstats), directed = T)
#        G = graph.adjacency(adj2)
      }
      
    if (verbose)
      cat(sprintf('Running on graph with %s aberrant junctions with dist %s and cn.dist %s with interval %s and chunksize %s\n', nrow(kag$ab.edges), dist, cn.dist, interval, chunksize))

    ## define edge source to edge sink distance
    ## this is minimum between (1) sum of vertex width of path from e2 source to e1 sink (excluding source and sink)
    ## and (2) sum of vertex width of path from e1 sink to e2 source (including source and sink)
    
    tmp = get.edges(G, E(G))
    E(G)$from = tmp[,1]
    E(G)$to = tmp[,2]    
    E(G)$weights.source = width(kag$tile[E(G)$from])
    
    ab.edges = cbind(rbind(kag$ab.edges[nnab, c('from', 'to'), '+'], kag$ab.edges[nnab, c('from', 'to'), '-']), junc.id = rep(nnab,2))

    if (nrow(ab.edges)==0)
      return(list(cycles = NULL, paths = NULL))
    
    #    emap = c(1:nrow(kag$ab.edges), -(1:nrow(kag$ab.edges)))
    emap = c(nnab, -nnab)
    
    D1 = D2 = array(Inf, dim = rep(nrow(ab.edges),2))
    
    uix = unique(c(ab.edges[,1], ab.edges[,2]))
    uixmap1 = match(ab.edges[,1], uix)
    uixmap2 = match(ab.edges[,2], uix)
    
    tmp = shortest.paths(G, uix, uix, weights = E(G)$weights.source, mode = 'out')


    ## deletion bridge, or reciprocal
    if (reciprocal)
    {
        ## "deletion bridge", i.e. source to sink bridge 
        D1 = t(sweep(tmp[uixmap1, uixmap2], 
                     1, width(kag$tile[ab.edges[,1]]))) ## subtract width of first vertex from path length (second vertex already excluded)
        D1[do.call('rbind', lapply(ab.edges[,2], function(x) ab.edges[,1] %in% x))] = NA ## edge case where e1 sink = e2 source
    }
    
    ## "amplification bridge", i.e. sink to source bridge
    if (hijacked)
    {    
        D2 = sweep(tmp[uixmap2, uixmap1],
                   2, -width(kag$tile[ab.edges[,1]])) ## add width of last vertex to path (first vertex already included)
    }

             
    D = matrix(pmin(D1, D2, na.rm = T), nrow = nrow(D1), ncol = nrow(D2))
    D.which = matrix(ifelse(D1<D2, 1, 2), nrow = nrow(D1), ncol = nrow(D2))
    D.which[is.na(D.which)] = 2
    D[is.infinite(D)] = NA
    D[cbind(1:nrow(D), 1:nrow(D))] = NA
    
    ## quasi pairs are ab edge pairs within a certain distance of each other on the graph
    quasi.pairs = which(D<dist, arr.ind = T)
    quasi.pairs.which = D.which[quasi.pairs]
    
    ## now need to check .. depending on whether edge pair is deletion bridge or amp bridge or fully reciprocal
    ## whether associated vertices show a copy change "in the right direction"
    
    ## to do this, we need to examine the vertices "in between" for a deletion bridge and the source / sink vertices
    ## in an amplification bridge, and see if they show a copy change with respect their reference parents
    
    ## for reciprocal pairs, the source and sink will be the same

    adj.ref = kag$adj; adj.ref[ab.edges[, 1:2]] = 0
    
    del.bridge.candidate = which(quasi.pairs.which == 1)
    v1 = ab.edges[quasi.pairs[del.bridge.candidate, 1], 2]
    v1.parent = apply(adj.ref[, v1, drop = FALSE], 2, function(x) which(x != 0)[1])
    v1.child = apply(adj.ref[v1, , drop = FALSE], 1, function(x) which(x != 0)[1])    
    v2 = ab.edges[quasi.pairs[del.bridge.candidate, 2], 1]
    v2.child = apply(adj.ref[v2, , drop = FALSE], 1, function(x) which(x != 0)[1])
    v2.parent = apply(adj.ref[, v2, drop = FALSE], 2, function(x) which(x != 0)[1])
    
    recip = del.bridge.candidate[which(v2.child == v1)]
    nonrecip = which(v2.child != v1) ## these need to fulfill the "deletion bridge criterion"
    
    # test for deletion bridge criterion, i.e. does v1 have greater copy number than its
    # child, and does v2 have greater copy number than its parent?
    if (!is.null(sol))
      del.bridge = del.bridge.candidate[nonrecip[(sol$segstats$cn[v2.child[nonrecip]] < sol$segstats$cn[v2[nonrecip]] & sol$segstats$cn[v1.parent[nonrecip]] < sol$segstats$cn[v1[nonrecip]]) | D[quasi.pairs][del.bridge.candidate[nonrecip]] < cn.dist]]
    else
      del.bridge = del.bridge.candidate
          
    amp.bridge.candidate = which(quasi.pairs.which == 2)
    v1 = ab.edges[quasi.pairs[amp.bridge.candidate, 1], 2]
    v1.parent = apply(adj.ref[, v1, drop = FALSE], 2, function(x) which(x != 0)[1])    
    v2 = ab.edges[quasi.pairs[amp.bridge.candidate, 2], 1]
    v2.child = apply(adj.ref[v2,,drop = FALSE], 1, function(x) which(x != 0)[1])
    
    # test for amp bridge criterion, does v1 have higher copy number than its parent, does v2 have higher copy number than its child?

    if (!is.null(sol))
      amp.bridge = amp.bridge.candidate[(sol$segstats$cn[v1.parent] < sol$segstats$cn[v1] & sol$segstats$cn[v2.child] < sol$segstats$cn[v2])
        | D[quasi.pairs][amp.bridge.candidate] < cn.dist]
    else
      amp.bridge = amp.bridge.candidate
    
    ## now put together all surviving edges into a graph and try to find cycles

    ## store data frame of edge pairs for bp graph
    ## NOTE: every node in bp graph is an edge in the original karyograph, and thus edges in the bp graph represent ordered <edge pairs>    
    bp.df = data.frame(
      e1 = quasi.pairs[c(recip, del.bridge, amp.bridge), 1], e2 = quasi.pairs[c(recip, del.bridge, amp.bridge), 2],
      from = ab.edges[quasi.pairs[c(recip, del.bridge, amp.bridge), 1], 1],
      to = ab.edges[quasi.pairs[c(recip, del.bridge, amp.bridge), 2], 2],
      type = c(rep('recip', length(recip)), rep('del', length(del.bridge)), rep('amp', length(amp.bridge))), stringsAsFactors = F)

    bp.df = bp.df[!is.na(bp.df$e1) & !is.na(bp.df$e2), ]
        
    ## make adj matrix of breakpoints, basically by matching bp1 and bp2 if "to" field of bp1 = "from" field of bp2
    ## here we are looking for <exact> matches because we are now going to join an edge to another edge if the target
    ## of one edge is the source of the next
#    adj.bp = matrix(0, nrow = nrow(bp.df), ncol = nrow(bp.df))
#    for (i in 1:ncol(adj.bp))    
#      adj.bp[i,] = bp.df$from %in% bp.df$to[i] & !is.na(bp.df$to[i])

    ## breakpoint graph links every edge to every other edge via "quasi pair" connection
    ## we find cycles and paths in this graph
    adj.bp = sparseMatrix(i = bp.df$e1, j = bp.df$e2, x = 1, dims = rep(nrow(ab.edges), 2))

    if (verbose)
      if (paths)
        cat(sprintf('Running with paths on breakpoint graph with dim %s vertices and %s edges\n', nrow(adj.bp), sum(adj.bp)))
      else
        cat(sprintf('Running without paths on breakpoint graph with dim %s vertices and %s edges\n', nrow(adj.bp), sum(adj.bp)))
    
    if (prod(dim(adj.bp))>0)
      {
        ## want to exclude any paths involving breaks and their pairs
#        tmp = split(1:nrow(ab.edges), ab.edges[,'junc.id'])
#        exclude.ij = cbind(ab.edges[,3], unlist(tmp))
#        exclude = sparseMatrix(exclude.ij[,1], exclude.ij[,2], x = 1)
        exclude = NULL
        ## dedup junctions with equiv connections (TODO: why are there dups?) to prevent blowup .. in this case only choosing one "path" for each dup edge
        dt = data.table(i = 1:nrow(adj.bp), j = mmatch(adj.bp, adj.bp[!duplicated(as.matrix(adj.bp)), , drop = F]), key = 'j')
        dtu = dt[!duplicated(j), ]
        pc = all.paths(adj.bp[dtu$i, dtu$i, drop = FALSE], all = paths, verbose = verbose, interval = interval, chunksize = chunksize, exclude = exclude)
        if (length(pc$cycles)>0)
            pc$cycles = lapply(pc$cycles, function(x) dtu[x, ]$i)
        if (length(pc$paths)>0)
            pc$paths = lapply(pc$paths, function(x) dtu[x, ]$i)
      }
    else
      return(list(paths = c(), cycles = c()))

    ## hack: this removes some "inconsistent" paths or cycles whose del bridge / amp bridge links
    ## involve aberrant edges that are contained in the respective cycle
    ## TODO: Can try to repair this more elegantly above, though prob no great way to prevent this,
    ## TODO: also this post hoc pruning may (rarely) result in a non exhaustive enumeration of cycles
    ## if there are other possible "bridge links" between members of a cycle that do not involve
    ## members of the cycle.  Fix: Best way to fix this would be actually recompute shortest paths after removing
    ## edges cresponding to edges in the path.  
    .check.pc = function(x, is.cycle = F)
      {
        if (is.cycle)
          tmp.edges = cbind(x, c(x[-1], x[1]))
        else
          tmp.edges = cbind(x[-length(x)], x[-1])        
        tmp.D.which = D.which[tmp.edges]  ## D.which keeps track of whether we linked these edges via D1 or D2
        if (any(ix <- tmp.D.which==1)) ## if 1 then we are looking for path from col 2 to col 1, so flip
          tmp.edges[ix,] = tmp.edges[ix, c(2:1)]
        tmp.ab.edges = cbind(ab.edges[cbind(tmp.edges[,1], tmp.D.which)], ab.edges[cbind(tmp.edges[,2], ifelse(tmp.D.which == 1, 2, 1))])
        tmp.sp = lapply(1:nrow(tmp.ab.edges), function(i)
          get.shortest.paths(G, tmp.ab.edges[i,1], tmp.ab.edges[i,2], weights = E(G)$weights.source, mode = 'out')$vpath[[1]])
        if (any(ix <- tmp.D.which==1))
          tmp.sp[ix] = lapply(tmp.sp[ix], function(x) x[-c(1, length(x))])
        bp.id = unique(unlist(lapply(tmp.sp, function(x) E(G, path = x)$bp.id)))
        return(any(x %in% bp.id))            
        # test distance to make sure
        ## sapply(tmp.sp, function(x, w) sum(w[x[-length(x)]]), as.numeric(width(kag$tile))) - width(kag$tile)[tmp.ab.edges[cbind(1:length(tmp.D.which), tmp.D.which)]]*ifelse(tmp.D.which == 1, 1, -1)
      }
    
    if (length(pc$cycles)>0)
      {
        pc$cycles = pc$cycles[!sapply(pc$cycles, .check.pc, is.cycle = T)]        
        pc$cycles = lapply(pc$cycles, function(x) sign(emap[x])*edge.ix[abs(emap[x])])
        pc$cycles = pc$cycles[!duplicated(sapply(pc$cycles, function(x) paste(unique(sort(x)), collapse = ' ')))]
        pc$cycles = pc$cycles[order(-sapply(pc$cycles, length))]
      }
    
    if (length(pc$paths)>0)
      {
        pc$paths = pc$paths[!sapply(pc$paths, .check.pc, is.cycle = F)]
        pc$paths = lapply(pc$paths, function(x) sign(emap[x])*edge.ix[abs(emap[x])])
        pc$paths = pc$paths[!duplicated(sapply(pc$paths, function(x) paste(unique(sort(x)), collapse = ' ')))]
        pc$paths = pc$paths[order(-sapply(pc$paths, length))]
      }

    
    return(pc)
    
#    return(list(
#                paths = lapply(pc$paths[sapply(pc$paths, length)>1], function(x) cbind(bp.df$e1[x], bp.df$e2[x[length(x)]]), NA),
#                cycles = lapply(pc$cycles, function(x) c(bp.df$e1[x], bp.df$e2[x[length(x)]])),
#                paths = lapply(pc$paths[sapply(pc$paths, length)>1], function(x) bp.df[x, ]),
#                cycles = lapply(pc$cycles, function(x) bp.df[x, ])
#                ))
  }
  


###########################
#' proximity
#' 
#' Takes a set of n "query" elements (GRanges object, e.g. genes) and determines their proximity to m "subject" elements
#' (GRanges object, e.g. regulatory elements) subject to set of rearrangement adjacencies (GRangesList with width 1 range pairs)
#' 
#' This analysis makes the (pretty liberal) assumption that all pairs of adjacencies that can be linked on a karyograph path are in
#' cis (i.e. share a chromosome) in the tumor genome.  
#'
#' @param query GRanges of "intervals of interest" eg regulatory elements
#' @param subject GRanges of "intervals of interest" eg genes
#' @param ra GRangesList of junctions (each a length 2 GRanges, similar to input to karyograph)
#' @param jab existing JaBbA object (overrides ra input)
#' @param verbose logical flag 
#' @param mc.cores how many cores (default 1)
#' @param max.dist maximum genomic distance to store and compute (1MB by default) should the maximum distance at which biological interactions may occur
#' @return
#' list of n x m sparse distance matrices:
#' $ra = subject-query distance in the rearranged genome for all loci < max.dist in tumor genome
#' $wt = subject-query distance in the reference genome for all loci < max.dist in tumor genome
#' $rel = subject-query distance in ra relative to wild type for above loci
#' NOTE: values x_ij in these matrices should be interpreted with a 1e-9 offset to yield the actual value y_ij
#' i.e. y_ij = x_ij-1e-9, x_ij>0, y_ij = NA otherwise (allows for sparse encoding of giant matrices)
#' @export
############################################
proximity = function(query, subject, ra = GRangesList(), jab = NULL, verbose = F, mc.cores = 1, 
  max.dist = 1e6 ## max distance to store / compute in the output matrix.cores
  )
  {
    if (!is.null(jab))
    {
        nnab = which(!ifelse(is.na(jab$ab.edges[,3,1]), TRUE, jab$ab.edges[,3,1]==0))
        ix = nnab[which(jab$adj[jab$ab.edges[nnab,1:2,1]]>0)]
        if (length(ix)>0)
          {
            ra1 = gr.flipstrand(gr.end(jab$segstats[jab$ab.edges[ix,1,1]], 1, ignore.strand = F))
            ra2 = gr.start(jab$segstats[jab$ab.edges[ix,2,1]], 1, ignore.strand = F)
            ra1 = GenomicRanges::shift(ra1, ifelse(as.logical(strand(ra1)=='+'), -1, 0))
            ra2 = GenomicRanges::shift(ra2, ifelse(as.logical(strand(ra2)=='+'), -1, 0))
            ra = grl.pivot(GRangesList(ra1,ra2))
          }
      }
    
    if (length(ra)==0)
      return(list())

    if (length(query)==0 | length(subject)==0)
        return(list())
    
    if (is.null(names(query)))
        names(query) = 1:length(query)

    if (is.null(names(subject)))
        names(subject) = 1:length(subject)
    
    query.nm = names(query);
    subject.nm = names(subject);

    query = query[, c()]
    subject = subject[, c()]

    query$id = 1:length(query)
    subject$id = 1:length(subject)
    
    qix.filt = gr.in(query, unlist(ra)+max.dist) ## to save time, filter only query ranges that are "close" to RA's    
    query = query[qix.filt]

    six.filt = gr.in(subject, unlist(ra)+max.dist) ## to save time, filter only query ranges that are "close" to RA's    
    subject = subject[six.filt]    

    if (length(query)==0 | length(subject)==0)
        return(list())
    
    query$type = 'query'
    subject$type = 'subject'
    
    gr = gr.fix(c(query, subject))
  
    kg = karyograph(ra, gr)              

    ## node.start and node.end delinate the nodes corresponding to the interval start and end
    ## on both positive and negative tiles of the karyograph    
    gr$node.start = gr$node.end = gr$node.start.n = gr$node.end.n = NA;     

    ## start and end indices of nodes
    tip = which(strand(kg$tile)=='+')
    tin = which(strand(kg$tile)=='-')
    gr$node.start = tip[gr.match(gr.start(gr,2), gr.start(kg$tile[tip]))]        
    gr$node.end = tip[gr.match(GenomicRanges::shift(gr.end(gr,2),1), gr.end(kg$tile[tip]))]
    gr$node.start.n = tin[gr.match(GenomicRanges::shift(gr.end(gr,2),1), gr.end(kg$tile[tin]))]
    gr$node.end.n = tin[gr.match(gr.start(gr,2), gr.start(kg$tile[tin]))]

    if (any(nix <<- is.na(gr$node.start)))
        gr$node.start[nix] = tip[nearest(gr.start(gr[nix],2), gr.start(kg$tile[tip]))]
    
    if (any(nix <<- is.na(gr$node.end)))
        gr$node.end[nix] = tip[nearest(GenomicRanges::shift(gr.end(gr[nix],2),1), gr.end(kg$tile[tip]))]
    
    if (any(nix <<- is.na(gr$node.end.n)))
      gr$node.end.n[nix] = tin[nearest(gr.start(gr[nix],2), gr.start(kg$tile[tin]))]

    if (any(nix <<- is.na(gr$node.start.n)))
        gr$node.start.n[nix] = tin[nearest(GenomicRanges::shift(gr.end(gr[nix],2),1), gr.end(kg$tile[tin]))]
        
    
#    gr$node.start = gr.match(gr.start(gr-1,2), gr.start(kg$tile))
#    gr$node.end = suppressWarnings(gr.match(gr.end(gr+1,2), gr.end(kg$tile)))

    ## so now we build distance matrices from query ends to subject starts
    ## and subject ends to query starts

    ## so for each query end we will find the shortest path to all subject starts
    ## and for each query start we will find the shortest.path from all subject ends
    ix.query = which(gr$type == 'query')
    ix.subj = which(gr$type == 'subject')

    node.start = gr$node.start
    node.end = gr$node.end
    node.start.n = gr$node.start.n
    node.end.n = gr$node.end.n
    
    w = width(kg$tile)

    E(kg$G)$weight = width(kg$tile)[E(kg$G)$to]

    ## ix.query and ix.subj give the indices of query / subject in gr
    ## node.start, node.end map gr to graph node ids
    ##
    ## these matrices are in dimensions of query and subject, and will hold the pairwise distances between
    ##
    D.rel = D.ra = D.ref = D.which = Matrix(data = 0, nrow = length(ix.query), ncol = length(ix.subj))

    ## "reference" graph (missing aberrant edges)
    G.ref = subgraph.edges(kg$G, which(E(kg$G)$type == 'reference'), delete.vertices = F)

    EPS = 1e-9

#    for (i in ix.query)
    tmp = mclapply(ix.query, function(i)
      {
        if (verbose)
          cat('starting interval', i, 'of', length(ix.query), '\n')
        
        ## D1 = shortest query to subject path, D2 = shortest subject to query path, then take shortest of D1 and D2
        ## for each path, the edge weights correspond to the interval width of the target node, and to compute the path
        ## length we remove the final node since we are measuring the distance from the end of the first vertex in the path
        ## to the beginning of the final vertex

        u.node.start = unique(node.start[ix.subj]) ## gets around annoying igraph::shortest.path issue (no dups allowed)
        u.node.end = unique(node.end[ix.subj])

        uix.start = match(node.start[ix.subj], u.node.start)
        uix.end = match(node.end[ix.subj], u.node.end)
        
        tmp.D1 = (shortest.paths(kg$G, node.end[i], u.node.start, weights = E(kg$G)$weight, mode = 'out') - w[u.node.start])[uix.start]
        tmp.D2 = (shortest.paths(kg$G, node.start[i], u.node.end, weights = E(kg$G)$weight, mode = 'in') - w[node.start[i]])[uix.end]
        tmp.D3 = (shortest.paths(kg$G, node.end.n[i], u.node.start, weights = E(kg$G)$weight, mode = 'out') - w[u.node.start])[uix.start]
        tmp.D4 = (shortest.paths(kg$G, node.start.n[i], u.node.end, weights = E(kg$G)$weight, mode = 'in') - w[node.start.n[i]])[uix.end]
        tmp.D = pmin(tmp.D1, tmp.D2, tmp.D3, tmp.D4)
        ix = which(tmp.D<max.dist)        
        D.ra[i, ix] = tmp.D[ix]+EPS
        D.which[i, ix] = apply(cbind(tmp.D1[ix], tmp.D2[ix], tmp.D3[ix], tmp.D4[ix]), 1, which.min)

        u.node.start = unique(node.start[ix.subj][ix]) ## gets around annoying igraph::shortest.path issue (no dups allowed)
        u.node.end = unique(node.end[ix.subj][ix])

        uix.start = match(node.start[ix.subj][ix], u.node.start)
        uix.end = match(node.end[ix.subj][ix], u.node.end)
        
        tmp.D1 = (shortest.paths(G.ref, node.end[i], u.node.start, weights = E(G.ref)$weight, mode = 'out') - w[u.node.start])[uix.start]
        tmp.D2 = (shortest.paths(G.ref, node.start[i], u.node.end, weights = E(G.ref)$weight, mode = 'in') - w[node.start[i]])[uix.end]
        tmp.D3 = (shortest.paths(G.ref, node.end.n[i], u.node.start, weights = E(G.ref)$weight, mode = 'out') - w[u.node.start])[uix.start]
        tmp.D4 = (shortest.paths(G.ref, node.start.n[i], u.node.end, weights = E(G.ref)$weight, mode = 'in') - w[node.start.n[i]])[uix.end]
        tmp.D = pmin(tmp.D1, tmp.D2, tmp.D3, tmp.D4)
        D.ref[i, ix] = tmp.D+EPS

        ## if subject and query intersect (on the reference) then we count both RA and Ref distance as 0
        ## (easier to do a simple range query here)
        ix.zero = gr.in(subject[ix], query[i])
        if (any(ix.zero))
          {
            D.ra[i, ix[ix.zero]] = 0
            D.ref[i, ix[ix.zero]] = 0
          }

        D.rel[i, ix] = ((D.ra[i, ix]-EPS) / (D.ref[i, ix]-EPS)) + EPS        
        
        if (verbose)
           cat('finishing interval', i, 'of', length(ix.query), ':', paste(round(D.rel[i, ix],2), collapse = ', '), '\n')

        return(list(D.rel = D.rel, D.ref = D.ref, D.ra = D.ra, D.which = D.which))        
       }, mc.cores = mc.cores)

    for (i in 1:length(tmp))
    {
        if (class(tmp[[i]]) != 'list')
            warning(sprintf('Query %s failed', ix.query[i]))
        else
        {
            D.rel = D.rel + tmp[[i]]$D.rel
            D.ra = D.ra + tmp[[i]]$D.ra
            D.ref = D.ref + tmp[[i]]$D.ref
            D.which = D.which + tmp[[i]]$D.which
        }
    }
    
    ## "full" size matrix
    rel = ra = ref = ra.which =
      Matrix(data = 0, nrow = length(qix.filt), ncol = length(six.filt), dimnames = list(dedup(query.nm), dedup(names(subject.nm))))
    rel[qix.filt, six.filt] = D.rel
    ra[qix.filt, six.filt] = D.ra
    ref[qix.filt, six.filt] = D.ref    
    ra.which[qix.filt, six.filt] = D.which    
    
    ## summary is data frame that has one row for each query x subject pair, relative distance, ra distance, and absolute distance
    tmp = which(rel!=0, arr.ind = T)
    colnames(tmp) = c('i', 'j');
    sum = as.data.frame(tmp)

    if (!is.null(query.nm))
      sum$query.nm = query.nm[sum$i]

    if (!is.null(subject.nm))
      sum$subject.nm = subject.nm[sum$j]
    
    sum$rel = rel[tmp]
    sum$ra = ra[tmp]
    sum$wt = ref[tmp]
    
    sum = sum[order(sum$rel), ] 
    sum = sum[sum$rel<1, ] ## exclude those with rel == 1
    
    ## reconstruct paths
    vix.query = matrix(NA, nrow = length(qix.filt), ncol = 4, dimnames = list(NULL, c('start', 'end', 'start.n', 'end.n')))
    vix.subject = matrix(NA, nrow = length(six.filt), ncol = 4, dimnames = list(NULL, c('start', 'end', 'start.n', 'end.n')))
    vix.query[qix.filt, ] = cbind(values(gr)[ix.query, c('node.start')], values(gr)[ix.query, c('node.start')], values(gr)[ix.query, c('node.start.n')], values(gr)[ix.query, c('node.end.n')])
    vix.subject[six.filt] = cbind(values(gr)[ix.subj, c('node.start')], values(gr)[ix.subj, c('node.start')], values(gr)[ix.subj, c('node.start.n')], values(gr)[ix.subj, c('node.end.n')])


    if (nrow(sum)==0)
        return(NULL)
        
    sum.paths = mcmapply(function(x, y, i)
    {
        if (verbose)
        message('path ', i, ' of ', nrow(sum), '\n')
        if ((ra.which[x, y]) == 1)
            get.shortest.paths(kg$G, vix.query[x, 'end'], vix.subject[y, 'start'], weights = E(kg$G)$weight, mode = 'out')$vpath[[1]]
        else if ((ra.which[x, y]) == 2)
            rev(get.shortest.paths(kg$G, vix.query[x, 'start'], vix.subject[y, 'end'], weights = E(kg$G)$weight, mode = 'in')$vpath[[1]])
        else if ((ra.which[x, y]) == 3)
            get.shortest.paths(kg$G, vix.query[x, 'end.n'], vix.subject[y, 'start'], weights = E(kg$G)$weight, mode = 'out')$vpath[[1]]
        else if ((ra.which[x, y]) == 4)
            rev(get.shortest.paths(kg$G, vix.query[x, 'start.n'], vix.subject[y, 'end'], weights = E(kg$G)$weight, mode = 'in')$vpath[[1]])
    }, sum$i, sum$j, 1:nrow(sum), SIMPLIFY = F, mc.cores = mc.cores)

#    sum$paths = lapply(sum.paths, function(x) x[-c(1, length(x))])
    sum$paths = lapply(sum$paths, as.numeric)
    sum$ab.edges = lapply(sum.paths, function(p) setdiff(E(kg$G, path = p)$bp.id, NA))
    
    return(list(sum = as.data.table(sum), rel = rel, ra = ra, wt = ref, G = kg$G, G.ref = G.ref, tile = kg$tile, vix.query = vix.query, vix.subject = vix.subject))
  }





########################################
#' karyoSim
#' 
#' Simulate (random) evolution of rearrangements according to input junctions, which are provided as a GRangesList, and
#' grouped into "events" by events list (list of numeric vectors or of lists of numeric vectors indexing "junctions")
#' 
#' Goal of the simulation is to instantiate a collection of junctions (+/- approach some copy number profile)
#' through a sequence of rearrangements and whole-chromosome copy changes
#' 
#' Junctions are sequences of signed reference intervals that are contiguous on the
#' on the tumor genome (usually pairs)
#' 
#' Each event consists of either
#' (1) a "quasi reciprocal sequence" (QRS) of junctions, implemented during a single "cell cycle", and are specified by vectors of junction indices,
##    where no indices are repeated (save the last and first, which specifies a cycle)
#' (2) a set of sets of junctions, specified as a list of list of junction indices, again without repetition, corresponding to complex events
#'     spanning multiple QRS's or "cell cycles" e.g. a BFB, which involve a replication step in between each QRS.  The subsequent QRS's
#'     (attempted to be) instantiated in cis to the last item in the previous QRS
#' (3) a GRanges object specifying a reference locus and meta data field $type = "loss" or "gain" specifying one or more pieces of reference genome that should
#'     be lost or gained at a given step.
#' 
#' - Events are interpreted as strings of one or more "quasi reciprocal sequence" (QRS) of junctions
#'   which may be closed / cyclic (if they begin and end with the same junction index) or open. in which case they will result
#'   in at least some interval loss.  We restrict QRS's to contain at most one repeated junction, and this has to be the first
#'   and the last item in the sequence.  "Quasi" reciprocal means that we allow some sequence to be lost or gained in between breaks.
#' - Every QRS is instantiated in the current genome, by mapping junctions, which are specified in haploid reference
#'   coordinates to intervals on the current genome.  By default, instantiation is chosen so that the source interval of 
#'   every junction in the QRS is on the same chromosome as the target interval on the previous junction in the QRS.  If this is the
#'   case, then we say that the current junction  "follows" the previous one in this WQRS instantiation.  Quasi-reciprocality is then applied
#'   by possibly adding intervals at the site of a break (i.e. if the target interval of the previous junction is upstream of the
#'   source interval of the next junction).   In situations where an instance of a subsequent junction cannot be found to followthe current
#'   junction, then the chain is either (1) interpreted as "broken", i.e. equiv to an unbalanced rearrangement (if strict.chain = F or (2)
#'   the event is discarded (if strict.chain = T)
#' - Junctions / events can have many possible instantiations at a given round of evolution.
#'   This is because a given haploid interval on the reference can be associated with many loci on the tumor chromosome
#'   (in the simplest case, two homologues of the same chromosome)
#'   By default, the following preferences are exercised for choosing junction instantiations:
#'   (1) if a chromosome strand can be found that contains all the intervals in the junction (2) a chromosome whose both strands
#'   contain all the intervals in the junction (3) a set of chromosome that instantiates the event as a chain of junctions
#'   these prefereences can be over-ridden by specifying instant.local and instant.chain flags
#' - After every cycle we do a "clean up" which involves (1) rejoining any pairs of broken ends that were partners at the previous
#'   iteration (2) removing any fragments that lack a telomere (if req.tel = T) or lack other req.gr (3) replacing reverse complements
#'   of chromosomes from previous iteration that were rearranged in the previous iteration with the reverse complements of their alteration
#'   products in the current iteration.
#' - Every junction is implemented <exactly once> during the evolutionary history, i.e. lightning does not strike twice, infinite
#'   sites model
#'   
#' p.chrom = prob of chrom event at each simulation step
#' p.chromloss = probability of chromosomal loss | chrom event (default 0.5)
#' p.chromgain = probability of chromosomal gain | chrom event (default 0.5)
#' lambda.chrom = poisson lambda of number of different chromosomes gained or lost at a chromosomal event
#' lambda.chromgain = poisson lambda of number of chromosomes gained at each "gain" event (default lambda.chrom)
#' labmda.chromloss = poisson lambda of number of chromosomes lostd at each "loss" event (default lambda.chrom)
#' 
#' p.wgd = prob of whole genome doubling at each simulation step
#' 
#' Optionally can provide a copy profile cn.profile (GRanges tiling genome with $cn meta data field) and heuristic will be applied
#' to attempt to "evolve" the simulation towards the observed copy profile (to be implemented)
#' 
#' Output is provided as
#' - (if full = F) list with fields
#'                 $chroms = Named GRangesList of final tumor chromosomes in terms of reference intervals
#'                 $gChain = gChain mapping reference to tumor genome
#'                 $cn = gRanges in reference genome of copy counts of reference intervals
#'                 $events = data frame of event indices with field $id (for event id), $desc (see below for description)
#' - (if full = T) List of lists, each item k corresponding to each stage k of evolution and 
#'                 containing the following fields:
#'                 $chroms = Named GRangesList of tumor chromosomes at that stage in terms of reference intervals
#'                 $gChain = gChain mapping reference to current genome k
#'                 $gChain.last = gChain mapping last evolution step to current (from reference in first item of history)
#'                 $cn = gRanges in reference genome of copy counts of reference intervals
#'                 $event = list with $id, $desc that gives the id and description of event
#'                          for chromosomal loss / gain $id = 'chromgain', or 'chromloss', $desc = indices chromosomes
#'                          for ra event, $id event id, desc = junctions involved
#' 
########################################
karyoSim = function(junctions, # GRangesList specifying junctions, NOTE: currently only allows "simple junctions", i.e. locus pairs, eventually will
                               # allow multi (i.e. two or more) range pairs
  events = NULL, # list of integer vectors or lists of integer vetcors corresponding to "events", list item can be GRanges with meta data field $type
                 # with values "loss" or "gain"
  p.chrom = 0, ## probability of chromosomal event at each time step
  p.wgd = 0,  ## conditional prob of wgd given chromosomal event
  p.chromgain = 0.5, ## conditional probability of chromosome gain given not WGD, chromosomal event
  cn = NULL, ## GRanges with 'cn' property
  req.gr = NULL, ## GRanges that every chromosome needs to overlap in order to make it to the next evolution time step
                  ## e.g. centromeres
  req.tel = TRUE, ## logical flag whether to require every chromosome to have telomeres at both ends at every evolution time step
  neo.tel = NULL, ## GRanges specifying intervals that qualify as neo-telomeres, these will only be applied if req.tel = TRUE
  haploid = T, ## tells us whether input genome is haploid, in which case we will begin the simulation with a "genome doubling"
  local.junction = T, ## if T, we prefer to instantiate intervals of a junction "locally" on the same chromosome,
                     # if F we allow all instantiations to be equally likely
  local.qrs = T, ## if T, we prefer QRS instantiations that operate on a single chromosome, if F then all are equally good
                 ## by "prefer", we mean that we score each instantiation, and then choose the best scoring (or "a best", if there
                 ## are ties)
  force.event = T, ## if T, will attempt to implement QRS / event even if QRS can only be instantiated partially or in fragments
  lambda.chrom = 0, lambda.chromgain = lambda.chrom, lambda.chromloss = lambda.chrom,
  full = F, ## full output?
  random.event = T, ## if not random.event, then provided event order will be followed, and blank events will trigger a (random) chromosomal event
  precedence = NULL, ## length(events) x length(events) binary matrix of DAG entries ij representing whether event i occurs before event j
  dist = 1000, ## distance at which to allow deletion breaks
  verbose = T,
  ... # other optional input to chromoplexy()  
  )
  {                  
    kag = karyograph(junctions)
    
    ## check events to make sure kosher
      
    if (!all(sapply(1:length(events), function(x) {
      if (is(events[[x]], 'GRanges'))
        {
          ev = events[[x]]
          if (is.null(values(ev)$type))
            stop('GRanges specifying copy number events must have $type field set to "gain" or "loss"')
          else if (any(!(values(ev)$type %in% c('gain', 'loss'))))
            stop('GRanges specifying copy number events must have $type field set to "gain" or "loss"')
          else
            T
        }      
      else if (is(events[[x]], 'vector'))
        {
          ev = as.numeric(events[[x]])
          if (any(!(abs(ev) %in% 1:length(junctions))))
            stop(sprintf('Event %s has junction index out of bounds', x))
          else
            T
        }
      else if (is.list(events[[x]]))
        all(sapply(1:length(events[[x]], function(y)
               {
                 ev = events[[x]][[y]]
                 if (any(!abs(ev)) %in% 1:length(junctions))
                   stop(sprintf('Event %s, subevent %s has junction index out of bounds', x, y))
                 else
                   T
               })))
    })))
      stop('Some events are of the wrong type')
                         
    ## junctions in terms of graph nodes
    junctions.kg = kag$ab.edges[, c(1:2), ]
        
    if (is.null(events))
      {
        pc = chromoplexy(kag, dist = dist, ...)
        events = c(pc$paths, lapply(pc$cycles, function(x) c(x, x[1])))
        singletons = setdiff(1:nrow(kag$ab.edges), unlist(events)) ## all events that are not part of a path or cycle
        if (length(singletons)>0)
          events = c(events, split(singletons, 1:length(singletons)))
      }
    
    ## this helps us keep track of how many junctions we have accounted for
    ## while choosing events
    cn.event = sapply(events, is, 'GRanges')    
    event2junction = sparseMatrix(i = 1, j = 1, x = 0, dims = c(length(events), nrow(junctions.kg)))
    if (any(!cn.event))
        {        
          ix = cbind(unlist(mapply(function(x,y) rep(x, y), which(!cn.event), sapply(events[!cn.event], length))), unlist(events[!cn.event]))
          ix = ix[!duplicated(ix), ]
          event2junction[ix] = 1
        }
    
    ## events are "done" if we have already used a junction / ra that belongs to that event
    done.events = rep(F, nrow(event2junction))
    done.junctions = rep(F, ncol(event2junction))
    k = 0 ## evolution time step
    
    #### some local utility functions
    ####
    ####
    
    ## makes list mapping reference signed intervals to chromosomal interval coordinates
    .rid2cid = function(intervals)
      {
        out = split(c(1:nrow(intervals), -(1:nrow(intervals))), c(intervals[, '+'], intervals[, '-']))
        out = out[1:max(c(intervals[, '+'], intervals[, '-']))]
        return(out)
      }
    
    ## updates state with chrom gain or loss
    .chrom_change = function(state, gain = NULL, loss = NULL)
      {
        if (is.null(gain))
          gain = c()

        if (is.null(loss))
          loss = c()
        
        keep = !(state$intervals[, 'i'] %in% loss)
        gain = state$intervals[, 'i'] %in% gain;

        gain.intervals = state$intervals[gain, , drop = F]
        gain.intervals[, 'i'] = gain.intervals[, 'i'] + 0.01 ## give these new intervals a unique new chrom name

        ## conatenate and rename
        tmp.intervals = rbind(state$intervals[keep, ], gain.intervals)        
        tmp.intervals = .munlist(lapply(split(1:nrow(tmp.intervals), tmp.intervals[,'i']), function(x) tmp.intervals[x, 3:ncol(tmp.intervals), drop = F]))
        
        out = list(
          intervals = tmp.intervals,
          rid2cid = .rid2cid(tmp.intervals),
          cid2prev = c(which(keep), which(gain))
          )
        
        return(out)
      }

    ## unlists and cbinds matrices (if dim = 2) or rbinds matrices (if dim = 1)
    ## whose first column specifies the list item index of the entry
    ## and second column specifies the sublist item index of the entry
    ## and the remaining columns specifies the value(s) of the vector
    ## or matrices.
    .munlist = function(x, dim = 1)
      {
        if (dim == 2)
          return(t(rbind(i = unlist(lapply(1:length(x), function(y) rep(y, ncol(x[[y]])))),
                         j = unlist(lapply(1:length(x), function(y) 1:ncol(x[[y]]))),
                         do.call('cbind', x))))
        else
          return(cbind(i = unlist(lapply(1:length(x), function(y) rep(y, nrow(x[[y]])))),
                       j = unlist(lapply(1:length(x), function(y) 1:nrow(x[[y]]))),
                       do.call('rbind', x)))
      }
    
    
    ## takes k vectors of length n_1 , ... , n_k and outputs a matrix
    ## of dimension (n_1 x ... x n_k) x k representing their cartesian product
    .cartprod = function(...)
      {
        vecs = list(...)
        if (length(vecs)==0)
          return(NULL)
        out = matrix(vecs[[1]], ncol = 1)
        if (length(vecs)==1)
          return(out)
        if (length(vecs)>1)
          for (i in 2:length(vecs))
            {
              y = vecs[[i]]
              ix = cbind(rep(1:nrow(out), length(y)), rep(1:length(y), each = nrow(out)))
              out = cbind(out[ix[,1], ], y[ix[,2]])
            }
        return(out)
      }
        
    ## main data structure to keep track of current state of chromosomal evolution
    current.state = list(
      intervals = list(),     ## matrix of n signed reference intervals on k chromosomes of current reference genome
                              ## "i" maps current chromosome, and "j" maps position in that chromosome
                              ## '+' col has rids on positive strand and '-' rids on negative strand of current genome
      rid2cid = list(), ## list mapping reference interval ids to signed current ids 
      cid2prev = c()   ## vector mapping current ids (cids) to signed cids of previous genome (prev)
      )
    
    ix = order(as.character(seqnames(kag$tile)), as.character(strand(kag$tile)))
    ix.pos = ix[which(as.logical( strand(kag$tile)[ix]=='+'))]
    ix.neg = ix[which(as.logical( strand(kag$tile)[ix]=='-'))]      
    current.state$intervals = .munlist(mapply(function(x, y) cbind("+"=x, "-"=y),
      split(ix.pos, as.character(seqnames(kag$tile)[ix.pos])),
      split(ix.neg, as.character(seqnames(kag$tile)[ix.neg]))))      
    current.state$rid2cid = .rid2cid(current.state$intervals)
    current.state$cid2prev = 1:nrow(current.state$intervals)
    
    ## map reference intervals to their rev comp
    int2rc = suppressWarnings(match(kag$tile, gr.flipstrand(kag$tile)))

    ## keep track of telomeric reference intervals (todo: specify centromeric intervals as input or other customizable characteristics
    ## that will specify chromosomes that are kept from timepoint to timepoint in the simulation)    
    is.tel = kag$tile$is.tel
    
    if (!is.null(neo.tel))
      is.tel = is.tel || gr.in(kag$tile, gr.stripstrand(neo.tel))
    
    if (!is.null(req.gr))
      is.req = gr.in(kag$tile, gr.stripstrand(req.gr))
    
    ## if we are in haploid land, first step is a "whole genome doubling" that will give us homologues   
    if (haploid)
      current.state = .chrom_change(current.state, gain = unique(current.state$intervals[, 'i']))


    ## history is a list of current.states
    history = list()
    
    step = 0;  ## step in evolution
    
#'    done.events = rowSums(event2junction)==0
    done.events = rep(F, length(events))
    done.this.round = rep(F, length(done.events)) ## also keep track of events that have been done this round
    
    while (!all(done.events | done.this.round)) ## only finish when we have done all events or at least tried them once in this round
      {
        if (verbose)
          cat(sprintf('Evolution step %s\n', step))
            
        history = c(history, list(current.state))
        step = step +1
        
        ## rearrangement event triggered either randomly based on p.chrom
        ## or non-randomly if our event[[step]] is non empty

        if (random.event)
          k = sample(which(!done.events), 1)
        else
          k = step

        force.cn = is(events[[k]], 'GRanges')
                
        ## if random.event = F and event is empty, then will trigger (random) chromosomal loss / gain (see else statement below)
        if (!force.cn & runif(1)>=p.chrom)
          {
            if (verbose)
              cat(sprintf('trying event %s: %s\n', k, sapply(done.events, function(x) paste('(', x, ')', collapse = ', '), collapse = ', ')))
            
            done.this.round[k] = T
            done.events = done.events | rowSums(event2junction[, events[[k]], drop = FALSE])!=0
            last.qrs = NULL ## will store last.qrs in current coordinates if we have a multi-qrs event
            
            this.event = events[[k]]
            
            if (!is.list(this.event))
              this.event = list(this.event)

            qrs.i = 0
            abort = F
            while (qrs.i < length(this.event) & !abort) ## iterate through qrs's in this event
              {
                if (verbose)
                  cat(sprintf('QRS %s of event %s\n', k, qrs.i))
                
                qrs.i = qrs.i+1
                this.qrs = this.event[[i]] ## vector of junctions
                is.cycle = this.qrs[length(this.qrs)] == this.qrs[1]
                
                ## instantiate this QRS
                ## i.e. assign the reference-centric junctions.kg with pairs of intervals in current genome
                
                ## first enumerate all paths involving instantiations of junctions in qrs 
                ## qrs.paths k x m matrix of k QRS instantiations, each consisting of a sequence of
                ## (signed) cids m across n junction id (n < m)
                ##
                ## qrs.juncid maps the columns of qrs paths to junction id in the sequence
                ##
                qrs.paths = array()
                qrs.juncid = c()
                qrs.fid = c() ## keeps track of "fragments" (in case force.event = F, and we are not

                j = 0
                while (j < length(this.qrs) & !abort)
                  {                    
                    j = j+1
                    
                    if (j == 1)
                      {
                        if (verbose)
                          cat(sprintf('junction  %s of qrs %s event %s in step %s\n', j, qrs.i, k, step))

                        ## take cartesian product of all instantiations of node 1 and node 2 in junction
                        qrs.paths = do.call(.matcart, current.state$rid2cid[junctions.kg[abs(this.qrs[j]), ,ifelse(this.qrs[j]>0, 1, 2)]])
                        qrs.juncid = rep(j, dim(qrs.paths)[3]) ## this keeps
                        qrs.fid = rep(1, dim(qrs.paths)[3])

                        ## if we have multi qrs event and last.qrs is defined, then we constrain
                        ## the first event to be in cis (i.e. on the same chromosome) as the previous
                        if (!is.null(last.qrs))
                          qrs.paths = qrs.paths[qrs.paths[, 1] %in% current.state$intervals[abs(last.qrs), 'i'], ]

                        if (verbose)
                          cat(sprintf('\t dim of qrs.paths: (%s, %s)\n', nrow(qrs.paths), ncol(qrs.paths)))

                        if (nrow(qrs.paths)==0)
                          {
                            abort = T
                            break
                          }
                      }
                    else
                      {
                        if (verbose)
                          cat(sprintf('junction %s of qrs %s event %s in step %s\n', j, qrs.i, k, step))
                        
                        qrs.paths.old = qrs.paths
                        tmp = do.call(.cartprod, c(list(1:nrow(qrs.paths)), current.state$rid2cid[junctions.kg[this.qrs[j],]]))
                        qrs.paths = cbind(qrs.paths[tmp[,1], ], tmp[, 2:ncol(tmp)])

                        if (verbose)
                          cat(sprintf('\t dim of qrs.paths: (%s, %s)\n', nrow(qrs.paths), ncol(qrs.paths)))
                        
                        ## ensure that instantiations of source of current .cid in qrs is in cis with previous (including strand)
                        ## we check if the instantation of the first interval in this cid
                        ## is on the same chromosome as the instantiation of the last interval in the
                        ## last cid.  If not, then we throw them out
                        ## TODO: will check whether instantiations are within some threshold distance of each other
                        ## 
                        ## if we run out of instantiations (i.e. the qrs cannot be fully instantiated) then we can 
                        ## either implement a fragmented qrs (if force.event = T) or abort and try a different event
                        ##
                        
                        keep = current.state$intervals[abs(qrs.paths[,length(qrs.juncid)+1]), 'i'] ==
                          current.state$intervals[abs(qrs.paths[,length(qrs.juncid)]), 'i'] ## check which have same chromosome
                        keep = keep & sign(qrs.paths[,length(qrs.juncid)+1]) == sign(qrs.paths[,length(qrs.juncid)]) ## check which have same strand

                        qrs.juncid = c(qrs.juncid, rep(j, length(junctions.kg[this.qrs[j],, '+'])))
                        if (!any(keep)) 
                          {
                            if (!force.event)                              
                              {
                                if (verbose)
                                  cat(sprintf('Aborting event %s at qrs %s\n', k, qrs.i))
                                abort = T
                                break
                              }
                            else ## if can't keep any, then we keep all and just start a new "fragment"
                              qrs.fid = c(qrs.fid, rep(qrs.fid[length(qrs.fid)]+1, length(junctions.kg[this.qrs[j],, '+'])))
                          }
                        else
                          {
                            qrs.paths = qrs.paths[keep, ]
                            qrs.fid = c(qrs.fid, rep(qrs.fid[length(qrs.fid)], length(junctions.kg[this.qrs[j],, '+'])))
                          }                        
                      }
                  }
                
                if (abort)
                  break                    
                
                qrs.score = junction.score = rep(0, nrow(qrs.paths))
                
                if (local.junction | local.qrs) 
                  {
                    ## score how many pairs are on same chromosome
                    tmp = do.call('cbind',
                      lapply(split(1:length(qrs.juncid), qrs.juncid),
                             function(x) current.state$intervals[abs(qrs.paths[,x[1]]), 'i'] == current.state$intervals[abs(qrs.paths[,x[length(x)]]), 'i']))
                    
                    if (local.junction)                      
                      junction.score = junction.score + rowSums(tmp)
                    
                    ## check if entire event is on a single chromosome
                    if (local.qrs)
                      qrs.score = qrs.score + apply(tmp, 1, prod)
                  }

                if (verbose)
                  cat(sprintf('\t final dim of qrs.paths: (%s, %s)\n', nrow(qrs.paths), ncol(qrs.paths)))
                
                ## sort with respect to score, keep best, sample uniformly from best
                ord.ix = order(qrs.score, junction.score, decreasing = T)
                keep = sample(ord.ix[qrs.score[ord.ix] == qrs.score[ord.ix[1]] & junction.score[ord.ix] == junction.score[ord.ix[1]]], 1)

                ## this.qrs.path is vector of signed cid's
                this.qrs.path = qrs.paths[keep,]
                
                ## apply junctions
                
                ## first "check out" the relevant chromosomes, we will replace these in the output genome
                ## these are only chromosomes on which junction instantiations (instantations of intervals
                ## at ends of junctions) lie, and not any internal junction interval
                tmp.ix = which(diff(c(0, qrs.juncid, 0))!=0)
                internal = !c(1:length(qrs.juncid) %in% c(tmp.ix, tmp.ix-1))
                qrs.adj = !(1:length(current.state$intervals) %in% c(tmp.ix, tmp.ix + 1))
                                
                ## make new "fragments" data structure initially representing
                ## all strands of "checked out" chromosomes, which we will break and join
                
                ## chroms and strands to "check out" of the current genome
                check.out = cbind(chr = current.state$intervals[abs(this.qrs.path[!internal]), 'i'], str = sign(this.qrs.path[!internal]))
                check.out = check.out[!duplicated(check.out), ]

                if (verbose)
                  cat(sprintf('\t checked out chroms: (%s)\n', paste(sign(this.qrs.path)[!internal]*current.state$intervals[abs(this.qrs.path[!internal]), 'i'], collapse = ',')))

                ## cook up list of fragment ids
                tmp.fid = split(cumsum(unlist(lapply(current.state$intervals[check.out[,1]], function(x) rep(1, length(x))))),
                  unlist(lapply(1:nrow(check.out), function(x) rep(x, ncol(current.state$intervals[[check.out[x,1]]])))))

                ## fragments is n x 4 matrix representing n "checked out" single strand DNA intervals across k fragments
                ## and their fragment number (i), fragment pos (j), signed current genome interval id (cid),
                ## flags specifying whether it is an amp bridge (is.bridge) and to be broken (to.break)
                ##                
                fragments = .munlist(mapply(function(chr, str) {
                  cid = which(current.state$intervals[, 'i'] == chr)
                  if (str==1)
                    cbind(cid = cid, to.break = F, is.bridge = F)
                  else ## rev prev strands
                    cbind(cid = -cid, to.break = F, is.bridge = F)
                  }, check.out[, 'chr'], to.break = F, check.out[, 'str']), dim = 1)

                if (verbose)
                  cat(sprintf('\t checked out %s fragments comprising %s intervals on %s chromosomes\n', length(unique(fragments[, 'i'])), nrow(fragments)), length(unique(check.out[, 'chr'])))
               
                ## this k x 2 matrix will keep track of left and right (current genome) neighbors of fragments
                ## "left" and "right" store either NA or the fragment number of the neighbor
                fragment.partners = array(NA, dim = c(unique(fragments[, 'i']), 2), dimnames = list(NULL, c('left', 'right')))

                ## current2frag = n x 2 matrix mapping signed cid --> fid
                ## we also keep mapping of current genome signed intervals to (non bridge) fragment interval,
                ## within a single QRS, this mapping will be one to one, since the only signed intervals
                ## that get duplicated are within amp bridges.
                ## 
                ## note: BFB's are implemented by connecting signed intervals to their reciprocal, which will result
                ## in duplication only after we "strand complete" the fragments produced by the QRS,
                ## however, during the implementation of the QRS, there will be a one to one mapping
                ## between current genome signed intervals and single stranded fragments
                ##
                current2frag = matrix(NA, nrow = length(current.state$cid2prev), ncol = 2,
                  dimnames = list(NULL, c('+', '-')))
                current2frag[cbind(abs(fragments[, 'cid']), ifelse(fragments[, 'cid']>0, 1, 2))] = 1:nrow(fragments)
                current2frag.nna = !is.na(current2frag) ## will be useful for updating this
                
                is.start = c(T, diff(qrs.fid)!=0) ## vector of qrs fragment starts (will only be one T entry if force.event = F)
                qrs.iter = split(1:length(qrs.juncid), qrs.juncid)
                                
                k = 0
                ## iterate through all the adjacencies / junctions in this qrs
                for (k in 1:(length(qrs.iter) - !is.cycle)) ## stop at next to last if is.cycle
                  {                
                    ## if first in qrs (or first in qrs fragment), then apply both breaks, and no amp bridge
                    ## if last in qrs cycle, then apply no breaks, and possibly two amp bridges (from previous, to first)
                    ## otherwise apply one break, and possibly one amp bridge (from previous)

                    if (verbose)
                      cat(sprintf('current fragments: \n%s\n',
                                  paste(paste('[', sapply(split(fragments[, 'cid'], fragments[, 'i']), paste, collapse = ' '), ']', sep = ''), collapse = '\n')))
                    
                    fragments[, 'to.break'] = F
                    
                    ## junction.cid = signed cid of instantiations of current junction intervals on genome
                    this.junction = this.qrs.path[qrs.iter[[k]]]

                    ## we locate their fragment locations (flocs) fragment_id fragment_pos
                    this.junction.fids = current2frag[cbind(abs(this.junction), ifelse(this.junction>0, 1, 2))]
                                        
                    if (is.start[k])
                      {
                        fragments[this.junction.fids, 'to.break'] = T

                        if (verbose)
                          {
                            frag.ix = fragments[, 'i'] == fragments[this.junction.fids, 'i']
                            frag1 = fragments[frag.ix[fragments[frag.ix, 'j'] <= fragments[this.junction.fids, 'j']], 'cid']
                            frag2 = fragments[frag.ix[fragments[frag.ix, 'j'] > fragments[this.junction.fids, 'j']], 'cid']
                            cat(sprintf('[%s] --> [%s], [%s]\n',
                                        paste(fragments[frag.ix, 'cid'], collapse = ' '),
                                        paste(frag1, collapse = ' '),
                                        paste(frag2, collapse = ' ')))
                          }                                                
                      }
                    else  
                      {
                        ## if not start, check for backward amp bridge between target
                        ## adj of previous junction and source adj of current
                        ##
                        ## amp bridge occurs only if last adj targ is upstrand of this edge source on current genome
                        last.adj.targ.qix = qrs.iter[[k-1]][length(qrs.iter[[k-1]])]
                        last.adj.targ = this.qrs.path[, last.adj.targ.qix]
                        this.adj.source = this.qrs.path[, qrs.iter[[k]][1]]

                        ## sanity check: last adj targ and this edge source should be on the same
                        ## chromosome and strand in the current genome
                        if (!(last.adj.targ[1] == this.adj.source[1] & last.adj.targ[3] == this.adj.source[3]))
                          stop('something wrong: last.adj.targ and current.adj are not on same chr and strand')

                        ## amp.bridge only if targ at or before source, which is left on pos strand and right on neg 
                        is.amp.bridge = ((last.adj.targ[3] == 1 & last.adj.targ[2] <= this.adj.source[2]) |
                                         (last.adj.targ[3] == 2 & last.adj.targ[2] >= this.adj.source[3]))
                        
                        ## backward amp bridge will be added to fragment opposite last target
                        ## amp bridge consists of intervals between last target and current source (inclusive)
                        if (is.amp.bridge)
                          {                                                                                                                
                            ## find interval where to add amp.bridge
                            ## this interval is opposite last target in current genome
                            ## which is left of last.adj.targ for neg strand and right of last.adj.targ for positive
                            this.chr = which(current.state$intervals[, 'i'] == last.adj.targ[1])
                            to.add.fid = current2frag[this.chr[last.adj.targ[2]], last.adj.targ[3]]
                            
                            ## make amp.bridge frag (i j cid is.bridge to.break)
                            amp.bridge.frag = cbind(cid = c(1, -1)[last.adj.targ[3]] * this.chr[last.adj.targ[2]:this.adj.source[2]])
                          
                            amp.bridge.frag = cbind(
                              i = fragments[to.add.fid, 'i'],
                              j = fragments[to.add.fid, 'j'] + 1:nrow(amp.bridge.frag), ## we are adding to the right so j will be the new index
                              amp.bridge.frag, is.bridge = T, to.break = F)
                          
                            if (verbose)
                              cat(sprintf('Making amp bridge %s\n', amp.bridge.frag))
                            
                            ## sanity check: is there something wrong, i.e. the interval that we are adding to
                            ## is not at the right end of its fragment
                            if (to.add.fid != nrow(fragments))
                              if (fragments[to.add.fid, 'i'] == fragments[to.add.fid+1, 'i'])
                                stop('something is wrong: right amp bridge to be added to internal fragment')
                            
                            ## add amp.bridge to right of to.add.fid in fragments
                            aft.ix = c()                            
                            bef.ix = 1:to.add.fid
                            if (to.add.fid != nrow(fragments))
                              aft.ix = (to.add.fid+1):nrow(fragments)
                            fragments = cbind(fragments[bef.ix, ], amp.bridge, fragments[aft.ix, ])
                            
                            ## update current2frag
                            pix = c(bef.ix, rep(NA, nrow(amp.bridge)), aft.ix)
                            current2frag[current2frag.nna] = match(current2frag[current2frag.nna], pix)
                            
                            ## remove any reference neighbors from the right side of this frag
                            ## (i.e. it can only be attached through a subsequent adjacency)
                            fragment.partners[flocs[to.add.fid,'i'], 'right'] = NA 
                          }                                                                        
                      }
                   
                    if (!(is.cycle & k == (length(qrs.iter)-1))) ## unless next to last in cycle, break target adjacency in junction
                      {
                        tmp.fid = this.junction.fids[length(this.junction.fids)]+1
                        if (fragments[tmp.fid, 'i'] != fragments[tmp.fid-1, 'i'])
                          stop('Something is wrong, we are breaking at the beginning of a fragment')                                                    
                        fragments[tmp.fid, ] = T
                      }
                    else 
                      {
                        ## if we are next to last iteration of cycle QRS,
                        ## check for forward amp bridge between target
                        ## adj of this junction and source adj of next (i.e. first)
                        ##
                        ## amp bridge occurs only if this adj target is upstrand of next edge source
                        ## on current genome
                        this.adj.targ = this.qrs.path[, qrs.iter[[k]][1]]
                        next.adj.source = this.qrs.path[, qrs.iter[[k+1]][1]]
                        
                        ## sanity check: this adj targ and next edge source should be on the same
                        ## chromosome and strand in the current genome
                        if (!(next.adj.source[1] == this.adj.targ[1] & next.adj.source[3] == this.adj.targ[3]))
                          stop('something wrong: next.adj.source and this.adj.targ are not on same chr and strand')

                        ## amp.bridge only if source at or before targ, which is left on pos strand and right on neg 
                        is.amp.bridge = ((this.adj.targ[3] == 1 & this.adj.targ[2] <= next.adj.source[2]) |
                                         (this.adj.targ[3] == 2 & this.adj.targ[2] >= next.adj.source[3]))
                                                
                        ## forwards amp bridge extends fragment opposite <next> source with intervals from 
                        ## <this target> until <next source>
                        if (is.amp.bridge)
                          {
                            ## find interval where to add amp.bridge
                            ## this interval is opposite last target in current genome
                            ## which is right of next.adj.source for neg strand and left of next.adj.source for positive
                            this.chr = which(current.state$intervals[, 'i'] == next.adj.source[1])
                            tmp.cid = ifelse(next.adj.source[3] == 1, this.chr[next.adj.source[2]-1], this.chr[next.adj.source[2]+1])
                            to.add.fid = current2frag[tmp.cid, next.adj.source[3]]
                                                        
                            ## make amp.bridge frag (rid fid is.bridge to.break)
                            amp.bridge.frag = cbind(cid = c(1, -1)[this.adj.targ[3] * this.adj.targ[3]] * this.chr[this.adj.targ[2]:next.adj.source[2]])

                            amp.bridge.frag = cbind(
                              i = fragments[to.add.fid, 'i'],
                              j = 1:nrow(amp.bridge.frag), ## we are adding to the left 
                              amp.bridge.frag, is.bridge = T, to.break = F)
                                                        
                            ## check to see if there is something wrong
                            if (fragments[to.add.fid, 'j'] != 1)
                              warning('something is wrong: left amp bridge to be added to internal fragment')

                            ## shift j on the current fragments
                            tmp.ix = fragments[, 'i'] == fragments[to.add.fid, 'i']
                            fragments[tmp.ix, 'j'] = fragments[tmp.ix, 'j'] + nrow(amp.bridge.frag) 
                            
                            ## add amp.bridge to <left> of to.add.fid in fragments                                                        
                            bef.ix = c()
                            if (to.add.fid != 1 )
                              bef.ix = 1:(to.add.fid-1)
                            aft.ix = to.add.fid:nrow(fragments)
                            
                            fragments = cbind(fragments[bef.ix, ], amp.bridge, fragments[aft.ix, ])
                            
                            pix = c(bef.ix, rep(NA, nrow(amp.bridge)), aft.ix)
                            current2frag[current2frag.nna] = match(current2frag[current2frag.nna], pix)
                            
                            ## remove reference partners from left side of this frag
                            fragment.partners[flocs[to.add.fid, 'i'], 'left'] = NA
                          }
                      }

                    ### so far, we have applied no breaks in <this> QRS iteration, we have only selected fragments to break
                    ### (and appended amp bridges to fragments broken in previous iterations)

                    ## to apply a break, we split and re-.munlist the old fragments list using to.break field and update fragment.partners

                    ## obtaining new fragments only requires relabeling the 'i' fields
                    old.fragments = fragments;
                    fragments = .munlist(lapply(split(1:ncol(fragments), fragments[, 'i'] + fragments[, 'to.break']/10),
                      function(x) fragments[x, ]), dim = 1)

                    ## updating fragment.partners requires letting separated partners
                    ## remember who they were just joined to, so they can be later rejoined
                    ## if they don't find a new partner
                    ## (unless they are amp bridges or their partner was re-fused)
                    new2old.frag = rep(NA, nrow(fragment.partners))
                    new2old.frag[fragments[,'i']] = old.fragments[, 'i']                    
                    fragment.partners = fragment.partners[new2old.frag, ]
                    fragment.partners[fragments[fragments[, 'to.break'], 'i'], 'right'] = fragments[fragments[, 'to.break'], 'i']+1
                    tmp.ix = which(fragments[, 'to.break'])+1
                    fragment.partners[fragments[tmp.ix, 'i'], 'left'] = fragments[tmp.ix, 'i'] - 1
                             
                    ## to apply fusion, we only have to move rows around in the fragments matrix
                    ## moving the left fragment in the fusion immediately before the right fragment in the fusion

                    ## update this junction.fids to current fragments matrix
                    this.junction.fids = current2frag[this.junction[, c('cid','str')]] 

                    ## junction connects fragment ending in interval junction.fids[1] to the fragment
                    ## beginning with junction.fids[length(junction.fids)] (with any other intervals connected in between)
                    ## (NOTE: current implementation treats junction only as interval pair)
                    left.ix = range(which(fragments[, 'i'] == fragments[junction.fids[1], 'i']))
                    right.ix = range(which(fragments[, 'i'] == fragments[junction.fids[length(junction.fids)], 'i']))

                    old.fragments = fragments;
                    fragments[c(left.ix, right.ix), ] =  fragments[left.ix[1], 'i'] ## rename fused fragment according to the left partner
                    pix = c(c(left.ix, right.ix), (1:nrow(fragments))[-c(left.ix, right.ix)]) ## save permutation

                    ## apply permutation to fragments and current2frag
                    fragments = fragments[pix, ]
                    current2frag[current2frag.nna] = match(current2frag[current2frag.nna], pix)

                    ## rename fragments so they are numbered in order (to make compatible with fragment partners)
                    fragments[, 'i'] = 1 + cumsum(c(0, diff(fragments[, 'i'])))

                    ## let the new fragment inheriting the name of the left partner inherit the partner of the right partner
                    fragment.partners[old.fragments[left.ix[1], 'i'], 'right'] = fragment.partners[old.fragments[right.ix[1], 'i'], 'right'] 

                    ## rewire fragment partners so that new fragment is first
                    fragment.partners = rbind(fragment.partners[old.fragments[left.ix[1], 'i'], ],
                      fragment.partners[-(old.fragments[c(left.ix[1], right.ix[1]), 'i']), ])
                  }

                ## now we have finished processing QRS
                
                ############
                ## Rejoin
                ############
                    
                ## match up and rejoin "left" and "right" partners
                ## to.rejoin is adj.matrix of directed graph that should only have paths
                to.rejoin = matrix(0, nrow = nrow(fragment.partners), ncol = nrow(fragment.partners))
                to.rejoin[cbind(fragment.partners[, 'left'],
                                fragment.partners[, 'right'][match(fragment.partners[, 'left'], fragment.partners[, 'right'])])] = 1            
                
                sources = which(colSums(to.rejoin,1)==0 & rowSums(to.rejoin,1)>0)
                rejoined.paths = list();
                for (s in 1:length(sources))
                  {                        
                    rejoined.paths[[s]] = sources[s]
                    next.v = which(to.rejoin[rejoined.paths[[s]][length(rejoined.paths[[s]])], ]!=0)
                    while (length(next.v)>0)
                      {                            
                        rejoined.paths[[s]] = c(rejoined.paths[[s]], next.v)
                        next.v = which(to.rejoin[rejoined.paths[[s]][length(rejoined.paths[[s]])], ]!=0)
                      }
                  }
                
                ## all other fragments will be in single fragment paths
                tmp.ix = setdiff(1:nrow(fragments.partners), unlist(rejoined.paths))
                non.rejoined.paths = split(tmp.ix, 1:length(tmp.ix))
                
                ## now permute and relabel fragment matrix according to paths
                all.paths = c(rejoined.paths, non.rejoined.paths)

                new.frag.ix = mapply(function(x, y) unlist(x[y]), split(1:nrow(fragments), fragments[, 'i']), all.paths, SIMPLIFY = F)
                fragments = .munlist(lapply(new.frag.ix, function(x) fragments[x, ]), dim = 1)

                ## update current2frag
                current2frag[current2frag.nna] = match(current2frag[current2frag.nna], unlist(new.frag.ix))
                
                ############
                ## Clean up
                ############
                
                ## now remove fragments that (1) do not have a "legal" telomeric interval at both ends and (2) (if is.req is not NULL)
                ## do not contain and do not contain an is.req interval 
                frag.rid = split(fragments[, 'rid'], fragments[, 'i'])
                
                keep = sapply(frag.rid, function(x) is.tel[x[1]] & is.tel(x[length(x)]))                
                if (!is.null(is.req))
                  keep = keep & sapply(frag.rid, function(x) any(is.req[x]))

                new.ix = fragments[,'i'] %in% which(keep)
                
                fragments = fragments[new.ix, ]
                current2frag[current2frag.nna] = match(current2frag[current2frag.nna], new.ix)

                if (verbose)
                  cat(sprintf('current fragments: \n%s\n',
                              paste(paste('[', sapply(split(fragments[, 'cid'], fragments[, 'i']), paste, collapse = ' '), ']', sep = ''), collapse = '\n')))
                
                ############
                ## update current.state
                ############
                
                ## we essentially are doing "DNA replication", since we are pulling both strands of fragments
                ## that we have joined here 

                ## double stranded intervals are format i j + -
                remove = current.state$intervals[, 'i'] %in% check.out[, 'chr'] ## mark current state chromosomes to remove

                ## construct intervals to check in 
                intervals.check.in = cbind(
                  i = fragments[, 'i']+0.1, ## provide fake chromosome label
                  j = fragments[, 'j'],
                  '+'  = sapply(fragments[, 'cid'], function(x) if (x>0) current.state$intervals[abs(x), '+'] else current.state$intervals[abs(x), '-']),
                  '-'  = sapply(fragments[, 'cid'], function(x) if (x>0) current.state$intervals[abs(x), '-'] else current.state$intervals[abs(x), '+'])
                  )

                ## renumber chromosomes
                tmp.intervals = rbind(intervals.check.in, current.state$intervals[!remove,])
                new.intervals = .munlist(split(1:nrow(tmp.intervals), tmp.intervals[, 'i']), function(x) tmp.intervals[x, c('+', '-')])

                if (verbose)
                  cat(sprintf('----------------------------\nChecking in chromosomes: \n%s\n-----------------------------\n',
                              paste(sapply(split(1:nrow(intervals.check.in), fragments[, 'i']),
                                           function(x) paste(new.intervals[x, '+'], new.intervals[x, '-'], sep = '\t', collapse = '\n'),
                                           collapse = '\n==========================\n'))))
                
                ## construct new state, including mapping to previous id's
                current.state = list(                  
                  intervals = new.intervals,
                  rid2cid = .rid2cid(new.intervals)
                  )
                
                ## if this is the second or later qrs in a multi-qrs event, then maintain "cid2prev" to map to the
                ## last current state in the history
                if (qrs.i>1) 
                  current.state$cid2prev = current.state$cid2prev[c(fragments[, 'cid'], which(!remove))]
                else
                  current.state$cid2prev = c(fragments[, 'cid'], which(!remove))

                ## map last element in qrs.paths to new current.state cid
                ## this is where the next event in a multi qrs event will attempt to be instantiated
                last.qrs = current2frag[qrs.paths[length(qrs.paths)], ifelse(qrs.paths[length(qrs.paths)]>0, 1, 2)]                
              }
            
            ## if we have aborted, then roll history back to last current.state
            ## since done.this.round[k] = T, we won't redo this event in this round
            if (abort)
              {
                current.state = history[[length(history)]]
                history = history[[-length(history)]]
                step = step - 1
              }
            else ## otherwise, update done.events with all events that intersect with an junction covered by this event
              {
                done.events = done.events | rowSums(event2junction[, which(event2junction[k, ])]) != 0
                done.this.round = rep(F, length(done.events)) ## reset set of events that have been done this round
              }            
          }
        else ## do a chromosomal event
          {
            if (force.cn) ## then need to instantiate GRanges event
              {
                this.event = events[[k]]

                if (verbose)
                  cat(sprintf('Implementing event %s (CN event) on chroms %s\n', k, unique(as.character(seqnames(this.event)))))
                
                if (is.null(this.event$type))
                  stop('All GRanges specifying copy number events must have $type field set to "gain" or "loss"')

                ev2tile = gr.findoverlaps(this.event, kag$tile, pintersect = T)
                
                ## match each row of this event to possible instantiations
                ev2cid = lapply(split(current.state$rid2cid[ev2tile$subject.id], ev2tile$query.id), function(x) do.call('c', x))
                ev2cid.type = this.event$type[unique(ev2tile$query.id)] ## gain or loss                

                ## now pick a random instantiation per event
                ev.cid = lapply(ev2cid, function(x) if (length(x)>0) x[sample(length(x))] else c())

                gain.cid = loss.cid = c()
                ev.gain = which(ev2cid.type == 'gain')
                ev.loss= which(ev2cid.type == 'loss')

                ## choose randomly among instantiations to choose chromosomes to lose or gain
                ## (choosing with replacement,so if a cid exists in multiple rows in this.event then they
                ## may be both chosen)
                if (length(ev.gain)>0)
                  gain.cid = unlist(lapply(ev.cid[ev.gain], function(x) if (length(x)>0) abs(x[sample(length(x),1)]) else c()))
                
                if (length(ev.loss)>0)
                  loss.cid = unlist(lapply(ev.cid[ev.loss], function(x) if (length(x)>0) abs(x[sample(length(x),1)]) else c()))

                gain  = current.state$intervals[gain.cid, 'i']
                loss  = current.state$intervals[loss.cid, 'i']
                
                current.state = .chrom_change(current.state, gain = gain, loss = loss)

                done.events[k] = T
              }
            else ## otherwise, random event 
              {            
                if (runif(1)<p.wgd & !force.cn)
                  current.state = .chrom_change(current.state, gain = 1:length(current.state$intervals))
                else
                  {
                    if (runif(1)<p.chromgain)
                      {
                        num.chrom = pmin(length(current.state$intervals), rpois(1, lambda = lambda.chromgain))
                        prob = sapply(current.state$intervals, function(x, w) sum(w[x[1, ]]), w = as.numeric(width(kag$tile)))
                        gain = sample(1:length(current.state$intervals), replace = F, prob = prob, size = num.chrom)
                        current.state = .chrom_change(current.state, gain = gain)
                      }
                    else
                      {
                        num.chrom = pmin(length(current.state$intervals), rpois(1, lambda = lambda.chromloss))
                        prob = sapply(current.state$intervals, function(x, w) sum(w[x[1, ]]), w = as.numeric(width(kag$tile)))
                        loss = sample(1:length(current.state$intervals), replace = F, prob = prob, size = num.chrom)
                        current.state = .chrom_change(current.state, loss = loss)
                      }
                  }
              }
          }
      }
    if (full)
      return(history)
    else
      return(current.state)
  }
          


#################
#################


####################################################################
#' ppgrid
#' 
#' least squares grid search for purity and ploidy modes
#'
#' @param segstats GRanges object of intervals with meta data fields "mean" and "sd" (i.e. output of segstats function)
#' @param allelic logical flag, if TRUE will also look for mean_high, sd_high, mean_low, sd_low variables and choose among top solutions from top copy number according to the best allelic fit
#' @param mc.cores integer number of cores to use (default 1)
#' @return data.frame with top purity and ploidy solutions and associated gamma and beta values, for use in downstream jbaMIP
#' @export
############################################
ppgrid = function(
  segstats, # n x 1 GRanges object with "mean" and "sd" value fields, optional field $ncn for "normal tissue" cn (default = 2)

  ########### optional args to describe the "valid modes"
    allelic = FALSE, ## if TRUE will also look for mean_high, sd_high, mean_low, sd_low variables and choose among top solutions from top copy number according to the best allelic fit
    purity.min = 0.01, 
    purity.max = 1.0,
    ploidy.step = 0.01,
    purity.step = 0.01,
    ploidy.min = 1.2, # ploidy bounds (can be generous)
    ploidy.max = 6,
    plot = F,
    verbose = F,
    mc.cores = 10
  )
{
    
  if (verbose)
      cat('Setting up matrices .. \n')

    if (is.na(ploidy.min)) ploidy.min = 1.2
    if (is.na(ploidy.max)) ploidy.max = 6
    if (is.na(purity.min)) purity.min = 0.01
    if (is.na(purity.max)) purity.max = 1
    
    ##  purity.guesses = seq(0, 1, purity.step)
    purity.guesses = seq(pmax(0, purity.min), pmin(1.00, purity.max), purity.step)
    ## ploidy.guesses = seq(pmin(0.5, ploidy.min), pmax(10, ploidy.max), ploidy.step)
    ploidy.guesses = seq(pmax(0.5, ploidy.min), pmax(0.5, ploidy.max), ploidy.step)
    
    if (allelic)
        if (!all(c('mean_high', 'mean_low', 'sd_high', 'sd_low') %in% names(values(segstats))))
        {              
            warning('If allelic = TRUE then must have meta data fields mean_high, mean_low, sd_high, sd_low in input segstats')
            allelic = FALSE
        }
    
    if (is.null(segstats$mean))
        stop('segstats must have field $mean')
    
    segstats = segstats[!is.na(segstats$mean) & !is.na(segstats$sd)]

    if (!is.null(segstats$ncn))
        segstats = segstats[segstats$ncn==2, ]
    
    ## if (is.null(segstats$ncn))
    ##     ncn = rep(2, length(mu))
    ## else
    ##     ncn = segstats$ncn
    
    mu = segstats$mean
    w = as.numeric(width(segstats))
    Sw = sum(as.numeric(width(segstats)))
    sd = segstats$sd
    m0 = sum(as.numeric(mu*w))/Sw
    
    if (verbose)
        cat(paste(c(rep('.', length(purity.guesses)), '\n'), collapse = ''))
    
    NLL = matrix(unlist(mclapply(1:length(purity.guesses), function(i)
    {
        if (verbose)
            cat('.')
        nll = rep(NA, length(ploidy.guesses))
        for (j in 1:length(ploidy.guesses))
        {
            alpha = purity.guesses[i]
            tau = ploidy.guesses[j]
            gamma = 2/alpha - 2
            beta = (tau + gamma)/m0 ## replaced with below 9/10/14
                                        #          beta = ( tau + tau_normal * gamma /2 ) / m0
                                        #          v = pmax(0, round(beta*mu-ncn*gamma/2))
            v = pmax(0, round(beta*mu-gamma))
            
            ## using normal cn
                                        #          nll[j] = sum((v-beta*mu+ncn*gamma/2)^2/((sd)^2))

                                        # REVERT
            nll[j] = sum((v-beta*mu+gamma)^2/((sd)^2))

                                        #  mv = pmax(20, pmin(20, max(v, na.rm = TRUE)))
                                        #  mv = 2
                                        # log likelihood matrix across approximately "all" integer copy states
                                        # we are obviously ignoring very high states in this estimate
                                        # but they are likely to have high sd and thus contribute less to the overall log likelihood
            ##           ll = -sapply(0:mv, function(x) (x-beta*mu+gamma)^2/((sd)^2)) ## OG VERSION
                                        #    ll = -sapply(0:mv, function(x) ((x+gamma)/beta-mu)^2 / sd^2)  ## NEWFANGLED VERSION
            ##        ml = apply(ll, 1, max) ##  get maximum likelihood
            ##       probs  = 1/rowSums(exp(ll-ml)) ## normalize to get posterior probabilities (assuming uniform prior)
            ##      nll[j] = sum(-log(probs))
        }
        return(nll)
    }, mc.cores = mc.cores)), nrow = length(purity.guesses), byrow = T)
    
    dimnames(NLL) = list(as.character(purity.guesses), as.character(ploidy.guesses))

    cat('\n')

    ## rix = as.numeric(rownames(NLL))>=purity.min & as.numeric(rownames(NLL))<=purity.max
    ## cix = as.numeric(colnames(NLL))>=ploidy.min & as.numeric(colnames(NLL))<=ploidy.max  
    ## NLL = NLL[rix, cix, drop = FALSE]
    
    a = rep(NA, nrow(NLL));
    b = rep(NA, ncol(NLL)+2)
    b.inf = rep(Inf, ncol(NLL)+2)
    #'  a = rep(Inf, nrow(NLL));
    #'  b = rep(Inf, ncol(NLL)+2)
    NLLc = rbind(b, cbind(a, NLL, a), b) ## padded NLL and all of its shifts
    NLLul = rbind(cbind(NLL, a, a), b.inf, b)
    NLLuc = rbind(cbind(a, NLL, a), b.inf, b)
    NLLur = rbind(cbind(a, a, NLL), b.inf, b)
    NLLcl = rbind(b, cbind(NLL, a, a), b)
    NLLcr = rbind(b, cbind(a, a, NLL), b)
    NLLll = rbind(b, b, cbind(NLL, a, a))
    NLLlc = rbind(b, b, cbind(a, NLL, a))
    NLLlr = rbind(b, b, cbind(a, a, NLL))

    if (min(c(ncol(NLL), nrow(NLL)))>1) ## up up down down left right left right ba ba start
        M = (NLLc < NLLul & NLLc < NLLuc & NLLc < NLLur & NLLc < NLLcl & NLLc < NLLcr & NLLc < NLLll & NLLc < NLLlc & NLLc < NLLlr)[-c(1, nrow(NLLc)), -c(1, ncol(NLLc)), drop = FALSE]
    else if (ncol(NLL)==1) ## one column, only go up and down
        M = (NLLc < NLLuc & NLLc < NLLlc)[-c(1, nrow(NLLc)), -c(1, ncol(NLLc)), drop = FALSE]
    else ## only row, only go left right
        M = (NLLc < NLLcl & NLLc < NLLcr)[-c(1, nrow(NLLc)), -c(1, ncol(NLLc)), drop = FALSE]

    if (length(M)>1)
    {
        ix = which(M, arr.ind= T);
        if (nrow(ix)>1)
        {
            C = hclust(d = dist(ix), method = 'single')
            cl = cutree(C, h = min(c(nrow(NLL), ncol(NLL), 2)))
            minima = ix[vaggregate(1:nrow(ix), by = list(cl), function(x) x[which.min(NLL[ix[x, drop = FALSE]])]), , drop = FALSE]
        }
        else
            minima = ix[1,, drop = FALSE]
    }
    else
        minima = cbind(1,1)

    out = data.frame(purity = as.numeric(rownames(NLL)[minima[,1]]), ploidy = as.numeric(colnames(NLL)[minima[,2]]), NLL = NLL[minima],
                     i = minima[,1], j = minima[,2])

    out = out[order(out$NLL), , drop = FALSE]
    rownames(out) = 1:nrow(out)
    ## Saturday, Sep 02, 2017 10:33:26 PM
    ## Noted floating point error, use the epsilon trick to replace '>='
    ## out = out[out$purity>=purity.min & out$purity<=purity.max & out$ploidy>=ploidy.min & out$ploidy<=ploidy.max, ]
    eps = 1e9
    out = out[out$purity - purity.min >= -eps &
              out$purity - purity.max <= eps &
              out$ploidy - ploidy.min >= -eps &
              out$ploidy - ploidy.max <= eps, ]
    out$gamma = 2/out$purity -2
    out$beta = (out$ploidy + out$gamma)/m0
    out$mincn = mapply(function(gamma, beta) min(round(beta*mu-gamma)), out$gamma, out$beta)
    out$maxcn = mapply(function(gamma, beta) max(round(beta*mu-gamma)), out$gamma, out$beta)
    
    ## group solutions with (nearly the same) slope (i.e. 1/beta), these should have almost identical
    ## NLL (also take into account in-list distance just be safe)
    if (nrow(out)>1)
        out$group = cutree(hclust(d = dist(cbind(100/out$beta, 1:nrow(out)), method = 'manhattan'), method = 'single'), h = 2)
    else
        out$group = 1
    out = out[out$group<=3, ,drop = FALSE] ## only pick top 3 groups

    if (allelic) ## if allelic then use allelic distance to rank best solution in group
    {          
        ## remove all NA allelic samples
        segstats = segstats[!is.na(segstats$mean_high) & !is.na(segstats$sd_high) & !is.na(segstats$mean_low) & !is.na(segstats$sd_low)]
        out$NLL.allelic = NA          
        mu = cbind(segstats$mean_high, segstats$mean_low)
        w = matrix(rep(as.numeric(width(segstats)), 2), ncol = 2, byrow = TRUE)
        Sw = sum(as.numeric(width(segstats)))*2
        sd = cbind(segstats$sd_high, segstats$sd_low)
        m0 = sum(as.numeric(mu*w))/Sw
        
        if (verbose)
            cat(paste(c(rep('.', length(purity.guesses)), '\n'), collapse = ''))
        
        for (i in 1:nrow(out))
        {
            cat(sprintf('Evaluating alleles for solution %s of %s\n', i, nrow(out)))
            alpha = out$purity[i]
            tau = out$ploidy[i]
                                        #                  gamma = 2/alpha - 2
            gamma = 1/alpha - 1 ## 1 since we are looking at hets
            beta = (tau + gamma)/m0 ## replaced with below 9/10/14
                                        #          beta = ( tau + tau_normal * gamma /2 ) / m0
                                        #          v = pmax(0, round(beta*mu-ncn*gamma/2))
            v = pmax(0, round(beta*mu-gamma))
            
            vtot = round(out$beta[i]*segstats$mean-out$gamma[i])
            vlow.mle = rep(NA, length(vtot))

            for (j in 1:length(vlow.mle))
            {
                if (vtot[j]==0)
                    vlow.mle[j] = 0
                else
                {
                    vlow = 0:floor(vtot[j]/2)
                    vhigh = vtot[j]-vlow
                    tmp.nll = cbind((vlow-beta*mu[j,2]+gamma)^2/(sd[j,2])^2, (vhigh-beta*mu[j, 1]+gamma)^2/((sd[j,1])^2))
                    vlow.mle[j] = vlow[which.min(rowSums(tmp.nll))]
                }
            }
            
            vlow.mle = apply(cbind(mu, sd, vtot), 1, function(x) {
                tot = x[5]
                if (tot == 0)
                    return(0)
                else
                {
                    vlow = 0:floor(tot/2)
                    vhigh = tot-vlow
                    muh = x[1]
                    mul = x[2]
                    sdh = x[3]
                    sdl = x[4]
                    tmp.nll = cbind((vlow-beta*mul+gamma)^2/(sdl)^2, (vhigh-beta*muh+gamma)^2/((sdh)^2))
                    return(vlow[which.min(rowSums(tmp.nll))])
                }                      
            })                  
            
            out$NLL.allelic[i] = sum((cbind(vtot-vlow.mle, vlow.mle)-beta*mu+gamma)^2/sd^2)
        }
        
        out$NLL.tot = out$NLL
        out$NLL = out$NLL.tot + out$NLL.allelic
        out.all = out
        ix = vaggregate(1:nrow(out), by = list(out$group), FUN = function(x) x[order(abs(out$NLL[x]))][1])
    }
    else ## otherwise choose the one that gives the lowest magnitude copy number
    {
        out.all = out
        ix = vaggregate(1:nrow(out), by = list(out$group), FUN = function(x) x[order(abs(out$mincn[x]), out$mincn[x]<0)][1])
    }
    
    out = out[ix, , drop = FALSE]
    out$NLL = vaggregate(out$NLL, by = list(out$group), FUN = min)
    
    out.all$keep = 1:nrow(out.all) %in% ix ## keep track of other ploidy group peaks for drawing purposes
    out.all = out.all[out.all$group %in% out$group, ] ## only draw the groups in the top solution
    
    if (plot)
    {
        library(ellipse)
        library(numDeriv)

        ## xval = as.numeric(rownames(NLL))
        ## yval = as.numeric(colnames(NLL))
        ## f = function(x) {              
        ##     i = x[1]; ## interpolate between closest values of NLL
        ##     im = which(i<=xval)[1]
        ##     ip = (i-xval[im-1])/diff(xval[c(im-1, im)]);  ## proportion of lower value to consider
        ##     j = x[2];
        ##     jm = which(j<=yval)[1]
        ##     jp = (j-yval[jm-1])/diff(yval[c(jm-1, jm)]);  ## proportion of lower value to consider
        ##     nllm = NLL[c(im-1, im), c(jm-1, jm)] ## piece of NLL matrix containing the low and high i and j matches
        ##     nllp = cbind(c(ip, 1-ip)) %*% rbind(c(jp, 1-jp)) ## proportion of values to input into interpolation
        ##     return(sum(-nllm*nllp))
        ## }
        
        plot(out.all$purity, out.all$ploidy, pch = 19,
             xlim = c(purity.min, purity.max), ylim = c(ploidy.min, ploidy.max), xlab = 'purity', ylab = 'ploidy', cex = 0.2, col = alpha('white', out.all$intensity))
        
        
        f = function(x) -NLL[((x[1]-1) %% nrow(NLL))+1, ((x[2]-1) %% ncol(NLL))+1]
          ir = range(as.numeric(rownames(NLL)))
          jr = range(as.numeric(colnames(NLL)))
          txf = function(z) cbind(affine.map(z[,1], ir, c(1, nrow(NLL))), affine.map(z[,2], jr, c(1, ncol(NLL))))

          levs = c(0.95, 0.99)
                                        #          levs = c(1-1e-7) 
          tmp.out = out.all
          tmp.out$NLL = as.data.table(tmp.out)[, NLL := min(NLL), by = group][, NLL]
          tmp.out$intensity = affine.map(tmp.out$NLL, c(1, 0.1))
          tmp.out = tmp.out[rev(1:nrow(tmp.out)), ]
                                        #          tmp.out$col = brewer.master(length(levels(tmp.out$group)), 'PuRd')[as.integer(tmp.out$group)]
          tmp.out$col = brewer.pal(length(unique(out$group))+2, 'PuRd')[match(tmp.out$group, unique(tmp.out$group))]
          tmp.out$rank = ''; ## hacky stuff to just plot ranks for the top per group solutions
          tmp.out$rank[tmp.out$keep] = match(tmp.out$group[tmp.out$keep], out$group)

          require(DiceKriging)
          bla = mapply(function(x, y, c, a)
          {
              ## grab k square from computed NLL values to krig around
              k = 4
              ir = pmin(nrow(NLL), pmax(1, (x-k):(x+k)))
              jr = pmin(ncol(NLL), pmax(1, (y-k):(y+k)))
              ij = expand.grid(ir, jr)
              xy = expand.grid(purity.guesses[ir], ploidy.guesses[jr])
              m = tryCatch(km(design = xy[, 1:2], response = -NLL[as.matrix(ij)]), error = function(e) NULL)                     
              ## custom function gives best krigged interpolation in
              ## region will be used for hessian computaiton
              f = function(x) predict(m, newdata = data.frame(Var1 = x[1], Var2 = x[2]), type = 'UK')$mean                                          
              for (lev in levs) ## TOFIX: MESS
              {
                  if (!is.null(m))
                      F = tryCatch(-hessian(f, c(purity.guesses[x], ploidy.guesses[y])), error = function(e) matrix(NA, ncol = 2, nrow = 2))
                  else
                      F = NA
                  
                  if (all(!is.na(F)))
                      V = solve(F)
                  else
                      V = NA                          
                  if (any(is.na(V)))
                      V = diag(rep(1,2))
                  M = cov2cor(V)
                  XY = ellipse(M, scale = .01*c(diff(par('usr')[1:2]),diff(par('usr')[3:4])), centre = c(purity.guesses[x], ploidy.guesses[y]), level = lev)
                  polygon(XY[,1], XY[,2], border = NA, col = alpha(c, a*affine.map(lev, c(1, 0.3))));
              }                                                                  
          }, tmp.out$i, tmp.out$j, tmp.out$col, tmp.out$intensity, SIMPLIFY = FALSE)
          
          tmp.out = tmp.out[tmp.out$keep, ]
          text(tmp.out$purity, tmp.out$ploidy, tmp.out$rank, col = alpha('black', 0.5))
          
          tmp.out = tmp.out[rev(1:nrow(tmp.out)), ]
          legend(0, par('usr')[4]*0.98, legend = paste(tmp.out$rank, ') ', sapply(tmp.out$purity, function(x) sprintf('%0.2f',x)), ', ',
                                                       sapply(tmp.out$ploidy, function(x) sprintf('%0.2f',x)),
                                                       ' (gamma = ',sapply(tmp.out$gamma, function(x) sprintf('%0.3f',x)),
                                                       ', beta = ',sapply(tmp.out$beta, function(x) sprintf('%0.3f',x)),
                                                       ')', sep = ''), fill = tmp.out$col, cex = 0.8, yjust = 1, ncol = 1)
      }

    out = out.all;
    out = out[order(out$group, !out$keep, out$NLL), ]
    out$rank = NA
    out$rank[out$keep] = 1:sum(out$keep)
    out$keep = out$i = out$j = NULL
    rownames(out) = NULL
    return(out)
}

############################
#' rel2abs
#' 
#' rescales CN values from relative to "absolute" (i.e. per cancer cell copy) scale given purity and ploidy
#' 
#' takes in gr with signal field "field"
#'
#' @param gr GRanges input with meta data field corresponding to mean relative copy "mean" in that interval
#' @param purity purity of sample
#' @param ploidy ploidy of sample
#' @param gamma gamma fit of solution (over-rides purity and ploidy)
#' @param beta beta fit of solution (over-rides purity and ploidy)
#' @param field meta data field in "gr" variable from which to extract signal, default "mean"
#' @param field.ncn meta data field in "gr" variable from which to extract germline integer copy number, default "ncn", if doesn't exist, germline copy number is assumed to be zero
#' @return
#' numeric vector of integer copy numbers
#' @export
############################################
rel2abs = function(gr, purity = NA, ploidy = NA, gamma = NA, beta = NA, field = 'ratio', field.ncn = 'ncn')
{        
    mu = values(gr)[, field]
    mu[is.infinite(mu)] = NA
    w = as.numeric(width(gr))
    w[is.na(mu)] = NA
    sw = sum(w, na.rm = T)
    mutl = sum(mu * w, na.rm = T)

    ncn = rep(2, length(mu))
    if (!is.null(field.ncn))
        if (field.ncn %in% names(values(gr)))
            ncn = values(gr)[, field.ncn]

    ploidy_normal = sum(w * ncn, na.rm = T) / sw  ## this will be = 2 if ncn is trivially 2
    
    if (is.na(gamma))
        gamma = 2*(1-purity)/purity

    if (is.na(beta))
        beta = ((1-purity)*ploidy_normal + purity*ploidy) * sw / (purity * mutl)
                                        #      beta = (2*(1-purity)*sw + purity*ploidy*sw) / (purity * mutl)
    

                                        # return(beta * mu - gamma)
    return(beta * mu - ncn * gamma / 2)
}


############################
#' pp2gb
#' 
#' converts purity / ploidy to gamma / beta (or reverse)
#' 
#' takes in gr with signal field "field"
#'
#' @param purity value between 0 and 1
#' @param ploidy value nonnegative
#' @param mu vector of n segment averages
#' @param w vector of n segment widths 
#' @param gamma non-negative value
#' @param beta non-negative value
#' @return
#' list with purity / ploidy / gamma / beta entries
#' @export
############################
pp2gb = function(purity = NA, ploidy = NA, mu = NA, w = NA, gamma = NA, beta = NA)
{                
    if (all(is.na(mu)) & all(is.na(w)))
        stop('mu and w should be non-empty')
    
    if (length(mu) != length(w))
        stop('mu and w should match up in length')

    w = as.numeric(w)
    mu[is.infinite(mu)] = NA
    w[is.na(mu)] = NA
    sw = sum(w, na.rm = T)
    mutl = sum(mu * w, na.rm = T)                
    ncn = rep(2, length(mu))
    ploidy_normal = sum(w * ncn, na.rm = T) / sw  ## this will be = 2 if ncn is trivially 2                
    
    if (is.na(gamma) & is.na(beta) & !is.na(purity) & !is.na(ploidy))
    {
        gamma = 2*(1-purity)/purity
        beta = ((1-purity)*ploidy_normal + purity*ploidy) * sw / (purity * mutl)
    }
    else if (!is.na(gamma) & !is.na(beta) & is.na(purity) & is.na(ploidy))
    {
        purity = 2/(gamma+2)
        v = beta * mu - ncn * gamma / 2
        ploidy = sum(v*w, na.rm = TRUE)/sum(w, na.rm = TRUE)
    }
    else
        stop('Either gamma and beta are empty OR purity and ploidy are empty')
    return(list(purity = purity, ploidy = ploidy, gamma = gamma, beta = beta))
}


############################
#' abs2rel
#' 
#' rescales CN values from relative to "absolute" (i.e. per cancer cell copy) scale given purity and ploidy
#' By default, output is normalized to 1 (i.e. assumes that the total relative copy number signal mass over the genome is 1)
#' 
#' takes in gr with signal field "field"
#' @param gr GRanges input with meta data field corresponding to mean relative copy "mean" in that interval
#' @param purity purity of sample
#' @param ploidy ploidy of sample
#' @param gamma gamma fit of solution (over-rides purity and ploidy)
#' @param beta beta fit of solution (over-rides purity and ploidy)
#' @param field meta character specifying meta data field in "gr" variable from which to extract signal, default "mean"
#' @param field.ncn character specifying meta data field in "gr" variable from which to extract germline integer copy number, default "ncn", if doesn't exist, germline copy number is assumed to be zero
#' @return
#' numeric vector of integer copy numbers
#' @export
############################
abs2rel = function(gr, purity = NA, ploidy = NA, gamma = NA, beta = NA, field = 'cn', field.ncn = 'ncn', total = 1)
{        
    abs = values(gr)[, field]
    w = width(gr)
    sw = sum(as.numeric(w))

    ncn = rep(2, length(gr))
    if (!is.null(field.ncn))
        if (field.ncn %in% names(values(gr)))
            ncn = values(gr)[, field.ncn]
    
    if (is.na(gamma))
        gamma = 2*(1-purity)/purity

    ploidy_normal = sum(w * ncn, na.rm = T) / sw  ## this will be = 2 if ncn is trivially 2
    
    if (is.na(beta))
        beta = ((1-purity)*ploidy_normal + purity*ploidy) * sw / (purity * mutl)
                                        #  beta = (2*(1-purity)*sw + purity*ploidy*sw) / (purity * total)    
    
                                        #    return((abs + gamma) / beta)
    return((abs + ncn*gamma/2) / beta)
}


##############################
#' pp.nll
#' 
#' computes neg log likelihood ($nll) for purity, ploidy combo and mle abs copy numbers ($v), returns as list
#' 
#' v can be over-ridden to compute NLL for other (eg non MLE) values
#'
#' @param segstats granges of segstats with field $mean and $sd corresponding to mean and standard deviation on estimates of fragment density
#' @param purity numeric purity value
#' @param ploidy numeric ploidy value
#' @param gamma numeric gamma value (over-rides purity, ploidy)
#' @param beta numeric beta value (over-rides purity, ploidy)
#' @param field.ncn character specifying meta data field specifying germline copy number, default value is "ncn", if absent then will assume 2
#' @param solution over-ride if not interested in MLE / weighted least squares fit, this should be an integer vector of length(segstats)
#' @return
#' negative log likelihood of MLE i.e. least squares model (or supplied solution) fit
#' @export
##############################
pp.nll = function(segstats, purity = NA, ploidy = NA, gamma = NA, beta = NA, field = 'mean', field.ncn = 'ncn', v = NULL)
{
    mu = segstats$mean
    w = width(segstats)
    Sw = sum(as.numeric(width(segstats)))
    sd = segstats$sd
    m0 = sum(as.numeric(mu*w))/Sw
    alpha = purity
    tau = ploidy
    ncn = rep(2, length(mu))
    
    if (!is.null(field.ncn))
        if (field.ncn %in% names(values(segstats)))
            ncn = values(segstats)[, field.ncn]

    ploidy_normal = sum(w * ncn, na.rm = T) / Sw  ## this will be = 2 if ncn is trivially 2
    
    if (is.na(gamma))
        gamma = 2/alpha - 2
    
    if (is.na(beta))
        beta = ( tau + ploidy_normal * gamma / 2 ) / m0
                                        # beta = (tau + gamma)/m0

    if (is.null(v))
        v = round(beta*mu-gamma)

    return(list(NLL = sum((v-beta*mu+ncn*gamma/2)^2/((sd)^2)), v = v))
                                        #    return(list(NLL = sum((v-beta*mu+gamma)^2/sd^2), v = v))
}



############# Utilities for preprocessing eg segmentation etc
#############

#############################################################
#' munlist
#' 
#' unlists a list of vectors, matrices, data frames into a n x k matrix
#' whose first column specifies the list item index of the entry
#' and second column specifies the sublist item index of the entry
#' and the remaining columns specifies the value(s) of the vector
#' or matrices.
#' 
#' force.cbind = T will force concatenation via 'cbind'
#' force.rbind = T will force concatenation via 'rbind'
#############################################################
munlist = function(x, force.rbind = F, force.cbind = F, force.list = F)
{
    x = x[!sapply(x, is.null)]
    
    if (length(x)==0)
        return(NULL)
    
    if (!any(c(force.list, force.cbind, force.rbind)))
    {
        if (any(sapply(x, function(y) is.null(dim(y)))))
            force.list = T
        if (length(unique(sapply(x, function(y) dim(y)[2]))) == 1)
            force.rbind = T
        if ((length(unique(sapply(x, function(y) dim(y)[1]))) == 1))
            force.cbind = T
    }        
    else
        force.list = T    
    
    if (force.list)
        return(cbind(ix = unlist(lapply(1:length(x), function(y) rep(y, length(x[[y]])))),
                     iix = unlist(lapply(1:length(x), function(y) 1:length(x[[y]]))), 
                     unlist(x)))
    else if (force.rbind)
        return(cbind(ix = unlist(lapply(1:length(x), function(y) rep(y, nrow(x[[y]])))), 
                     iix = unlist(lapply(1:length(x), function(y) 1:nrow(x[[y]]))),
                     do.call('rbind', x)))
    else if (force.cbind)
        return(t(rbind(ix = unlist(lapply(1:length(x), function(y) rep(y, ncol(x[[y]])))),
                       iix = unlist(lapply(1:length(x), function(y) 1:ncol(x[[y]]))),
                       do.call('cbind', x))))
}


############################
#' over-ride Rcplex
#' 
#' 
############################
Rcplex2 = function (cvec, Amat, bvec, Qmat = NULL, lb = 0, ub = Inf, control = list(), 
    objsense = c("min", "max"), sense = "L", vtype = NULL, n = 1) 
{
    stopifnot((is(Amat, "matrix") && is.double(Amat)) || is(Amat, 
        "dsparseMatrix") || is(Amat, "simple_triplet_matrix"))
    numrows <- nrow(Amat)
    numcols <- ncol(Amat)
    if (!is.null(Qmat)) {
        stopifnot((is(Qmat, "matrix") && is.double(Qmat)) || 
            is(Qmat, "dsparseMatrix") || is(Qmat, "simple_triplet_matrix"), 
            nrow(Qmat) == numcols, ncol(Qmat) == numcols)
    }
    stopifnot(length(cvec) == numcols, is.double(cvec), length(bvec) == 
        numrows, is.double(bvec))
    if (length(lb) == 1L) 
        lb <- rep(lb, numcols)
    if (length(ub) == 1L) 
        ub <- rep(ub, numcols)
    stopifnot(length(lb) == numcols, is.double(lb), length(ub) == 
        numcols, is.double(ub))
    if (missing(objsense)) 
        objsense <- "min"
    stopifnot(objsense %in% c("min", "max"))
    if (objsense == "min") {
        objsensei <- 1L
    }
    else {
        objsensei <- -1L
    }
    if (length(sense) == 1L) 
        sense <- rep(sense, numrows)
    stopifnot(length(sense) == numrows, is.character(sense), 
        all(sense %in% c("L", "G", "E")))
    if (!is.null(vtype)) {
        if (length(vtype) == 1L) 
            vtype <- rep(vtype, numcols)
        stopifnot(length(vtype) == numcols, is.character(vtype), 
            all(vtype %in% c("C", "I", "B")))
        isMIP <- ifelse(any(vtype != "C"), 1L, 0L)
    }
    else isMIP <- 0L
    n <- as.integer(n)
    if (is.na(n)) 
        n <- 2100000000L
    Acpx <- Rcplex:::toCPXMatrix(Amat)
    Qcpx <- Rcplex:::toCPXMatrix(Qmat)
    isQP <- ifelse(is.null(Qcpx), 0L, 1L)
    control <- check.Rcplex.control(control, isQP)
    control <- Rcplex:::split.control.list(control)
    on.exit(.C("Rcplex_free"))
    res <- .Call("Rcplex", as.integer(numcols), as.integer(numrows), 
        as.integer(objsensei), as.double(cvec), as.double(bvec), 
        Acpx, Qcpx, as.double(lb), as.double(ub), as.character(sense), 
        as.character(vtype), as.integer(isQP), as.integer(isMIP), 
        as.integer(n), control$C, as.integer(control$R$maxcalls))
    if (isMIP) {
        intvars <- which(vtype != "C")
        .canonicalize <- function(x) {
            names(x) <- c("xopt", "obj", "status", "extra")
            names(x$extra) <- c("nodecnt", "slack")
            if (control$R$round) 
                x$xopt[intvars] <- round(x$xopt[intvars])
            x
        }
        res <- if (n > 1L) 
            lapply(res, .canonicalize)
        else .canonicalize(res)
    }
    else {
        names(res) <- c("xopt", "obj", "status", "extra")
        names(res$extra) <- c("lambda", "slack")
    }
    return(res)
}


### over-ride of Cplex function
check.Rcplex.control = function (control, isQP) 
{
    con <- list(trace = 1L, method = 0L, preind = 1L, aggind = -1L, 
        itlim = 2100000000L, epagap = 0, epgap = 1e-04, tilim = 1e+74, 
        disjcuts = 0L, mipemphasis = 0L, cliques = 0L, nodesel = 1L, 
        probe = 0L, varsel = 0L, flowcovers = 0L, solnpoolagap = 0, 
        solnpoolgap = 0, solnpoolintensity = 0L, maxcalls = 500L, 
        round = 0L)
    con[names(control)] <- control
    if (!is.null(con$trace)) {
        con$trace <- as.integer(con$trace)
        if (!con$trace %in% c(0L, 1L)) {
            warning("Improper value for trace parameter. Using default.")
            con$trace <- 1L
        }
    }
    if (!is.null(con$method)) {
        con$method <- as.integer(con$method)
        if (!con$method %in% 0L:4L) {
            warning("Improper value for method parameter. Using default.")
            con$method <- 0L
        }
#        if (isQP && (!con$method %in% c(0L, 4L))) {
#            warning("Improper value for method parameter. Using default.")
#            con$method <- 0L
#        }
    }
    if (!is.null(con$preind)) {
        con$preind <- as.integer(con$preind)
        if (!con$preind %in% c(0L, 1L)) {
            warning("Improper value for preind parameter. Using default.")
            con$preind <- 1L
        }
    }
    if (!is.null(con$aggind)) {
        con$aggind <- as.integer(con$aggind)
        if (con$aggind < 0L && con$aggind != -1L) {
            warning("Improper value for aggind parameter. Using default.")
            con$aggind <- -1L
        }
    }
    if (!is.null(con$itlim)) {
        con$itlim <- as.integer(con$itlim)
        if ((con$itlim < 0L) || (con$itlim > 2100000000L)) {
            warning("Improper value for itlim parameter. Using default.")
            con$itlim <- 2100000000L
        }
    }
    if (!is.null(con$epagap)) {
        con$epagap <- as.numeric(con$epagap)
        if (con$epagap < 0) {
            warning("Improper value for epagap parameter. Using default.")
            con$epagap <- 0
        }
    }
    if (!is.null(con$epgap)) {
        con$epgap <- as.numeric(con$epgap)
        if (con$epgap < 0 || con$epgap > 1) {
            warning("Improper value for epgap parameter. Using default.")
            con$epgap <- 1e-04
        }
    }
    if (!is.null(con$tilim)) {
        con$tilim <- as.numeric(con$tilim)
        if (con$tilim < 0) {
            warning("Improper value for tilim parameter. Using default")
            con$tilim <- 1e+75
        }
    }
    if (!is.null(con$disjcuts)) {
        con$disjcuts <- as.integer(con$disjcuts)
        if (!con$disjcuts %in% -1L:3L) {
            warning("Improper value for disjcuts parameter: Using default.")
            con$disjcuts <- 0L
        }
    }
    if (!is.null(con$mipemphasis)) {
        con$mipemphasis <- as.integer(con$mipemphasis)
        if (!con$mipemphasis %in% 0L:4L) {
            warning("Improper value for mipemphasis parameter: Using default.")
            con$mipemphasis <- 0L
        }
    }
    if (!is.null(con$cliques)) {
        con$cliques <- as.integer(con$cliques)
        if (!con$cliques %in% -1L:2L) {
            warning("Improper value for cliques parameter: Using default.")
            con$cliques <- 0L
        }
    }
    if (!is.null(con$nodesel)) {
        con$nodesel <- as.integer(con$nodesel)
        if (!con$nodesel %in% 0L:3L) {
            warning("Improper value for nodesel parameter: Using default.")
            con$nodesel <- 1L
        }
    }
    if (!is.null(con$probe)) {
        if (!con$probe %in% -1L:3L) {
            warning("Improper value for probe parameter: Using default.")
            con$probe <- 0L
        }
    }
    if (!is.null(con$varsel)) {
        con$varsel <- as.integer(con$varsel)
        if (!con$varsel %in% -1L:4L) {
            warning("Improper value for varsel parameter: Using default.")
            con$varsel <- 0L
        }
    }
    if (!is.null(con$flowcovers)) {
        con$flowcovers <- as.integer(con$flowcovers)
        if (!con$flowcovers %in% -1L:2L) {
            warning("Improper value for flowcovers parameter: Using default")
            con$flowcovers <- 0L
        }
    }
    if (!is.null(con$solnpoolagap)) {
        con$solnpoolagap <- as.numeric(con$solnpoolagap)
        if (!con$solnpoolagap >= 0) {
            warning("Improper value for solnpoolagap parameter: Using default")
            con$solnpoolagap <- 0
        }
    }
    if (!is.null(con$solnpoolgap)) {
        con$solnpoolgap <- as.numeric(con$solnpoolgap)
        if (!con$solnpoolgap >= 0) {
            warning("Improper value for solnpoolagap parameter: Using default")
            con$solnpoolgap <- 0
        }
    }
    if (!is.null(con$solnpoolintensity)) {
        con$solnpoolintensity <- as.integer(con$solnpoolintensity)
        if (!con$solnpoolintensity %in% 0L:4L) {
            warning("Improper value for solnpoolintensity parameter: Using default.")
            con$solnpoolagap <- 0L
        }
    }
    if (!is.null(con$maxcalls)) {
        if (con$maxcalls <= 0L) {
            warning("Improper value for maxcalls parameter. Using default.")
            con$maxcalls <- 500L
        }
    }
    if (!is.null(con$round)) {
        if (!con$round %in% c(0, 1)) {
            warning("Improper value for round option: Using default")
            con$round <- 0
        }
    }
    return(con)
}


##################################
#' segment
#' 
#' Wrapper around cumSeg to segment numeric data in an input GRanges with signal meta data field (e.g. $signal)
#' Returns a GRAnges of piecewise constant regions with their associated value
#' 
##################################
cumseg = function(gr, field = 'signal', log = T, type = 'bic', alg = 'stepwise', S = 1, verbose = F, mc.cores = 1, ...)
  {
    require(cumSeg)

    if (!(field %in% names(values(gr))))
      stop(sprintf('Field "%s" not a meta data columnin the input GRanges', field))

    if (log)
      values(gr)[, field] = log(values(gr)[, field])

    good.ix = !(is.infinite(values(gr)[, field]) | values(gr)[, field]==0)
    good.ix[is.na(values(gr)[, field])] = F
    gr = gr[good.ix]

    if (length(gr) == 0)
      return(gr[c(), field])
    
    args = list(...)

    if (!("sel.control" %in% names(args)))
      args$sel.control = sel.control(type = 'bic', alg = alg, S = S)
            
    grl = split(gr, as.character(seqnames(gr)))
    
    out = do.call('c', mclapply(names(grl), function(chr)
      {
        if (verbose)
          cat('Starting ', chr, '\n')
        this.gr = grl[[chr]]

        if (length(this.gr)<=1)
          return(si2gr(this.gr)[chr])
        
        tmp = do.call(jumpoints, c(list(values(this.gr)[, field], start(this.gr)), args))
        out = GRanges(chr, IRanges(c(1, tmp$psi), c(tmp$psi, seqlengths(this.gr)[[chr]])), seqlengths = seqlengths(this.gr))
        values(out)[, field] = NA
        values(out)[, field] = tmp$est.means
        
        if (verbose)
          cat('\t... generated ', length(out), ' segs\n')
        
        return(out)
      }, mc.cores = mc.cores))

    if (log)
      values(out)[, field] = exp(values(out)[, field])

    return(out)        
  }


######################################
#' multicoco
#' 
#' multi-scale coverage correction
#' 
#' Given gc and mappability coverage correction at k "nested" scales finds the coverage
#' assignment at the finest scale that yields the best correction at every scale
#' 
#' FUN = is a function that takes in a data frame / granges with
#' $reads and other covariate functions $gc, $map and outputs a vector of corrected read counts
#' 
#' cov is a constant with GRanges of coverage samples with (by default) fields $reads, $map, $gc 
#' 
#' base = is the multiple with which to perform "numlevs" additional scales of correction
#' 
#####################################
multicoco = function(cov,
    numlevs = 1, ## numbers of scales at which to correct
    base = 100, ## scale multipler
    fields = c("gc", "map"), # fields of gc which to use as covariates
    iterative = TRUE,
    presegment = TRUE, ## whether to presegment
    min.segwidth = 5e6, ## when presegmenting, min seg width
    mono = FALSE, ## just do single iteration at finest scale
    verbose = T,
    imageroot = NULL, ## optional file root to which to dump images of correction
    FUN = NULL, ## function with which to correct coverage (by default loess correction modified from HMMcopy that takes in granges with fields $reads and other fields specified in "fields"
    ..., ## additional args to FUN
    mc.cores = 1)
    {
        require(Rcplex)
        if (verbose)
           cat('Converting to data.table\n')

        WID = max(width(cov))
        
        cov.dt = as.data.table(as.data.frame(cov))
        
        sl = structure(as.numeric(1:length(seqlevels(cov))), names = seqlevels(cov))       

        if (verbose)
            cat('Grouping intervals\n')

        ## compute level means
        ## lev 0 = raw data
        ## lev 1 = base-fold collapsed
        ## lev 2 = base^2-fold collapsed
        ## etc
        parentmap= list() ## data.tables that map lev labels at level k  to parent lev labels
        cov.dt[, lev0:=as.character(1:length(seqnames))]
        for (k in 1:numlevs)
            {
                if (verbose)
                    cat('Defining', base^k, 'fold collapsed ranges\n')
                cov.dt[, eval(paste("lev", k, sep = '')) := as.character(sl[seqnames] + as.numeric(Rle(as.character(1:length(start)), rep(base^k, length(start)))[1:length(start)])/length(start)), by = seqnames]
                parentmap[[k]] = data.table(parent = cov.dt[, get(paste("lev", k, sep = ''))], child = cov.dt[, get(paste("lev", k-1, sep = ''))], key = 'child')[!duplicated(child), ]
            }

        if (presegment) ## perform rough segmentation at highest level
            {
                sl = seqlengths(cov)
                if (verbose)
                    cat('Presegmenting at ', as.integer(WID*base^(numlevs)), ' bp scale \n')
                require(DNAcopy)
                tmp.cov = seg2gr(cov.dt[,list(chr = seqnames[1], start = min(start), end = max(end), strand = strand[1], reads = mean(reads, na.rm = T)), by = get(paste("lev", numlevs, sep = ''))], seqlengths = sl)
                ix = which(!is.na(values(tmp.cov)[, 'reads']))
                cna = CNA(log(values(tmp.cov)[, 'reads'])[ix], as.character(seqnames(tmp.cov))[ix], start(tmp.cov)[ix], data.type = 'logratio')
                tmp = print(segment(smooth.CNA(cna), alpha = 1e-5, verbose = T))
                tmp = tmp[!is.na(tmp$loc.start) & !is.na(tmp$chrom) & !is.na(tmp$loc.end), ]
                seg = sort(seg2gr(tmp, seqlengths = sl))
                seg = seg[width(seg)>min.segwidth] ## remove small segments
                seg.dt = as.data.table(as.data.frame(seg))
                seg = seg2gr(seg.dt[, list(seqnames = seqnames,
                    start = ifelse(c(FALSE, seqnames[-length(seqnames)]==seqnames[-1]), c(1, start[-1]), 1),
                    end = ifelse(c(seqnames[-length(seqnames)]==seqnames[-1], FALSE), c(start[-1]-1, Inf), seqlengths(seg)[seqnames]))], seqlengths = sl)
                seg = gr.val(seg, tmp.cov, 'reads') ## populate with mean coverage
                seg$reads = seg$reads/sum(as.numeric(seg$reads*width(seg))/sum(as.numeric(width(seg)))) ## normalize segs by weigthed mean (so these become a correction factor)
            }
        else
            seg = NULL
        
        if (verbose)
            cat('Aggregating coverage within levels \n')
        
        ## list of data frames showing scales of increasing collapse

        cov.dt[, ix := 1:nrow(cov.dt)]
        
        cmd1 = paste('list(ix.start = ix[1], ix.end = ix[length(ix)], reads = mean(reads, na.rm = T),', paste(sapply(fields, function(f) sprintf("%s = mean(%s, na.rm = T)", f, f)), collapse = ','), ')', sep = '')
        
        cmd2 = paste('list(lab = lev0, reads,', paste(fields, collapse = ','), ', seqnames, start, end)', sep = '')

        if (mono)
            {
                if (verbose)
                    cat('Mono scale correction \n')
                 grs = list(cov.dt[, eval(parse(text=cmd2))])
                 numlevs = 1
             }
        else
            {
                grs = c( list(cov.dt[, eval(parse(text=cmd2))]),
                    lapply(1:numlevs, function(x)
                        {
                            if (verbose)
                                cat('Aggregating coverage in level', x,  '\n')
                            out = cov.dt[, eval(parse(text=cmd1)), keyby = list(lab = get(paste('lev', x, sep = '')))]
                            out[, ":="(seqnames = cov.dt$seqnames[ix.start], end = cov.dt$end[ix.start], start = cov.dt$start[ix.start])]
                            out[, ":="(ix.start= NULL, ix.end = NULL)]
                            return(out)
                        }))
            }
        
        setkey(grs[[1]], 'lab')
                               
        ## modified from HMMCopy to 
        ## (1) take arbitrary set of covariates, specified by fields vector
        ## (2) employ as input an optional preliminary (coarse) segmentation with which to adjust signal immediately prior to loess
        ## NOTE: (this only impacts the loess fitting, does not impose any segmentation on thed ata)
        ##
        if (is.null(FUN))            
            FUN = function(x, fields = fields, samplesize = 5e4, seg = NULL, ## seg is a Granges with meta data field $reads
                verbose = T, doutlier = 0.001, routlier = 0.01)
                {                    
                    if (!all(fields %in% names(x)))
                        stop(paste('Missing columns:', paste(fields[!(fields %in% names(x))], collapse = ',')))
                    
                    x$valid <- TRUE
                    for (f in fields)
                        {
                            x$valid[is.na(x[, f])] = FALSE
                            x$valid[which(is.infinite(x[, f]))] = FALSE
                        }
                   
                    if (verbose)
                        cat('Quantile filtering response and covariates\n')

                    range <- quantile(x$reads[x$valid], prob = c(routlier, 1 - routlier), na.rm = TRUE)

                    if (verbose)
                        cat(sprintf("Response min quantile: %s max quantile: %s\n", round(range[1],2), round(range[2],2)))
                    
                    domains = lapply(fields, function(f) quantile(x[x$valid, f], prob = c(doutlier, 1 - doutlier), na.rm = TRUE))
                    names(domains) = fields
                    
                    x$ideal <- x$valid
                    x$ideal[x$reads<=range[1] | x$reads>range[2]] = FALSE
                                        
                    for (f in fields)
                        x$ideal[x[, f] < domains[[f]][1] | x[, f] > domains[[f]][2]] = FALSE

                    if (verbose)
                        cat(sprintf('Nominated %s of %s data points for loess fitting\n', sum(x$ideal), nrow(x)))

                    set <- which(x$ideal)

                    if (length(set)<=10)
                        {
                            warning("Not enough samples for loess fitting - check to see if missing or truncated data?")
                            return(x$reads)
                        }
                    
                    for (f in fields)
                        {                            
                            if (verbose)
                                message(sprintf("Correcting for %s bias...", f))
                            
                            select <- sample(set, min(length(set), samplesize))

                            x2 = x[, c('reads', f)]
                            x2$covariate = x2[, f]

                            x2s = x2[select, ]
                            
                            if (!is.null(seg)) ## here we apply a prelmiinary segmentation to correct for large scale copy variation allow more power to reveal the covariate signal
                                {
                                    if (verbose)
                                        message('Applying preliminary segmentation prior to loess fitting')
                                    
                                    x.grs = gr.val(seg2gr(x[select, ], seqlengths = NULL), seg, 'reads')
                                    x2s$reads = x2s$reads/x.grs$reads
                                }
                            
                            fit = loess(reads ~ covariate, data = x2s, span = 0.3)
                            if (is.na(fit$s))
                                {
                                    warning("Using all points since initial loess failed")
                                    fit = loess(reads ~ covariate, data = x2[select, ], span = 1)
                                }

                            tryCatch(
                                {
                                    if (!is.na(fit$s))
                                        {
                                            domain = domains[[f]]

                                            yrange <- quantile(x2s$reads, prob = c(routlier, 1 - routlier), na.rm = TRUE)                                            
                                            df = data.frame(covariate = seq(domain[1], domain[2], 0.001))
       
                                            if (!is.null(imageroot))
                                                {
                                                    out.pdf = paste(imageroot, ifelse(grepl("/$", imageroot), '', '.'), f,'_correction.pdf', sep = '')
                                                    if (verbose) {
                                                        cat("Dumping figure to", out.pdf, "\n")
                                                    }

                                                    pdf(out.pdf, height = 6, width = 6) 
                                                    plot(x2s$covariate, x2s$reads, col = alpha('black', 0.1), pch = 19, cex = 0.4, xlim = domain, ylim = yrange, ylab = sprintf('signal before %s correction', f), xlab = f);
                                                    lines(df$covariate, predict(fit, df), col = 'red', lwd = 2)
                                                    dev.off()
                                                }                               
                                            x$reads = x2$reads/predict(fit, x2) ## apply correction
                                        }
                                    else
                                        print("Loess failed, yielding NA loess object, continuing without transforming data")
                                }, error = function(e) print("Unspecified loess or figure output error"))
                        }
                    return(x$reads)
                }

        if (verbose)
            cat('Correcting coverage at individual scales\n')
        
        ## level 1,2, ..., k corrections
        ## these are the computed corrected values that will be input into the objective function

        if (iterative) ## iteratively correct
            {
                correction = NULL
                for (i in rev(1:length(grs)))
                    {
                        cat('Correcting coverage at ', WID*base^(i-1), 'bp scale, with', nrow(grs[[i]]), 'intervals\n')
                        if (i != length(grs))                            
                            grs[[i]]$reads = grs[[i]]$reads/correction[parentmap[[i]][grs[[i]]$lab, parent], cor]

                        if (WID*base^(i-1) > 1e5) ## for very large intervals do not quantile trim, only remove NA
                            grs[[i]]$reads.corrected = FUN(as.data.frame(grs[[i]]), fields, doutlier = 0, seg = seg)
                        else
                            grs[[i]]$reads.corrected = FUN(as.data.frame(grs[[i]]), fields, seg = seg);  


                        if (is.null(correction))
                            correction = data.table(lab = grs[[i]]$lab, cor = grs[[i]]$reads / grs[[i]]$reads.corrected, key = 'lab')
                        else
                            {
                                ## multiply new correction and old correction
                                old.cor = correction[parentmap[[i]][grs[[i]]$lab, parent], cor]
                                new.cor = grs[[i]]$reads / grs[[i]]$reads.corrected                                     
                                correction = data.table(lab = grs[[i]]$lab,  cor = old.cor * new.cor, key = 'lab') ## relabel with new nodes
                            }                        
                    }

                cov.dt$reads.corrected = grs[[1]][cov.dt$lev0, ]$reads.corrected
                rm(grs)
                gc()
#                cov.dt[, reads.corrected := grs[[1]][lev0, reads.corrected]]
            }
        else ## parallel mode
            {

                ## compute marginal values / ratios for reads and covariates
                if (verbose)
                    cat('Computing marginal values of read coverage and covariates\n')
                
                ## now compute marginal ratio relative next-level mean for all levels except for top level
                for (k in 1:numlevs)
                    {
                        if (verbose)
                            cat('Defining marginal coverage for', base^k, 'fold collapsed ranges\n')
                        ix = parentmap[[k]][grs[[k]]$lab, parent]
                        grs[[k]]$reads = grs[[k]]$reads / grs[[k+1]][ix, reads]
                        print('fucken bitch yeah yeah')

                        for (f in fields)                    
                            grs[[k]][, eval(parse(text = sprintf("%s := grs[[k]]$%s / grs[[k+1]][ix, %s]",f, f, f)))]
                    }
                                                
                grs = mclapply(grs, function(x) {
                    if (verbose)
                        cat('Correcting coverage at ', base^(k-1), 'fold collapse, with', nrow(grs[[k]]), 'intervals\n')
                    x$cor = FUN(as.data.frame(x), fields);
                    return(x)
                }, mc.cores = mc.cores)
                
                gc()
                
                for (k in 1:length(grs))
                    {
                        cov.dt[, eval(paste('cor', k-1, sep = '')) := as.numeric(NA)]
                        cov.dt[, eval(paste('cor', k-1, sep = ''))] = grs[[k]][cov.dt[, get(paste('lev', k-1, sep = ''))], ]$cor
                    }
                
                ulev = unique(cov.dt[, get(paste('lev', numlevs, sep = ''))])
                
                cov.dt$lid = factor(cov.dt[, get(paste('lev', numlevs, sep = ''))])
                cov.dt[, gid := 1:length(start)]                
                cov.dt[, eval(parse(text = paste("reads.corrected := ", paste('cor', 0:numlevs, sep = '', collapse = "*"))))]
                
                ##        out = rbindlist(mclapply(ulev, function(x) .optimize_subtree(cov.dt[get(paste('lev', numlevs,sep = '')) == x], numlevs), mc.cores = mc.cores))                
                ## .optimize_subtree = function(this.chunk, numlevs)   {
                ##     browser()
                ##     this.chunk.og = this.chunk
                ##     this.chunk.og[, reads.corrected := NA] 
                ##     this.chunk.og$id = as.character(1:nrow(this.chunk.og))
                ##     setkey(this.chunk.og, 'id')
                ##     this.chunk = this.chunk.og[!is.na(cor0), ]

                ##     if (verbose)                
                ##         cat('chunk', as.integer(this.chunk.og$lid)[1], 'of', length(levels(this.chunk.og$lid)), '\n')
                    
                    
                ##     if (nrow(this.chunk)>0)
                ##         {
                ##             y_lev = rbindlist(
                ##                 c(lapply(0:(numlevs-1), function(x) this.chunk[, list(lev = x, cor = get(paste('cor', x, sep = ''))[1], parent_lab = get(paste('lev', x+1, sep = ''))[1]), by = list(lab = get(paste('lev', x, sep = '')))]), list(this.chunk[, list(lev = numlevs, cor = get(paste('cor', numlevs, sep = ''))[1], parent_lab = NA), by = list(lab = get(paste('lev', numlevs, sep = '')))])))

                ##                         #                    y_lev = y_lev[!is.na(cor), ]
                            
                ##             ## build tree structure of rows of y_lev i.e. the variables in our optimization
                ##             ## map nodes to their parents 
                ##             y_lev[, parent.ix := match(paste(lev+1, parent_lab), paste(lev, lab))]
                ##             y_lev[, id := 1:length(parent.ix)]

                ##             if (any(!is.na(y_lev$parent.ix)))
                ##                 {
                ##                     ## make constraints matrix - one constraint per unique parent
                ##                     ##
                ##                     ## this defines parents in terms of their children (mean function)
                ##                     A = y_lev[!is.na(parent.ix), sparseMatrix(as.integer(factor(c(parent.ix, parent.ix))), c(parent.ix, id), x = c(rep(-1, length(parent.ix)), rep(1, length(id))))]
                ##                     A = cBind(A, A*0) ## add room for residual variables
                ##                     Arhs = rep(0, nrow(A)) ## right hand side of A is 0

                ##                     ## residual constraints relate nodes, their parents and the fits contained in the "cor" columns
                ##                     ## if q_ki is the fit for the ith node in the kth level
                ##                     ## then x_ki - q_ki*x_p(ki) + r_ki = 0

                ##                     ## this encodes basic residual matrix whose rows have the form:  x_ki - q_ki*x_p(ki) = residual
                ##                     ## except for the top level, which is a single parentless node, and has the form:  x_ki - q_ki = residual
                ##                     ##
                ##                     ## the variables are indexed 1:nrow(y_lev), and the residuals are the next nrow(y_lev) indices
                ##                     R = sparseMatrix(rep(1:nrow(y_lev),2), c(1:nrow(y_lev), nrow(y_lev) + 1:nrow(y_lev)), x = rep(c(1,-1), each = nrow(y_lev)))
                ##                     rownames(R) = as.character(y_lev$id)
                ##                     R[y_lev[!is.na(parent.ix), cbind(id, parent.ix)]] = -y_lev[!is.na(parent.ix), cor]
                ##                     Rrhs = y_lev[, ifelse(is.na(parent.ix), cor, 0)]
                ##                     R = R[!is.na(y_lev$cor), ]
                ##                     Rrhs = Rrhs[!is.na(y_lev$cor)]
                                    
                ##                     ## make objective function
                ##                     ##
                ##                     ## problem:
                ##                     ## find x, r minimizing ||r||
                ##                     ##
                ##                     ## Ax = 0 (encodes node-parent relationships where parent = mean(children)
                ##                     ## x_ki - q_ki*x_p(ki) - r_ki = 0 (for k<numlevs)
                ##                     ## x_ki - r_ki = q_ki (for k = numlevs, a single node)
                ##                     ## 
                ##                     y_lev[, obj.weight := 1/length(id), by = lev]

                ##                     Q = diag(c(y_lev$obj.weight*0, 1 + 0*y_lev$obj.weight)) ## only put 1 weights on the matrix entries corresponding to the residual

                ##                     vtype = rep('C', ncol(A)) ## all variables are continuous
                ##                     sense = rep('E', nrow(A) + nrow(R)) ## all constraints specify equality
                ##                     tilim = 10;
                ##                     cvec = rep(0, ncol(A)) ## all variables are continuous

                ##                     ## now we need to encode the sum relationships between x0 etc.
                ##                     sol = Rcplex(cvec = cvec, Amat = rBind(A, R), bvec = c(Arhs, Rrhs), sense = sense, Qmat = Q, lb = rep(c(0, -Inf), each = nrow(y_lev)), ub = Inf, n = 1, objsense = "min", vtype = vtype, control = list(tilim = tilim))
                ##                     y_lev$opt = sol$xopt[1:nrow(y_lev)]
                ##                     setkey(this.chunk, 'lev0')
                ##                     this.chunk[y_lev[lev==0, lab], reads.corrected := y_lev[lev == 0, ]$opt]
                ##                     this.chunk[, reads.corrected := y_lev[lev == 0, ]$opt]
                ##                     this.chunk.og$reads.corrected = this.chunk[this.chunk.og$id, ]$reads.corrected
                ##                 }
                ##         }            
                ##     return(this.chunk.og)
                ## }
            }
        
        if (verbose)
            cat('Converting to GRanges\n')

        gc()
        
        out = seg2gr(as.data.frame(cov.dt), seqlengths = seqlengths(cov))

        if (verbose)
            cat('Made GRanges\n')
        
        gc()
        return(out)
    }

#####
#####

############################################
#' ra_breaks
#'
#' takes in either file or data frame from dranger or snowman or path to BND / SV type vcf file
#' and returns junctions in VCF format.
#' 
#' The default output is GRangesList each with a length two GRanges whose strands point AWAY from the break.  If get.loose = TRUE (only relevant for VCF)
#'
#' @name ra_breaks
#' @import VariantAnnotation
#' @export
############################################
ra_breaks = function(rafile, keep.features = T, seqlengths = hg_seqlengths(), chr.convert = T, snowman = FALSE, swap.header = NULL,  breakpointer = FALSE, seqlevels = NULL, force.bnd = FALSE, skip = NA, 
    get.loose = FALSE ## if TRUE will return a list with fields $junctions and $loose.ends
  )
  {
      if (is.character(rafile))
          {
              if (grepl('(.bedpe$)', rafile))                  
              {
                      ra.path = rafile                      
                      cols = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'name', 'score', 'str1', 'str2')

                      ln = readLines(ra.path)
                      if (is.na(skip))
                      {
                              nh = min(c(Inf, which(!grepl('^((#)|(chrom))', ln))))-1
                              if (is.infinite(nh))
                                  nh = 1
                          }
                      else
                          nh = skip

                                       
                      if ((length(ln)-nh)==0)
                          if (get.loose)
                              return(list(junctions = GRangesList(GRanges(seqlengths = seqlengths))[c()], loose.ends = GRanges(seqlengths = seqlengths)))
                          else                              
                              return(GRangesList(GRanges(seqlengths = seqlengths))[c()])
#                          return(GRangesList())
                      
                          
                      if (nh ==0)
                          rafile = fread(rafile, header = FALSE)
                      else
                          {
                              
                              rafile = tryCatch(fread(ra.path, header = FALSE, skip = nh), error = function(e) NULL)
                              if (is.null(rafile))
                                  rafile = tryCatch(fread(ra.path, header = FALSE, skip = nh, sep = '\t'), error = function(e) NULL)

                              if (is.null(rafile))
                                  rafile = tryCatch(fread(ra.path, header = FALSE, skip = nh, sep = ','), error = function(e) NULL)

                              if (is.null(rafile))
                                  stop('Error reading bedpe')                                  
                          }
                      setnames(rafile, 1:length(cols), cols)
                      rafile[, str1 := ifelse(str1 %in% c('+', '-'), str1, '*')]
                      rafile[, str2 := ifelse(str2 %in% c('+', '-'), str2, '*')]
#                      rafile[, str1 := ifelse(str1=='+', '-', '+')]
                                        #                      rafile[, str2 := ifelse(str2=='+', '-', '+')]
                      
                  }
              else if (grepl('(vcf$)|(vcf.gz$)', rafile))
                  {
                    library(VariantAnnotation)

                    if (grepl('gz$', rafile))
                      p = pipe(paste('zcat', rafile, ' | grep -P "\\#\\#contig"'))
                    else
                      p = pipe(paste('grep -P "\\#\\#contig"', rafile))
                    header = readLines(p)
                    vcf_sl= fread(paste(gsub("[<>]", "", gsub("\\,length\\=", " ", gsub('.*contig=<ID=', '', header))), collapse = '\n'))[ , V2 := max(V2), by = V1][, structure(V2, names = V1)]
                    vcf = readVcf(rafile, Seqinfo(seqnames = names(vcf_sl), seqlengths = vcf_sl))
                    
                    #if (!is.null(seqlengths))
                    #  vcf = gr.fix(vcf, seqlengths, drop = TRUE)
                    
                    if (!('SVTYPE' %in% names(info(vcf)))) {
                      warning('Vcf not in proper format.  Is this a rearrangement vcf?')
                      return(GRangesList());
                    }
                    
                      ## vgr = rowData(vcf) ## parse BND format                      
                      vgr = read_vcf(rafile, swap.header = swap.header)
                      
                      ## no events
                      if (length(vgr) == 0)
                        return (GRangesList())

                      ## fix mateids if not included
                      if (!"MATEID"%in%colnames(mcols(vgr))) {
                        nm <- vgr$MATEID <- names(vgr)
                        ix <- grepl("1$",nm)
                        vgr$MATEID[ix] = gsub("(.*?)(1)$", "\\12", nm[ix])
                        vgr$MATEID[!ix] = gsub("(.*?)(2)$", "\\11", nm[!ix])
                        vgr$SVTYPE="BND"
                      }
                      
                      if (!any(c("MATEID", "SVTYPE") %in% colnames(mcols(vgr))))
                        stop("MATEID or SVTYPE not included. Required")
                      
                      vgr$mateid = vgr$MATEID
                      if (is.null(vgr$SVTYPE))
                          vgr$svtype = vgr$SVTYPE
                      else
                          vgr$svtype = vgr$SVTYPE

                      if (!is.null(info(vcf)$SCTG))
                          vgr$SCTG = info(vcf)$SCTG

                      if (force.bnd)
                          vgr$svtype = "BND"
                      
                      if (sum(vgr$svtype == 'BND')==0)
                          warning('Vcf not in proper format.  Will treat rearrangements as if in BND format')

                      if (!all(vgr$svtype == 'BND'))
                          warning(sprintf('%s rows of vcf do not have svtype BND, ignoring these', sum(vgr$svtype != 'BND')))

                      bix = which(vgr$svtype == "BND")
                      vgr = vgr[bix]
                      alt <- sapply(vgr$ALT, function(x) x[1])
                      vgr$first = !grepl('^(\\]|\\[)', alt) ## ? is this row the "first breakend" in the ALT string (i.e. does the ALT string not begin with a bracket)
                      vgr$right = grepl('\\[', alt) ## ? are the (sharp ends) of the brackets facing right or left
                      vgr$coord = as.character(paste(seqnames(vgr), ':', start(vgr), sep = ''))
                      vgr$mcoord = as.character(gsub('.*(\\[|\\])(.*\\:.*)(\\[|\\]).*', '\\2', alt))
                      vgr$mcoord = gsub('chr', '', vgr$mcoord)

                      if (all(is.na(vgr$mateid)))
                          if (!is.null(names(vgr)) & !any(duplicated(names(vgr))))
                              {
                                  warning('MATEID tag missing, guessing BND partner by parsing names of vgr')
                                  vgr$mateid = paste(gsub('::\\d$', '', names(vgr)), (sapply(strsplit(names(vgr), '\\:\\:'), function(x) as.numeric(x[length(x)])))%%2 + 1, sep = '::')
                              }
                          else if (!is.null(vgr$SCTG))
                              {
                                  warning('MATEID tag missing, guessing BND partner from coordinates and SCTG')
                                  require(igraph)
                                  ucoord = unique(c(vgr$coord, vgr$mcoord))
                                  vgr$mateid = paste(vgr$SCTG, vgr$mcoord, sep = '_')

                                  if (any(duplicated(vgr$mateid)))
                                      {
                                          warning('DOUBLE WARNING! inferred mateids not unique, check VCF')
                                          bix = bix[!duplicated(vgr$mateid)]
                                          vgr = vgr[!duplicated(vgr$mateid)]
                                      }
                              }
                          else
                              stop('MATEID tag missing')

                      vgr$mix = as.numeric(match(vgr$mateid, names(vgr)))

                      pix = which(!is.na(vgr$mix))

                      vgr.pair = vgr[pix]

                      if (length(vgr.pair)==0)
                        stop('No mates found despite nonzero number of BND rows in VCF')

                      vgr.pair$mix = match(vgr.pair$mix, pix)
                      vix = which(1:length(vgr.pair)<vgr.pair$mix )
                      vgr.pair1 = vgr.pair[vix]
                      vgr.pair2 = vgr.pair[vgr.pair1$mix]

                      if (length(vgr.pair1)==0 | length(vgr.pair2)==0)
                        stop('No mates found despite nonzero number of BND rows in VCF')

                      ## now need to reorient pairs so that the breakend strands are pointing away from the breakpoint
                      
                      ## if "first" and "right" then we set this entry "-" and the second entry "+"
                      tmpix = vgr.pair1$first & vgr.pair1$right
                      if (any(tmpix))
                          {
                              strand(vgr.pair1)[tmpix] = '-'
                              strand(vgr.pair2)[tmpix] = '+'
                          }

                      ## if "first" and "left" then "-", "-"
                      tmpix = vgr.pair1$first & !vgr.pair1$right
                      if (any(tmpix))
                          {
                              strand(vgr.pair1)[tmpix] = '-'
                              strand(vgr.pair2)[tmpix] = '-'
                          }

                      ## if "second" and "left" then "+", "-"
                      tmpix = !vgr.pair1$first & !vgr.pair1$right
                      if (any(tmpix))
                          {
                              strand(vgr.pair1)[tmpix] = '+'
                              strand(vgr.pair2)[tmpix] = '-'
                          }

                      ## if "second" and "right" then "+", "+"
                      tmpix = !vgr.pair1$first & vgr.pair1$right
                      if (any(tmpix))
                          {
                              strand(vgr.pair1)[tmpix] = '+'
                              strand(vgr.pair2)[tmpix] = '+'
                          }

                      pos1 = as.logical(strand(vgr.pair1)=='+') ## positive strand junctions shift left by one (i.e. so that they refer to the base preceding the break for these junctions
                      if (any(pos1))
                          {
                              start(vgr.pair1)[pos1] = start(vgr.pair1)[pos1]-1
                              end(vgr.pair1)[pos1] = end(vgr.pair1)[pos1]-1
                          }

                      pos2 = as.logical(strand(vgr.pair2)=='+') ## positive strand junctions shift left by one (i.e. so that they refer to the base preceding the break for these junctions
                      if (any(pos2))
                          {
                              start(vgr.pair2)[pos2] = start(vgr.pair2)[pos2]-1
                              end(vgr.pair2)[pos2] = end(vgr.pair2)[pos2]-1
                          }
                      ra = grl.pivot(GRangesList(vgr.pair1[, c()], vgr.pair2[, c()]))

                      this.inf = values(vgr)[bix[pix[vix]], ]

                      if (is.null(this.inf$POS))
                          this.inf = cbind(data.frame(POS = ''), this.inf)
                      if (is.null(this.inf$CHROM))
                          this.inf = cbind(data.frame(CHROM = ''), this.inf)

                      if (is.null(this.inf$MATL))
                          this.inf = cbind(data.frame(MALT = ''), this.inf)
                      
                      this.inf$CHROM = seqnames(vgr.pair1)
                      this.inf$POS = start(vgr.pair1)
                      this.inf$MATECHROM = seqnames(vgr.pair2)
                      this.inf$MATEPOS = start(vgr.pair2)
                      this.inf$MALT = vgr.pair2$AL

                      ## NOT SURE WHY BROKEN
                      ## tmp = tryCatch(cbind(values(vgr)[bix[pix[vix]],], this.inf), error = function(e) NULL)
                      ## if (!is.null(tmp))
                      ##     values(ra) = tmp
                      ## else
                      ##     values(ra) = cbind(vcf@fixed[bix[pix[vix]],], this.inf)
                      
                      values(ra) = this.inf
                      
                      if (is.null(values(ra)$TIER))
                          values(ra)$tier = ifelse(values(ra)$FILTER == "PASS", 2, 3) ## baseline tiering of PASS vs non PASS variants
                      else
                          values(ra)$tier = values(ra)$TIER

                      if (!get.loose)
                          return(ra)
                      else
                          {
                              npix = is.na(vgr$mix)
                              vgr.loose = vgr[npix, c()] ## these are possible "loose ends" that we will add to the segmentation

                              ## NOT SURE WHY BROKEN                              
                              tmp =  tryCatch( values(vgr)[bix[npix], ],
                                  error = function(e) NULL)                             
                              if (!is.null(tmp))
                                  values(vgr.loose) = tmp
                              else
                                  values(vgr.loose) = cbind(vcf@fixed[bix[npix], ], info(vcf)[bix[npix], ])
                              
                              return(list(junctions = ra, loose.ends = vgr.loose))
                          }
                  }
              else
                  rafile = read.delim(rafile)
          }
            
     if (is.data.table(rafile))
         {
             require(data.table)
             rafile = as.data.frame(rafile)
         }

    if (nrow(rafile)==0)
        {
            out = GRangesList()
            values(out) = rafile
            return(out)
        }
  
    if (snowman) ## flip breaks so that they are pointing away from junction
      {
        rafile$str1 = ifelse(rafile$strand1 == '+', '-', '+')
        rafile$str2 = ifelse(rafile$strand2 == '+', '-', '+')
      }
      
    if (!is.null(seqlevels)) ## convert seqlevels from notation in tab delim file to actual
      {
        rafile$chr1 = seqlevels[rafile$chr1]
        rafile$chr2 = seqlevels[rafile$chr2]        
      }

     
    if (is.null(rafile$str1))
      rafile$str1 = rafile$strand1

    if (is.null(rafile$str2))
      rafile$str2 = rafile$strand2
     if (!is.null(rafile$pos1) & !is.null(rafile$pos2))
         {
             if (breakpointer)
                 {
                     rafile$pos1 = rafile$T_BPpos1
                     rafile$pos2 = rafile$T_BPpos2
                 }
             
             if (!is.numeric(rafile$pos1))
                 rafile$pos1 = as.numeric(rafile$pos1)

             if (!is.numeric(rafile$pos2))
                 rafile$pos2 = as.numeric(rafile$pos2)

             ## clean the parenthesis from the string

             rafile$str1 <- gsub('[()]', '', rafile$str1)
             rafile$str2 <- gsub('[()]', '', rafile$str2)

             ## goal is to make the ends point <away> from the junction where - is left and + is right
             if (is.character(rafile$str1) | is.factor(rafile$str1))
                 rafile$str1 = gsub('0', '-', gsub('1', '+', gsub('\\-', '1', gsub('\\+', '0', rafile$str1))))
             
             if (is.character(rafile$str2) | is.factor(rafile$str2))
                 rafile$str2 = gsub('0', '-', gsub('1', '+', gsub('\\-', '1', gsub('\\+', '0', rafile$str2))))

             
             if (is.numeric(rafile$str1))
                 rafile$str1 = ifelse(rafile$str1>0, '+', '-')

             if (is.numeric(rafile$str2))
                 rafile$str2 = ifelse(rafile$str2>0, '+', '-')
             
             rafile$rowid = 1:nrow(rafile)

             bad.ix = is.na(rafile$chr1) | is.na(rafile$chr2) | is.na(rafile$pos1) | is.na(rafile$pos2) | is.na(rafile$str1) | is.na(rafile$str2) | rafile$str1 == '*'| rafile$str2 == '*' | rafile$pos1<0 | rafile$pos2<0
             
             rafile = rafile[which(!bad.ix), ]
             
             if (nrow(rafile)==0)
                 return(GRanges())
             
             seg = rbind(data.frame(chr = rafile$chr1, pos1 = rafile$pos1, pos2 = rafile$pos1, strand = rafile$str1, ra.index = rafile$rowid, ra.which = 1, stringsAsFactors = F),
                 data.frame(chr = rafile$chr2, pos1 = rafile$pos2, pos2 = rafile$pos2, strand = rafile$str2, ra.index = rafile$rowid, ra.which = 2, stringsAsFactors = F))

             if (chr.convert)
                 seg$chr = gsub('chr', '', gsub('25', 'M', gsub('24', 'Y', gsub('23', 'X', seg$chr))))
             
             out = seg2gr(seg, seqlengths = seqlengths)[, c('ra.index', 'ra.which')];
             out = split(out, out$ra.index)
         }
     else if (!is.null(rafile$start1) & !is.null(rafile$start2) & !is.null(rafile$end1) & !is.null(rafile$end2))
         {
             ra1 = gr.flipstrand(GRanges(rafile$chr1, IRanges(rafile$start1, rafile$end1), strand = rafile$str1))
             ra2 = gr.flipstrand(GRanges(rafile$chr2, IRanges(rafile$start2, rafile$end2), strand = rafile$str2))
             out = grl.pivot(GRangesList(ra1, ra2))             
         }
     
     
     if (keep.features)
         values(out) = rafile[, ]

      if (!get.loose)
          return(out)
      else                              
          return(list(junctions = out, loose.ends = GRanges()))

     return(out)
 }


## 
##
.correct.slack = function(ra.sol)
  {
    adj = ra.sol$adj
    eslack.out = ra.sol$segstats$eslack.out
    eslack.in = ra.sol$segstats$eslack.in
    cn = ra.sol$segstats$cn
    str = as.logical(strand(ra.sol$segstats) == '+')
    sn = as.character(seqnames(ra.sol$segstats))
    for (i in (1:length(cn))[c(-1, -length(cn))])
      {
        if (str[i])
          {
            if (sn[i] == sn[i-1] & eslack.out[i-1]>0 & eslack.in[i]>0)
              {
                cat('+')
                a = min(eslack.out[i-1], eslack.in[i])
                eslack.out[i-1] = eslack.out[i-1]-a
                eslack.in[i] = eslack.in[i]-a
                adj[i-1, i] = adj[i-1, i] + a
              }
          }
        else if (!str[i] & sn[i] == sn[i+1] & eslack.out[i+1]>0 & eslack.in[i]>0)
          {
            cat('-')
            a = min(eslack.out[i+1], eslack.in[i])
            eslack.out[i+1] = eslack.out[i+1]-a
            eslack.in[i] = eslack.in[i]-a
            adj[i+1, i] = adj[i+1, i] + a
          } 
      }
    ra.sol$adj = adj;
    ra.sol$segstats$eslack.out = eslack.out
    ra.sol$segstats$eslack.in = eslack.in
    ra.sol$eslack.out = eslack.out
    ra.sol$eslack.in = eslack.in
    cat('\n')
    return(ra.sol)
  }


## cplex set max threads (warning can only do once globally per machine, so be wary of multiple hosts running on same machine)
##
.cplex_customparams = function(out.file, numthreads = 0, nodefileind = NA, treememlim = NA)
{
  param_lines = "CPLEX Parameter File Version 12.6.0.0"    
  
  param_lines = c(param_lines, paste("CPX_PARAM_THREADS", numthreads, sep = '\t'))
  
  if (!is.na(nodefileind))
    param_lines = c(param_lines, paste("CPX_PARAM_NODEFILEIND", nodefileind, sep = '\t'))
  
  if (!is.na(treememlim))
    {
#'      param_lines = c(param_lines, paste("CPX_PARAM_WORKDIR", getwd(), sep = '\t'))
      param_lines = c(param_lines, paste("CPX_PARAM_TRELIM", treememlim, sep = '\t'))
    }
  
  writeLines(param_lines, out.file)
  Sys.setenv(ILOG_CPLEX_PARAMETER_FILE=out.file)
}
  
######################## utilities
########################

## prints Ax = b equations encoded by matrix Ab whose last column is "b"
##
print_eq = function(Ab, xlabels = NULL)
  {
    Ab = as.matrix(Ab)
        
    A = Ab[, -ncol(Ab), drop = F]
    b = Ab[, ncol(Ab)]
    
    if (!is.null(xlabels))
      colnames(A) = xlabels
    else if (is.null(colnames(A)))
      colnames(A) = paste('x', 1:ncol(A), sep = '')
    else if (length(ix <- which(nchar(colnames(A)) == 0)) != 0)
      colnames(A)[ix] = paste('x', ix, sep = '')
       
    cat('', paste(sapply(1:nrow(A), function(x,y)
                        {
                          ix = which(A[x, ]!=0)
                          if (length(ix)>0)
                            paste(signif(as.numeric(A[x, ix]),2), y[ix], sep = ' ', collapse = ' +\t')
                          else
                            '0'
                        },
                        colnames(A)), '\t=', b, '\n'))
  }


#' jgraph
#'
#' takes in a jabba object and threshold for clusters and "quasi-reciprocal"
#' junctions
#'
#' @name jgraph
#' @export
jgraph = function(jab, thresh_cl = 1e6, all = FALSE, thresh_r = 1e3, clusters = FALSE)
{   
    ## identify sinks and sources
    jix = which(values(jab$junctions)$cn>0)    
    so = gr.end(jab$segstats[as.integer(jab$ab.edges[jix,1,1:2])], ignore.strand = FALSE)
    si = gr.start(jab$segstats[as.integer(jab$ab.edges[jix,2,1:2])], ignore.strand = FALSE)
    jid = rep(jix, 2)
    recip = rep(c(FALSE, TRUE), each = length(jix))


    ## compute directed and undirected distances from sinks to sources
    gundir = jabba.dist(jab, si, so, directed = FALSE)
    gdir = jabba.dist(jab, si, so, directed = TRUE)


    gedges = gdir<=thresh_cl
    redges = gundir<=thresh_r & gdir>thresh_cl ## reciprocal edges are indirect that are non-indirect connections


    ## pick only closest g edge leaving a node
    ## this will be asymmetric
    ## since a_j may be closest to b_i
    ## but aa_i may not be closest to bb_j
    gfilt = t(apply(gdir, 1, function(x) sign((1:length(x) %in% which.min(x)))))

    ## key mathematical property:
    ## (identical to stranded adjacency matrix in jabba)
    ## recip, recip is identical to t(!recip, !recip)
    ## BUT recip, !recip and !recip, recip are distinct and symmetric

    ## just try these checks all should be empty
    ## which(gedges[!recip, !recip] != t(gedges[recip, recip]), arr.ind = TRUE)
    ## which(gedges[!recip, recip] != t(gedges[!recip, recip]), arr.ind = TRUE)
    ## which(gedges[recip, !recip] != t(gedges[recip, !recip]), arr.ind = TRUE)


    Gg.dist = pmin(gdir[recip, recip], gdir[!recip, recip], gdir[recip, !recip], gdir[!recip, recip])
    Gr.dist = pmin(gundir[recip, recip], gundir[!recip, recip], gundir[recip, !recip], gundir[!recip, recip])


    ## Gg and Gr should be symmetric
    Gg = sign(gedges[recip, recip]) + sign(gedges[!recip, recip]) + sign(gedges[recip, !recip]) + sign(gedges[!recip, !recip])
    Gr = sign(redges[recip, recip]) + sign(redges[recip, !recip]) + sign(redges[!recip, recip]) + sign(redges[!recip, !recip])
    Gr[cbind(1:nrow(Gr), 1:nrow(Gr))] = 0 ## no self reciprocal edges, as tempting as that might be

    ## Ggf won't be symmetric but we will make it and also make sure we don't remove self loops
    Ggf = sign(gfilt[recip, recip]) + sign(gfilt[!recip, recip]) + sign(gfilt[recip, !recip]) + sign(gfilt[!recip, !recip])
    Ggf = Ggf+t(Ggf) + diag(rep(1, nrow(Ggf)))

    M = Matrix(0,
               nrow = length(jab$junctions),
               ncol = length(jab$junctions))
    
    out = list(Gg = M,
               Gr = M,
               Ggd = M,
               Grd = M)                      
    
    out$Gg[jix, jix] = sign(Gg)*sign(Ggf)
    out$Gr[jix, jix] = sign(Gr)

    ## some fake reciprocals leak through somehow - get rid!! (TOFIX)
                                        #    out$Gr[which(out$Gg!=0, arr.ind = TRUE)] = 0

    if (clusters)
    {
        tmp.out = tryCatch(split(1:nrow(out$Gr),
                             igraph::clusters(igraph::graph.adjacency(out$Gg + out$Gr))$membership), error = function(e) NULL)
        if (is.null(tmp.out))
            tmp.out = split(1:nrow(out$Gr),
                                 igraph::clusters(igraph::graph.adjacency(as(out$Gg + out$Gr, 'matrix')))$membership)
        out = tmp.out[rev(order(sapply(tmp.out, length)))]
        names(out) = 1:length(out)
    }

    return(out)
    
}

######################## wrappers
######################## todo: wrap these into a single "jabba" call

#' segment
seg_stub = function(cov.file, out.file, field = 'ratio', mc.cores = 4)
  {
    ## resegment from scratch, using cbs, and recompute karyoMIP w out edge cons,
    cat('reading', cov.file, '\n')
#'    this.cov = seg2gr(read.delim(cov.file, strings = F), seqlengths = hg_seqlengths())
    this.cov = readRDS(cov.file)
    cat('starting segmentation')
    out = cumseg(this.cov, field = field, verbose = T, mc.cores = mc.cores)
    cat('finished segmenting\n')
    cat(length(out), ' segments produced\n')
    saveRDS(out, out.file)    
  }

#' karyograph
karyograph_stub = function(seg.file, ## path to rds file of initial genome partition (values on segments are ignored)
    cov.file, ## path to rds file GRanges of coverage with meta data field "field"
    nseg.file = NULL, ## rds file of GRanges providing integer copy numbers for normal segments in the genome
    het.file = NULL,
    ra = NULL,
    dranger.file = NULL,
    out.file,
    ra.file = NULL,
    force.seqlengths = NULL,
    purity =  NA,
    ploidy = NA,    
    field = 'ratio', mc.cores = 1, max.chunk = 1e8, subsample = NULL)
    {
        loose.ends = GRanges()

    if (!is.null(ra))
        this.ra = ra
    else
    {
        
        if (!is.null(dranger.file))
        {
            this.dranger = read.delim(dranger.file, strings = F)
            
            if (ncol(this.dranger)<=1) ## wrong separator
                this.dranger = read.delim(dranger.file, sep = ',', strings = F)
            
            if (!is.null(this.dranger$strand1) & !is.null(this.dranger$strand2)) 
            {
                ## looks like snowman input flip breaks so that they are pointing away from junction
                this.dranger$str1 = ifelse(this.dranger$strand1 == '+', '-', '+')
                this.dranger$str2 = ifelse(this.dranger$strand2 == '+', '-', '+')
            }
            
            this.dranger$chr1 = gsub('23', 'X', gsub('24', 'Y', this.dranger$chr1))
            this.dranger$chr2 = gsub('23', 'X', gsub('24', 'Y', this.dranger$chr2))
            this.ra = ra_breaks(this.dranger, seqlengths = hg_seqlengths())
        }        
        else if (grepl('(\\.bedpe)|(\\.vcf$)|(\\.vcf\\.gz$)', ra.file))
        {        
            tmp.ra = ra_breaks(ra.file, seqlengths = hg_seqlengths(), get.loose = T)
            if (length(tmp.ra)==0){
                this.ra = gr.fix(GRangesList(), hg_seqlengths())
                loose.ends = GRanges(seqlengths = hg_seqlengths())
            } else {
                this.ra = tmp.ra$junctions
                loose.ends = tmp.ra$loose.ends
            }
        }
        else
            this.ra = readRDS(ra.file)
    }
    ## if we don't have normal segments then coverage file will be our bible for seqlengths


    if (is.character(cov.file))
    {
        if (grepl('\\.rds$', cov.file))        
            this.cov = readRDS(cov.file)
        else
        {
            this.cov = import.ucsc(cov.file)
            field = 'score';
        }
    }
    else
        this.cov = cov.file
    
    ## now make sure we have the "best" seqlengths
    .fixsl = function(sl, gr) {sl[seqlevels(gr)] = pmax(seqlengths(gr), sl[seqlevels(gr)]); return(sl)}

    if (is.null(force.seqlengths))
        sl = .fixsl(seqlengths(this.ra), this.cov)
    else
        sl = .fixsl(force.seqlengths, this.cov)
                                        #      sl = .fixsl(seqlengths(this.ra), this.cov)

    if (!is.null(nseg.file))
    {
        if (is.character(nseg.file))              
            nseg = readRDS(nseg.file)
        else
            nseg = nseg.file          
        sl = .fixsl(sl, nseg)
    }
    
    ## make sure all sl's are equiv
    if (is.character(seg.file))
        this.seg = gr.fix(readRDS(seg.file), sl, drop = T)[, c()]
    else
        this.seg = seg.file
    
    if (length(loose.ends>0))
    {
        cat('Adding loose ends from vcf file to seg file')
        this.seg = grbind(this.seg, gr.fix(loose.ends, sl, drop = T))              
    }
    
    this.ra = gr.fix(this.ra, sl, drop = T)

    this.kag = karyograph(this.ra, this.seg)

    if (is.null(nseg.file))
        this.kag$segstats$ncn = 2
    
    cat('Computing segstats\n')
    hets.gr = NULL

    if (!is.null(het.file))
    {
        hets = fread(het.file)
        cat('loaded hets\n')
        hets.gr = seg2gr(hets[pmin(ref.frac.n, 1-ref.frac.n) > 0.2 & (ref.count.n + alt.count.n)>20, ])
    }







    if (length(hets.gr)>0){
        ## pretend we don't have hets at all
        
        this.kag$segstats = segstats(this.kag$tile, this.cov, field = field, prior_weight = 1, max.chunk = max.chunk, subsample = subsample, asignal = hets.gr, afield = c('ref.count.t', 'alt.count.t'), mc.cores = mc.cores)
    }
    else
        this.kag$segstats = segstats(this.kag$tile, this.cov, field = field, prior_weight = 1, max.chunk = max.chunk, subsample = subsample, mc.cores = mc.cores)
    
    this.kag$segstats$ncn = 2
    
    if (!is.null(nseg.file))
        if (is.null(nseg$cn))        
            stop('Normal seg file does not have "cn" met data field')
        else
        {
            this.kag$segstats$ncn = round(gr.val(this.kag$segstats, nseg, 'cn')$cn)
            this.kag$segstats$mean[is.na(this.kag$segstats$ncn)] = NA ## remove segments for which we have no normal copy number
        }
    
    
    ## 6/15 temp fix for sd on short segments, which we overestimate for now
    cov.thresh = pmin(1e5, median(width(this.cov)))
    cat('!!!!!!!!!!! cov.thresh for fix.sd is', cov.thresh, '\n')
    fix.sd  = width(this.kag$segstats)<(3*cov.thresh)
                                        #    this.kag$segstats$mean[make.na] = NA
    this.kag$segstats$sd[fix.sd] = sqrt(this.kag$segstats$mean[fix.sd])
    
    #'      if (is.character(tryCatch(png(paste(out.file, '.ppgrid.png', sep = ''), height = 500, width = 500), error = function(e) 'bla')))

    ss.tmp = this.kag$segstats[width(this.kag$segstats)>1e4, ] ## don't use ultra short segments
    pdf(paste(out.file, '.ppgrid.pdf', sep = ''), height = 10, width = 10)

    if (!is.null(het.file))
        pp = ppgrid(ss.tmp, verbose = T, plot = F, mc.cores = mc.cores,
                    purity.min = ifelse(is.na(purity), 0, purity), purity.max = ifelse(is.na(purity),1, purity),
                    ploidy.min = ifelse(is.na(ploidy), 1.2, ploidy), ploidy.max = ifelse(is.na(ploidy), 6, ploidy), allelic = TRUE)
    else
        pp = ppgrid(ss.tmp, verbose = T, plot = F, mc.cores = mc.cores,
                    purity.min = ifelse(is.na(purity), 0, purity), purity.max = ifelse(is.na(purity),1, purity),
                    ploidy.min = ifelse(is.na(ploidy), 1.2, ploidy), ploidy.max = ifelse(is.na(ploidy), 6, ploidy), allelic = FALSE)
    
    dev.off()
    
    saveRDS(pp, paste(out.file, '.ppgrid.solutions.rds', sep = '')) ## save alternate solutions
    
    this.kag$purity = pp[1,]$purity
    this.kag$ploidy = pp[1,]$ploidy
    this.kag$beta = pp[1,]$beta
    this.kag$gamma = pp[1,]$gamma
    this.kag$segstats$cn = rel2abs(this.kag$segstats, beta = this.kag$beta, gamma = this.kag$gamma, field = 'mean')
    
    if (is.character(tryCatch(png(paste(out.file, '.ppfit.png', sep = ''), height = 1000, width = 1000), error = function(e) 'bla')))
        pdf(paste(out.file, '.ppfit.pdf', sep = ''), height = 10, width = 10)
    tmp.kag = this.kag

    tryres <- try( tmp.kag$segstats <- tmp.kag$segstats[ which( tmp.kag$segstats$ncn == 2 & !fix.sd & as.logical( strand(tmp.kag$segstats) =='+' ) ) ] ) 
    if ( class( tryres ) == "try-error" ) browser()
    
    if (length(tmp.kag$segstats)<10)
        warning('number of segments used for purity ploidy extremely low .. check coverage data')
    .plot_ppfit(tmp.kag)            
    dev.off()

    y1 = 10

    if (is.character(tryCatch(png(paste(out.file, '.inputdata.png', sep = ''), height = 1000, width = 1000), error = function(e) 'bla')))
        pdf(paste(out.file, '.inputdata.pdf', sep = ''), height = 10, width = 10)
    
    plot(c(gTrack(gr.fix(sample(this.cov, 5e4), this.kag$segstats), y.field = 'ratio', col = alpha('black', 0.3)),
           gTrack(this.kag$segstats, y.field = 'mean', angle = 0, col = 'gray10', border = alpha('black', 0.2))), links = this.kag$junctions, y1 = y1)
    dev.off()

    saveRDS(this.kag, out.file)
}

## diagnostic function used by karyograph_stub
.plot_ppfit = function(kag, xlim = c(-Inf, Inf))
{
    tmp = kag$segstats  ## only plot seg that we haven't fixed SD for and that have normal cn 1, to minimize confusion
    dupval = sort(table(tmp$mean), decreasing = TRUE)[1]
    if (!is.na(dupval))
        if (dupval>5)
            tmp = tmp[-which(as.character(tmp$mean) == names(dupval))]

    if (length(tmp)==0)
        return()
    segsamp = pmin(sample(tmp$mean, 1e6, replace = T, prob = width(tmp)), xlim[2])
    hist(pmax(xlim[1], pmin(xlim[2], segsamp)), 1000, xlab = 'Segment intensity', main = sprintf('Purity: %s Ploidy: %s Beta: %s Gamma: %s', kag$purity, kag$ploidy, round(kag$beta,2), round(kag$gamma,2)), xlim = c(pmax(0, xlim[1]), pmin(xlim[2], max(segsamp, na.rm = T))))
    abline(v = 1/kag$beta*(0:1000) + kag$gamma/kag$beta, col = alpha('red', 0.3), lty = c(4, rep(2, 1000)))
}

#' run cbs
cbs_stub = function(cov.file, out.file, field = 'ratio', alpha = 1e-5)
{
    require(DNAcopy)
    ## resegment from scratch, using cbs, and recompute karyoMIP w out edge cons,

    cat('reading', cov.file, '\n')
                                        #   this.cov = seg2gr(read.delim(cov.file, strings = F), seqlengths = hg_seqlengths())
    this.cov = readRDS(cov.file)
    ix = which(!is.na(values(this.cov)[, field]))
    cat('sending ', length(ix), ' segments\n')
    cna = CNA(log(values(this.cov)[, field])[ix], as.character(seqnames(this.cov))[ix], start(this.cov)[ix], data.type = 'logratio')
    cat('finished making cna\n')
    seg = segment(smooth.CNA(cna), alpha = 1e-5, verbose = T)
    cat('finished segmenting\n')
    out = seg2gr(print(seg))
    cat(length(out), ' segments produced\n')
    saveRDS(out, out.file)    
}

#' run ramip
ramip_stub = function(kag.file, out.file, mc.cores = 1, mem = 16, tilim = 1200, slack.prior = 0.001, gamma = NA, beta = NA, customparams = T,
                      purity.min = NA, purity.max = NA,
                      ploidy.min = NA, ploidy.max = NA, 
                      edge.nudge = 0,  ## can be scalar (equal nudge to all ab junctions) or vector of length readRDS(kag.file)$junctions
                      ab.force = NULL ## indices of aberrant junctions to force include into the solution
                      )
{
    cat('Starting', kag.file, '\n')

    this.kag = readRDS(kag.file)    

    if (is.null(this.kag$gamma) | is.null(this.kag$beta))
    {
        pp = ppgrid(this.kag$segstats, verbose = T, plot = T, purity.min = purity.min, purity.max = purity.max, ploidy.min = ploidy.min, ploidy.max = ploidy.max)
        this.kag$beta = pp[1,]$beta
        this.kag$gamma = pp[1,]$gamma
    }

    if (!is.na(gamma))
    {
        cat(sprintf('Overriding gamma with %s\n', gamma))
        this.kag$gamma = gamma
    }

    if (!is.na(beta))
    {
        cat(sprintf('Overriding beta with %s\n', beta))
        this.kag$beta = beta
    }
    
    if (customparams)
    {
        max.threads = Sys.getenv("LSB_DJOB_NUMPROC")
        if (nchar(max.threads) == 0)
            max.threads = Inf
        else
            max.threads = as.numeric(max.threads)
        max.threads = min(max.threads, mc.cores)
        if (is.infinite(max.threads))
            max.threads = 0 

        param.file = paste(out.file, '.prm', sep = '')
        .cplex_customparams(param.file, max.threads, treememlim = mem * 1e3)

        Sys.setenv(ILOG_CPLEX_PARAMETER_FILE = normalizePath(param.file))
        print(Sys.getenv('ILOG_CPLEX_PARAMETER_FILE'))
    }
    
    adj.nudge = this.kag$adj*0;
    adj.nudge[this.kag$ab.edges[,1:2,1]] = 1*edge.nudge ## if edge.nudge is length ab.edges, then corresponding edges will be nudged

    adj.lb = NULL
    if (!is.null(ab.force))
    {
        print(paste('Enforcing lower bounds on aberrant junctions:', paste(ab.force, collapse = ',')))
        adj.lb = this.kag$adj*0
        adj.lb[rbind(this.kag$ab.edges[ab.force, ,1])[, 1:2, drop = FALSE]] = 1
        adj.lb[rbind(this.kag$ab.edges[ab.force, ,2])[, 1:2, drop = FALSE]] = 1
    }

    ra.sol = jbaMIP(this.kag$adj, this.kag$segstats, beta.guess = this.kag$beta, gamma.guess = this.kag$gamma, tilim = tilim, slack.prior = slack.prior, cn.prior = NA, mipemphasis = 0, ignore.cons = T, adj.lb = adj.lb , adj.nudge = adj.nudge, cn.ub = rep(500, length(this.kag$segstats)), verbose = T)
    saveRDS(ra.sol, out.file)

    if (customparams)
    {
        system(paste('rm', param.file))
        Sys.setenv(ILOG_CPLEX_PARAMETER_FILE='')    
        cat('Finished', kag.file, '\n')
    }
}

#' run super-enhancer analysis
#' 
superenhancer_stub = function(rg.span.file, se.bed.file, this.ra.file, out.file, verbose = T, max.dist = 1e6)  
{
    library(RColorBrewer)
    rg.span = readRDS(rg.span.file)
    se.bed = readRDS(se.bed.file)
    this.ra = ra_breaks(read.delim(this.ra.file, strings = F))
    px = proximity(rg.span, se.bed, this.ra, verbose = T, max.dist=1e6)
    saveRDS(px, out.file)
}

#' run gene-gene px analysis
#' 
genepx_stub = function(rg.span.file, this.ra.file, out.file, verbose = T, max.dist = 1e6)  
{
    rg.span = readRDS(rg.span.file)
    this.ra = ra_breaks(read.delim(this.ra.file, strings = F))
    px = proximity(rg.span, rg.span, this.ra, verbose = T, max.dist=1e6)
    saveRDS(px, out.file)
}

#' run gene-chromhmm px analysis
#' 
chromhmmpx_stub = function(rg.span.file, chromhmm.file, this.ra.file, out.file, verbose = T, max.dist = 1e6)
{
    rg.span = readRDS(rg.span.file)
    chromhmm = readRDS(chromhmm.file)
    this.ra = ra_breaks(read.delim(this.ra.file, strings = F))
    px = proximity(rg.span, chromhmm, this.ra, verbose = T, max.dist=1e6)
    saveRDS(px, out.file)
}


#' gr.tile.map
#'
#' Given two tilings of the genome (eg at different resolution)
#' query and subject outputs a length(query) list whose items are integer vectors of indices in subject
#' overlapping that overlap that query (strand non-specific)
#'
#' @note Assumes that input query and subject have no gaps (including at end) or overlaps, i.e. ignores end()
#' coordinates and only uses "starts"
#' @param query Query
#' @param subject Subject
#' @param mc.cores number of cores
#' @param verbose Default FALSE
#' @export
############################################
gr.tile.map = function(query, subject, mc.cores = 1, verbose = FALSE)
{
    ix.q = order(query)
    ix.s = order(subject)

        q.chr = as.character(seqnames(query))[ix.q]
        s.chr = as.character(seqnames(subject))[ix.s]

        ql = split(ix.q, q.chr)
        sl = split(ix.s, s.chr)

        tmp = mcmapply(
            function(x,y)
	        {
                    if (length(y)==0)
                        return(NULL)
                    all.pos = c(start(query)[x], start(subject)[y])
                    is.q = c(rep(T, length(x)), rep(F, length(y)))
                    all.ix = c(x, y)
                    ord.ix = order(all.pos)
                    all.pos = all.pos[ord.ix]
                    is.q = is.q[ord.ix]
                    all.ix = all.ix[ord.ix]
                    out = matrix(NA, nrow = length(all.pos), ncol = 2)
                    last.x = last.y = NA
                    for (i in 1:length(all.pos))
                        {
                            if (verbose)
                                if (i %% 100000 == 0) cat('Iteration', i, 'of', length(all.pos), '\n')
                            if (is.q[i])
                                {
                                    out[i, ] = c(all.ix[i], last.y)

                                    if (i<length(all.pos)) ## edge case where subject and query intervals share a start point, leading to two consecutive all.pos
                                        if (all.pos[i] == all.pos[i+1])
                                            out[i, ] = NA

                                    last.x = all.ix[i]
                                }
                            else
                                {
                                    out[i, ] = c(last.x, all.ix[i])

                                    if (i<length(all.pos)) ## edge case where subject and query intervals share a start point, leading to two consecutive all.pos
                                        if (all.pos[i] == all.pos[i+1])
                                            out[i, ] = NA

                                    last.y = all.ix[i]
                                }
                        }
                    out = out[Matrix::rowSums(is.na(out))==0, ]
                    return(out)
                }, ql, sl[names(ql)], mc.cores = mc.cores, SIMPLIFY = FALSE)

        m = munlist(tmp)[, -c(1:2), drop = FALSE]
        out = split(m[,2], m[,1])[as.character(1:length(query))]
        names(out) = as.character(1:length(query))
        return(out)
    }


##################################
#' @name vaggregate
#' @title vaggregate
#'
#' @description
#' same as aggregate except returns named vector
#' with names as first column of output and values as second
#'
#' Note: there is no need to ever use aggregate or vaggregate, just switch to data.table
#' 
#' @param ... arguments to aggregate
#' @return named vector indexed by levels of "by"
#' @author Marcin Imielinski
##################################
vaggregate = function(...)
  {
    out = aggregate(...);
    return(structure(out[,ncol(out)], names = do.call(paste, lapply(names(out)[1:(ncol(out)-1)], function(x) out[,x]))))
  }

write.tab = function(x, ..., sep = "\t", quote = F, row.names = F)
  {
      if (!is.data.frame(x))
          x = as.data.frame(x)

      write.table(x, ..., sep = sep, quote = quote, row.names = row.names)
  }

alpha = function(col, alpha)
{
  col.rgb = col2rgb(col)
  out = rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha)
  names(out) = names(col)
  return(out)
}


#' @name jab2json
#' @title jab2json
#'
#' @description
#' 
#' Dumps JaBbA graph into json
#'
#'
#'
#' @param jab input jab object
#' @param file output json file
#' @author Marcin Imielinski
jab2json = function(jab, file, maxcn = 100, maxweight = 100)
{    

    #' ++ = RL
    #' +- = RR
    #' -+ = LL
    qw = function(x) paste0('"', x, '"')

    ymin = 0;
    ymax = maxcn;
    
    nodes = jab$segstats %Q% (strand == "+")
    id = rep(1:length(nodes), 2)
    id.type = ifelse(nodes$loose, 'loose_end', 'interval')
    str = ifelse(as.character(strand(jab$segstats))=='+', 1, -1)

    node.dt = data.table(
        iid = 1:length(nodes),
        chromosome = qw(as.character(seqnames(nodes))),
        startPoint = as.character(start(nodes)),
        strand = "*",
        endPoint = as.character(end(nodes)),    
        title = as.character(1:length(nodes)),
        type = ifelse(nodes$loose, "loose_end", "interval"),
        y = pmin(maxcn, nodes$cn))

    aadj = jab$adj*0
    rix = which(rowSums(is.na(jab$ab.edges[, 1:2, '+']))==0)
    aadj[rbind(jab$ab.edges[rix, 1:2, '+'], jab$ab.edges[rix, 1:2, '+'])] = 1
    ed = which(jab$adj!=0, arr.ind = TRUE)

    if (nrow(ed)>0)
        {
            ed.dt = data.table(
                so = id[ed[,1]],
                so.str = str[ed[,1]],
                si = id[ed[,2]],
                weight = jab$adj[ed],
                title = "",
                type = ifelse(aadj[ed], 'ALT', 'REF'),
                si.str = str[ed[,2]])[, sig := ifelse(so<si,
                                                      paste0(so * so.str, '_', -si*si.str),
                                                      paste0(-si * si.str, '_', so*so.str)
                                                      )][!duplicated(sig), ][, cid := 1:length(weight), ][,
                                                                                                          ":="(so = so*so.str, si = -si*si.str)]             
            connections.json = ed.dt[, paste0(
                c("connections: [", paste(
                                        "\t{",
                                        "cid: ", cid,
                                        ", source: ", so,
                                        ", sink:", si,
                                        ", title: ", qw(title),
                                        ", type: ", qw(type),
                                        ", weight: ", pmin(maxweight, weight),
                                        "}",
                                        sep = "",
                                        collapse = ',\n'),             
                  "]"),
                collapse = '\n')
                ]    
        }

    intervals.json = node.dt[, paste0(
        c("intervals: [", paste(
                              "\t{",
                              "iid: ", iid,
                              ", chromosome: ", chromosome,
                              ", startPoint: ", startPoint,
                              ", endPoint: ", endPoint,
                              ", y: ", y,
                              ", title: ", qw(title),
                              ", type: ", qw(type),
                              ", strand: ", qw(strand),
                              "}",
                              sep = "",
                              collapse = ',\n'),
          "]"),
        collapse = '\n')
        ]

    meta.json =
        paste('meta: {\n\t',
              paste(
                  c(paste('"ymin:"', ymin),
                  paste('"ymax:"', ymax)),
                  collapse = ',\n\t'),
              '\n}')

    out = paste(c("var json = {",
                  paste(
                      c(meta.json,
                      intervals.json,              
                      connections.json),
                      collapse = ',\n'
                  ),"}"),                                                          
                  sep = "")
    
    writeLines(out, file)
}


#' @name gr2json
#' @title gr2json
#'
#' @description
#' 
#' Dumps GRanges into JSON with metadata features as data points in  "intervals"
#'
#'
#'
#' @param GRange input jab object
#' @param file output json file
#' @author Marcin Imielinski
gr2json = function(intervals, file, y = rep("null", length(intervals)), labels = '', maxcn = 100, maxweight = 100)
{    

    #' ++ = RL
    #' +- = RR
    #' -+ = LL
    qw = function(x) paste0('"', x, '"')

    ymin = 0;
    ymax = maxcn;
    
    nodes = intervals
    id = rep(1:length(nodes), 2)

    node.dt = data.table(
        iid = 1:length(nodes),
        chromosome = qw(as.character(seqnames(nodes))),
        startPoint = as.character(start(nodes)),
        strand = as.character(strand(nodes)),
        endPoint = as.character(end(nodes)),
        y = y, 
        title = labels)

    oth.cols = setdiff(names(values(nodes)), colnames(node.dt))
    node.dt = as.data.table(cbind(node.dt, values(nodes)[, oth.cols]))

    oth.cols = union('type', oth.cols)
    if (is.null(node.dt$type))
        node.dt$type = 'interval'
    
    intervals.json = node.dt[, paste0(
        c("intervals: [", paste(
                              "\t{",
                              "iid: ", iid,
                              ", chromosome: ", chromosome,
                              ", startPoint: ", startPoint,
                              ", endPoint: ", endPoint,
                              ", y: ", y,
                              ", title: ", qw(title),
                               ", strand: ", qw(strand),
                              eval(parse(text = ## yes R code making R code making JSON .. sorry .. adding additional columns
                                             paste0("paste0(",
                                                    paste0('", ', oth.cols, ':", qw(', oth.cols, ')', collapse = ','),
                                                    ")", collapse = ''))),                                  
                              "}",                              
                              sep = "",
                              collapse = ',\n'),
          "]"),
        collapse = '\n')
        ]

    meta.json =
        paste('meta: {\n\t',
              paste(
                  c(paste('"ymin:"', ymin),
                  paste('"ymax:"', ymax)),
                  collapse = ',\n\t'),
              '\n}')

    out = paste(c("var json = {",
                  paste(
                      c(meta.json,
                        intervals.json          
                        ),
                      collapse = ',\n'
                  ),"}"),                                                          
                sep = "")

    writeLines(out, file)
}
