## Marcin Imielinski
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


#' 
#' @import igraph
#' @import Matrix
#' @importFrom gplots col2hex
#' @import GenomicRanges
#' @import parallel
#' @import data.table
#' @import DNAcopy
#' @import gUtils
#' @importFrom graphics abline hist title
#' @importFrom grDevices col2rgb dev.off pdf png rgb
#' @importFrom stats C aggregate dist loess median ppois predict runif setNames
#' @importFrom utils read.delim write.table
#' @importFrom methods as is
#' @import gTrack
#' @useDynLib JaBbA
#' 

#' @name JaBbA
#' @title JaBbA
#' @description
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
#' @param junctions  GRangesList of junctions  (i.e. bp pairs with strands oriented AWAY from break) OR path to junction VCF file (BND format), dRanger txt file or rds of GRangesList
#' @param coverage  GRanges of coverage OR path to tsv of cov file w GRanges style columns, rds of GRanges or .wig / .bed file of (+/- normalized, GC corrected) fragment density
#' @param field  field of coverage GRanges to use as fragment density signal (only relevant if coverage is GRanges rds file)
#' @param seg  optional path to existing segmentation, if missing then will segment coverage using DNACopy with standard settings
#' @param cfield  character, junction confidence meta data field in ra
#' @param tfield  character, tier confidence meta data field in ra
#' @param outdir  out directory to dump into, default ./
#' @param nseg  optional path to normal seg file with $cn meta data field
#' @param hets  optional path to hets.file which is tab delimited text file with fields seqnames, start, end, alt.count.t, ref.count.t, alt.count.n, ref.count.n
#' @param name  prefix for sample name to be output to seg file
#' @param cores  number of cores to use (default 1)
#' @param nseg  path to data.frame or GRanges rds of normal seg file with coordinates and $cn data field specifying germline integer copy number
#' @param subsample  numeric between 0 and 1 specifying fraction with which to  sub-sample high confidence coverage data
#' @param tilim  integer scalar timeout in seconds for jbaMIP computation (default 1200 seconds)
#' @param mem  numeric scalar max memory in GB for MIP portion of algorithm (default 16)
#' @param edgenudge  numeric hyper-parameter of how much to nudge or reward aberrant junction incorporation, default 0.1 (should be several orders of magnitude lower than average 1/sd on individual segments), a nonzero value encourages incorporation of perfectly balanced rearrangements which would be equivalently optimal with 0 copies or more copies.
#' @param slack.penalty  penalty to put on every loose.end copy, should be calibrated with respect to 1/(k*sd)^2 for each segment, i.e. that we are comfortable with junction balance constraints introducing k copy number deviation from a segments MLE copy number assignment (the assignment in the absence of junction balance constraints)
#' @param overwrite  logical flag whether to overwrite existing output directory contents or just continue with existing files.
#' @param use.gurobi  logical flag specifying whether to use gurobi (if TRUE) instead of CPLEX (if FALSE) .. up to user to make sure the respective package is already installed
#' @param reiterate  integer scalar specifying how many (re-)iterations of jabba to do, rescuing lower tier junctions that are near loose ends (requires junctions to be tiered via a grangeslist or VCF metadata field $tfield), tiers are 1 = must use, 2 = may use, 3 = use only in iteration>1 if near loose end
#' @param allin if TRUE, use all available junctions except for tier 3 INDELs in the first interation
#' @param indel if TRUE, force the small isolated tier 2 events into the model
#' @param rescue.window integer scalar bp window around which to rescue lower tier junctions
#' @param strict logical flag specifying whether to only include junctions that exactly overlap segs
#' @param mc.cores integer how many cores to use to fork subgraphs generation (default = 1)
#' @param init jabba object (list) or path to .rds file containing previous jabba object which to use to initialize solution, this object needs to have the identical aberrant junctions as the current jabba object (but may have different segments and loose ends, i.e. is from a previous iteration)
#' @return gGraph (gGnome package) of balanced rearrangement graph
#' 
#' @examples
#'\dontrun{
#' library(JaBbA)
#' junctions = system.file("extdata", "junctions.vcf", package = 'JaBbA')
#' coverage = system.file("extdata", "coverage.txt", package = 'JaBbA')
#' 
#' ## run analysis without hets 
#' jab = JaBbA(junctions = junctions, coverage = coverage)
#'}
#' @import DNAcopy
#' @export
JaBbA = function(junctions, # path to junction VCF file, dRanger txt file or rds of GRangesList of junctions (with strands oriented pointing AWAY from breakpoint)
                 junctions.unfiltered = NULL, 
                 coverage, # path to cov file, rds of GRanges
                 seg = NULL, # path to seg file, rds of GRanges
                 outdir = './JaBbA', # out directory to dump into
                 cfield = NULL, # character, junction confidence meta data field in ra
                 tfield = "tier", # character, tier confidence meta data field in ra
                 nudge.balanced = FALSE,  ## if TRUE nudge chains of balanced (or quasi balanced junctions)
                 thresh.balanced = 500, ## threshold for balanced junctions
                 nseg = NULL, # path to normal seg file with $cn meta data field
                 hets = NULL, # path to hets.file which is tab delimited text file with fields seqnames, start, end, alt.count.t, ref.count.t, alt.count.n, ref.count.n
                 name = 'tumor', ## prefix for sample name to be output to seg file
                 purity = NA,
                 ploidy = NA,
                 field = 'ratio', ## character, meta data field to use from coverage object to indicate numeric coveragendance, coverage,
                 subsample = NULL, ## numeric scalar between 0 and 1, how much to subsample coverage per segment
                 tilim = 2400, ## timeout for MIP portion: 40 min per subgraph
                 ## mem = 32, ## max memory for MIP portion
                 reiterate = 0, ## how many (additional) times to iterate beyond the first iteration
                 rescue.window = 1e5, ## new 1e5 ## window around loose ends at which to rescue low tier junctions
                 init = NULL, ## previous JaBbA object to use as a solution
                 edgenudge = 0.1, ## hyper-parameter of how much to "nudge" or reward edge use, will be combined with cfield information if provided
                 use.gurobi = FALSE, ## use gurobi instead of CPLEX
                 slack.penalty = 1e2, ## nll penalty for each loose end copy
                 overwrite = FALSE, ## whether to overwrite existing output in outdir
                 mc.cores = 1,
                 strict = FALSE,
                 max.threads = Inf,
                 max.mem = 16,
                 indel = TRUE, ## default TRUE ## whether to force the small isolated tier 2 events into the model
                 all.in = FALSE, ## default FALSE ## whether to use all available junctions in the first interation
                 verbose = TRUE ## whether to provide verbose output
                 )
{
    system(paste('mkdir -p', outdir))
    jmessage('Starting analysis in ', normalizePath(outdir))
    if (overwrite)
        jmessage('Overwriting previous analysis contents')
    ra = junctions
    reiterate = 1 + as.integer(reiterate)

    if (is.character(ra))
    {
        if (!file.exists(ra))
        {
            stop(paste('Junction path', ra, 'does not exist'))
        }

        if (grepl("rds$", ra)){
            ra.all = readRDS(ra)
        } else {
            ra.all = read.junctions(ra)
        }
    } else {
        ra.all = ra
    }

    if (verbose)
    {
        jmessage("Read in ", length(ra.all), " total junctions")
    }

    if (length(ra.all)==0){
        if (verbose)
        {
            jmessage("Empty junction input. Will do just one round of optimization.")
        }
        reiterate = 0
    }

    ## xtYao Tuesday, Jun 19, 2018 04:52:17 PM
    ## Only when tier exists or unfiltered junctions provided, do we do the iterations
    ## if unfiltered set is given first parse it
    if (!is.null(junctions.unfiltered)){
        if (inherits(junctions.unfiltered, "character") & file.exists(junctions.unfiltered)){
            if (grepl(".rds$", junctions.unfiltered)){
                ra.uf = readRDS(junctions.unfiltered)
            } else {
                ra.uf = read.junctions(junctions.unfiltered)
            }
        } else if (inherits(junctions.unfiltered, "GRangesList")){
            ra.uf = junctions.unfiltered
        }
    }
    
    if (is.null(tfield)){
        tfield = 'tier'
    }

    ## if no tier field in junctions, set all of them to 2
    if (!(tfield %in% names(values(ra.all))))
    {
        warning("Tier field", tfield, "missing: giving every junction the same tier, i.e. all have the potential to be incorporated")
        values(ra.all)$tier = 2
    }

    if (exists("ra.uf")){
        ## merge ra.all with ra.uf
        ## junctions from ra.all will always have tier 2
        ra.all.uf = ra.merge(ra.all, ra.uf, pad=0, ind=TRUE) ## hard merge
        ## those match a record in junction, will be assigned to the tier in junction       
        values(ra.all.uf)$tier[which(!is.na(values(ra.all.uf)$seen.by.ra1))] =
                            values(ra.all)[, tfield][values(ra.all.uf)$seen.by.ra1]
        ## the rest will be tier 3
        values(ra.all.uf)$tier[which(is.na(values(ra.all.uf)$seen.by.ra1))] = 3
        ra.all = ra.all.uf
    }

    if (length(unique(values(ra.all)[, tfield]))==1) {
        jmessage("Only one tier of junctions found, cancel iteration if requested")
        reiterate = 1
    }

    ## final sanity check before running
    if (!all(unique(values(ra.all)[, tfield]) %in% 1:3)){
        stop(sprintf('Tiers in tfield can only have values 1,2,or 3'))
    }

    ## if we are iterating more than once
    if (reiterate>1){
        continue = TRUE
        this.iter = 1;

        values(ra.all)$id = 1:length(ra.all)
        saveRDS(ra.all, paste(outdir, '/junctions.all.rds', sep = ''))

        ## big change tonight, I'm gonna start with all of the tiers in the first round
        ## and then in each of following iterations keep the ones incorporated
        ## plus the ones that didn't but fall inside the range of a lo
        if (all.in & length(ra.all)>0){
            t3 = values(ra.all)[, tfield]==3

            if (any(t3)){
                ## save every t3 except small indel
                t3.indel = which.indel(ra.all[which(t3)])
                t3.non.indel = which(t3)[setdiff(seq_along(which(t3)), t3.indel)]
                values(ra.all)[t3.non.indel, tfield] = 2
                jmessage('All-in mode: ', length(t3.non.indel),
                         ' tier 3 junctions being included yielding ',
                         sum(values(ra.all)[, tfield]==2), ' total junctions\n')
            }
        }
        
        last.ra = ra.all[values(ra.all)[, tfield]<3]
        
        ## #    if (verbose)
        ## #    {
        ## jmessage('Starting JaBbA iterative with ', length(last.ra), ' upper tier (tier 1 and 2) junctions')
        jmessage('Starting JaBbA iterative with ', length(last.ra), ' junctions')
        ## jmessage('Will progressively add tier 3 junctions within ', rescue.window, 'bp of a loose end in prev iter')
        jmessage('Will progressively add junctions within ', rescue.window, 'bp of a loose end in prev iter')
        jmessage('Iterating for max ', reiterate, ' iterations or until convergence (i.e. no new junctions added)')
        ## jmessage('(note: adjust iterations and rescue window via reiterate= and rescue.window= parameters)')
        ##   }

        while (continue) {
            gc()

            this.iter.dir = paste(outdir, '/iteration', this.iter, sep = '')
            system(paste('mkdir -p', this.iter.dir))

            jmessage('Starting iteration ', this.iter, ' in ', this.iter.dir, ' using ', length(last.ra), ' junctions')

            this.ra.file = paste(this.iter.dir, '/junctions.rds', sep = '')
            saveRDS(last.ra, this.ra.file)

            jab = jabba_stub(
                junctions = this.ra.file,
                seg = seg,
                coverage = coverage,
                hets = hets,
                nseg = nseg,
                cfield = cfield,
                tfield = tfield,
                nudge.balanced = as.logical(nudge.balanced),
                outdir = this.iter.dir,
                mc.cores = as.numeric(mc.cores),
                max.threads = as.numeric(max.threads),
                max.mem = as.numeric(max.mem),
                edgenudge = as.numeric(edgenudge),
                tilim = as.numeric(tilim),
                strict = strict,
                name = name,
                use.gurobi = as.logical(use.gurobi),
                field = field,
                subsample = subsample,
                slack.penalty = as.numeric(slack.penalty),
                mipstart = init, 
                ploidy = as.numeric(ploidy),
                purity = as.numeric(purity),
                indel = as.logical(indel),
                overwrite = as.logical(overwrite),
                verbose = as.numeric(verbose)
            )
            gc()

            jab = readRDS(paste(this.iter.dir, '/jabba.simple.rds', sep = ''))
            jabr = readRDS(paste(this.iter.dir, '/jabba.raw.rds', sep = ''))
            le = jab$segstats[jab$segstats$loose]

            ## Annotate ra.all
            all.input = readRDS(paste0(outdir, "/junctions.all.rds"))
            all.ov = ra.overlaps(all.input, jab$junctions, pad=0, arr.ind=TRUE)
            if (ncol(all.ov)==2){
                all.ov = data.table(data.frame(all.ov))
                all.ov[, this.cn := values(jab$junctions)$cn[ra2.ix]]
                values(all.input)[, paste0("iteration", this.iter, ".cn")] =
                    all.ov[, setNames(this.cn, ra1.ix)][as.character(seq_along(all.input))]
            } else {
                values(all.input)[, paste0("iteration", this.iter, ".cn")] = NA
            }
            saveRDS(all.input, paste0(outdir, "/junctions.all.rds"))

            ## junction rescue
            ## rescues junctions that are within rescue.window bp of a loose end
            new.ra.id = union(values(jab$junctions)$id[which(values(jab$junctions)$cn>0)],## got used, stay there
                              values(ra.all)$id[which(grl.in(ra.all, le + rescue.window, some = T))]) ## near a loose ends, got another chance
            if (tfield %in% colnames(ra.all)){
                high.tier.id = values(ra.all)$id[which(as.numeric(values(ra.all)[, tfield])<3)]
                new.ra.id = union(new.ra.id, high.tier.id)
            }
            new.ra = ra.all[which(values(ra.all)$id %in% new.ra.id)]
            ## new.ra  = ra.all[union(values(last.ra)$id,
            ##                        values(ra.all)$id[grl.in(ra.all, le + rescue.window, some = T)])]
            ## new.junc.id = setdiff(new.ra.id, values(jab$junctions)$id[which(values(jab$junctions)$cn>0)])
            new.junc.id = setdiff(new.ra.id, values(jab$junctions)$id)
            ## num.new.junc = length(setdiff(values(new.ra)$id, values(last.ra)$id)==0)
            num.new.junc = length(new.junc.id)
            jcn = rep(0, nrow(jab$ab.edges))
            abix = rowSums(is.na(rbind(jab$ab.edges[, 1:2, 1])))==0
            if (any(abix)){
                jcn[abix] = jab$adj[rbind(jab$ab.edges[abix, 1:2, 1])]
            }
            num.used.junc = length(which(jcn>0))

            t3 = values(new.ra)[, tfield]==3
            jmessage('Adding ', num.new.junc,
                     ' new junctions, including ', sum(t3),
                     ' tier 3 junctions, yielding ', num.used.junc,
                     ' used junctions and ', length(new.ra), ' total junctions\n')

            if (any(t3)){
                values(new.ra)[t3, tfield] = 2
            }

            if (num.new.junc==0 | this.iter >= reiterate)
                continue = FALSE
            else {
                last.ra = new.ra
                this.iter = this.iter + 1
            }


            pp1 = readRDS(paste0(outdir, '/iteration1/karyograph.rds.ppgrid.solutions.rds'))
            purity = pp1$purity[1]
            ploidy = pp1$ploidy[1]

            seg = readRDS(paste0(outdir,'/iteration1/seg.rds')) ## why read from the first iteration??
            ## seg = jabr$segstats %Q% (strand=="+")

            if (verbose)
            {
                jmessage("Setting mipstart to previous iteration's jabba graph")
            }

            init = jab

            if (verbose)
            {
                jmessage('Using purity ', round(purity,2), ' and ploidy ', round(ploidy,2), ' across ', length(seg), ' segments used in iteration 1')
            }
        }

        this.iter.dir = paste(outdir, '/iteration', this.iter, sep = '')

        system(sprintf('cp %s/* %s', this.iter.dir, outdir))
        jmessage('Done Iterating')
    } else  {
        ## if all.in, convert all tier 3 to tier 2
        if (tfield %in% colnames(values(ra.all))){
            t3 = (values(ra.all)[, tfield] == 3)
            if (all.in & length(ra.all)>0){
                if (any(t3)){
                    ## save every t3 except small indel
                    t3.indel = which.indel(ra.all[which(t3)])
                    t3.non.indel = which(t3)[setdiff(seq_along(which(t3)), t3.indel)]
                    values(ra.all)[t3.non.indel, tfield] = 2
                }
            } else {
                ## if not all.in, only use t2 or t1
                ra.all = ra.all[setdiff(seq_along(ra.all), which(t3))]
            }
        }        
        jab = jabba_stub(
            junctions = ra,
            seg = seg,
            coverage = coverage,
            hets = hets,
            nseg = nseg,
            cfield = cfield,
            tfield = tfield,
            nudge.balanced = as.logical(nudge.balanced),
            outdir = outdir,
            mc.cores = as.numeric(mc.cores),
            max.threads = as.numeric(max.threads),
            max.mem = as.numeric(max.mem),
            edgenudge = as.numeric(edgenudge),
            tilim = as.numeric(tilim),
            strict = strict,
            name = name,
            use.gurobi = as.logical(use.gurobi),
            field = field,
            subsample = subsample,
            slack.penalty = as.numeric(slack.penalty),
            mipstart = init, 
            ploidy = as.numeric(ploidy),
            purity = as.numeric(purity),
            indel = as.logical(indel),
            overwrite = as.logical(overwrite),
            verbose = as.numeric(verbose)
        )
    }

    return(jab)
}


#' @name jabba_stub
#' @rdname internal
#' @title jabba_stub
#' @description
#' Internal function to run single iteration of JaBbA
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
#' @param junctions  GRangesList of junctions  (i.e. bp pairs with strands oriented AWAY from break) OR path to junction VCF file (BND format), dRanger txt file or rds of GRangesList
#' @param coverage  GRanges of coverage OR path to cov file, rds of GRanges or .wig / .bed file of (normalized, GC corrected) fragment density
#' @param field  field of coverage GRanges to use as fragment density signal (only relevant if coverage is GRanges rds file)
#' @param seg  optional path to existing segmentation, if missing then we will segment coverage using DNACopy with standard settings
#' @param cfield  character, junction confidence meta data field in ra
#' @param tfield  character, tier confidence meta data field in ra
#' @param outdir  out directory to dump into, default ./
#' @param nseg  optional path to normal seg file with $cn meta data field
#' @param hets  optional path to hets.file which is tab delimited text file with fields seqnames, start, end, alt.count.t, ref.count.t, alt.count.n, ref.count.n
#' @param name  prefix for sample name to be output to seg file
#' @param mc.cores  number of cores to use (default 1)
#' @param use.gurobi  logical flag whether to use gurobi vs CPLEX
#' @param nseg  path to data.frame or GRanges rds of normal seg file with coordinates and $cn data field specifying germline integer copy number
#' @param subsample  numeric between 0 and 1 specifying how much to sub-sample high confidence coverage data
#' @param tilim  timeout for jbaMIP computation (default 1200 seconds)
#' @param edgenudge  numeric hyper-parameter of how much to nudge or reward aberrant junction incorporation, default 0.1 (should be several orders of magnitude lower than average 1/sd on individual segments), a nonzero value encourages incorporation of perfectly balanced rearrangements which would be equivalently optimal with 0 copies or more copies.
#' @param slack.penalty  penalty to put on every loose.end copy, should be calibrated with respect to 1/(k*sd)^2 for each segment, i.e. that we are comfortable with junction balance constraints introducing k copy number deviation from a segments MLE copy number assignment (the assignment in the absence of junction balance constraints)
#' @param init jabba object (list) or path to .rds file containing previous jabba object which to use to initialize solution, this object needs to have the identical aberrant junctions as the current jabba object (but may have different segments and loose ends, i.e. is from a previous iteration)
#' @param overwrite  flag whether to overwrite existing output directory contents or just continue with existing files.
#' @import DNAcopy
jabba_stub = function(
                      junctions, # path to junction VCF file, dRanger txt file or rds of GRangesList of junctions (with strands oriented pointing AWAY from breakpoint)
                      coverage, # path to cov file, rds of GRanges
                      seg = NULL, # path to seg file, rds of GRanges
                      cfield = NULL, # character, junction confidence meta data field in ra
                      tfield = NULL, # character, tier confidence meta data field in ra
                      nudge.balanced = FALSE,  ## if TRUE nudge chains of balanced (or quasi balanced junctions)
                      thresh.balanced = 500, ## threshold for balanced junctions
                      outdir = './', # out directory to dump into
                      nseg = NULL, # path to normal seg file with $cn meta data field
                      hets = NULL, # path to hets.file which is tab delimited text file with fields seqnames, start, end, alt.count.t, ref.count.t, alt.count.n, ref.count.n
                      name = 'tumor', ## prefix for sample name to be output to seg file
                      mc.cores = 1, # default 1
                      max.threads = Inf,
                      max.mem = 16,
                      purity = NA,
                      ploidy = NA,
                      strict = FALSE,
                      mipstart = NULL, 
                      field = 'ratio', ## character, meta data field to use from coverage object to indicate numeric coveragendance, coverage,
                      subsample = NULL, ## numeric scalar between 0 and 1, how much to subsample coverage per segment
                      tilim = 1200, ## timeout for MIP portion
                      ## mem = 16, ## max memory for MIP portion
                      init = NULL, ## previous JaBbA object to use as a solution
                      edgenudge = 0.1, ## hyper-parameter of how much to "nudge" or reward edge use, will be combined with cfield information if provided
                      slack.penalty = 1e2, ## nll penalty for each loose end cop
                      use.gurobi = FALSE,
                      indel = TRUE,
                      overwrite = F, ## whether to overwrite existing output in outdir
                      verbose = TRUE
                      )
{
    kag.file = paste(outdir, 'karyograph.rds', sep = '/')
    hets.gr.rds.file = paste(outdir, 'hets.gr.rds', sep = '/')
    junctions.txt.file = paste(outdir, 'junctions.txt', sep = '/')
    junctions.rds.file = paste(outdir, 'junctions.rds', sep = '/')
    jabba.raw.rds.file = paste(outdir, 'jabba.raw.rds', sep = '/')
    jabba.rds.file = paste(outdir, 'jabba.rds', sep = '/')
    jabba.json.file = paste(outdir, 'jabba.json', sep = '/')
    jabba.gg.rds.file = paste(outdir, 'jabba.gg.rds', sep = '/')
    jabba.vcf.file = paste(outdir, 'jabba.vcf', sep = '/')
    jabba.cnv.vcf.file = paste(outdir, 'jabba.cnv.vcf', sep = '/')
    jabba.simple.rds.file = paste(outdir, 'jabba.simple.rds', sep = '/')
    jabba.simple.gg.rds.file = paste(outdir, 'jabba.simple.gg.rds', sep = '/')
    jabba.simple.vcf.file = paste(outdir, 'jabba.simple.vcf', sep = '/')
    jabba.simple.cnv.vcf.file = paste(outdir, 'jabba.simple.cnv.vcf', sep = '/')
    jabba.png.file = paste(outdir, 'jabba.png', sep = '/')
    jabba.simple.png.file = paste(outdir, 'jabba.simple.png', sep = '/')
    seg.tab.file = paste(outdir, 'jabba.seg.txt', sep = '/')
    seg.gr.file = paste(outdir, 'jabba.seg.rds', sep = '/')
    seg.adj.file = paste(outdir, 'jabba.adj.txt', sep = '/')
    nozzle.file = paste(outdir, 'nozzle', sep = '/')

    if (is.character(coverage))
    {
        if (!file.exists(coverage))
        {
            stop(paste('Coveraeg path ', coverage, 'does not exist'))
        }
        
        if (grepl('\\.rds$', coverage))
        {
            coverage = readRDS(coverage)
        }
        else if (grepl('(\\.txt$)|(\\.tsv$)|(\\.csv$)', coverage))
        {
            tmp = fread(coverage)
            coverage = dt2gr(tmp, seqlengths = tmp[, max(end), by = seqnames][, structure(V1, names = seqnames)])
        }
        else
        {
            jmessage('Importing seg from UCSC format')
            coverage = import.ucsc(coverage)
            field = 'score';
            coverage = gr.fix(coverage)
        }
    }
    else
        coverage = coverage

    if (!inherits(coverage, 'GRanges'))
        coverage = dt2gr(coverage)


    if (verbose)
    {
        jmessage("Read in coverage data across ", prettyNum(length(coverage), big.mark = ','), " bins and ", length(unique(seqnames(coverage))), ' chromosomes')
    }

    if (!(field %in% names(values(coverage))))
    {
        new.field = names(values(coverage))[1]
        warning(paste0('Field ', field, ' not found in coverage GRanges metadata so using ', new.field, ' instead'))
        field = new.field
    }

    seg.fn = paste0(outdir, '/seg.rds')
    if (!overwrite & file.exists(seg.fn))
    {
        if (verbose)
        {
            jmessage('Using previous segmentation found in jabba directory')
        }
        seg = readRDS(seg.fn)
    } else {
        if (is.null(seg))
        {
            if (verbose)
            {
                jmessage('No segmentation provided, so performing segmentation using CBS')
            }
            set.seed(42)
            vals = values(coverage)[, field]
            new.sl = seqlengths(coverage)
            ix = which(!is.na(vals))
            cna = DNAcopy::CNA(log(vals[ix]), as.character(seqnames(coverage))[ix], start(coverage)[ix], data.type = 'logratio')
            seg = DNAcopy::segment(DNAcopy::smooth.CNA(cna), alpha = 1e-5, verbose = FALSE)

            if (verbose)
            {
                jmessage('Segmentation finished')
            }
            seg = seg2gr(seg$out, new.sl) ## remove seqlengths that have not been segmented
            seg = gr.fix(seg, seqlengths(coverage), drop = T)

            if (verbose)
            {
                jmessage(length(seg), ' segments produced')
            }
            names(seg) = NULL
        }
        else
        {
            if (is.character(seg))
            {
                if (!file.exists(seg))
                {
                    stop(paste('Seg path ', seg, 'does not exist'))
                }
                
                if (grepl('\\.rds$', seg))
                {
                    seg = readRDS(seg)
                }
                else
                {
                    jmessage('Importing seg from UCSC format')
                    seg = import.ucsc(seg)
                    field = 'score';
                }
            }
        }
    }

    saveRDS(seg, seg.fn)

    if (!inherits(seg, 'GRanges'))
        seg = dt2gr(seg, seqlengths(coverage), drop = TRUE)

    if (!is.null(hets))
        if (!file.exists(hets))
        {
            warning(sprintf('hets file "%s" not found, ignoring hets\n', hets))
            hets = NULL
        }

    if (strict)
    {
        ends = c(gr.start(seg),
                 gr.end(gr.end(seg)))
        nog = length(ra)
        ra = ra[grl.in(ra, ends, logical, logical = FALSE) == 2]
        jmessage('Applying strict junction filtering of junctions to only those that land at segment ends')
        jmessage('Leaving ', length(ra), ' junctions from an initial set of ', nog)
    }
    
    ra = junctions

    if (!is.character(ra))
    {
        if (overwrite | !file.exists(kag.file))
            karyograph_stub(seg, coverage, ra = ra, out.file = kag.file, nseg.file = nseg, field = field, purity = purity, ploidy = ploidy, subsample = subsample, het.file = hets, verbose = verbose)
        else
            warning("Skipping over karyograph creation because file already exists and overwrite = FALSE")
    }
    else
    {
        if (grepl('rds$', ra) |
            grepl('vcf$', ra) |
            grepl('vcf\\.gz$', ra) |
            grepl('bedpe$', ra))
        {
            if (overwrite | !file.exists(kag.file))
                karyograph_stub(seg, coverage, ra.file = ra, out.file = kag.file, nseg.file = nseg, field = field, subsample = subsample, purity = purity, ploidy = ploidy, het.file = hets, mc.cores = mc.cores, verbose = verbose)
            else
                warning("Skipping over karyograph creation because file already exists and overwrite = FALSE")
        } else  {
            if (overwrite | !file.exists(kag.file))
                karyograph_stub(seg, coverage, junction.file = gsub('all.mat$', 'somatic.txt', ra), nseg.file = nseg, out.file = kag.file, field = field, purity = purity, ploidy = ploidy, subsample = subsample, mc.cores = mc.cores, het.file = hets, verbose = verbose)
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
    ab.exclude = NULL
    
    if (!is.null(tfield))
    {
        if (tfield %in% names(values(kag$junctions)))
        {
            ## default forcing tier 1
            ab.force = which(gsub('tier', '', as.character(values(kag$junctions)[, tfield]))=='1')

            ## when the switch is on:
            ## small DUP/DEL junctions in tier 2, bump them up to tier 1
            if (indel){
                tier2.ix = which(
                    gsub('tier', '', as.character(values(kag$junctions)[, tfield]))=='2')
                if (length(tier2.ix)>0){
                    ## max size hardcoded for now
                    which.like.indel = tier2.ix[
                        which.indel(kag$junctions[tier2.ix], max.size = 1e4) 
                    ]
                } else {
                    which.like.indel = numeric(0)
                }

                ab.force = union(ab.force, which.like.indel)
            }
            
            if (verbose)
            {
                jmessage('Found tier field enforcing >=1 CN at ', length(ab.force), ' junctions')
                if (exists("which.like.indel")){
                    jmessage(length(which.like.indel), ' of them are INDEL like isolated events')
                }                
            }

            ab.exclude = which(gsub('tier', '', as.character(values(kag$junctions)[, tfield]))=='3')

            if (verbose)
            {
                jmessage('Found tier field enforcing 0 CN at ', length(ab.exclude), ' junctions')
            }

        }
    } else { ## no tfield given, assume everything is tier 2
        if (indel){
            which.like.indel = which.indel(kag$junctions, max.size = 1e4) 
            ab.force =union(ab.force, which.like.indel)
        }
    }

    gc()

    juncs = kag$junctions
    bpss = grl.unlist(juncs)

    ## s.penalty = 1/slack.penalty
    ## e.penalty = edgenudge * 1.1
    ## if (e.penalty > s.penalty){
    ##     edgenudge = s.penalty * 0.9 ## don't let 
    ## }
    if (nudge.balanced) {
        balanced.jix = c()
        if (length(juncs)>0) {
            balanced.jix = chromoplexy(kag, filt.jab=F, verbose=T, junc.only=T, dist=thresh.balanced)
            dp.jix = which(gUtils::ra.duplicated(juncs, pad=1500))
            balanced.jix = setdiff(balanced.jix, dp.jix)
        }
        edgenudge = edgenudge * as.numeric(1:length(juncs) %in% balanced.jix) ## only adds edge nudge to the balanced junctions
    } else { ## nudge everything ..
        if (length(edgenudge)==1) edgenudge = rep(edgenudge, length(juncs))
        if (length(juncs)>0){   ## hot fix for preventing nudging of NA segments
            bps.cov = gr.val(bpss, coverage, val = 'ratio')
                                        #        bps.cov = bpss %$% coverage
            na.jix = unique(bps.cov$grl.ix[is.na(bps.cov$ratio)])
            if (length(na.jix)>0){
                ## if (verbose)
                ## {
                ##   jmessage("Cancel edge nudge for ", length(na.jix), " junctions.")
                ## }
                edgenudge[na.jix] = 0
            }
        }
    }

    ## some edges should be excluded:
    ## completely "dark" reference contigs
    nothing.contig = gr2dt(kag$segstats)[, .(nothing = all(is.na(mean))), by=seqnames][nothing==TRUE, seqnames]
    ## both breakpoints in NA regions
    if (length(juncs)>0){
        junc.dt = data.table(data.frame(values(juncs)))
        junc.dt = cbind(junc.dt,
                        data.table(matrix(kag$ab.edges[,,1],
                                          nrow=nrow(kag$ab.edges),
                                          dimnames=dimnames(kag$ab.edges)[1:2])))
        junc.dt[, ":="(mean.a=kag$segstats$mean[from],
                       mean.b=kag$segstats$mean[to],
                       chr.a = as.character(seqnames(kag$segstats[from])),
                       chr.b = as.character(seqnames(kag$segstats[to])))]
        junc.dt[, both.na := is.na(mean.a) & is.na(mean.b)]
        both.na.ix = junc.dt[, which(both.na==TRUE)] ## both breakpoint in NA

        no.man.land = junc.dt[, which(chr.a %in% nothing.contig | chr.b %in% nothing.contig)]
        ## either breakpoint in a contig that's completely NA
        ## excluding those whose both bp in NA regions or mapped to completely NA contigs
        ab.exclude = union(ab.exclude, union(both.na.ix, no.man.land))
        ab.force = setdiff(ab.force, ab.exclude)
        edgenudge[ab.exclude] = 0
        ## furthermore, some extra edges should not be nudged
        either.na.ix = junc.dt[, which(both.na==FALSE & (is.na(mean.a) | is.na(mean.b)))]
        edgenudge[either.na.ix] = 0

        if (verbose){
            jmessage("Excluding ", length(both.na.ix), " aberrant junctions whose both breakpoints are in NA coverage regions")
            jmessage("Cancel nudge for ", length(either.na.ix), " aberrant junctions where one of the 2 breakpoint is in NA coverage regions")
        }
    }
    
    if (verbose){
        jmessage("In sum, we are forcing ", length(ab.force),
                 " junctions, excluding ", length(ab.exclude),
                 " junctions, and nudging ", sum(edgenudge>0), " junctions")
    }

    saveRDS(ab.exclude, paste0(outdir, "/ab.exclude.rds"))
    saveRDS(ab.force, paste0(outdir, "/ab.force.rds"))
    saveRDS(edgenudge, paste0(outdir, "/edge.nudge.rds"))

    if (!is.null(init))
    {
        if (is.character(init))
            init = readRDS(init)
    }

    if (overwrite | !file.exists(jabba.raw.rds.file))
    { 
        ramip_stub(kag.file,
                   jabba.raw.rds.file,
                   mc.cores = mc.cores,
                   max.threads = max.threads,
                   ## mem = mem,
                   mem = max.mem,
                   tilim = tilim,
                   edge.nudge = edgenudge,
                   use.gurobi = use.gurobi,
                   ab.force = ab.force,
                   ab.exclude = ab.exclude,
                   init = init,
                   verbose = verbose,
                   purity.min = purity,
                   mipstart = mipstart,
                   purity.max = purity,
                   ploidy.min = ploidy,
                   ploidy.max = ploidy,
                   slack.prior = 1/slack.penalty)
    }

    kag = readRDS(kag.file)
    jab = readRDS(jabba.raw.rds.file)

    if (overwrite | !file.exists(jabba.rds.file))
    {
        jabd = JaBbA.digest(jab, kag)
    }
    else
    {
        jabd = readRDS(jabba.rds.file)
    }

    jabd$purity = jab$purity
    jabd$ploidy = jab$ploidy

    if (overwrite | !file.exists(jabba.simple.rds.file))
    {
        if (verbose)
        {
            jmessage('simplifying segments in JaBbA graph but keeping all edges (including copy 0), dumping to jabba.rds')
        }

        jabd.simple = JaBbA.digest(jab, kag, keep.all = F) ## simplified
    }
    else {
        if (verbose)
        {
            jmessage('simplifying segments in JaBbA graph but removing all copy 0 aberrant edges, dumping to jabba.simple.rds')
        }
        jabd.simple = readRDS(jabba.simple.rds.file)
    }

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

    if (verbose)
    {
        jmessage('Checking for hets')
    }

    if (!is.null(hets))
        if (file.exists(hets.gr.rds.file))
            tryCatch(
            {
                jmessage('Loading hets')
                hets.gr = readRDS(hets.gr.rds.file)
                jmessage('Computing alleles for jabd ')
                jabd = c(jabd, jabba.alleles(jabd, hets.gr, verbose = TRUE, uncoupled=TRUE)[c('asegstats', 'aadj', 'agtrack')])
                jmessage('Computing alleles for jabd simple ')
                jabd.simple = c(jabd.simple, jabba.alleles(jabd.simple, hets.gr, verbose = TRUE, uncoupled=TRUE)[c('asegstats', 'aadj', 'agtrack')])
                jmessage('Done computing alleles')
            },
            error = function(e) print("Jabba allelic generation failed"))

    jab$segstats = gr.fix(jab$segstats)
    jabd$segstats = gr.fix(jabd$segstats)
    jabd.simple$segstats = gr.fix(jabd.simple$segstats)

    if (overwrite | !file.exists(jabba.simple.rds.file))
    {
        saveRDS(jabd$segstats, seg.gr.file)
        saveRDS(jab, jabba.raw.rds.file)
        saveRDS(jabd, jabba.rds.file)
        saveRDS(jabd.simple, jabba.simple.rds.file)

        jab.gg = gGnome::gGraph$new(jab = jabd)
        jab.simple.gg = gGnome::gGraph$new(jab = jabd.simple)

        saveRDS(jab.simple.gg, jabba.simple.gg.rds.file)
        saveRDS(jab.gg, jabba.gg.rds.file)
    }

    tryCatch(
    {
        jabba2vcf(jabd, jabba.vcf.file);
        jabba2vcf(jabd, jabba.cnv.vcf.file, cnv = TRUE)
        jabba2vcf(jabd.simple, jabba.simple.vcf.file)
        jabba2vcf(jabd.simple, jabba.simple.cnv.vcf.file, cnv = TRUE)
    }, error = function(e) print("Jabba VCF generation failed"))

    if (nrow(jabd$edges)>0){
        seg.adj = cbind(data.frame(sample = rep(name, nrow(jabd$edges))), jabd$edges[, c('from', 'to', 'cn', 'type')])
        write.tab(seg.adj, seg.adj.file)
    }

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

    tmp.cov = sample(coverage, pmin(length(coverage), 5e5))
    tmp.cov = gr.fix(tmp.cov, jabd$segstats)

    y1 = pmax(5, max(jabd$segstats$cn)*1.1)
    jabd$gtrack$y1 = y1
    jabd.simple$gtrack$y1 = y1

    td.cov = gTrack(tmp.cov, y.field = field, col = alpha('black', 0.2), name = 'Cov', y1 = (y1 + jab$gamma)/jab$beta)

    if (verbose)
    {
        jmessage('Generating figures')
    }

    if (overwrite | !file.exists(jabba.png.file))
    {
        if (is.character(tryCatch(png(jabba.png.file, width = 2000, height = 1000), error = function(e) 'bla')))
            pdf(gsub('png$', 'pdf', jabba.png.file), width = 10, height = 10)

        jun = jabd$junctions
        values(jun)$col = ifelse(values(jun)$cn>0, 'red', alpha('gray', 0.2))

        if (is.null(jabd$agtrack)){
            plotted = tryCatch(plot(c(td.cov, jabd$gtrack), links = jun), error = function(e) return(NULL))
        } else {
            plotted = tryCatch(plot(c(jabd$agtrack, td.cov, jabd$gtrack), links = jun), error = function(e) return(NULL))
        }
        
        if (is.null(plotted)){
            if (verbose){
                jmessage("Something wrong with plotting JaBbA results. Please try it later.")
            }
        }

        dev.off()
    }

    if (overwrite | !file.exists(jabba.simple.png.file))
    {

        jun = jabd.simple$junctions
        values(jun)$col = ifelse(values(jun)$cn>0, 'red', alpha('gray', 0.2))

        if (is.character(tryCatch(png(jabba.simple.png.file, width = 2000, height = 1000), error = function(e) 'bla')))
            pdf(gsub("png$", "pdf", jabba.simple.png.file), width = 20, height = 10)

        ## if (is.null(jabd.simple$agtrack))
        ##     plot(c(td.cov, jabd.simple$gtrack), links = jun)
        ## else
        ##     plot(c(jabd.simple$agtrack, td.cov, jabd.simple$gtrack), links = jun)

        if (is.null(jabd$agtrack)){
            plotted = tryCatch(plot(c(td.cov, jabd.simple$gtrack), links = jun),
                               error = function(e) return(NULL))
        } else {
            plotted = tryCatch(plot(c(jabd.simple$agtrack, td.cov, jabd.simple$gtrack), links = jun),
                               error = function(e) return(NULL))
        }
        
        if (is.null(plotted)){
            if (verbose){
                jmessage("Something wrong with plotting JaBbA simplified results. Please try it later.")
            }
        }
        
        dev.off()
    }

    jmessage('Done .. job output in: ', normalizePath(outdir))

    return(readRDS(jabba.simple.gg.rds.file))
}


###############################
## NOT EXPORTED AND INTERNAL TO JABBA
##
##############################

#' karyograph
#' @name karyograph
#' @rdname internal
karyograph_stub = function(seg.file, ## path to rds file of initial genome partition (values on segments are ignored)
                           cov.file, ## path to rds file GRanges of coverage with meta data field "field"
                           nseg.file = NULL, ## rds file of GRanges providing integer copy numbers for normal segments in the genome
                           het.file = NULL,
                           ra = NULL,
                           junction.file = NULL,
                           out.file,
                           use.ppurple = TRUE,
                           ra.file = NULL,
                           verbose = FALSE,
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
        if (!is.null(junction.file))
        {
            these.junctions = read.delim(junction.file, strings = F)

            if (ncol(these.junctions)<=1) ## wrong separator
                these.junctions = read.delim(junction.file, sep = ',', strings = F)

            if (!is.null(these.junctions$strand1) & !is.null(these.junctions$strand2))
            {
                ## looks like snowman input flip breaks so that they are pointing away from junction
                these.junctions$str1 = ifelse(these.junctions$strand1 == '+', '-', '+')
                these.junctions$str2 = ifelse(these.junctions$strand2 == '+', '-', '+')
            }

            these.junctions$chr1 = gsub('23', 'X', gsub('24', 'Y', these.junctions$chr1))
            these.junctions$chr2 = gsub('23', 'X', gsub('24', 'Y', these.junctions$chr2))
            this.ra = read.junctions(these.junctions, seqlengths = hg_seqlengths())
        }
        else if (grepl('(\\.bedpe)|(\\.vcf$)|(\\.vcf\\.gz$)', ra.file))
        {
            tmp.ra = read.junctions(ra.file, seqlengths = hg_seqlengths(), get.loose = T)
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
        {
            if (file.exists(nseg.file))
            {
                if (grepl('\\.rds$', nseg.file, ignore.case = TRUE))
                {
                    nseg = readRDS(nseg.file)
                }
                else
                {
                    nseg = dt2gr(fread(nseg.file))
                }
            } else {
                stop('Did not find nseg file!')
            }

        }
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
        if (verbose)
        {
            jmessage('Adding loose ends from vcf file to seg file')
        }
        this.seg = grbind(this.seg, gr.fix(loose.ends, sl, drop = T))
    }

    this.ra = gr.fix(this.ra, sl, drop = T)

    this.kag = karyograph(this.ra, this.seg)

    if (is.null(nseg.file))
        this.kag$segstats$ncn = 2

    hets.gr = NULL

    if (!is.null(het.file))
    {
        if (grepl(".rds$", het.file)){
            hets = readRDS(het.file)
        } else {
            hets = fread(het.file)
        }
        
        if (verbose)
        {
            jmessage('loaded hets')
        }

        if (inherits(hets, "data.frame")){
            if (!is.null(hets$alt.count.n) & !is.null(hets$ref.count.n))
                ## old format, apply het filter ourselves
            {
                hets$ref.frac.n = hets$alt.count.n / (hets$alt.count.n + hets$ref.count.n)
                ##      hets.gr = dt2gr(hets[pmin(ref.frac.n, 1-ref.frac.n) > 0.2 & (ref.count.n + alt.count.n)>20, ])
                hets.gr = dt2gr(hets[pmin(ref.frac.n, 1-ref.frac.n) > 0.2 & (ref.count.n + alt.count.n)>=2, ])
                hets.gr$alt = hets.gr$alt.count.t
                hets.gr$ref = hets.gr$ref.count.t
            }
            else ## new, standard format, with $alt and $ref field
            {
                hets.gr = dt2gr(hets)
                if (all(c("alt", "ref") %in% colnames(hets))){
                    jmessage("Valid hets already")
                    ## hets.gr$alt.count.t = hets.gr$alt
                    ## hets.gr$ref.count.t = hets.gr$ref
                } else if (all(c("alt.count.t", "ref.count.t") %in% colnames(hets))){
                    hets.gr$alt = hets.gr$alt.count.t
                    hets.gr$ref = hets.gr$ref.count.t
                    hets.gr = hets.gr %Q% (!is.na(alt)) %Q% (!is.na(ref))
                } else {
                    jmessage("hets is not in valid format, ignore")
                    hets.gr = NULL
                }            
            }
        } else if (inherits(hets, "GRanges")){
            if (all(c("alt.count.t", "ref.count.t") %in% colnames(values(hets)))){
                hets.gr = hets
                hets.gr$alt = hets$alt.count.t
                hets.gr$ref = hets$ref.count.t
                hets.gr = hets.gr %Q% (!is.na(alt)) %Q% (!is.na(ref))
            }
        } else {
            jmessage("hets is neither data.table nor GRanges, ignore.")            
        }

        if (!is.null(hets.gr)){
            ## save hets object for later
            saveRDS(hets.gr, paste(dirname(out.file), 'hets.gr.rds', sep = '/'))
        }        
    }

    if (length(hets.gr)>0){
        ## pretend we don't have hets at all
        this.kag$segstats = segstats(this.kag$tile, this.cov, field = field, prior_weight = 1, max.chunk = max.chunk, subsample = subsample, asignal = hets.gr, afield = c('ref', 'alt'), mc.cores = mc.cores)
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
                                        #    jmessage('!!!!!!!!!!! cov.thresh for fix.sd is', cov.thresh, '\n')
    fix.sd  = width(this.kag$segstats)<(3*cov.thresh)
                                        #    this.kag$segstats$mean[make.na] = NA
    this.kag$segstats$sd[fix.sd] = sqrt(this.kag$segstats$mean[fix.sd])

                                        #      if (is.character(tryCatch(png(paste(out.file, '.ppgrid.png', sep = ''), height = 500, width = 500), error = function(e) 'bla')))
    ss.tmp = this.kag$segstats[width(this.kag$segstats)>1e4, ] ## don't use ultra short segments
    pdf(paste(out.file, '.ppgrid.pdf', sep = ''), height = 10, width = 10)

    purity = as.numeric(purity)
    ploidy = as.numeric(ploidy)
    if (!is.na(purity) & !is.na(ploidy)) ## purity and ploidy are completely set
    {
        pp = data.table(purity = purity, ploidy = ploidy)
    }
    else if (use.ppurple)
    {

        if (is.na(purity))
        {
            purity = seq(0, 1, 0.1)
        }

        if (is.na(ploidy))
        {
            ploidy = seq(1, 6, 0.2)
        }

        this.cov$y = values(this.cov)[, field]

        if (verbose)
        {
            jmessage('Computing purity and ploidy with Ppurple')
        }

        max.chunk = 1e3
        numchunks = ceiling(length(ss.tmp)/max.chunk)
        if (numchunks>length(purity)*length(ploidy)){
            pp = ppurple(cov = this.cov, hets = hets.gr, seg = ss.tmp,
                         purities = purity, ploidies = ploidy,
                         verbose = verbose,
                         mc.cores = mc.cores,
                         ## numchunks = numchunks,
                         ignore.sex = TRUE)
        } else {
            pp = ppurple(cov = this.cov, hets = hets.gr, seg = ss.tmp,
                         purities = purity, ploidies = ploidy,
                         verbose = verbose,
                         mc.cores = mc.cores,
                         numchunks = numchunks,
                         ignore.sex = TRUE)
        }        
    }
    ## else
    ##   {
    ##     if (!is.null(het.file))
    ##       {
    ##         pp = ppgrid(ss.tmp, verbose = verbose, plot = F, mc.cores = mc.cores,
    ##                     purity.min = ifelse(is.na(purity), 0, purity), purity.max = ifelse(is.na(purity),1, purity),
    ##                     ploidy.min = ifelse(is.na(ploidy), 1.2, ploidy), ploidy.max = ifelse(is.na(ploidy), 6, ploidy), allelic = TRUE)
    ##       }
    ##     else
    ##       {
    ##         pp = ppgrid(ss.tmp, verbose = verbose, plot = F, mc.cores = mc.cores,
    ##                     purity.min = ifelse(is.na(purity), 0, purity), purity.max = ifelse(is.na(purity),1, purity),
    ##                     ploidy.min = ifelse(is.na(ploidy), 1.2, ploidy), ploidy.max = ifelse(is.na(ploidy), 6, ploidy), allelic = FALSE)
    ##       }
    ##   }

    mu = this.kag$segstats$mean
    mu[is.infinite(mu)] = NA
    w = as.numeric(width(this.kag$segstats))
    w[is.na(mu)] = NA
    sw = sum(w, na.rm = T)
    ncn = this.kag$segstats$ncn
    ploidy_normal = sum(w * ncn, na.rm = T) / sw  ## this will be = 2 if ncn is trivially 2
    mutl = sum(mu * w, na.rm = T)
    pp$beta = ((1-pp$purity)*ploidy_normal + pp$purity*pp$ploidy) * sw / (pp$purity * mutl)
    pp$gamma = 2*(1-pp$purity)/pp$purity

    saveRDS(pp, paste(out.file, '.ppgrid.solutions.rds', sep = '')) ## save alternate solutions

    this.kag$purity = pp[1,]$purity
    this.kag$ploidy = pp[1,]$ploidy
    this.kag$beta = pp[1,]$beta
    this.kag$gamma = pp[1,]$gamma
    this.kag$segstats$cn = rel2abs(this.kag$segstats, purity = this.kag$purity, ploidy = this.kag$ploidy, field = 'mean') ## cn is the copy number b4 rounding
    saveRDS(this.kag, out.file)

    if (is.character(tryCatch(png(paste(out.file, '.ppfit.png', sep = ''), height = 1000, width = 1000), error = function(e) 'bla')))
        pdf(paste(out.file, '.ppfit.pdf', sep = ''), height = 10, width = 10)
    tmp.kag = this.kag

    if (length(tmp.kag$segstats)<10)
        warning('number of segments used for purity ploidy extremely low .. check coverage data')
    .plot_ppfit(tmp.kag)
    dev.off()

    if (verbose)
    {
        jmessage('Built gGraph with ', length(this.kag$tile), ' nodes, ', sum(this.kag$adj!=0), ' edges, purity ', round(this.kag$purity,2), ', and ploidy ', round(this.kag$ploidy,2))
    }

    y1 = 10

    if (is.character(tryCatch(png(paste(out.file, '.inputdata.png', sep = ''), height = 1000, width = 1000), error = function(e) 'bla'))){
        pdf(paste(out.file, '.inputdata.pdf', sep = ''), height = 10, width = 10)
        plot(c(gTrack(gr.fix(sample(this.cov, pmin(length(this.cov), 5e4)), this.kag$segstats), y.field = 'ratio', col = alpha('black', 0.3)),
               gTrack(this.kag$segstats, y.field = 'mean', angle = 0, col = 'gray10', border = alpha('black', 0.2))), links = this.kag$junctions, y1 = y1)
        dev.off()
    }
}

## diagnostic function used by karyograph_stub
#' @name .plot_ppfit
#' @rdname internal
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

#' run ramip
#' @name .ramip_stub
#' @rdname internal
ramip_stub = function(kag.file,
                      out.file,
                      mc.cores = 1,
                      max.threads = Inf,
                      mem = 16,
                      tilim = 1200,
                      slack.prior = 0.001,
                      gamma = NA,
                      beta = NA,
                      customparams = T,
                      purity.min = NA, purity.max = NA,
                      ploidy.min = NA, ploidy.max = NA,
                      init = NULL,
                      mipstart = NULL,
                      use.gurobi = FALSE,
                      verbose = FALSE,
                      edge.nudge = 0,  ## can be scalar (equal nudge to all ab junctions) or vector of length readRDS(kag.file)$junctions
                      ab.force = NULL, ## indices of aberrant junctions to force include into the solution
                      ab.exclude = NULL ## indices of aberrant junctions to force exclude from the solution
                      )
{
    outdir = normalizePath(dirname(kag.file))
    this.kag = readRDS(kag.file)

    ## if (is.null(this.kag$gamma) | is.null(this.kag$beta))
    ## {
    ##     pp = ppgrid(this.kag$segstats, verbose = verbose, plot = T, purity.min = purity.min, purity.max = purity.max, ploidy.min = ploidy.min, ploidy.max = ploidy.max)
    ##     this.kag$beta = pp[1,]$beta
    ##     this.kag$gamma = pp[1,]$gamma
    ## }

    if (!is.na(gamma))
    {
        if (verbose)
        {
            jmessage(sprintf('Overriding gamma with %s\n', gamma))
        }
        this.kag$gamma = gamma
    }

    if (!is.na(beta))
    {
        if (verbose)
        {
            jmessage(sprintf('Overriding beta with %s\n', beta))
        }
        this.kag$beta = beta
    }

    if (customparams)
    {
        MAX.THREADS = Sys.getenv("LSB_DJOB_NUMPROC")
        if (nchar(MAX.THREADS) == 0)
            MAX.THREADS = Inf
        else
            MAX.THREADS = as.numeric(MAX.THREADS)
        max.threads = min(max.threads, MAX.THREADS)
        if (is.infinite(max.threads))
            max.threads = 0

        param.file = paste(out.file, '.prm', sep = '')
        .cplex_customparams(param.file, max.threads, treememlim = mem * 1e3)

        Sys.setenv(ILOG_CPLEX_PARAMETER_FILE = normalizePath(param.file))
        if (verbose)
        {
            jmessage('Creating ILOG CPLEX PARAMETER FILE in ', Sys.getenv('ILOG_CPLEX_PARAMETER_FILE'))
        }
    }

    adj.nudge = this.kag$adj*0;
    adj.nudge[this.kag$ab.edges[,1:2,1]] = 1*edge.nudge ## if edge.nudge is length ab.edges, then corresponding edges will be nudged

    adj.lb = NULL
    if (!is.null(ab.force))
    {
        if (verbose)
        {
            if (length(ab.force)>0)
            {
                jmessage(paste('Enforcing lower bounds on aberrant junctions:', paste(ab.force, collapse = ',')))
            }
        }
        adj.lb = this.kag$adj*0
        adj.lb[rbind(this.kag$ab.edges[ab.force, ,1])[, 1:2, drop = FALSE]] = 1
        adj.lb[rbind(this.kag$ab.edges[ab.force, ,2])[, 1:2, drop = FALSE]] = 1
        saveRDS(adj.lb, paste0(outdir, "/adj.lb.rds"))
    }

    if (!is.null(ab.exclude))
    {
        if (verbose)
        {
            if (length(ab.exclude)>0)
            {
                jmessage(paste('Excluding aberrant junctions:', paste(ab.exclude, collapse = ',')))
            }
        }
        adj.ub = this.kag$adj*0
        adj.ub[rbind(this.kag$ab.edges[ab.exclude, ,1])[, 1:2, drop = FALSE]] = -1
        adj.ub[rbind(this.kag$ab.edges[ab.exclude, ,2])[, 1:2, drop = FALSE]] = -1
        saveRDS(adj.ub, paste0(outdir, "/adj.ub.rds"))
    } else {
        adj.ub = NULL
    }

    ## if mipstart is not given, construct the naive solution
    ## if mipstart is given (a gGnome or JaBBA) object
    ## here we create an mipstart "adj" matrix for the new graph
    ## by looking up the junctions in the current graph in the old object
    if (is.null(mipstart)){
        if (file.exists("mipstart.gg.rds")){
            jmessage("Using existing mipstart in the current directory")
            mipstart = readRDS("mipstart.gg.rds")
        } else {
            jmessage("Adjusting the kag (naive solution) as mipstart (initial solution).")
            mipstart = gGnome::gread(kag.file)
            segs = mipstart$segstats
            segs$cn = pmax(round(segs$cn), 0) ## negative given 0
            segs$cn[is.na(segs$cn)] = 0 ## NA given 0

            es = mipstart$edges
            es[, ":="(from.cn = segs$cn[from], to.cn = segs$cn[to])]
            es = es[type != "loose"]
            es[type=="reference", cn := 0]
            es[type=="aberrant", cn := 0]
            
            if (!is.null(ab.force)){
                jmatch = data.table(ra.overlaps(mipstart$junctions, this.kag$junctions))
                ab.force.mipstart = jmatch[, setNames(ra1.ix, ra2.ix)][as.character(ab.force)]
                jdt = data.table(data.frame(values(mipstart$junctions)))
                jdt[ab.force.mipstart, force.in := TRUE]
                jdt[, eid := paste(from1, to1)]
                jdt[, reid := paste(from2, to2)]
                es[eid %in% jdt[force.in==TRUE, c(eid, reid)] |
                   eid %in% jdt[force.in==TRUE, c(eid, reid)],
                   cn := 1]
            }

            es[, from.remain := from.cn - sum(cn), by=from]
            es[, to.remain := to.cn - sum(cn), by=to]

            if (any(es[, from.remain<0 | to.remain<0])){
                segs$cn[es[from.remain<0, from]] = es[from.remain<0, abs(from.remain)]
                segs$cn[es[to.remain<0, to]] = es[to.remain<0, abs(to.remain)]
                ## re-annotate the es
                es[, ":="(from.cn = segs$cn[from], to.cn = segs$cn[to])]
                es[, from.remain := from.cn - sum(cn), by=from]
                es[, to.remain := to.cn - sum(cn), by=to]
            }

            es[type=="reference", cn := pmin(from.remain, to.remain)]

            mipstart = gGraph$new(segs = segs, es = es)$make.balance()
            mipstart = gGraph$new(segs = segs, es = es)$fillin()
            saveRDS(mipstart, "mipstart.gg.rds")
        }        
    }
    
    if (!is.null(mipstart)) 
    {
        if (verbose)
            jmessage('Applying mipstarts from previous jabba solution')
        ## for mipstart graph
        mgre = suppressWarnings(gr.end(mipstart$segstats,1, ignore.strand = FALSE))
        mgrs = suppressWarnings(gr.start(mipstart$segstats,1, ignore.strand = FALSE))
        mgend = gr.string(mgre)
        mgstart = gr.string(mgrs)
        mij = Matrix::which(mipstart$adj!=0, arr.ind = TRUE)
        mijs = data.table(mstr = paste(mgend[mij[,1]], mgstart[mij[,2]]),
                          cn = mipstart$adj[mij])
        
        setkey(mijs, mstr)

        ## for new (i.e. this) graph
        gre = suppressWarnings(gr.end(this.kag$segstats,1, ignore.strand = FALSE))
        grs = suppressWarnings(gr.start(this.kag$segstats,1, ignore.strand = FALSE))
        gend = gr.string(gre)
        gstart = gr.string(grs)
        ij = Matrix::which(this.kag$adj!=0, arr.ind = TRUE)
        ijs = data.table(i = ij[,1], j = ij[,2], gstr = paste(gend[ij[,1]], gstart[ij[,2]]), mipstart = as.numeric(NA))
        ijs$is.ref = GenomicRanges::shift(gre[ij[,1]], ifelse(as.logical(strand(gre)[ij[,1]]=="+"), 1, -1)) == grs[ij[,2]]
        setkey(ijs, gstr)

        mcn = mipstart$segstats$cn[gr.match(this.kag$segstats, mipstart$segstats)]
        ## for all remaining ref edges just pick the cn as the cn of the segment that it matches in mipstart
        ijs[gstr %in% mijs$mstr, mipstart := mijs[.(gstr), cn]]
        ijs[is.na(mipstart) & is.ref==TRUE, mipstart := pmin(mcn[i], mcn[j])] 
        ijs[is.na(mipstart), mipstart := 0] ## all remaining are 0

        mipstart = sparseMatrix(ijs$i, ijs$j, x = ijs$mipstart,
                                dims = dim(this.kag$adj))
    }

    ## too hardcoded!!
    ## nothing.contig = gr2dt(this.kag$segstats)[
    ##   , .(nothing = all(is.na(mean))), by=seqnames][
    ##     nothing==TRUE, seqnames]
    ## if (verbose){
    ##     jmessage("Finally ignoring ",
    ##              length(nothing.contig),
    ##              " contigs in the reference genome completely not covered.")
    ## }

    ra.sol = jbaMIP(this.kag$adj,
                    this.kag$segstats,
                    beta.guess = this.kag$beta,
                    gamma.guess = this.kag$gamma,
                    tilim = tilim,
                    slack.prior = slack.prior,
                    cn.prior = NA,## why not used?
                    mipemphasis = 0,
                    ignore.cons = T,
                    mipstart = mipstart, ## make mipstart if not provided
                    adj.lb = adj.lb,
                    adj.ub = adj.ub,
                    use.gurobi = use.gurobi,
                    mc.cores = mc.cores,
                    adj.nudge = adj.nudge,
                    cn.ub = rep(500, length(this.kag$segstats)),
                    ## cn.ub = ifelse(as.character(seqnames(this.kag$segstats)) %in% nothing.contig,
                    ##                0, 500),
                    verbose = verbose)
    saveRDS(ra.sol, out.file)

    ## ## report the optimization status
    opt.report = do.call(`rbind`,
                         lapply(seq_along(ra.sol$sols),
                                function(cl){
                                    x = ra.sol$sols[[cl]]
                                    if (inherits(x$nll.cn, "Matrix") |
                                        inherits(x$nll.cn, "matrix")){
                                        nll.cn = x$nll.cn[1, 1]
                                    } else {
                                        nll.cn = NA
                                    }
                                    data.table(cl = cl,
                                               obj = ifelse(is.null(x$obj), NA, x$obj),
                                               status = ifelse(is.null(x$status), NA, x$status),
                                               nll.cn = nll.cn,
                                               nll.opt = x$nll.opt,
                                               gap.cn = x$gap.cn)
                                }))
    saveRDS(opt.report, paste0(outdir, "/opt.report.rds"))
    
    if (customparams)
    {
        system(paste('rm', param.file))
        Sys.setenv(ILOG_CPLEX_PARAMETER_FILE='')
    }
}

##############################
#' @name .ramip_stub
#' @rdname internal
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
###########################################
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
                                        #    pc.na = sapply(vall, function(x) sum(is.na(x))/length(x))

                                        #   mu[which(pc.na>na.thresh)] = NA
                                        #   sd[which(pc.na>na.thresh)] = NA

        target$mean = NA;
        target$sd = sqrt(prior_beta / (prior_alpha + 1)); ## These are prior parameters

        mu = sapply(vall, mean, na.rm = TRUE)
        ix = !is.na(mu)
        target$mean[ix] = mu[ix]

        ## ## loess var estimation
        ## ## i.e. we fit loess function to map segment mean to variance across the sample
        ## ## the assumption is that such a function exists 
        ## var = sapply(vall, var, na.rm = TRUE)       
        ## ##        target$nbins = sapply(map, length)[as.character(abs(as.numeric(names(target))))]
        ## target$nbins = sapply(vall, function(x) sum(!is.na(x)))[as.character(abs(as.numeric(names(target))))]
        ## loe = loess(var ~ mu, weights = target$nbins)
        ## target$var = pmax(predict(loe, target$mean), min(var, na.rm = TRUE))
        
        ## loess var estimation
        ## i.e. we fit loess function to map segment mean to variance across the sample
        ## the assumption is that such a function exists 
        sample.var = sapply(vall, var, na.rm = TRUE)        ## computing sample variance for each segment
        ##        target$nbins = sapply(map, length)[as.character(abs(as.numeric(names(target))))]
        target$nbins = sapply(vall, function(x) sum(!is.na(x)))[as.character(abs(as.numeric(names(target))))]
        target$nbins.tot = sapply(map, length)[as.character(abs(as.numeric(names(target))))]
        target$nbins.nafrac = 1-target$nbins/target$nbins.tot

        tmp = data.table(var = sample.var, mean = target$mean, nbins = target$nbins, na.frac = target$nbins.nafrac)
        loe = tmp[nbins>2 & na.frac<0.5, loess(var ~ mean, weights = nbins)]

        ## inferring segment specific variance using loess fit of mean to sample variance across dataset
        target$var = pmax(predict(loe, target$mean), min(sample.var, na.rm = TRUE))

        ## clean up NA values which are below or above the domain of the loess function which maps mean -> variance
        ## basically assign all means below the left domain bound of the function the variance of the left domain bound
        ## and analogously for all means above the right domain bound
        na.var = is.na(target$var)
        rrm = range(target$mean[!na.var])
        rrv = predict(loe, rrm)
        target$var[target$mean<=rrm[1]] = rrv[1]
        target$var[target$mean>=rrm[2]] = rrv[2]    

        ## computing sd / sem for each target 
        target$sd = sqrt((2*target$var)/target$nbins)

        ## final clean up
        good.bin = signal[which(!is.na(values(signal)[, field]) & !is.infinite(values(signal)[, field]))]
        ## target$good.prop = (target+5e5) %O% good.bin
        target$good.prop = (target+1e5) %O% good.bin
        ## table(target$good.prop < 0.85)
        ## t.td = gTrack(target %Q% (good.prop<0.75))
        ## c.td = gTrack(signal, y.field = field, circles=T, lwd.border=0.2)
        ## ppdf(plot(c(c.td, t.td), streduce((target %Q% (good.prop < 0.75))+1e6)[1:20]), width=25)
        ## bad.nodes = which(target$good.prop < 0.75)
        bad.nodes = which(target$good.prop < 0.9)
        target$raw.mean = target$mean
        target$raw.sd = target$sd
        target$mean[bad.nodes] = NA
        target$sd[bad.nodes] = NA
        target$bad.nodes = seq_along(target) %in% bad.nodes
        jmessage("Definining coverage good quality nodes as 90% bases covered by non-NA and non-Inf values in +/-100KB region")
        jmessage("Hard setting ", sum(width(target[bad.nodes]))/1e6, " Mb of the genome to NA that didn't pass our quality threshold")      
    }

    return(target)
}

#' @name jmessage
#' @rdname internal
jmessage = function(..., pre = 'JaBbA')
    message(pre, ' ', paste0(as.character(Sys.time()), ': '), ...)

#' @name jbaMIP
#' @rdname internal
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
############################################
jbaMIP = function(adj, # binary n x n adjacency matrix ($adj output of karyograph)
                  segstats, # n x 1 GRanges object with "mean" and "sd" value fields
                  mipstart = NULL, ## sparse adjacency matrix of mipstarts (0 = NA, 0+eps + 0, k>=1 = k)
########### optional args
                  beta = NA, # beta guess (optional)
                  gamma = NA, # gamma guess (optional)
                  field.ncn = 'ncn', # will use this field to take into account normal copy number in transformation of relative to integer copy number
                  tilim = 20, mipemphasis = 0, epgap = 0.01, # MIP params
                  ploidy.min = 0.1, # ploidy bounds (can be generous)
                  ploidy.max = 20,
                  ploidy.normal = NULL, ## usually inferred from ncn field but can be entered for subgraph analysis
                  ## #  purity.guess = NA,
                  ## #  ploidy.guess = NA,
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
                  adj.ub = NULL,
                  adj.nudge = 0*adj, # linear objective function coefficients for edges (only which(adj!=0) components considered)
                  na.node.nudge = TRUE,
                  ecn.out.ub = rep(NA, length(segstats)), ## upper bound for cn of edges leaving nodes
                  ecn.in.ub = rep(NA, length(segstats)),  ## upper bound for cn of edges entering nodes
                  use.gurobi = FALSE, # otherwise will use cplex
                  nsolutions = 1,
                  verbose = F,
                  debug = F,
                  mc.cores = 1, ## only matters if partition = T
                  ignore.edge = FALSE, ignore.cons = TRUE, edge.slack = TRUE,
                  slack.prior = 1,
                  ... # passed to optimizer
                  )
{
    if (length(segstats) != nrow(adj))
        stop('length(segstats) !=  nrow(adj)')

    if (is.null(adj.lb))
        adj.lb = 0*adj

    ## save the naive solutions
    segstats$kag.cn = segstats$cn
    
    ## wrapper that calls jbaMIP recursively on subgraphs after "fixing"
    if (partition & !is.na(gamma.guess) & !is.na(beta.guess))
    {

        ##      m = segstats$mean*beta.guess - gamma.guess

        ## transform means from data space into copy number space
        m = rel2abs(segstats, gamma = gamma.guess, beta = beta.guess, field = 'mean', field.ncn = field.ncn)

        ## transform sds from data space into copy number space (only need to multiply by beta)
        segstats$sd = segstats$sd * beta.guess

        cnmle = round(m) ## MLE estimate for CN
        residual.min = ((m-cnmle)/(segstats$sd))^2
        residual.other =
            apply(cbind(
        (m-cnmle-1)/segstats$sd,
        (m-cnmle+1)/segstats$sd
        )^2,
        1, min)
        residual.diff = residual.other - residual.min ## penalty for moving to closest adjacent copy state

        ## we fix nodes for which the penalty for moving to non (locally) optimal copy state
        ## is greater than k / slack.prior penalty (where k is some copy difference
        ## that we would never imagine a "reasonable" slack to have to over-rule
        ## fix = as.integer(which(residual.diff>(8/slack.prior)))
        ## 8 is a constant that is conservative, but basically assumes that no node will have more than 4 neighbors (todo: make adjustable per node)
        ## fix = as.integer(which(residual.diff>(4/slack.prior)))
        ## 8 is a constant that is conservative, let's try 4
        ## let's not fix to zero
        fix = as.integer(which(residual.diff>(4/slack.prior) &
                               cnmle != 0))

        ## If we have too few fixed nodes, we will have too few subgraphs to optimize,
        ## each bigger and harder to solve
        if (verbose)
        {
            jmessage('Fixing ', length(fix), ' nodes that are unmovable by slack ')
        }

        ##
        ## now we will create a graph of unfixed nodes and fixed node "halves"
        ## i.e. we split each fixed node to a node that is receiving edges
        ## and a node that is sending edges
        ##
        unfix = as.numeric(setdiff(1:length(segstats), fix))
        G = graph(as.numeric(t(Matrix::which(adj!=0, arr.ind = T))), n = length(segstats), directed = T)
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
            tofix = Matrix::which(adj[unfix, fix]!=0, arr.ind = T)
            fromfix = Matrix::which(adj[fix, unfix]!=0, arr.ind = T)
        }
        else
        {
            tofix = c()
            fromfix = c()
        }

        if (length(fix)>0)
            fixtofix = Matrix::which(adj[fix, fix]!=0, arr.ind = T)
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
        {
            jmessage('Partitioned graph into ', length(cll), ' connected components with the size of the highest 10 components being:\n',
                     paste(sapply(cll[1:min(10, length(cll))], length), collapse = ','), '')
        }

        cn.fix = ifelse(1:length(segstats) %in% fix, cnmle, NA)

        ## force "non lazy" evaluation of args in order to avoid weird R ghosts (WTF) downstream in do.call
        args = as.list(match.call())[-1]
        args = structure(lapply(names(args), function(x) eval(parse(text = x))), names = names(args))


        if (is.null(ploidy.normal))
        {
            if (field.ncn %in% names(values(segstats)))
            {
                args$ploidy.normal = as.data.table(segstats)[, sum(ncn*as.numeric(width), na.rm = TRUE)/sum(ncn*0+1*as.numeric(width), na.rm = TRUE)]
                
            }
        }

        sols = mclapply(1:length(cll), function(k, args)
        {
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

            if (!is.null(mipstart))
                args$mipstart = mipstart[uix, uix]

            args$adj.nudge = adj.nudge[uix, uix, drop = F]
            args$na.node.nudge = na.node.nudge
            args$adj.lb = adj.lb[uix, uix, drop = F]
            if (!is.null(adj.ub)){
                args$adj.ub = adj.ub[uix, uix, drop = F] ## xt added 5/4
            }
            args$segstats = segstats[uix]
            args$cn.fix = cn.fix[uix]
            args$cn.lb = cn.lb[uix]
            args$cn.ub = cn.ub[uix]
            args$partition = F
            args$nsolutions = 1
            args$ploidy.min = 0 ## no ploidy constraints
            ##          args$ploidy.max = max(c(100, cnmle[ix]), na.rm = T)*1.5
            args$ploidy.max = Inf


            if (verbose)
                jmessage('Junction balancing subgraph ', k, ' of ', length(cll), ' which has ', length(uix), ' nodes comprising ',
                         round(sum(as.numeric(width(segstats[uix])))/2/1e6, 2), ' MB and ', length(unique(seqnames((segstats[uix])))),
                         ' chromosomes, including chrs ', paste(names(sort(-table(as.character(seqnames((segstats[uix])))))[1:min(4,
                                                                                                                                  length(unique(seqnames((segstats[uix])))))]), collapse = ', '))

            if (k==1){
                saveRDS(args, "first.args.rds")
            }
            
            out = do.call('jbaMIP', args)

            ## ## timed out??!
            ## if (out$status==107){
            ##     jmessage("Subgraph ", k, " did not reach optima within time limit: ", args$tilim, "s")
            ## }

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
                                       function(x) {ix = Matrix::which(adj[,x]!=0); paste(ix, '(', out$adj[ix,x], ')', '->', sep = '', collapse = ',')})
        out$segstats$edges.out = sapply(1:length(out$segstats),
                                        function(x) {ix = Matrix::which(adj[x, ]!=0); paste('->', ix, '(', out$adj[x,ix], ')', sep = '', collapse = ',')})

        ncn = rep(2, length(segstats))
        if (!is.null(field.ncn))
            if (field.ncn %in% names(values(segstats)))
                ncn = values(segstats)[, field.ncn]


        nnix = !is.na(out$segstats$mean) & !is.na(out$segstats$sd) & !is.na(out$segstats$cn)

        ## new obj allowing variable normal copy number
        out$obj = 1/4*sum(((out$segstats$cn[nnix] + ncn[nnix]/2*out$gamma - out$beta*out$segstats$mean[nnix])/out$segstats$sd[nnix])^2) +
            1/slack.prior * (sum(out$segstats$eslack.in + out$segstats$eslack.out, na.rm = T)) ## 1/4 because our original objective is 1/2 for pos strand intervals only
        out$nll.cn = (1/2*sum(((out$segstats$cn[nnix] + out$gamma - out$beta*out$segstats$mean[nnix])/out$segstats$sd[nnix])^2))
        out$nll.opt = (1/2*sum(((cnmle[nnix] + out$gamma - out$beta*out$segstats$mean[nnix])/out$segstats$sd[nnix])^2))
        out$gap.cn = as.numeric(1 - out$nll.opt / out$nll.cn)
        out$sols = sols

        return(out)
    }
    
    ## map intervals to their reverse complement to couple their copy number (and edge variables)
    pos.ix = which( as.logical( strand(segstats)=='+') )
    neg.ix = which( as.logical( strand(segstats)=='-') )

    ## "original vertices"
    og.ix = pos.ix

    ## map flipping positive to negative vertices
    rev.ix = match(segstats, gr.flipstrand(segstats))

    ## "duplicates" of og.ix i.e. revcomp vertices
    dup.ix = suppressWarnings(neg.ix[match(segstats[og.ix], gr.flipstrand(segstats[neg.ix]))])

    if (!identical(segstats$mean[og.ix] , segstats$mean[dup.ix]) & !identical(segstats$sd[og.ix] , segstats$sd[dup.ix]))
        stop('Segstats mean or sd not identical for all pos / neg strand interval pairs: check segstats computation')

    edges = Matrix::which(adj!=0, arr.ind = T)

    varmeta = data.frame() ## store meta data about variables to keep track
    consmeta = data.frame() ## store meta data about constraints to keep track
    n = 2*nrow(adj) + nrow(edges) + 2;  # number of vertices + slack variables, number of edges, and beta, gamma parameters.
    v.ix = 1:nrow(adj)
    s.ix = length(v.ix) + v.ix

    varmeta = data.frame(id = v.ix, subid = 1:length(v.ix), label = paste('interval', 1:length(v.ix), sep = ''), type = 'interval',
                         stringsAsFactors = F)
    varmeta = rbind(varmeta, data.frame(id = s.ix, subid = 1:length(s.ix), label = paste('residual', 1:length(s.ix), sep = ''),
                                        type = 'residual', stringsAsFactors = F))


    
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

    vtype = rep('C', n); vtype[c(v.ix, e.ix)] = 'I' ## edge and vertex copy numbers are integers
    lb = rep(0, n); lb[s.ix] = -Inf;

    if (nrow(edges)>0){
        lb[e.ix] = adj.lb[edges] ## lower bound on edges, some are forced in
    }

    if (any(ix <<- !is.na(cn.fix)))
        if (verbose>1)
        {
            jmessage('Fixing copy states on ', sum(ix), ' vertices')
        }

    ## implement lower bounds and fixes
    if (any(!is.na(cn.lb))){
        lb[v.ix[!is.na(cn.lb)]] = cn.lb[!is.na(cn.lb)]
    }

    ub = rep(Inf, n);

    if (any(!is.na(cn.ub)))
        ub[v.ix[!is.na(cn.ub)]] = cn.ub[!is.na(cn.ub)]

    ## set upper bound for some edges
    if (!is.null(adj.ub)){
        ub[e.ix] = ifelse(adj.ub[edges]<0, 0, Inf) ## upper bound on edges, some are forced out
        ## for debug only
        ## saveRDS(ub, "ub.rds")
    }

    varmeta$vtype = vtype
    varmeta$lb = lb
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

    ## vertices that will actually have constraints (i.e. those that have non NA segstats )
    ## XT: if some node's sd is 0, evaluate it to NA
    if (length(zero.sd.ix <- which(segstats$sd==0))>0){
        segstats$sd[zero.sd.ix] <- NA
    }
    v.ix.c = setdiff(v.ix[!is.na(segstats$mean) & !is.na(segstats$sd)], dup.ix)

    v.ix.na = which(is.na(segstats$mean) | is.na(segstats$sd))

                                        # weighted mean across vertices contributing to mean
    if (length(v.ix.c))
        mu.all = (width(segstats)[v.ix.c] %*% segstats$mean[v.ix.c]) / sum(as.numeric(width(segstats)[v.ix.c]))
    else
        mu.all = NA

    ## if (verbose)
    ## {
    ##                                     #      jmessage('cn constraints ..')
    ## }

    if (length(v.ix.c)>0)
    {
        ## take into account (variable) normal cn
        ncn = rep(2, length(segstats))
        if (!is.null(field.ncn))
            if (field.ncn %in% names(values(segstats)))
                ncn = values(segstats)[, field.ncn]

        if (is.null(ploidy.normal))
            ploidy.normal = sum(width(segstats)[v.ix.c]*ncn[v.ix.c]) / sum(as.numeric(width(segstats))[v.ix.c])

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
        Acn[length(v.ix.c)+1, gamma.ix] = ploidy.normal/2; ## taking into account (normal) variable cn
        Acn[length(v.ix.c)+1, beta.ix] = -mu.all;
        bcn = rep(0, nrow(Acn)) ## 
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

    ## dup constraints on vertices
    ## constrain every vertex to get the same copy number as its reverse complement
    Dcn = Zero[rep(1, length(dup.ix)),, drop = F];
    Dcn[cbind(1:nrow(Dcn), dup.ix)] = 1
    Dcn[cbind(1:nrow(Dcn), og.ix)] = -1
    dcn = rep(0, nrow(Dcn))
    sensedcn = rep("E", nrow(Dcn))

    consmeta = rbind(consmeta, data.frame(type = 'Dup', label = paste('Dup', 1:nrow(Dcn)), sense = 'E', b = dcn, stringsAsFactors = F))

    if (edge.slack)
    {
        ## dup constraints on (reverse complement) edge.slack
        ## (these make sure that reverse complement edge.slacks are
        ## given the same solution as their reverse complement)
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

    ## if (!is.na(cn.prior))
    ## {
    ##     if (verbose)
    ##         cat('cn prior .. \n')
    ##     Pcn = Zero[rep(1, length(v.ix)+1), ]
    ##     Pcn[cbind(v.ix, v.ix)] = 1
    ##     Pcn[v.ix, ploidy.ix] = -1
    ##     Pcn[cbind(v.ix, d.ix)] = -1
    ##     Pcn[length(v.ix)+1, v.ix] = width(segstats)/sum(as.numeric(width(segstats)))
    ##     Pcn[length(v.ix)+1, ploidy.ix] = -1;
    ##     bpcn = rep(0, nrow(Pcn))
    ##     Acn = rBind(Acn, Pcn)
    ##     bcn = c(bcn, bpcn)
    ##     sensecn = c(sensecn, rep("E", length(bpcn)))
    ##     lb[d.ix] = -Inf;

    ##     consmeta = rbind(consmeta, data.frame(type = 'CNPrior', label = paste('CNPrior', 1:nrow(Pcn)), sense = 'E', b = bpcn, stringsAsFactors = F))
    ## }

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

    if (is.na(beta.guess) | is.na(gamma.guess))
    {
        ## ploidy constraints
        Aineq = Zero[rep(1, 2), , drop = F];
        Aineq[1:2 , v.ix] = rbind(width(segstats)/sum(as.numeric(width(segstats))), width(segstats)/sum(as.numeric(width(segstats))))
        bineq = c(pmax(0, ploidy.min), pmin(Inf, ploidy.max))
        senseineq = c("G", "L")
        
        ## only constrain ploidy if beta / gamma not specified (ie via non NA guess)
        consmeta = rbind(consmeta, data.frame(type = 'PloidyConstraint', label = paste('PloidyConstraint', 1:nrow(Aineq)), sense = senseineq, b = bineq, stringsAsFactors = F))
    }
    else
    {
        Aineq = NULL
        bineq = NULL
        senseineq = NULL
    }


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

    if (!ignore.edge & nrow(edges)>0)
    {
        ## add edge consistency criteria
        ## for every node that is source of an edge
        ## ensure that sum of weights on outgoing edges
        ## = node weight
        ## do the same for nodes that are targets of edges

        v.ix.s = unique(edges[,1])
        v.ix.t = unique(edges[,2])
        Bs = Zero[v.ix.s, , drop = F]
        Bt = Zero[v.ix.t, , drop = F]

        Bs[cbind(1:nrow(Bs), v.ix.s)] = 1
        Bt[cbind(1:nrow(Bt), v.ix.t)] = 1
        Bs[cbind(match(edges[,1], v.ix.s), e.ix)] = -1
        Bt[cbind(match(edges[,2], v.ix.t), e.ix)] = -1

        if (edge.slack)
        {
            Bs[cbind(1:nrow(Bs), es.s.ix[v.ix.s])] = -1  # provide "fake" edges bringing flux into and out of vertex
            Bt[cbind(1:nrow(Bt), es.t.ix[v.ix.t])] = -1
        }

        B = rBind(Bs, Bt)

        consmeta =
            rbind(consmeta,
                  data.frame(type = 'EdgeSource', label = paste('EdgeSource', 1:nrow(Bs)),
                             sense = 'E', b = 0, stringsAsFactors = F),
                  data.frame(type = 'EdgeTarget', label = paste('EdgeSource', 1:nrow(Bt)),
                             sense = 'E', b = 0, stringsAsFactors = F))
        
        ## populate linear constraints
        Aed = B;
        bed = rep(0, nrow(B))
        senseed = rep("E", length(bed))

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


    ##  EPS = 0.000000000001
    ## if (na.node.nudge) ## added Monday, Oct 23, 2017 05:17:25 PM to prevent unconstrained nodes from blowing up
    ##   if (length(v.ix.na)>0)
    ##   {
    ##     if (verbose)
    ##       jmessage('NA nudging node!!')
    ##     Qobj[cbind(v.ix.na, v.ix.na)] = EPS
    ##   }

    ## #  if (!ignore.cons)
    ## #    Qobj[s.ix[length(s.ix)], s.ix[length(s.ix)]] = 1

    ## # linear portion of objective function
    cvec = Zero[,1]

    if (edge.slack)
    {
        cvec[c(es.s.ix, es.t.ix)] = 1/slack.prior

        ## let any specified "loose ends" have unpenalized slack
        if (length(loose.ends)>0)
        {
            cvec[c(es.s.ix[loose.ends], es.t.ix[loose.ends])] = 0
        }

        if (verbose>1)
        {
            jmessage(sprintf('Total mass on cn portion of objective function: %s. Total mass on edge slack: %s', sum(Qobj[cbind(s.ix, s.ix)]), sum(cvec[cbind(es.s.ix, es.t.ix)])))
        }
        if (is.infinite(sum(Qobj[cbind(s.ix, s.ix)]))){
            jmessage("Things are gonna fall apart. Brace yourself. There is some node with zero sd.")
        }

    }
    
    if (nrow(edges)>0)
    {
        if (verbose)
        {
                                        #        jmessage('Adding ', sum(adj.nudge[edges]), " of edge nudge across", sum(adj.nudge[edges]!=0), "edges")
        }

        ## Future to do: weigh the edges in objective functions
        ## what is the conversion from supporting reads to copy number space?
        en = max(abs(adj.nudge), na.rm=T)
        if (!is.na(en)){
            if (abs(en)>0){
                e.penalty = abs(en * 1.1)
            } else {
                e.penalty = 0.01
            }
        } else {
            e.penalty = 0.01
        }

        cvec[e.ix] = e.penalty-adj.nudge[edges] ### reward each edge use in proportion to position in edge nudge
    }


    ## if (!is.na(cn.prior))
    ##     Qobj[cbind(d.ix, d.ix)] = 1/cn.prior^2

    ## the slack prior will determine the degree of "coupling" enforced between neighboring copy states
    ## this should be high if we think that our rearrangement annotation is quite complete
    ## in the end, there will be tension between enforcing edge consistency and consistency with means / sd
    ## abundances at intervals


    if (!any(is.na(purity.prior)))
        Qobj[cbind(pd.ix, pd.ix)] = 1/purity.prior[2]^2

    if (!is.null(mipstart))
    {
        if (is.na(gamma.guess) | is.na(beta.guess))
        {
            warning("Can't do JaBbA mipstart without setting purity and ploidy ... ignoring mipstart")
        } else
        {

            mips.dt = as.data.table(Matrix::which(mipstart>0, arr.ind = TRUE))
            setnames(mips.dt, c("row", "col"))
            mips.dt[, cn := mipstart[cbind(row, col)]]
            setkeyv(mips.dt, c("row", "col"))

            ## convert everything to data.tables
            varmeta = as.data.table(varmeta)
            consmeta = as.data.table(consmeta)
            consmeta[, id := 1:.N]
            varmeta[, id := 1:.N]
            varmeta[, mipstart := as.numeric(0)]

            ## first mipstart the node copy number n_hat by rounding
            ## mu_hat = ifelse(!is.na(segstats$mean), segstats$mean*beta.mipstart-ncn/2*gamma.mipstart, 0)
            ix = varmeta[type == 'interval', id]
            ## n_hat = pmin(pmax(round(mu_hat), lb[ix]), ub[ix])
            cnr = mips.dt[, .(cn = sum(cn, na.rm = TRUE)),  keyby = 'row']
            cnc = mips.dt[, .(cn = sum(cn, na.rm = TRUE)),  keyby = 'col']
            varmeta[type == 'interval', mipstart := cnc[.(subid), cn]]
            varmeta[type == 'interval', mipstart := pmax(mipstart, cnr[.(subid), cn], na.rm = TRUE)]
            varmeta[type == 'interval' & is.na(mipstart), mipstart := 0]


            varmeta[type == 'edge', mipstart := mips.dt[.(as.data.table(edges[subid, ])), cn]]
            varmeta[type == 'edge' & is.na(mipstart), mipstart := 0]
            varmeta[, mipstart := pmax(pmin(mipstart, ub), lb)]
            ## ## mipstart each edge copy number by iterating through edges and peeling off the min copy number of source
            ## ## and sink
            ## 
            ## e_done = rep(FALSE, length(e_hat))
            ## rev.eix = mmatch(edges, cbind(rev.ix[edges[,2]], rev.ix[edges[,1]])) ## match edges to their reverse complements
            ## jmessage('Computing initial mipstart on edge weight')
            ## for (i in 1:nrow(edges))
            ## {
            ##   if (!e_done[i]) ## specify both edge and its reverse complement (so if e_hat[i] is non NA, it has already been filled)
            ##   {
            ##     if (edges[i, 1] == edges[i, 2]) ## be careful overdrawing fold back edges, otherwise can generate strand imbalance
            ##       e_hat[rev.eix[i]] = e_hat[i] = max(e_hat[i], floor(min(n_hat[edges[i, 1]], n_hat[edges[i, 2]], na.rm = TRUE)/2))
            ##     else
            ##       e_hat[rev.eix[i]] = e_hat[i] = max(e_hat[i], min(n_hat[edges[i, 1]], n_hat[edges[i, 2]], na.rm = TRUE))
            ##     ## make sure to also change reverse complement
            ##     n_hat[rev.ix[edges[i,1]]] = n_hat[edges[i,1]] = n_hat[edges[i,1]] - e_hat[i]
            ##     n_hat[rev.ix[edges[i,2]]] = n_hat[edges[i,2]] = n_hat[edges[i,2]] - e_hat[i]
            ##     e_done[i] = e_done[rev.eix[i]] = TRUE
            ##   }
            ## }

            ## now mipstart each target and source slack
            ## which is just the difference between the copy number at
            ## c_i and the incoming / outgoing edges

            e_hat = varmeta[type == 'edge', mipstart]
            n_hat = varmeta[type == 'interval', mipstart]
            ## Bs stores the constraints c_i - - slack_s_i - sum_j \in Es(i) e_j for all i in six
            Bs = Amat[consmeta[type == 'EdgeSource', id],]

            ## sometimes it's too big!
            tmp.Bs.interval = Bs[, varmeta[type == "interval", id]]
            if (prod(dim(tmp.Bs.interval))>.Machine$integer.max){
                chunk.num = ceiling(ncol(tmp.Bs.interval)/floor(.Machine$integer.max/nrow(tmp.Bs.interval)))
                chunk.ix = cut(seq_len(ncol(tmp.Bs.interval)), chunk.num, labels=FALSE)
                six = lapply(seq_len(chunk.num),
                             function(chunk){
                                 jmessage("Processing chunk ", chunk)
                                 apply(tmp.Bs.interval[, which(chunk.ix==chunk), drop=FALSE],
                                       1, function(x) which(x!=0))
                             })
                six = unlist(six)
            } else {
                six = apply(tmp.Bs.interval, 1, function(x) which(x!=0))
            }
            rm("tmp.Bs.interval"); gc()

            s_slack_hat = rep(0, length(n_hat))
            if (length(varmeta[type == "edge", id])>0)
                s_slack_hat[six] = Bs[, varmeta[type == "edge", id]] %*% e_hat + n_hat[six]
            s_slack_hat[is.na(s_slack_hat)] = 0
            varmeta[type == 'source.slack', mipstart := s_slack_hat]

            ## Bt stores the constraints c_i - - slack_t_i - sum_j \in Et(i) e_j for all i in tix
            Bt = Amat[consmeta[type == 'EdgeTarget', id],]
            ## tix = apply(Bt[, varmeta[type == "interval", id]], 1, function(x) which(x!=0))
            ## sometimes it's too big!
            tmp.Bt.interval = Bt[, varmeta[type == "interval", id]]
            if (prod(dim(tmp.Bt.interval))>.Machine$integer.max){
                chunk.num = ceiling(ncol(tmp.Bt.interval)/floor(.Machine$integer.max/nrow(tmp.Bt.interval)))
                chunk.ix = cut(seq_len(ncol(tmp.Bt.interval)), chunk.num, labels=FALSE)
                tix = lapply(seq_len(chunk.num),
                             function(chunk){
                                 jmessage("Processing chunk ", chunk)
                                 apply(tmp.Bt.interval[, which(chunk.ix==chunk), drop=FALSE],
                                       1, function(x) which(x!=0))
                             })
                tix = unlist(tix)
            } else {
                tix = apply(tmp.Bt.interval, 1, function(x) which(x!=0))
            }
            rm("tmp.Bt.interval"); gc()
            
            t_slack_hat = rep(0, length(n_hat))
            if (length(varmeta[type == "edge", id])>0)
                t_slack_hat[tix] = Bt[, varmeta[type == "edge", id]] %*% e_hat + n_hat[tix]
            t_slack_hat[is.na(t_slack_hat)] = 0

            varmeta[type == 'target.slack', mipstart := t_slack_hat]
            varmeta[label == 'beta', mipstart := beta.guess]
            varmeta[label == 'gamma', mipstart := gamma.guess]

            ## ## fix negative loose ends, by adding copy number to the node and then slack at the other end
            ## ## (these occur from nonzero lower bounds on junctions)
            ## ## i.e. for negative source slack add copy number to node, and positive target slack
            ## if (any(neg.source <- varmeta[type == 'source.slack', which(mipstart<0)]))
            ## {
            ##   tmp = varmeta[type == 'source.slack', mipstart[neg.source]]
            ##   varmeta[type == 'source.slack' & subid %in% neg.source, mipstart := 0]
            ##   varmeta[type == 'interval' & subid %in% neg.source, mipstart := mipstart - tmp]
            ##   varmeta[type == 'target.slack' & subid %in% neg.source , mipstart := mipstart - tmp]
            ## }

            ## if (any(neg.target <- varmeta[type == 'target.slack', which(mipstart<0)]))
            ## {
            ##   tmp = varmeta[type == 'target.slack', mipstart[neg.target]]
            ##   varmeta[type == 'target.slack' & subid %in% neg.target, mipstart := 0]
            ##   varmeta[type == 'interval' & subid %in% neg.target, mipstart := mipstart - tmp]
            ##   varmeta[type == 'source.slack' & subid %in% neg.target , mipstart := mipstart - tmp]
            ## }

            ## now we have a lot of slacks
            ## so we need to "straighten" out all high variance nodes
            ## with respect to their nearest low variance node
            ## we draw an arbitrary distinction between low and high variance using

            ## m = rel2abs(segstats, gamma = gamma.mipstart, beta = beta.mipstart, field = 'mean', field.ncn = field.ncn)
            ## cnmle = round(m) ## MLE estimate for CN
            ## residual.min = ((m-cnmle)/(segstats$sd))^2
            ## residual.other = apply(cbind((m-cnmle-1)/segstats$sd, (m-cnmle+1)/segstats$sd)^2, 1, min)
            ## residual.diff = residual.other - residual.min ## penalty for moving to closest adjacent copy state
            ## fix = as.integer(which(residual.diff>(1/(10*slack.prior)))) ## 8 is a constant that is conservative, but basically assumes that no node will have more than 4 neighbors (todo: make adjustable per node)
            ## branch = Matrix::rowSums(adj!=0)>1 | Matrix::colSums(adj!=0)>1

            ## ## path endpoints are either fix or branch
            ## endpoints = union(fix, branch)

            ## compute residual as difference between rounded and "mean" value
            n_hat = varmeta[type == 'interval', mipstart]
            mu_hat = ifelse(!is.na(segstats$mean), segstats$mean*beta.guess-ncn/2*gamma.guess, 0)
            eps_hat = mu_hat - n_hat
            varmeta[type == 'residual', mipstart := eps_hat]
        }
    }

    ## cap astronomical Qobj values so that CPLEX / gurobi does not freak out about large near-infinite numbers
    ## astronomical = value that is 1e8 higher than lowest value
    qr = range(setdiff(diag(Qobj), 0))
    CPLEX.INFIN = 1e9
    Qobj[cbind(1:nrow(Qobj), 1:nrow(Qobj))] = pmin(CPLEX.INFIN*qr[1], Qobj[cbind(1:nrow(Qobj), 1:nrow(Qobj))])

    ## setup MIP
    if (use.gurobi) # translate into gurobi
    {
        if (verbose)
        {
            jmessage('Running gurobi!')
        }

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

        sol = gurobi::gurobi(model, params = c(list(TimeLimit=tilim), list(...)));
        sol$xopt = sol$x;
    }
    else
    {
        control = c(list(...), list(trace = ifelse(verbose>=2, 1, 0), tilim = tilim, epgap = epgap ,mipemphasis = mipemphasis))
        if (!is.null(mipstart)) ## apply mipstart if provided
            control$mipstart = varmeta$mipstart

        ## slightly adjust the mipstart
        
        sol = Rcplex2(cvec = cvec, Amat = Amat, bvec = b, sense = sense, Qmat = Qobj, lb = lb, ub = ub, n = nsolutions, objsense = "min", vtype = vtype, control = control)
    }
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


    sol.l = lapply(sol.l, function(sol)
    {
        vcn = round(sol$xopt[v.ix])
        ecn = round(sol$xopt[e.ix])
        sol$residual = round(sol$xopt[s.ix])
        sol$beta = sol$xopt[beta.ix]
        sol$gamma = sol$xopt[gamma.ix]
        sol$purity = 2/(2+sol$gamma)
        sol$ploidy = (vcn%*%width(segstats))/sum(as.numeric(width(segstats)))
        sol$adj = adj*0;
        if (length(v.ix.c)>0)
            sol$nll.cn = (sol$xopt[s.ix[v.ix.c]]%*%Qobj[s.ix[v.ix.c], s.ix[v.ix.c]])%*%sol$xopt[s.ix[v.ix.c]]
        else
            sol$nll.cn = NA

        if (length(v.ix.c)>0)
            sol$nll.opt = pp.nll(segstats[v.ix.c], gamma = sol$gamma, beta = sol$beta, field = 'mean', field.ncn = field.ncn)$NLL
        else
            sol$nll.opt = NA
        
        ## supposed to be how far away from naive MLE is the optima
        sol$gap.cn = as.numeric(1 - sol$nll.opt / sol$nll.cn) 
        sol$adj[edges] = ecn;
        sol$segstats = segstats
        sol$segstats$cn = round(vcn)
        sol$segstats$ecn.in = round(Matrix::colSums(sol$adj))
        sol$segstats$ecn.out = round(Matrix::rowSums(sol$adj))
        sol$segstats$edges.in = sapply(1:length(sol$segstats),
                                       function(x) {ix = Matrix::which(adj[,x]!=0); paste(ix, '(', sol$adj[ix,x], ')', '->', sep = '', collapse = ',')})
        sol$segstats$edges.out = sapply(1:length(sol$segstats),
                                        function(x) {ix = Matrix::which(adj[x, ]!=0); paste('->', ix, '(', sol$adj[x,ix], ')', sep = '', collapse = ',')})

        if (edge.slack)
        {
            sol$segstats$eslack.in = round(sol$xopt[es.t.ix])
            sol$segstats$eslack.out = round(sol$xopt[es.s.ix])
            sol$eslack.in = round(sol$xopt[es.t.ix])
            sol$eslack.out = round(sol$xopt[es.s.ix])
        }

        sol$ploidy.constraints = c(ploidy.min, ploidy.max)
        sol$beta.constraints = c(beta.min, beta.max)
        sol$cn.prior = cn.prior
        sol$slack.prior = slack.prior

        return(sol)
    });

    sol.l = sol.l[order(sapply(sol.l, function(x) x$obj))]

    if (length(sol.l)==1)
        sol.l = sol.l[[1]]

    return(sol.l)
}



####################################################################
#' @name JaBbA.digest
#' @rdname internal
#' JaBbA.digest
#'
#' @details
#' processes JaBbA object
#' (1) collapsing segments with same copy number that lack loose ends
#' (2) (optional) +/- adds segments correponding to loose ends
#' (3) outputting edges data frame with colors, and other formatting information
#' (4) outputting junctions GRangesList with copy number, color, lty and other plotting components
#'
#'
#' @param jab JaBbA object "undigested"
#' @param kag karyograph (original karyograph input to JaBbA), if NULL then will "redigest" JaBbA object
#' @param verbose logical flag
#' @param keep.all keep.all (default TRUE) whether to keep 0 copy junctions or collapse segments across these as well
############################################
JaBbA.digest = function(jab, kag = NULL, verbose = T, keep.all = T)
{
    if (is.null(kag))
    {
                                        #              if (verbose)
                                        #                 cat('Redigesting..\n')
        ## redigesting
        ix = rep(TRUE, length(jab$segstats))
        if (!is.null(jab$segstats$loose))
            ix = !jab$segstats$loose

        jab = list(segstats = jab$segstats[ix], adj = jab$adj[ix, ix], ab.edges = jab$ab.edges, purity = jab$purity, ploidy = jab$ploidy)

        #' mimielinski Friday, Jan 26, 2018 06:45:21 PM
        #' fixing eslack
        segstats = gr.fix(jab$segstats)
        segstats$is.left.end = start(segstats)==1
        segstats$is.right.end = end(segstats) == seqlengths(segstats)[as.character(seqnames(segstats))]
        segstats$eslack.in = segstats$cn-colSums(jab$adj)
        segstats$eslack.out = segstats$cn-rowSums(jab$adj)
        segstats$eslack.in[(segstats$is.left.end & as.logical(strand(segstats)=='+')) |
                           (segstats$is.right.end & as.logical(strand(segstats)=='-'))] = 0
        segstats$eslack.out[(segstats$is.right.end & as.logical(strand(segstats)=='+')) |
                            (segstats$is.left.end & as.logical(strand(segstats)=='-'))] = 0
        jab$segstats = segstats

        tmp.adj = jab$adj*0
        jab$segstats$ref.child = gr.match(gr.flipstrand(flank(gr.flipstrand(jab$segstats),1)), jab$segstats, ignore.strand = F)
        jab$segstats$ref.parent = gr.match(flank(jab$segstats,1), jab$segstats, ignore.strand = F)
        ix = suppressWarnings(cbind(1:length(jab$segstats), jab$segstats$ref.child))
        nnaix = Matrix::rowSums(is.na(ix))==0
        if (any(nnaix))
            tmp.adj[ix[nnaix,, drop = F]] = 1 ## fill in all reference edges

        ## mark all non-NA aberrant edges
        nnab = !ifelse(is.na(jab$ab.edges[,3,1]), TRUE, jab$ab.edges[,3,1]==0)  ## TRUE only if ab.edges has non zero and non NA edge id in ab.edges[,3]
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

    #' mimielinski Friday, Jan 26, 2018 07:43:18 PM
    #' rewriting to get rid of strange edge cases
    #' resulting from unnecessarily having to use coordinates
    #' to match up loose ends with their nodes

    
    if (any(jab$segstats$eslack.out>0 | jab$segstats$eslack.in>0, na.rm = TRUE))
    {
        sink.ix = which(jab$segstats$eslack.out>0)
        sinks = gr.end(jab$segstats[sink.ix],ignore.strand = FALSE)
        sinks$cn = jab$segstats$eslack.out[sink.ix]
        sinks$partner.id = sink.ix
        sinks$id = nrow(adj) + 1:length(sinks)
        sinks$loose = TRUE
        sinks$right = as.logical(strand(sinks)=='+')

        source.ix = which(jab$segstats$eslack.in>0)
        sources = gr.start(jab$segstats[source.ix],ignore.strand = FALSE)
        sources$cn = jab$segstats$eslack.in[source.ix]
        sources$partner.id = source.ix
        sources$id = nrow(adj) + length(sinks) + 1:length(sources)
        sources$loose = TRUE
        sources$right = as.logical(strand(sources)=='+')

        nlends = length(sources) + length(sinks)

        ## pad original matrix with new nodes
        adj.plus = rBind(cBind(adj, sparseMatrix(1,1,x = 0, dims = c(nrow(adj), nlends))),
                         cBind(sparseMatrix(1,1,x = 0, dims = c(nlends, ncol(adj))), sparseMatrix(1,1,x = 0, dims = c(nlends, nlends))))

        ## add new edges
        adj.plus[cbind(sinks$partner.id, sinks$id)] = sinks$cn+0.01
        adj.plus[cbind(sources$id, sources$partner.id)] = sources$cn+0.01
        adj = adj.plus
        segstats = grbind(jab$segstats, sinks, sources)
        segstats$loose = F
        values(segstats) = rrbind(values(jab$segstats), values(sinks), values(sources))
    }
    else
    {
        segstats = jab$segstats
        segstats$loose = FALSE
    }


    ##  lends = loose.ends(jab, kag)
    ##  if (!is.null(lends))
    ##  lends = lends[lends$type != 'type4'] ## don't include type 4 loose ends (i.e. unused rearrangements)

    ## if (length(lends)>0)
    ##   {
    ##     strand(lends) = '+'
    ##     lends = c(lends, gr.flipstrand(lends))
    ##     lends$partner.id = gr.match((lends), jab$segstats, ignore.strand = F)
    ##     lends$id = nrow(adj) + c(1:length(lends))
    ##     lends$right = end(lends) == end(jab$segstats)[lends$partner.id]
    ##     adj.plus = rBind(cBind(adj, sparseMatrix(1,1,x = 0, dims = c(nrow(adj), length(lends)))),
    ##       cBind(sparseMatrix(1,1,x = 0, dims = c(length(lends), ncol(adj))), sparseMatrix(1,1,x = 0, dims = c(length(lends), length(lends)))))

    ##     ## right side ends of '+' and left side ends of '-' are sinks
    ##     sink.ix = as.logical((lends$right & as.logical( strand(lends)=='+') )| (!lends$right & as.logical( strand(lends)=='-')) )
    ##     adj.plus[cbind(lends$partner.id, lends$id)[sink.ix, , drop = F]] = lends$cn[sink.ix]+0.01
    ##     adj.plus[cbind(lends$id, lends$partner.id)[!sink.ix, , drop = F]] = lends$cn[!sink.ix]+0.01
    ##     adj = adj.plus
    ##     lends$loose = T
    ##     segstats$loose = F
    ##     segstats = grbind(jab$segstats, lends)
    ##     values(segstats) = rrbind(values(jab$segstats), values(lends))
    ##   }
    ## else
    ##   segstats$loose = F

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
    edge.ix = Matrix::which(collapsed$adj!=0, arr.ind = T)

    if (nrow(edge.ix)>0)
    {
        out$adj[edge.ix] = round(adj[cbind(out$segstats$end.ix[edge.ix[,1]], out$segstats$start.ix[edge.ix[,2]])])
        adj.new.ix[edge.ix] = 1:nrow(edge.ix)
    }

    out$ab.edges = array(NA, dim = c(nrow(kag$ab.edges), 3, 2), dimnames = list(NULL, c('from', 'to', 'edge.ix'), c('+', '-')))

    ## match ab edges to new graph, excluding any edges that aren't included in graph (i.e. not given >0 copy number)
    ## (tmp.ix may map some ab.edges to "internal" vertices, so need to weed these out via keep)
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

    out$td = out$gtrack = gTrack(ss, y.field = 'cn', edges = out$edges[order(out$edges$cn), ], name ='JaBbA', angle = 0)
    out$purity = jab$purity
    out$ploidy = jab$ploidy

    return(out)
}

####################
#' @name jabba.alleles
#' @rdname internal
#' jabba.alleles
#'
#' Populates allelic value s for JaBbA object.  This does not explicitly impose junction balance constraints on alleles, but rather just computes
#' the maximum likelihood estimate given allelic counts and the inferred total copy number on a given segment according to JaBbA
#'
#' @param jab JaBbA object
#' @param het.sites GRanges with meta data fields (see below) for alt and rref count
#' @param alt.count.field character specifying alt.count meta data field in input het.sites (default $alt)
#' @param ref.count.field character specifying ref.count meta data field in input het.sites (default $ref)
#' @param split.ab logical flag whether to split aberrant segmetns (segmentss with ab edge entering or leaving prior to computing allelic states (default FALSE)
#' @param uncoupled logical flag whether to not collapse segments after inferring MLE estimate (default FALSE), if FALSE will try to merge adjacent segments and populate allele-specific junctions with copy numbers on the basis of the MLE fit on individual allelic segments
#' @param conservative if TRUE then will leave certain allelic segments "unphased" if one cannot sync the high / low interval state with the incoming and / or outgoing junction state
#' @param verbose logical flag
#' @return
#' list with following fields:
#' $segstats = GRanges of input segments with $cn.high and $cn.low segments populated
#' $asegstats = GRanges of allelic segments (length is 2*length(segstats)) with high and low segments each having $cn, this is a "melted" segstats GRAnges
#' $agtrack = gTrack of allelic segments and supporting input het.sites
#' $aadj = allelic adjacency matrix of allele specific junctions
#' $ab.ix = indices of aberrant edges in $aadj
#' $ref.ix = indices of reference edges in $aadj
############################################
jabba.alleles = function(
                         jab,
                         het.sites, ## granges with meta data fields for alt.count and
                         alt.count.field = 'alt',
                         ref.count.field = 'ref',
                         baf.field = 'baf.t',
                         split.ab = F, ## if split.ab == T, then will split across any "aberrant" segment (i.e. segment with ab edge entering or leaving prior to computing allelic states (note: this might create gaps)
                         uncoupled = FALSE, ## if uncoupled, we just assign each high low allele the MLE conditioning on the total copy number
                         conservative = FALSE, ## if TRUE then will leave certain allelic segments "unphased" if one cannot sync the high / low interval state with the incoming and / or outgoing junction state
                         verbose = F
                         )
{
    if (!all(c(alt.count.field, ref.count.field) %in% names(values(het.sites)))){
        warning('count fields not found in meta data of het.sites input, trying BAF...')
        if (!(baf.field %in% names(values(het.sites))))
            stop('BAF field not found in meta data of het.sites input either!')
        else{
            ## outputs are re.seg$low and re.seg$high
            ## test deviations of observed BAF from expected by beta distribution
            if (verbose)
                jmessage('Processing', length(het.sites),
                         'het sites using fields', baf.field, '\n')

        }
    } else {
        ## stop('count fields not found in meta data of het.sites input')

        if (verbose)
        {
            jmessage('Processing ', length(het.sites), ' het sites using fields ', alt.count.field, ' and ', ref.count.field)
        }
        
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
        {
            jmessage('Computed high / low counts and matched to segs')
        }


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
        {
            jmessage('Computed SNP ploidy and allelic copy centers')
        }

        ## now test deviation from each absolute copy combo using poisson model
        ## i.e. counts ~ poisson(expected mean)
        ##
        re.seg$low = sapply(1:length(re.seg), function(i)
        {
            ##        if (verbose)
            ##          cat('.')
            x = lows[[i]]
            if (length(x)==0)
                return(NA)
            y = highs[[i]]
            tot.cn = cn[i]
            ll = sapply(0:(floor(tot.cn/2)), function(j) sum(ppois(x,centers[j+1], log.p = T) + ppois(y,centers[tot.cn-j+1],log.p = T)))
            ll = ll - min(ll)
            return(which.max(ll)-1)
        })

        re.seg$high = re.seg$cn-re.seg$low
    }
    ## #########################################################################
    ## borderline, below are common to both methods
    jab$segstats$cn.low = round(gr.val(jab$segstats, re.seg, 'low', na.rm = TRUE)$low)
    jab$segstats$cn.high = round(gr.val(jab$segstats, re.seg, 'high', na.rm = TRUE)$high)
    na.ix = (!gr.val(jab$segstats, re.seg, 'low', FUN = function(x,w,na.rm) any(!is.na(x)))$low) |
        (!gr.val(jab$segstats, re.seg, 'high', FUN = function(x,w,na.rm) any(!is.na(x)))$high)
    jab$segstats$cn.low[na.ix] = jab$segstats$cn.high[na.ix] = NA

    ## ###########
    ## phasing
    ## ###########

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
        jmessage('Starting phasing ')
    ##

    
    ##
    

    for (k in 1:nrow(ref.jun))
    {
        i = ref.jun[k, 1]
        j = ref.jun[k, 2]
        a = acn[ref.jun[k,1],]
        b = acn[ref.jun[k,2],]

        phased.out[amap[i, ]] = FALSE
        phased.in[amap[j, ]] = FALSE

        pairs.ij = cbind(rep(c(1:2), 2), rep(c(1:2), each = 2)) ## 4 possible matches
        m = setdiff(which(a[pairs.ij[,1]] == b[pairs.ij[,2]]), NA)

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

            if (length(a.ab <- Matrix::which(adj.ab[i,]!=0))>0)
            {
                ## if a.ab (partner) is already phased then unpopulate the non-ab allelic junction, otherwise populate both alleles of partner
                ## BUG: a.ab is length 2????
                ## hack: replace a.ab with a.ab[1]
                if (any(ph <- aadj[amap[i, fm.ij[1]], amap[a.ab[1], ]] !=0))
                {
                    aadj[amap[i, fm.ij[1]], amap[a.ab[1], ph]] = adj.ab[i, a.ab[1]]
                    aadj[amap[i, m.ij[1]], amap[a.ab[1], ph]] = 0
                }
                else
                    ## otherwise diffuse copy into both alleles of the partner (will be resolved when we resolve phase for the partner interval)
                    ## or collapse unphased nodes back
                    aadj[amap[i, fm.ij[1]], amap[a.ab[1], ]] = adj.ab[i, a.ab[1]]/2

                if (!conservative)
                    if (a[fm.ij[1]] < adj.ab[i, a.ab]) # if the allelic node can't handle the outgoing allelic edge flux, so unphase
                        phased.out[amap[i, ]] = FALSE
            }

            if (length(b.ab <- Matrix::which(adj.ab[,j]!=0))>0)
            {
                ## if b.ab (partner) is already phased then concentrate all of the junction copy into the aberrant allele of this interval
                ## BUG: why b.ab is length 2???? I thought we resolved this long ago
                ## hack: replace a.ab with a.ab[1]
                if (any(ph <- aadj[amap[b.ab[1], ], amap[j, fm.ij[2]]] !=0))
                {
                    aadj[amap[b.ab[1], ph], amap[j, fm.ij[2]]] = adj.ab[b.ab[1], j]
                    aadj[amap[b.ab[1], ph], amap[j, m.ij[2]]] = 0
                }
                else
                    ## otherwise diffuse copy into both alleles of the partner (will be resolved when we resolve phase for the partner interval)
                    ## or collapse unphased nodes back
                    aadj[amap[b.ab[1],], amap[j, fm.ij[2]]] = adj.ab[b.ab[1], j]/2

                if (!conservative)
                    if (b[fm.ij[2]] < adj.ab[b.ab, j]) # the allelic node cn can't handle the incoming allelic edge flux, so unphase
                        phased.in[amap[j, ]] = FALSE
            }
        }
    }

    if (verbose)
        jmessage('Finished phasing, finalizing.')

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
        jmessage('Annotating allelic vertices')

    tmp.string = gr.string(asegstats, mb = F, other.col = 'type'); tmp.string2 = gr.string(gr.flipstrand(asegstats), mb = F, other.col = 'type')
    asegstats$flip.ix = match(tmp.string, tmp.string2)
    asegstats$phased = !unphased

    asegstats.final$edges.in = sapply(1:length(asegstats.final),
                                      function(x) {ix = Matrix::which(aadj.final[,x]!=0); paste(ix, '(', aadj.final[ix,x], ')', '->', sep = '', collapse = ',')})
    asegstats.final$edges.out = sapply(1:length(asegstats.final),
                                       function(x) {ix = Matrix::which(aadj.final[x, ]!=0); paste('->', ix, '(', aadj.final[x,ix], ')', sep = '', collapse = ',')})

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
        agtrack = atd,
        aadj = aadj.final,
        ab.ix = Matrix::which((m %*% adj.ab %*% t(m))!=0, arr.ind = T),
        ref.ix = Matrix::which((m %*% adj.ref %*% t(m))!=0, arr.ind = T)
    )

    return(out)
}



##############################
#' @name pp.nll
#' @rdname internal
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


#############################################################
#' @name munlist
#' @rdname internal
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
                                        #      param_lines = c(param_lines, paste("CPX_PARAM_WORKDIR", getwd(), sep = '\t'))
        param_lines = c(param_lines, paste("CPX_PARAM_TRELIM", treememlim, sep = '\t'))
    }

    writeLines(param_lines, out.file)
    Sys.setenv(ILOG_CPLEX_PARAMETER_FILE=out.file)
}



#' @name gr.tile.map
#' @rdname internal
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
#' @rdname internal
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

#' @name write.tab
#' @rdname internal
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

############################
#' @name rel2abs
#' @rdname internal
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

#####################################################
#' @name all.paths
#' @rdname internal
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
#####################################################
all.paths = function(A, all = F, ALL = F, sources = c(), sinks = c(), source.vertices = sources, sink.vertices = sinks,
                     exclude = NULL, ## specifies illegal subpaths, all such paths / cycles and
                     ## their supersets will be excluded, specified as k x nrow(A) matrix of vertex sets
                     verbose = FALSE,...)
{
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

###############################################
#' @name collapse.paths
#' @rdname internal
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

                                        #  if (verbose)
                                        #     cat('graph size:', nrow(out), 'nodes\n')

    ## first identify all nodes with exactly one parent and child to do initial collapsing of graph
    singletons = which(Matrix::rowSums(out)==1 & Matrix::colSums(out)==1)

                                        #  if (verbose)
                                        #      cat('Collapsing simple paths..\n')

    sets = split(1:nrow(G), 1:nrow(G))
    if (length(singletons)>0)
    {
        tmp = out[singletons, singletons]
        cl = clusters(graph(as.numeric(t(Matrix::which(tmp, arr.ind = TRUE))), n = nrow(tmp)), 'weak')$membership
        dix = unique(cl)
        if (length(dix)>0)
        {
            for (j in dix)
            {
                                        #                          if (verbose)
                                        #                              cat('.')

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
    {
                                        #    jmessage('done\nnow fixing branches\n')
    }

    todo = rep(FALSE, nrow(G))
    todo[Matrix::rowSums(out)==1 | Matrix::colSums(out)==1] = TRUE

    while (sum(todo)>0)
    {
        sets.last = sets
        out.last = out

                                        #        if (verbose)
                                        #            if ((sum(todo) %% 200)==0)
                                        #                cat('todo:', sum(todo), 'num sets:', sum(!sapply(sets, is.null)), '\n')

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
#' @name sparse_subset
#' @rdname internal
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

##################################
#' @name convex.basis
#' @rdname internal
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
        if (nrow(K_i)==0 | ncol(K_i)==0)
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


############################################
#' @name read.junctions
#' @rdname internal
#' @title read.jucntions
#' @description
#' read.junctions
#'
#' takes in either file or data frame from various formats including BND VCF, bedpe, and others, and returns GRangesList of junctions
#' and returns junctions in VCF format.
#'
#' The default output is GRangesList each with a length two GRanges whose strands point AWAY from the break.  If get.loose = TRUE (only relevant for VCF)
#'
#' @import VariantAnnotation
#' @export
############################################
read.junctions = function(rafile, keep.features = T, seqlengths = hg_seqlengths(), chr.convert = T, snowman = FALSE, swap.header = NULL,  breakpointer = FALSE, seqlevels = NULL, force.bnd = FALSE, skip = NA,
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
            vcf = suppressWarnings(VariantAnnotation::readVcf(rafile, Seqinfo(seqnames = names(seqlengths), seqlengths = seqlengths)))
            if (!('SVTYPE' %in% names(VariantAnnotation::info(vcf)))) {
                warning('Vcf not in proper format.  Is this a rearrangement vcf?')
                return(GRangesList());
            }

            ## vgr = rowData(vcf) ## parse BND format
            vgr = suppressWarnings(read_vcf(rafile, swap.header = swap.header))

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

            if (!is.character(vgr$mateid))
            {
                vgr$mateid = unstrsplit(vgr$MATEID)
                if (any(nix<-nchar(vgr$mateid)==0))
                    vgr$mateid[nix] = NA
            }

            if (is.null(vgr$SVTYPE))
                vgr$svtype = vgr$SVTYPE
            else
                vgr$svtype = vgr$SVTYPE

            if (!is.null(VariantAnnotation::info(vcf)$SCTG))
                vgr$SCTG = VariantAnnotation::info(vcf)$SCTG

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
                    values(vgr.loose) = cbind(vcf@fixed[bix[npix], ],
                                              VariantAnnotation::info(vcf)[bix[npix], ])

                return(list(junctions = ra, loose.ends = vgr.loose))
            }
        }
        else
            rafile = read.delim(rafile)
    }

    if (is.data.table(rafile))
    {
        rafile = as.data.frame(rafile)
    }

    out = GRangesList()

    if (nrow(rafile)==0)
    {
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


#' @name karyograph
#' @rdname internal
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
############################################
karyograph = function(junctions, ## this is a grl of breakpoint pairs (eg output of read.junctions(dranger.df) where dranger is df of dranger output)
                      tile = NULL, # pre-existing set of intervals on top of which to build a graph (eg endpoints from a copy number based segmentation)
                      label.edges = FALSE
                      )
{
    if (length(junctions)>0)
    {
        bp.p = grl.pivot(junctions)
        bp1 = suppressWarnings(gr.end(gr.fix(bp.p[[1]]), 1, ignore.strand = F))
        bp2 = suppressWarnings(gr.start(gr.fix(bp.p[[2]]), 1, ignore.strand = F))


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
                     G = igraph::graph.adjacency(A), ab.adj = A != 0, ab.edges = NULL, junctions = junctions))
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
    tmp.nm = as.character(c(1:length(tile), -(1:length(tile))))
    tile = c(tile, gr.flipstrand(tile))
    names(tile) = tmp.nm

    ## apply ix to adj.ref and adj.ab, and create "adj" which has union of reference and aberrant junctions
    ## and adj.source which remembers whether edge ij was reference (value = 1) or aberrant (value = 2)
    adj.source = sign(adj.ref)+2*sign(adj.ab)
    adj = sign(adj.ref)+sign(adj.ab)
    tryres <- try( edges <- Matrix::which(adj!=0, arr.ind=T), silent=T ) ## num edge x 2 matrix of vertex pairs

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
                                  function(x) {ix = Matrix::which(adj[,x]!=0); paste(ix, '->', sep = '', collapse = ',')})
        tile$edges.out[ix] = sapply(ix,
                                    function(x) {ix = Matrix::which(adj[x, ]!=0); paste('->', ix,  sep = '', collapse = ',')})
    }

    return(list(tile = tile, adj = adj, G = G, ab.adj = adj.ab != 0, ab.edges = ab.edges, junctions = junctions))
}




########################################
#' @name jabba2vcf
#' @rdname internal
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
#########################################
jabba2vcf = function(jab, fn = NULL, sampleid = 'sample', hg = NULL, cnv = FALSE)
{
    if (is.null(hg))
        hg = tryCatch(skidb::read_hg(), error = function(e) NULL)

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
            REF = tryCatch(as.character(ffTrack::get_seq(hg, gr.stripstrand(gr.start(ss,1)))), error = function(e) NULL)
            if (is.null(REF))
            {
                REF = 'N'
            }
            body = data.frame("CHROM" = as.character(seqnames(ss)), POS = start(ss),
                              ID = paste(sampleid, '_seg', six, sep = ''),
                              REF = REF,
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

#################################################
#' @name chromoplexy
#' @rdname internal
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
#####################################################
chromoplexy = function(kag = NULL, # output of karyograph
                       jab = NULL, ## optional alternate input, if NOT null then this will be used in place of kag
                       sol = NULL, ## if sol is null, then copy state is ignored when determining amp or del bridges
                       all = F, ## if TRUE, will try to enumerate all possible cycles, otherwise will return (an arbitrary) minimal decomposition into the shortest "chains" of balanced rearrangements
                       ref.only = F, ## if T will only compute distance criteria on reference (i.e. won't use any subsequent rearrangements)
                       filt.jab = T, ## filter out 0 copy edges if input is a jabba object
                       reciprocal = TRUE, ## aka deletion bridge
                       hijacked = TRUE,  ## aka amplification bridge
                       paths = F,
                       dist = 1e3,
                       cn.dist = dist,
                       verbose = F,
                       interval = 400,
                       junc.only=TRUE,
                       mc.cores = 1,
                       chunksize = 5000)
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
        G = graph(as.numeric(t(Matrix::which(adj2!=0, arr.ind = T))), n = length(kag$segstats), directed = T)
    }

    ## define edge source to edge sink distance
    ## this is minimum between (1) sum of vertex width of path from e2 source to e1 sink (excluding source and sink)
    ## and (2) sum of vertex width of path from e1 sink to e2 source (including source and sink)

    tmp = igraph::get.edges(G, E(G))
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

    ## quasi pairvvs are ab edge pairs within a certain distance of each other on the graph
    quasi.pairs = which(D<dist, arr.ind = T)
    quasi.pairs.which = D.which[quasi.pairs]


    ## now need to check .. depending on whether edge pair is deletion bridge or amp bridge or fully reciprocal
    ## whether associated vertices show a copy change "in the right direction"

    ## to do this, we need to examine the vertices "in between" for a deletion bridge and the source / sink vertices
    ## in an amplification bridge, and see if they show a copy change with respect their reference parents

    ## for reciprocal pairs, the source and sink will be the same

    adj.ref = kag$adj; adj.ref[ab.edges[, 1:2]] = 0
    if (nrow(quasi.pairs) * nrow(adj.ref) > .Machine$integer.max){
        warning("Exceeding size limit. Empty integer will be returned. We will fix it later.")
        return(integer(0))
    }

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

    if (junc.only){
        G.bp = igraph::graph_from_adjacency_matrix(adj.bp)
        comp = components(G.bp, "strong")
        good.comp = which(comp$csize>1)
        good.ix = which(comp$membership %in% good.comp)
        return(unique(ab.edges[good.ix, 3]))
    } else {
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
        }


        ## xtYao modified: mclapply to replace lapply and sapply
        if (length(pc$cycles)>0)
        {
            pc$cycles = pc$cycles[!unlist(mclapply(pc$cycles, .check.pc, is.cycle = T, mc.cores=mc.cores))]
            pc$cycles = mclapply(pc$cycles, function(x) sign(emap[x])*edge.ix[abs(emap[x])], mc.cores=mc.cores)
            pc$cycles = pc$cycles[!duplicated(unlist(mclapply(pc$cycles, function(x) paste(unique(sort(x)), collapse = ' '), mc.cores=mc.cores)))]
            pc$cycles = pc$cycles[order(-unlist(mclapply(pc$cycles, length, mc.cores = mc.cores)))]
        }

        if (length(pc$paths)>0)
        {
            pc$paths = pc$paths[!unlist(mclapply(pc$paths, .check.pc, is.cycle = F, mc.cores = mc.cores))]
            pc$paths = mclapply(pc$paths, function(x) sign(emap[x])*edge.ix[abs(emap[x])], mc.cores = mc.cores)
            pc$paths = pc$paths[!duplicated(unlist(mclapply(pc$paths, function(x) paste(unique(sort(x)), collapse = ' '), mc.cores = mc.cores)))]
            pc$paths = pc$paths[order(-unlist(mclapply(pc$paths, length, mc.cores = mc.cores)))]
        }

        return(pc)
    }
}






#' @name read_vcf
#' @rdname internal
#' @title read_vcf
#'
#' @description
#'
#' wrapper around variantAnnotation reads VCF into granges or data.table format
#'
#' @author Marcin Imielinski
read_vcf = function(fn, gr = NULL, hg = 'hg19', geno = NULL, swap.header = NULL, verbose = FALSE, add.path = FALSE, tmp.dir = '~/temp/.tmpvcf', ...)
{
    in.fn = fn

    if (verbose)
        cat('Loading', fn, '\n')

    if (!is.null(gr))
    {
        tmp.slice.fn = paste(tmp.dir, '/vcf_tmp', gsub('0\\.', '', as.character(runif(1))), '.vcf', sep = '')
        cmd = sprintf('bcftools view %s %s > %s', fn,  paste(gr.string(gr.stripstrand(gr)), collapse = ' '), tmp.slice.fn)
        if (verbose)
            cat('Running', cmd, '\n')
        system(cmd)
        fn = tmp.slice.fn
    }

    if (!is.null(swap.header))
    {
        if (!file.exists(swap.header))
            stop(sprintf('Swap header file %s does not exist\n', swap.header))

        system(paste('mkdir -p', tmp.dir))
        tmp.name = paste(tmp.dir, '/vcf_tmp', gsub('0\\.', '', as.character(runif(1))), '.vcf', sep = '')
        if (grepl('gz$', fn))
            system(sprintf("zcat %s | grep '^[^#]' > %s.body", fn, tmp.name))
        else
            system(sprintf("grep '^[^#]' %s > %s.body", fn, tmp.name))

    system(sprintf("cat %s.header %s.body > %s", tmp.name, tmp.name, tmp.name))
    vcf = VariantAnnotation::readVcf(tmp.name, hg, ...)
    system(sprintf("rm %s %s.body %s.header", tmp.name, tmp.name, tmp.name))
  }
  else
    vcf = VariantAnnotation::readVcf(fn, hg, ...)

    out = granges(vcf)

  if (!is.null(values(out)))
    values(out) = cbind(values(out), VariantAnnotation::info(vcf))
  else
    values(out) = VariantAnnotation::info(vcf)

    if (!is.null(geno))
    {

        if (!is.logical(geno))
            geno = TRUE


        if (geno)
            for (g in  names(geno(vcf)))
            {
                geno = names(geno(vcf))
                warning(sprintf('Loading all geno fields:\n\t%s', paste(geno, collapse = ',')))
            }

        gt = NULL
        if (length(g)>0)
        {
            for (g in geno)
            {
                m = as.data.frame(geno(vcf)[[g]])
                names(m) = paste(g, names(m), sep = '_')
                if (is.null(gt))
                    gt = m
                else
                    gt = cbind(gt, m)
            }
            values(out) = cbind(values(out), as(gt, 'DataFrame'))
        }
    }

    if (!is.null(gr))
        system(paste('rm', tmp.slice.fn))

    if (add.path)
        values(out)$path = in.fn

    return(out)
}

#' @name levapply
#' @rdname internal
#' @title levapply
#'
#' @description
#' Applies FUN locally to levels of x and returns vector of length()
#' (eg can do a "local" order within levels)
#'
#' @param x input vector of data
#' @param by length(x) vector of categorical labels
#' @param FUN function that takes a length k vector and outputs a length k vector, used for processing each "level" of by
#' @return length(x) vector of outputs, the results of applying FUN to each "by" defined level of x
#' @author Marcin Imielinski
levapply = function(x, by, FUN = 'order')
{
    if (!is.list(by))
        by = list(by)

    f = factor(do.call('paste', c(list(sep = '|'), by)))
    ixl = split(1:length(x), f);
    ixv = lapply(ixl, function(y) x[y])
    res = structure(unlist(lapply(ixv, FUN)), names = unlist(ixl))
    out = rep(NA, length(x))
    out[as.numeric(names(res))] = res;
    return(out)
}

chr2num = function(x, xy = FALSE)
{
    if (inherits(x, 'factor') | inherits(x, 'Rle'))
        x = as.character(x)

    out = gsub('chr', '', x);

    if (!xy)
        out = as.numeric(gsub('M', '25', gsub('Y', '24', gsub('X', '23', out))))

    return(out)
}

#' @name which.indel
#' @rdname internal
#' @title which.indel
#'
#' @description
#' Among a GRangesList of junction set, find the indices of isolated, small scale tDup or DEL
#' They are in the grey area from SV to INDEL.
#'
#' @param juncs GRangesList of junctions
#' @param max.size the size cutoff in bp, any pair of breakpoints below this
#' with the correct orintation wil be called
#' @return indices of the identified junctions
which.indel = function(juncs,
                       max.size = 1e4){
    bps = grl.unlist(juncs)
    sort.grl.ix = rle((bps %Q% (order(seqnames, start)))$grl.ix)
    ## criterion 1: they are non-overlapping with others
    iso.ix = sort.grl.ix$values[which(sort.grl.ix$lengths==2)]
    out.ix = iso.ix
    if (length(out.ix)==0){
        return(out.ix)
    }
    juncs.iso = juncs[out.ix]
    ## criterion 2: they are smaller than max.size
    iso.sizes = sv.size(juncs.iso, ignore.strand = TRUE)
    small.ix = which(iso.sizes <= max.size)
    out.ix = iso.ix[small.ix]
    if (length(out.ix)==0){
        return(out.ix)
    }
    ## criterion 3: they need to have opposite directions
    out.ix = gr2dt(bps %Q% (grl.ix %in% out.ix))[
      , oppo := all(c("+", "-") %in% strand), by=grl.ix][oppo==TRUE, unique(grl.ix)]
    return(out.ix)
}

#' @name sv.size
#' @rdname internal
#' @description
#' Simply the distance between pairs of breakpoints
#' @param juncs GRangesList of junctions
#' @param mc.cores parallel
#' @param ignore.strand usually TRUE
#' @return numerical vector of the same length, Inf means they r not facing each other
sv.size = function(juncs,
                   ...){
    bps = gUtils::grl.pivot(juncs)
    return(IRanges::distance(bps[[1]], bps[[2]], ...))
}

## #' @name capply
## #' @description
## #' Wrapper around \{code}apply function
## capply = function(X, MARGIN, FUN){
##     if (prod(dim(X))<.Machine$integer.max){
##         ## no need
##         return(apply(X, margin, FUN))
##     }
##     ## if margin is 1, MARGIN.2 is 2, vise versa
##     MARGIN.2 = ifelse(MARGIN==1, 2, 1)
    
##     chunk.num = ceiling(dim(tmp.Bt.interval)[margin] / floor(.Machine$integer.max/dim(tmp.Bt.interval)[]))
##     chunk.ix = cut(seq_len(ncol(tmp.Bt.interval)), chunk.num, labels=FALSE)
##     tix = lapply(seq_len(chunk.num),
##                  function(chunk){
##                      jmessage("Processing chunk ", chunk)
##                      apply(tmp.Bt.interval[, which(chunk.ix==chunk), drop=FALSE],
##                            1, function(x) which(x!=0))
##                  })
##     tix = unlist(tix)
## }
