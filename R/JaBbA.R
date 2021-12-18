#' @importFrom gUtils streduce si2gr seg2gr rrbind ra.overlaps ra.duplicated parse.gr hg_seqlengths grl.unlist grl.pivot grl.in grl.eval grl.bind grbind gr2dt gr.val gr.tile.map gr.tile
#' @importFrom gUtils gr.stripstrand gr.sum gr.string gr.start gr.end gr.simplify gr.setdiff gr.sample gr.reduce gr.rand gr.quantile gr.nochr
#' @importFrom gUtils gr.match gr.in gr.flipstrand gr.fix gr.findoverlaps gr.duplicated gr.dist gr.disjoin gr.breaks dt2gr "%^%" "%Q%" "%&%" "%$%"
#' @importFrom gGnome gG balance
#' @importFrom GenomicRanges GRanges GRangesList values split match setdiff
#' @importFrom gTrack gTrack
#' @importFrom igraph graph induced.subgraph V E graph.adjacency clusters
#' @importFrom optparse make_option OptionParser parse_args print_help
#' @importFrom data.table data.table as.data.table setnames setkeyv fread setkey
#' @importFrom Matrix which rowSums colSums Matrix sparseMatrix t diag
#' @importFrom parallel mclapply
#' @importFrom gplots col2hex
#' @importFrom graphics plot abline hist title
#' @importFrom grDevices col2rgb dev.off pdf png rgb
#' @importFrom stats C aggregate dist loess median ppois predict runif setNames hclust cutree acf glm ks.test lag quantile
#' @importFrom utils read.delim write.table
#' @importFrom sequenza segment.breaks baf.model.fit get.ci
#' @importFrom rtracklayer import
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom DNAcopy CNA segment smooth.CNA
#' @importFrom methods as is
#' @useDynLib JaBbA

## appease R CMD CHECK misunderstanding of data.table syntax by declaring these global variables
low.count=high.count=seg=chromosome=alpha_high=alpha_low=beta_high=beta_low=predict.var=dup=psid=.N=es=res=esid=pid=ub=lb=esid=dup=lb=ub=dup=y=dup=pid=lb=es.s.ix=km=es.t.ix=adjusted.ratio=ref.count.t=alt.count.t=depth.normal=depth.tumor=good.reads=zygosity.normal=Bf=ALT=alt.count.n=bad=both.na=chr.a=chr.b=cn=eid=FILTER=force.in=FORMAT=from=from.cn=from.remain=from1=from2=GENO=grl.ix=gstr=i=id=ID=INFO=is.ref=j=label=mean.a=mean.b=mstr=nbins=nothing=oppo=ord=out=QUAL=ra=ra1.ix=ra2.ix=ref.count.n=ref.frac.n=reid=str1=str2=subid=this.cn=to=to.cn=to.remain=to1=to2=type=V1=var=NULL

.onLoad <- function(libname, pkgname) {
    op <- options()
    op.JaBbA <- list(
        JaBbA.verbose = TRUE,
        JaBbA.reference = "hg19"
        ## devtools.path = "~/R-dev",
        ## devtools.install.args = "",
        ## devtools.name = "Your name goes here",
        ## devtools.desc.author = "First Last <first.last@example.com> [aut, cre]",
        ## devtools.desc.license = "What license is it under?",
        ## devtools.desc.suggests = NULL,
        ## devtools.desc = list()
    )
    toset <- !(names(op.JaBbA) %in% names(op))
    if(any(toset)) options(op.JaBbA[toset])

    ## test for CPLEX environment
    cplex.dir = Sys.getenv("CPLEX_DIR")
    if (is.null(cplex.dir)){
        jerror("CPLEX_DIR environment variable not found!")
    } else if (!file.exists(paste0(cplex.dir, "/cplex"))) {
        jerror("${CPLEX_DIR}/cplex not found")
    } else if (!file.exists(paste0(cplex.dir, "/cplex/include")) ||
               !file.exists(paste0(cplex.dir, "/cplex/lib"))){
        jerror("${CPLEX_DIR}/cplex/[(include)|(lib)] do not both exist")
    } else {
        jmessage("Found CPLEX environment in: ", cplex.dir)
        ## jmessage("CPLEX version: ", cplex.version)
    }

    invisible()
}

#' @name JaBbA
#' @title JaBbA
#' @description
#' Main function to invoke junction balance analysis (JaBbA). Two input arguments are required: junctions and coverage. One output contains all information of the results, the gGraph saved in "jabba.simple.gg.rds".
#'
#' @param junctions rearrangement junctions  (i.e. breakpoint pairs with orientations). Supports BEDPE, BND VCF formats, Junction objects defined in gGnome, and GRangesList object. If providing GRangesList, the orientation must be "+" for a junction that fuses the side with larger coordinates and vice versa
#' @param coverage high-density read coverage data of constant-width genomic bins. Supports BED, BigWig, delimited text formats, and GRanges object
#' 
#' @param juncs.uf supplement junctions in the same format as \code{junctions}
#' @param blacklist.junctions rearrangement junctions to be excluded from consideration
#' @param whitelist.junctions rearrangement junctions to be forced to be incorporated
#' @param geno logical whether the input junctions have GENO field in the metadata. If so, will match the name argument to the corresponding sample and make positive junctions tier 2, the others tier 3.
#' @param indel character of the decision to "exclude" or "include" small(< min.nbins * coverage bin width) isolated INDEL-like events into the model. Default NULL, do nothing.
#' @param cfield  character, junction confidence meta data field in ra
#' @param tfield  character, tier confidence meta data field in ra. tiers are 1 = must use, 2 = may use, 3 = use only in iteration>1 if near loose end. Default "tier".
#' @param reiterate integer scalar specifying how many extra re-iterations allowed, to rescue lower confidence junctions that are near loose end. Default 0. This requires junctions to be tiered via a metadata field tfield.
#' @param rescue.window integer window size in bp within which to look for lower confidence junctions. Default 1000.
#' @param nudge.balanced logical whether to attempt to add a small incentive for chains of quasi-reciprocal junctions.
#' @param thresh.balanced numeric maximum distance between a pair of reciprocal junctions. Default 500.
#' ## TODO think whether this is the reason balanced junction pairs tend to have too many copies
#' @param edgenudge  numeric hyper-parameter of how much to nudge or reward aberrant junction incorporation. Default 0.1 (should be several orders of magnitude lower than average 1/sd on individual segments), a nonzero value encourages incorporation of perfectly balanced rearrangements which would be equivalently optimal with 0 copies or more copies.
#' @param strict logical flag specifying whether to only include junctions that exactly overlap segs
#' @param all.in whether to use all of the junctions but the tier 3 INDELs all at once
#' 
#' @param field name of the metadata column of coverage that contains the data. Default "ratio" (coverage ratio between tumor and normal). If using dryclean, usually it is "foreground".
#' @param seg  optional path to existing segmentation, if missing then will segment coverage using DNACopy with standard settings
#' @param nseg  optional path to normal seg file with $cn meta data field
#' @param hets  optional path to hets.file which is tab delimited text file with fields seqnames, start, end, alt.count.t, ref.count.t, alt.count.n, ref.count.n
#' @param purity cellularity value of the sample
#' @param ploidy ploidy value of the sample (segment length weighted copy number)
#' @param pp.method character select from "ppgrid" and "sequenza" to infer purity and ploidy if not both given. Default, "sequenza".
#' @param min.nbins integer minimum number of bins of coverage of a segment
#' @param max.na numeric between (0, 1), any vertex with more than this proportion missing coverage data is ignored
#' @param blacklist.coverage GRanges marking regions of the genome where coverage is unreliable
#' @param cn.signif numeric (0, 1), significance level of CN change point when seg is not given (the larger the more sensitive),
#' `alpha` parameter in DNAcopy::segment, default 1E-5
#' 
#' @param slack.penalty numeric penalty to put on every loose end (or copy number if loose.penalty.mode is "linear"). Default 100.
#' @param loose.penalty.mode either \code{"linear"} or \code{"boolean"}, for penalizing each copy or each count of a loose end
#' @param tilim integer time limit MIP solver on each subgraph. Default 2400 (seconds).
#' @param max.threads maximum thread number CPLEX/Gurobi is allowed to use
#' @param max.mem maximum memory CPLEX/Gurobi is allowed to use
#' @param epgap relative optimality gap tolerated by the solver
#' @param use.gurobi logical flag specifying whether to use gurobi (if TRUE) instead of CPLEX (if FALSE). Default FALSE.
#' @param outdir  out directory to dump into, default `./`
#' @param name  sample name
#' @param mc.cores integer how many cores to use to fork subgraphs generation. Default 1. CAUTION as more cores will multiply the number of max threads used by CPLEX.
#' @param overwrite  logical flag whether to overwrite existing output directory contents or just continue with existing files.
#' @param verbose logical whether to print out the most verbose version of progress log
#' @param init jabba object (list) or path to .rds file containing previous jabba object which to use to initialize solution, this object needs to have the identical aberrant junctions as the current jabba object (but may have different segments and loose ends, i.e. is from a previous iteration)
#' @param dyn.tuning logical whether to let JaBbA dynamically tune the convergence criteria, default TRUE
#' @param lp (logical) whether to run as linear program, default FALSE
#' @param ism (logical) wehther to add ism constraints. default FALSE. only used if lp = TRUE
#' @param fix.thres (numeric) freeze the CN of large nodes with cost penalty exceeding this multiple of lambda. default -1 (no nodes are fixed)
#' @param min.bins (numeric) minimum number of bins needed for a valid segment CN estimate (default 5)
#' @param filter_loose (logical) run loose end annotation? (default FALSE)
#' @export
JaBbA = function(## Two required inputs
                 junctions,
                 coverage,
                 ## options about junctions
                 juncs.uf = NULL,
                 blacklist.junctions = NULL,
                 whitelist.junctions = NULL,
                 geno = FALSE,
                 indel = NULL,
                 cfield = NULL, 
                 tfield = "tier",
                 reiterate = 0,
                 rescue.window = 1e3,
                 rescue.all = TRUE,
                 nudge.balanced = FALSE,
                 thresh.balanced = 500,
                 edgenudge = 0.1,
                 strict = FALSE,
                 all.in = FALSE,
                 ## options about coverage
                 field = 'ratio',
                 seg = NULL,
                 max.na = -1,
                 blacklist.coverage = NULL,
                 nseg = NULL,
                 hets = NULL,
                 purity = NA,
                 ploidy = NA,
                 pp.method = "sequenza",
                 min.nbins = 5,
                 cn.signif = 1e-5,
                 ## options about optimization
                 slack.penalty = 1e2,
                 loose.penalty.mode = "boolean",
                 tilim = 6000,
                 max.threads = Inf,
                 max.mem = 16,
                 epgap = 1e-4,
                 use.gurobi = FALSE,
                 ## options about general pipeline
                 outdir = './',
                 name = 'tumor',
                 mc.cores = 1,
                 overwrite = FALSE,
                 verbose = TRUE,
                 init = NULL,
                 dyn.tuning = TRUE,
                 lp = FALSE,
                 ism = FALSE,
                 fix.thres = -1,
                 min.bins = 5,
                 filter_loose = FALSE)
{
    system(paste('mkdir -p', outdir))
    jmessage('Starting analysis in ', outdir <- normalizePath(outdir))
    cdir = normalizePath(getwd())
    setwd(outdir)
    if (overwrite)
        jmessage('Overwriting previous analysis contents')
    ra = junctions
    reiterate = 1 + as.integer(reiterate)

    ## if (is.character(ra))
    ## {
    ##     if (!file.exists(ra))
    ##     {
    ##         jerror(paste('Junction path', ra, 'does not exist'))
    ##     }
    ra.all = read.junctions(ra, geno = geno) ## GRL

    if (is.null(ra.all)){
        jwarning("no junction file is given, will be running JaBbA without junctions!")
        ra.all = GRangesList()
    }

    ## } else if (is.null(ra)) {
    ##     ra.all = GRangesList()
    ## } else {
    ##     ra.all = ra
    ## }

    if (inherits(ra.all, "list") && all(sapply(ra.all, inherits, "GRangesList"))){
        ## this is a multisample junction input
        ids = names(ra.all)
        match.nm = which(grepl(name, ids))
        if (length(match.nm)==1){
            ra.all = ra.all[[match.nm]]            
        } else if (length(match.nm)==0){
            jerror("There's no junction matching this sample name: ", name)
        } else {
            jwarning("there are more than one sample id in the junciton input ",
                     ids[match.nm],
                     " match this run id ",
                     name,
                     " we are only using the first one.")
            ra.all = ra.all[[match.nm[1]]]
        }
    }
    
    if (!inherits(ra.all, "GRangesList")){
        jerror("The given input `ra` is not valid.")
    }

    ## temporary filter for any negative coords
    ## bad.bp = unname(grl.unlist(ra.all)) %Q% (start<0)
    ## if (length(bad.bp)>0){
    ##     jmessage("Warning!! ", length(bad.bp), " breakpoints in ",
    ##              length(bad.ix <- unique(bad.bp$grl.ix)), " junctions ",
    ##              "have negative coordinates, discard.")
    ##     ra.all = ra.all[setdiff(seq_along(ra.all), bad.ix)]
    ## }
    
    if (verbose)
    {
        jmessage("Read in ", length(ra.all), " total input junctions")
    }

    ## xtYao Tuesday, Jun 19, 2018 04:52:17 PM
    ## Only when tier exists or unfiltered junctions provided, do we do the iterations
    ## if unfiltered set is given first parse it
    ra.uf = NULL
    if (!is.null(juncs.uf)){
        jmessage("Loading supplementary junctions")
        ra.uf = read.junctions(juncs.uf, geno=FALSE)
    }

    if (is.null(tfield)){
        tfield = 'tier'
    }

    ## if no tier field in junctions, set all of them to 2
    if (!(tfield %in% names(values(ra.all))))
    {
        jwarning("Tier field ", tfield, " missing: giving every junction the same tier, i.e. all have the potential to be incorporated")
        values(ra.all)[, tfield] = rep(2, length.out = length(ra.all))
    }

    
    if (!is.null(ra.uf)){
        ## merge ra.all with ra.uf
        ## junctions from ra.all will always have tier 2
        ## junctions from ra.uf will always have tier 3
        ## at this point
        ## ra.all.uf = ra.merge(ra.all, ra.uf, pad=0, ind=TRUE) ## hard merge
        ra.all.uf = ra.merge(ra.all, ra.uf, pad=1000, ind=TRUE) ## soft merge
        ## those match a record in junction, will be assigned to the tier in junction
        if (tfield %in% names(values(ra.all.uf))) {
            values(ra.all.uf)[, tfield][which(!is.na(values(ra.all.uf)$seen.by.ra1))] =
                values(ra.all)[, tfield][values(ra.all.uf)$seen.by.ra1]
            values(ra.all.uf)[, tfield][which(is.na(values(ra.all.uf)$seen.by.ra1))] = 3
        } else {
            values(ra.all.uf)[, tfield] = 3
        }

        ## mark any NA tier junctions as tier 3
        values(ra.all.uf)[, tfield][which(is.na(values(ra.all.uf)[, tfield]))] = 3
        
        ## the rest will be tier 3
        ra.all = ra.all.uf
        ## FIXME: ra.merge still gives duplicates
        ## FIXME: substitute these ra.xx fuctions to Junction class in gGnome
        ## extra dedup
        ## ndup = which(!ra.duplicated(ra.all.uf))
        ## ra.all = ra.all.uf[ndup]
    }

    if (!is.null(blacklist.junctions)){2
        blacklist.junctions = read.junctions(blacklist.junctions)
        if (length(blacklist.junctions)>0){
            black.ix = which(gUtils::grl.in(ra.all, blacklist.junctions))
            if (length(black.ix)>0){
                jmessage("Removing ", length(black.ix), " junctions matched the blacklist")
                ra.all = ra.all[setdiff(seq_along(ra.all), black.ix)]
            }
        }
    }

    if (!is.null(whitelist.junctions)){
        whitelist.junctions = read.junctions(whitelist.junctions)
        if (length(whitelist.junctions)>0){
            white.ix = which(gUtils::grl.in(ra.all, whitelist.junctions))
            if (length(white.ix)>0){
                jmessage("Forcing incorporation of ", length(white.ix), " junctions matched the whitelist")
                values(ra.all)[white.ix, tfield] = 1
            }
        }
    }

    if (length(ra.all)>0){
        ## final sanity check before running
        if (!all(unique(values(ra.all)[, tfield]) %in% 1:3)){
            jerror('Tiers in tfield can only have values 1,2,or 3')
        }
        jmessage("In the input, There are ", sum(values(ra.all)[, tfield]==1), " tier 1 junctions; ",
                 sum(values(ra.all)[, tfield]==2), " tier 2 junctions; ",
                 sum(values(ra.all)[, tfield]==3), " tier 3 junctions.")
    }


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

        ## and then bump t2 to t1
        t2 = values(ra.all)[, tfield]==2
        if (any(t2)){
            values(ra.all)[t2, tfield] = 1
            jmessage("All-in mode: ", length(t2),
                     " tier 2 junctions forced into the model")
        }
    }

    ## if we are iterating more than once
    if (reiterate>1){

        ## important: rescue.all should always be TRUE if not running filter.loose
        if ((!rescue.all) & (!filter_loose)) {
            jmessage("Resetting rescue.all to TRUE as filter.loose is FALSE")
            rescue.all = TRUE
        }
        
        continue = TRUE
        this.iter = 1;

        values(ra.all)$id = seq_along(ra.all)
        saveRDS(ra.all, paste(outdir, '/junctions.all.rds', sep = ''))

        last.ra = ra.all[values(ra.all)[, tfield]<3]

        jmessage('Starting JaBbA iterative with ', length(last.ra), ' junctions')
        jmessage('Will progressively add junctions within ', rescue.window, 'bp of a loose end in prev iter')
        jmessage('Iterating for max ', reiterate, ' iterations or until convergence (i.e. no new junctions added)')

        while (continue) {
            gc()

            this.iter.dir = paste(outdir, '/iteration', this.iter, sep = '')
            system(paste('mkdir -p', this.iter.dir))

            jmessage('Starting iteration ', this.iter, ' in ', this.iter.dir, ' using ', length(last.ra), ' junctions')

            if (this.iter>1){
                kag1 = readRDS(paste0(outdir, '/iteration1/karyograph.rds'))
                ploidy = kag1$ploidy
                purity = kag1$purity
                jmessage("Using ploidy ", ploidy,
                         " and purity ", purity,
                         " consistent with the initial iteration")

                if (lp) {
                    jmessage("Using segments from JaBbA output of previous iteration")
                    loose.ends.fn = file.path(outdir,
                                              paste0("iteration", this.iter - 1),
                                              "loose.end.stats.rds")
                    seg.fn = file.path(outdir, paste0("iteration", this.iter - 1), "jabba.simple.rds")

                    seg = readRDS(seg.fn)$segstats[, c()]
                    seg = seg %Q% (strand(seg) == "+")
                    seg = gr.stripstrand(seg)
                }
                
            }

            this.ra.file = paste(this.iter.dir, '/junctions.rds', sep = '')
            saveRDS(last.ra, this.ra.file)

            jab = jabba_stub(
                junctions = this.ra.file,
                seg = seg,
                coverage = coverage,
                blacklist.coverage = blacklist.coverage,
                hets = hets,
                nseg = nseg,
                cfield = cfield,
                tfield = tfield,
                nudge.balanced = as.logical(nudge.balanced),
                outdir = this.iter.dir,
                mc.cores = as.numeric(mc.cores),
                max.threads = as.numeric(max.threads),
                max.mem = as.numeric(max.mem),
                max.na = max.na,
                edgenudge = as.numeric(edgenudge),
                tilim = as.numeric(tilim),
                strict = strict,
                name = name,
                use.gurobi = as.logical(use.gurobi),
                field = field,
                epgap = epgap,
                ## subsample = subsample,
                slack.penalty = as.numeric(slack.penalty),
                loose.penalty.mode = loose.penalty.mode,
                mipstart = init,
                ploidy = as.numeric(ploidy),
                purity = as.numeric(purity),
                pp.method = pp.method,
                indel = indel,
                min.nbins = min.nbins,
                overwrite = as.logical(overwrite),
                verbose = as.numeric(verbose),
                dyn.tuning = dyn.tuning,
                geno = geno,
                cn.signif = cn.signif,
                lp = lp,
                ism = ism,
                fix.thres = fix.thres,
                min.bins = min.bins,
                filter_loose = filter_loose)
            
            gc()

            jab = readRDS(paste(this.iter.dir, '/jabba.simple.rds', sep = ''))
            jabr = readRDS(paste(this.iter.dir, '/jabba.raw.rds', sep = ''))
            le = gr.stripstrand(jab$segstats %Q% (loose==TRUE & strand=="+"))
            if (length(le)==0){
                jmessage("No more loose ends to resolve, terminating.")
                break
            } else if (!rescue.all){
                le = le %Q% which(passed==TRUE)
                if (length(le)==0){
                    jmessage("No more plausible loose ends, terminating")
                    break
                }
            } else {
                jmessage("Rescuing all ", length(le), " loose ends, regardless of confidence.")
            }
            
            ## determine orientation of loose ends
            le.right = le %&% gr.start(jab$segstats %Q% (loose==FALSE))
            strand(le.right) = "+"
            le.left = le %&% gr.end(jab$segstats %Q% (loose==FALSE))
            strand(le.left) = "-"
            le = grbind(le.right, le.left)

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
            ## got used, stay there
            ## but not loose ends overlapping an exorbitant number of junctions
            le.keep = which((le %N% (stack(ra.all) + rescue.window)) < 20)
            tokeep = which(values(jab$junctions)$cn>0) 
            new.ra.id = unique(c(
                values(jab$junctions)$id[tokeep],
                ## near a loose ends, got another chance
                values(ra.all)$id[which(grl.in(ra.all,
                                               le[le.keep] + rescue.window,
                                               some = T,
                                               ignore.strand = FALSE))],
                ## tier 2 or higher must stay for all iterations
                values(ra.all)$id[which(values(ra.all)$tier==2)]
            )) 
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
            ## keep using the initial purity ploidy values
            pp1 = readRDS(paste0(
                outdir,
                '/iteration1/karyograph.rds.ppgrid.solutions.rds'))
            purity = pp1$purity[1]
            ploidy = pp1$ploidy[1]

            seg = readRDS(paste0(outdir,'/iteration1/seg.rds')) ## read from the first iteration

            
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
        jab = readRDS(paste0(outdir, "/jabba.simple.gg.rds"))
        jmessage('Done Iterating')
    } else {
        ## if all.in, convert all tier 3 to tier 2
        ## if (tfield %in% colnames(values(ra.all))){
        ##     t3 = (values(ra.all)[, tfield] == 3)
        ##     if (all.in & length(ra.all)>0){
        ##         if (any(t3)){
        ##             ## save every t3 except small indel
        ##             t3.indel = which.indel(ra.all[which(t3)])
        ##             t3.non.indel = which(t3)[setdiff(seq_along(which(t3)), t3.indel)]
        ##             values(ra.all)[t3.non.indel, tfield] = 2
        ##             t3 = values(ra.all)[, tfield] == 3
        ##         }
        ##     }
        ##     ## if not all.in, only use t2 or t1
        ##     ra.all = ra.all[setdiff(seq_along(ra.all), which(t3))]
        ## }
        jab = jabba_stub(
            junctions = ra.all,
            seg = seg,
            coverage = coverage,
            blacklist.coverage = blacklist.coverage,
            hets = hets,
            nseg = nseg,
            cfield = cfield,
            tfield = tfield,
            nudge.balanced = as.logical(nudge.balanced),
            outdir = outdir,
            mc.cores = as.numeric(mc.cores),
            max.threads = as.numeric(max.threads),
            max.mem = as.numeric(max.mem),
            max.na = max.na,
            edgenudge = as.numeric(edgenudge),
            tilim = as.numeric(tilim),
            strict = strict,
            epgap = epgap,
            name = name,
            use.gurobi = as.logical(use.gurobi),
            field = field,
            ## subsample = subsample,
            slack.penalty = as.numeric(slack.penalty),
            mipstart = init,
            ploidy = as.numeric(ploidy),
            purity = as.numeric(purity),
            pp.method = pp.method,
            indel = indel,
            min.nbins = min.nbins,
            loose.penalty.mode = loose.penalty.mode,
            overwrite = as.logical(overwrite),
            verbose = as.numeric(verbose),
            dyn.tuning = dyn.tuning,
            geno = geno,
            cn.signif = cn.signif,
            lp = lp,
            ism = ism,
            fix.thres = fix.thres,
            min.bins = min.bins,
            filter_loose = filter_loose)
    }
    setwd(cdir)
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
#' ## @param subsample  numeric between 0 and 1 specifying how much to sub-sample high confidence coverage data
#' @param tilim  timeout for jbaMIP computation (default 1200 seconds)
#' @param edgenudge  numeric hyper-parameter of how much to nudge or reward aberrant junction incorporation, default 0.1 (should be several orders of magnitude lower than average 1/sd on individual segments), a nonzero value encourages incorporation of perfectly balanced rearrangements which would be equivalently optimal with 0 copies or more copies.
#' @param slack.penalty  penalty to put on every loose.end copy, should be calibrated with respect to 1/(k*sd)^2 for each segment, i.e. that we are comfortable with junction balance constraints introducing k copy number deviation from a segments MLE copy number assignment (the assignment in the absence of junction balance constraints)
#' @param init jabba object (list) or path to .rds file containing previous jabba object which to use to initialize solution, this object needs to have the identical aberrant junctions as the current jabba object (but may have different segments and loose ends, i.e. is from a previous iteration)
#' @param overwrite  flag whether to overwrite existing output directory contents or just continue with existing files.
#' @param lp whether to run as linear program (default FALSE)
#' @param ism logical whether to add ISM constraints default FALSE
#' @param fix.thres (numeric) multiple of lambda above which to fix nodes
#' @param min.bins (numeric) min number of coverage bins for a valid CN estimate
#' @param filter_loose (logical) run loose end analysis?
jabba_stub = function(junctions, # path to junction VCF file, dRanger txt file or rds of GRangesList of junctions (with strands oriented pointing AWAY from breakpoint)
                      coverage, # path to cov file, rds of GRanges
                      blacklist.coverage = NULL,
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
                      max.na = -1,
                      purity = NA,
                      ploidy = NA,
                      pp.method = "sequenza",
                      strict = FALSE,
                      epgap = 1e-4,
                      mipstart = NULL,
                      field = 'ratio', ## character, meta data field to use from coverage object to indicate numeric coveragendance, coverage,
                      ## subsample = NULL, ## numeric scalar between 0 and 1, how much to subsample coverage per segment
                      tilim = 1200, ## timeout for MIP portion
                      ## mem = 16, ## max memory for MIP portion
                      init = NULL, ## previous JaBbA object to use as a solution
                      edgenudge = 0.1, ## hyper-parameter of how much to "nudge" or reward edge use, will be combined with cfield information if provided
                      slack.penalty = 1e2, ## nll penalty for each loose end cop
                      use.gurobi = FALSE,
                      loose.penalty.mode = "boolean",
                      indel = "exclude", ## default should be nothing
                      min.nbins = 5, ## default to 5, if larger bin can reduce
                      overwrite = F, ## whether to overwrite existing output in outdir
                      lp = FALSE, ## whether to run as linear program
                      ism = FALSE, ## add ISM constraints (only used if lp = TRUE)
                      fix.thres = -1,
                      min.bins = 5,
                      verbose = TRUE,
                      dyn.tuning = TRUE,
                      geno = FALSE,
                      filter_loose = FALSE,
                      outlier.thresh = 0.9999,
                      cn.signif = 1e-5)
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
    purity.ploidy.txt.file = paste(outdir, 'purity.ploidy.txt', sep = '/')
    jabba.png.file = paste(outdir, 'jabba.png', sep = '/')
    jabba.simple.png.file = paste(outdir, 'jabba.simple.png', sep = '/')
    seg.tab.file = paste(outdir, 'jabba.seg', sep = '/')
    seg.gr.file = paste(outdir, 'jabba.seg.rds', sep = '/')
    seg.adj.file = paste(outdir, 'jabba.adj.txt', sep = '/')
    le.class.file.rds = paste(outdir, 'loose.end.stats.rds', sep = '/')
    nozzle.file = paste(outdir, 'nozzle', sep = '/')

    if (is.character(coverage))
    {
        if (!file.exists(coverage))
        {
            jerror(paste('Coverage path ', coverage, 'does not exist'))
        }

        if (!file.size(coverage))
        {
            jerror(paste('Coverage file ', coverage, 'is empty'))
        }

        if (grepl('\\.rds$', coverage))
        {
            coverage = readRDS(coverage)
        }
        else if (grepl('(\\.txt$)|(\\.tsv$)|(\\.csv$)', coverage))
        {
            tmp = fread(coverage)
            coverage = try(dt2gr(tmp))
            if (inherits(tmp, "try-error")){
                jerror("Input coverage data ", coverage, " cannot be converted to genomic intervals.")
            }
            ## coverage = dt2gr(tmp, seqlengths = tmp[, max(end), by = seqnames][, structure(V1, names = seqnames)])
        }
        else
        {
            jmessage('Importing coverage from UCSC format')
            coverage = rtracklayer::import(coverage)
            names(values(coverage)) = field ## reset name of coverage field
            coverage = gr.fix(coverage)
        }
    }
    else
        coverage = coverage

    if (!inherits(coverage, 'GRanges'))
        coverage = dt2gr(coverage)

    ## the most frequent width in a sample of coverage points
    binwidth = as.numeric(names(sort(
        table(width(sample(coverage, 1000, replace=TRUE))), decreasing = TRUE
    )[1]))
    
    if (verbose)
    {
        jmessage(paste(
            "Read in",
            prettyNum(length(coverage), big.mark = ','),
            paste0(binwidth, "bp bins of coverage data across"),
            length(unique(seqnames(coverage))), 'chromosomes'))
    }

    if (!(field %in% names(values(coverage))))
    {
        new.field = names(values(coverage))[1]
        jwarning(paste0('Field ',
                        field,
                        ' not found in coverage GRanges metadata so using ',
                        new.field,
                        ' instead'))
        field = new.field
    }

    ## filter out the data in blacklist.coverage if any
    if (!is.null(blacklist.coverage)){
        if (is.character(blacklist.coverage)){
            if (file.exists(blacklist.coverage)){
                if (grepl("rds$", blacklist.coverage)){
                    blacklist.coverage = readRDS(blacklist.coverage)
                } else if (grepl("[(txt)|(tsv)|(csv)]$", blacklist.coverage)){
                    blacklist.coverage = try(gUtils::dt2gr(data.table::fread(blacklist.coverage)))
                } else if (grepl("bed$", blacklist.coverage)){
                    blacklist.coverage = rtracklayer::import.bed(blacklist.coverage)
                }                 
            }
        } else if (inherits(blacklist.coverage, "data.frame")){
            blacklist.coverage = try(gUtils::dt2gr(data.table(blacklist.coverage)))
        }
        ## if at this point we have a GRanges, then proceed
        if (inherits(blacklist.coverage, "GRanges")){
            bad.ix = which(coverage %^% blacklist.coverage)
            if (length(bad.ix)>0){
                values(coverage)[bad.ix, field] = NA
            }
        } else {
            jwarning("'--blacklist.coverage' cannot be parsed into a GRanges, please check")
        }
    }

    if (!is.null(outlier.thresh)) {
        field.thresh = quantile(values(coverage)[[field]], probs = outlier.thresh, na.rm = TRUE)
        bad.ix = which(values(coverage)[[field]] > field.thresh)
        values(coverage)[bad.ix, field] = NA
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
        if (is.null(seg) || (is.character(seg) && (!file.exists(seg) || !file.size(seg))))
        {
            if (verbose)
            {
                jmessage('No segmentation provided, so performing segmentation using CBS')
            }
            set.seed(42)
            binw = median(sample(width(coverage), 30), na.rm=T)
            vals = as.double(values(coverage)[, field])
            ## if `vals` contain minute values that are basically deemed as zero in R
            ## we add a tiny little bit of non-zero value to all of them before segmenting
            if (length(zero.ix <- which(vals==0))>0){
                tiny.val = .Machine$double.eps
                vals = vals + tiny.val
                jmessage(
                    length(zero.ix),
                    " coverage data points have zero value, adding a tiny value ",
                    tiny.val, " to prevent log error.")
            }
            new.sl = GenomeInfoDb::seqlengths(coverage)
            ix = which(!is.na(vals))

            cna = DNAcopy::CNA(
                log(vals[ix]),
                as.character(seqnames(coverage))[ix], start(coverage)[ix],
                data.type = 'logratio')
            ## Addy at some point suggested 1e-5 to prevent way too many unnecessary segmentations at 200bp resolution
            segs = DNAcopy::segment(DNAcopy::smooth.CNA(cna),
                                    alpha = cn.signif, ## old at 1e-5
                                    verbose = FALSE)
            if (verbose)
            {
                jmessage('Segmentation finished')
            }
            seg = gUtils::seg2gr(segs$out, new.sl) ## remove seqlengths that have not been segmented
            seg = gr.fix(seg, GenomeInfoDb::seqlengths(coverage), drop = T)
            names(seg) = NULL

            ## Filter out small gaps in CBS output (<=1e3)
            ## gap.seg = IRanges::gaps(seg)[which(as.character(strand(seg))=="*" & width(seg) > (5 * binw))]
            ## if (length(gap.seg)>0){
            ##     bps = c(gr.start(seg)[, c()], gr.start(gap.seg)[, c()])
            ## } else {
            ##     bps = gr.start(seg)[, c()]
            ## }
            
            ## ## keep the breakpoints of the big enough gaps (>10*binwidth), they may contain bad regions
            ## new.segs = gUtils::gr.stripstrand(gUtils::gr.breaks(bps, gUtils::si2gr(seqlengths(bps))))[, c()]
            ## names(new.segs) = NULL
            ## seg = gUtils::gr.fix(new.segs, GenomeInfoDb::seqlengths(coverage), drop = T)

            ## if (verbose)
            ## {
            ##     jmessage(length(seg), ' segments produced')
            ## }
        }
        else
        {
            if (is.character(seg))
            {                
                if (grepl('\\.rds$', seg))
                {
                    seg = readRDS(seg)
                } else if (grepl("\\.txt(.gz)?$", seg)){
                    seg = dt2gr(fread(seg))
                }
                else
                {
                    jmessage('Importing seg from UCSC format')
                    seg = rtracklayer::import(seg)
                    ## field = 'score';
                }
            }
        }
    }

    if (!inherits(seg, 'GRanges'))
        seg = dt2gr(seg, GenomeInfoDb::seqlengths(coverage))

    if (!is.null(hets))
    {
        if (is.character(hets))
        {
            if (!file.exists(hets))
            {
                jwarning(sprintf('hets file "%s" not found, ignoring hets\n', hets))
                hets = NULL
            }
        }
    }

    ra = junctions
    if (inherits(ra, "character")){
        ra = read.junctions(junctions, geno = geno)
    } else if (!inherits(ra, "GRangesList")){
        jerror("`ra` must be GRangesList here")
    }
    jmessage(paste("Loaded", length(ra), "junctions from the input."))

    if (strict)
    {
        ends = c(gr.start(seg),
                 gr.end(gr.end(seg)))
        nog = length(ra)
        ra = ra[grl.in(ra, ends, logical, logical = FALSE) == 2]
        jmessage('Applying strict junction filtering of junctions to only those that land at segment ends')
        jmessage('Leaving ', length(ra), ' junctions from an initial set of ', nog)
    }

    ## determine binwidth
    ## times the factor of minimum number of bins for a segment to be considered
    max.indel.size = min.nbins * binwidth
    ## know what junctions to include and exclude before karyogrpah_stub
    ab.force = NULL
    ab.exclude = NULL
    if (!is.null(tfield))
    {
        if (tfield %in% names(values(ra)))
        {
            ## default forcing tier 1
            ab.force = which(gsub('tier', '', as.character(values(ra)[, tfield]))=='1')

            ## when the switch is on:
            ## small DUP/DEL junctions in tier 2, bump them up to tier 1
            tier2.ix = which(gsub('tier', '', as.character(values(ra)[, tfield]))=='2')
            if (length(tier2.ix)>0){
                ## max size hardcoded for now
                which.like.indel = tier2.ix[
                    which.indel(ra[tier2.ix], max.size = max.indel.size)]
            } else {
                which.like.indel = numeric(0)
            }

            if (verbose)
            {
                jmessage('Found tier field enforcing >=1 CN at ', length(ab.force), ' junctions')                
            }

            ab.exclude = which(gsub('tier', '', as.character(values(ra)[, tfield]))=='3')
            if (verbose)
            {
                jmessage('Removing ', length(ab.exclude), ' tier 3 junctions')
            }
        }
    } else { ## no tfield given, assume everything is tier 2
        ## if (indel){
        which.like.indel = which.indel(ra, max.size = max.indel.size)
        ##     ab.force = union(ab.force, which.like.indel)
        ## }
    }

    if (exists("which.like.indel")){
        values(ra)$like.indel = is.element(seq_along(ra), which.like.indel)
        jmessage(length(which.like.indel), ' INDEL-like isolated events')
    }
    ## depending on the value of "indel", include or exclude indels or do nothing
    if (!is.null(indel)){
        if (indel=="include"){
            ab.force = union(ab.force, which.like.indel)
            if (verbose){
                jmessage(length(which.like.indel), ' INDEL-like isolated junctions mandatorily INCLUDED in the final model')
            }
        } else if (indel=="exclude"){
            ab.exclude = union(ab.exclude, which.like.indel)
            if (verbose){
                jmessage(length(which.like.indel), ' INDEL-like isolated junctions mandatorily EXCLUDED in the final model')
            }
        }
    } else {
        jwarning("doing nothing special to the small INDEL-like isolated junctions")
    }
    
    ## clean up the seqlevels before moving on
    seg.sl = seqlengths(seg)
    cov.sl = seqlengths(coverage)
    ra.sl = seqlengths(ra)
    union.sn = union(union(names(seg.sl), names(cov.sl)), names(ra.sl))
    union.sl = data.table(
        seqnames = union.sn,
        seg.sl = seg.sl[union.sn],
        cov.sl = cov.sl[union.sn],
        ra.sl = ra.sl[union.sn])
    ## keep the largest length so we don't lose any data
    ## it is more union for a certain data type to have shorter genome length than having longer
    union.sl = union.sl[
      , sl := pmax(seg.sl, cov.sl, ra.sl, na.rm = T)][
      , setNames(sl, seqnames)]
    new.seg = seg[which(is.element(as.character(seqnames(seg)), names(union.sl)))]
    if (length(new.seg)<length(seg)){
        jmessage(length(seg)-length(new.seg), " segments are discarded because they fall out of the ref genome.")
    }
    seg = new.seg
    seqlevels(seg) = seqlevels(seg)[which(is.element(seqlevels(seg), names(union.sl)))]
    new.coverage = coverage[which(is.element(as.character(seqnames(coverage)), names(union.sl)))]
    seqlevels(coverage) = seqlevels(coverage)[
        which(is.element(seqlevels(coverage), names(union.sl)))]
    if (length(new.coverage)<length(coverage)){
        jmessage(length(coverage)-length(new.coverage), " coverage points are discarded because they fall out of the ref genome.")
    }
    coverage = new.coverage
    tmp = grl.unlist(ra) 
    tmp.md = values(ra)
    nms = names(tmp)
    names(tmp) = NULL
    good.tmp = which(as.character(seqnames(tmp)) %in% names(union.sl))
    tmp = tmp[good.tmp]
    names(tmp) = nms[good.tmp]
    ra = split(tmp[, c()], tmp$grl.ix)
    values(ra) = tmp.md[as.numeric(names(ra)),]
    ## in case some junction lost one breakpoint
    intact.ix = which(elementNROWS(ra)==2)
    new.ra = ra[intact.ix]
    if (length(new.ra)<length(ra)){
        jmessage(length(ra)-length(new.ra), " rearrangements are discarded because they fall out of the ref genome.")
        ## at this point, there's got to be at least 1 coverage point to start the program
        if (length(coverage)==0){
            jerror("Empty coverage data. Please check if their reference chromsome name match the other inputs.")
        }
    }
    ra = new.ra
    seqlevels(ra) = seqlevels(ra)[which(is.element(seqlevels(ra), names(union.sl)))]
    jmessage("Conform the reference sequence length of: seg, coverage, and ra, to be: \n",
             paste0("\t", names(union.sl), ":", union.sl, collapse = "\n"))

    
    ## extend the end and start of the segmentation to the chromosome endpoints
    bp.dt = as.data.table(seg)[as.character(seqnames) %in% names(union.sl)]
    bp.dt[, first := (start == min(.SD$start, na.rm = TRUE)), by = seqnames]
    bp.dt[, last := (end == max(.SD$end, na.rm = TRUE)), by = seqnames]
    bp.dt[first == TRUE, start := 1]
    bp.dt[last == TRUE, end := union.sl[seqnames]]

    ## keep only valid junctions after this seqlength normalization process
    bp.dt = bp.dt[end > start,]

    ## convert back to GRanges
    seg = dt2gr(bp.dt[, .(seqnames, start, end, strand = "*")], seqlengths = union.sl)

    ## xtYao #' Wednesday, Jun 02, 2021 02:43:35 PM
    ## Move the gap filling step after the correction of seqlengths
    ## zchoo Friday, Apr 23, 2021 10:39:19 AM
    ## make gap filtering a general preprocessing step
    ## filter small gaps between segments containing less than ten bins
    binw = median(sample(width(coverage), 30), na.rm = TRUE)
    all.gaps = IRanges::gaps(seg)
    gap.seg = all.gaps[which(as.character(strand(all.gaps)) == "*" & width(all.gaps) > 5 * binw)]

    if (verbose) {
        n.gaps = sum(width(all.gaps) > 0, na.rm = TRUE)
        jmessage("Number of gaps with nonzero width: ", n.gaps)
        jmessage("Number of segments before gap filtering: ", n.gaps + length(seg))
    }

    if (length(gap.seg)>0){
        bps = c(gr.start(seg)[, c()], gr.start(gap.seg)[, c()])
    } else {
        bps = gr.start(seg)[, c()]
    }

    ## create new segments from the breakpoints of the old segments plus big gaps
    new.segs = gUtils::gr.stripstrand(gUtils::gr.breaks(bps, gUtils::si2gr(seqlengths(bps))))[, c()]
    names(new.segs) = NULL
    seg = gUtils::gr.fix(new.segs, union.sl, drop = T)

    if (verbose)
    {
        jmessage(length(seg), ' segments produced after gap filtering')
    }

    saveRDS(seg, seg.fn)

    
    if (overwrite | !file.exists(kag.file)){
        karyograph_stub(seg,
                        coverage,
                        ra = ra,
                        out.file = kag.file,
                        nseg.file = nseg,
                        field = field,
                        purity = purity,
                        ploidy = ploidy,
                        pp.method = pp.method,
                        ## subsample = subsample,
                        het.file = hets,
                        verbose = verbose,
                        ab.exclude = ab.exclude,
                        ab.force = ab.force,
                        max.na = max.na,
                        lp = lp)
    } else {
        jwarning("Skipping over karyograph creation because file already exists and overwrite = FALSE")
    }

    kag = readRDS(kag.file)
    ab.force = kag$ab.force
    ab.exclude = integer(0)

    if (!is.null(cfield))
    {
        if (cfield %in% names(values(kag$junctions)))
        {
            val = values(kag$junctions)[, cfield]
            val[is.na(val)] = 0
            edgenudge = val * edgenudge
        }
    }

    gc()

    juncs = kag$junctions ## already removed ab.exclude!!!
    bpss = grl.unlist(juncs)


    if (nudge.balanced) {
        balanced.jix = c()
        if (length(juncs)>0) {
            jmessage("Brand new function for reciprocal junctions calling.")
            balanced.jix = unlist(reciprocal.cycles(juncs, thresh = 1e3, mc.cores = mc.cores, verbose = verbose>1))
            dp.jix = which(gUtils::ra.duplicated(juncs, pad=1500))
            balanced.jix = setdiff(balanced.jix, dp.jix)
        }
        ## only adds edge nudge to the balanced junctions
        edgenudge = edgenudge * as.numeric(seq_along(juncs) %in% balanced.jix) 
    } else { ## nudge everything ..
        if (length(edgenudge)==1) edgenudge = rep(edgenudge, length(juncs))
        if (length(juncs)>0){   ## hot fix for preventing nudging of NA segments
            bps.cov = gr.val(bpss, coverage, val = 'ratio')
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
    nothing.contig = gr2dt(kag$segstats)[, list(nothing = all(is.na(mean))), by=seqnames][nothing==TRUE, seqnames]
    ## both breakpoints in NA regions
    if (length(juncs)>0){
        junc.dt = data.table(data.frame(values(juncs)))
        junc.dt[, ":="(from = NULL, to = NULL)]
        junc.dt = cbind(junc.dt,
                        data.table(matrix(kag$ab.edges[,,1],
                                          nrow=nrow(kag$ab.edges),
                                          dimnames=dimnames(kag$ab.edges)[1:2])))
        junc.dt[!is.na(from) & !is.na(to), ":="(mean.a=kag$segstats$mean[from],
                       mean.b=kag$segstats$mean[to],
                       chr.a = as.character(seqnames(kag$segstats[from])),
                       chr.b = as.character(seqnames(kag$segstats[to])))]
        junc.dt[, both.na := is.na(mean.a) & is.na(mean.b)]
        both.na.ix = junc.dt[, which(both.na==TRUE)] ## both breakpoint in NA

        no.man.land = junc.dt[, which(chr.a %in% nothing.contig | chr.b %in% nothing.contig)]
        ## either breakpoint in a contig that's completely NA
        ## excluding those whose both bp in NA regions or mapped to completely NA contigs
        ab.exclude = union(ab.exclude, union(both.na.ix, no.man.land))
        ## ab.exclude = union(ab.exclude, both.na.ix)
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


    ## save the final indices of selected edges
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
                   mem = max.mem,
                   tilim = tilim,
                   edge.nudge = edgenudge,
                   use.gurobi = use.gurobi,
                   ab.force = ab.force,
                   ab.exclude = ab.exclude, ## we now exclude things during karyograph_stub
                   ## ab.exclude = integer(0), 
                   init = init,
                   verbose = verbose,
                   purity.min = purity,
                   mipstart = mipstart,
                   epgap = epgap,
                   purity.max = purity,
                   ploidy.min = ploidy,
                   ploidy.max = ploidy,
                   slack.prior = 1/slack.penalty,
                   loose.penalty.mode = loose.penalty.mode,
                   dyn.tuning = dyn.tuning,
                   lp = lp,
                   ism = ism,
                   tfield = tfield,
                   fix.thres = fix.thres,
                   min.bins = min.bins)
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
    seg.out = cbind(
        sample = name,
        gr2dt(jabd$segstats)[
            loose==FALSE & strand=="+",
            .(chr = seqnames,
              start, end, width, cn)])
    names(seg.out)[1:4] = c('track.name', 'chrom', 'start', 'end')
    ## seg.out$seg.id = 1:nrow(seg.out)
    ## cols = c('track.name', 'chrom', 'start', 'end', 'cn', 'seg.id')
    ## seg.out = seg.out[, c(cols, setdiff(names(seg.out), cols))]
    write.tab(seg.out, seg.tab.file)
    jabd$segstats$seg.id = seq_along(jabd$segstats)
    
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
                jmessage('Done computing alleles'); file.remove(hets.gr.rds.file)
            },
            error = function(e) print("Jabba allelic generation failed"))

    ## annotate convergence status
    ## jabd.simple$segstats = jabd.simple$segstats %$% jab$segstats[, c("cl")]
    ## jabd.simple$segstats = gr.val(jabd.simple$segstats, jab$segstats[, c("epgap")], mean = FALSE, FUN = max, na.rm = TRUE, verbose = TRUE)
    jabd.simple$segstats = jabd.simple$segstats %$% jab$segstats[, c("cl", "epgap")]
    ## opti = readRDS(paste0(outdir, "/opt.report.rds"))
    ## sapply(strsplit(head(jabd.simple$segstats$cl), ","), as.numeric)

    jab$segstats = gr.fix(jab$segstats)
    jabd$segstats = gr.fix(jabd$segstats)
    jabd.simple$segstats = gr.fix(jabd.simple$segstats)    

    ## dependency function: dflm
    .dflm = function(x, last = FALSE, nm = '')
    {
        if (is.null(x))
            out = data.frame(name = nm, method = as.character(NA), p = as.numeric(NA), estimate = as.numeric(NA), ci.lower = as.numeric(NA),  ci.upper = as.numeric(NA), effect = as.character(NA))
        else if (any(c('lm', 'betareg') %in% class(x)))
        {

            coef = as.data.frame(summary(x)$coefficients)
            colnames(coef) = c('estimate', 'se', 'stat', 'p')
            if (last)
                coef = coef[nrow(coef), ]
            coef$ci.lower = coef$estimate - 1.96*coef$se
            coef$ci.upper = coef$estimate + 1.96*coef$se
            if (!is.null(summary(x)$family))
            {
                fam = summary(x)$family$family
                if (summary(x)$family$link %in% c('log', 'logit'))
                {
                    coef$estimate = exp(coef$estimate)
                    coef$ci.upper= exp(coef$ci.upper)
                    coef$ci.lower= exp(coef$ci.lower)
                }
            }
            else
                fam = 'Unknown'

            if (!last)
                nm = paste(nm, rownames(coef))
            out = data.frame(name = nm, method = fam, p = signif(coef$p, 3), estimate = coef$estimate, ci.lower = coef$ci.lower, ci.upper = coef$ci.upper, effect = paste(signif(coef$estimate, 3), ' [',  signif(coef$ci.lower,3),'-', signif(coef$ci.upper, 3), ']', sep = ''))
        }
        else
        {
            ci.lower = ifelse(is.null(x$conf.int[1]), NA, x$conf.int[1])
            ci.upper = ifelse(is.null(x$conf.int[2]), NA, x$conf.int[2])
            out = data.table(name = nm,
                             method = x$method,
                             p = signif(x$p.value, 3),
                             estimate = ifelse(is.null(x$estimate), NA, x$estimate),
                             ci.lower,
                             ci.upper)
            ## FIXME: some model doesn't have `estimate` field
            if (!is.null(x$estimate)){
                out[, effect := paste0(signif(x$estimate, 3), 
                                       ' [',  signif(x$conf.int[1],3),
                                       '-', signif(x$conf.int[2], 3), ']')]
            } else {
                out[, effect := "error"]
            }
        }
        out$effect = as.character(out$effect)
        out$name = as.character(out$name)
        out$method = as.character(out$method)
        rownames(out) = NULL
        return(as.data.table(out))
    }

    if (filter_loose) {

        jmessage("Starting loose end annotation")

        ## start building the model
        ## gather loose ends from sample
        gg = gG(jabba = jabd)
        ll = gr2dt(gr.start(gg$nodes[!is.na(cn) & loose.cn.left>0]$gr))[, ":="(lcn = loose.cn.left, strand = "+")]
        lr = gr2dt(gr.end(gg$nodes[!is.na(cn) & loose.cn.right>0]$gr))[, ":="(lcn = loose.cn.right, strand = "-")]

        ## TODO
        ## arbitrarily defined based on Julie's paper and prelim TGCT 1kb dryclean data
        if (binwidth<1000){
            PTHRESH = 3.4e-7
        } else {
            PTHRESH = 2e-6
        }       

        if ((nrow(ll)+nrow(lr))>0){
            l = rbind(ll, lr)[, ":="(sample = name)] ## FIXME
            l[, leix := 1:.N]
            l = dt2gr(l)
            le.class = filter.loose(gg, cov = coverage, field = field, l = l, PTHRESH = PTHRESH, verbose = TRUE, max.epgap = epgap)
            saveRDS(le.class, le.class.file.rds)
            n.le = dt2gr(le.class)
            jabd.simple$segstats =
                grbind(
                    jabd.simple$segstats %Q% (loose==FALSE),
                    jabd.simple$segstats %Q% (loose==TRUE) %$% n.le[, c("passed", "true.pos")])
            jabd$segstats =
                grbind(
                    jabd$segstats %Q% (loose==FALSE),
                    jabd$segstats %Q% (loose==TRUE) %$% n.le[, c("passed", "true.pos")])
            jmessage("Loose end quality annotated")
        }
    } else {
        jmessage("Skipping loose end annotation")
    }


    
    if (overwrite | !file.exists(jabba.simple.rds.file))
    {
        jmessage("Saving results")
        saveRDS(jabd$segstats, seg.gr.file)
        saveRDS(jab, jabba.raw.rds.file)
        saveRDS(jabd, jabba.rds.file)
        saveRDS(jabd.simple, jabba.simple.rds.file)

        purity.ploidy.dt = data.table(purity = jab$purity,
                                      ploidy = jab$ploidy)
        fwrite(purity.ploidy.dt, purity.ploidy.txt.file, sep = '\t')

        jab.gg = gGnome::gGraph$new(jab = jabd)
        tmp.jabd.simple = jabd.simple
        values(tmp.jabd.simple$junctions)$cn = NULL
        jab.simple.gg = gGnome::gGraph$new(jab = tmp.jabd.simple)

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
        names(tmp[[1]]) = seq_along(kag$junctions)
        names(tmp[[2]]) = seq_along(kag$junctions)
        ra1 = as.data.frame(tmp[[1]])
        ra2 = as.data.frame(tmp[[2]])
        names(ra1) = paste('bp1_', names(ra1), sep = '')
        names(ra2) = paste('bp2_', names(ra2), sep = '')
        junc.txt = as.data.frame(values(kag$junctions))
        write.tab(cbind(ra1, ra2, junc.txt), junctions.txt.file)
    }
    else
        writeLines(c("\t"), junctions.txt.file)

    ## all junctions are based on what's left in the karyograph
    saveRDS(kag$junctions, junctions.rds.file)

    tmp.cov = sample(coverage, pmin(length(coverage), 5e5))
    tmp.cov = gr.fix(tmp.cov, jabd$segstats)
    # add ncn values
    tmp.cov = tmp.cov %$% kag$segstats[,'ncn']
    # transform using rel2abs
    tmp.cov$cn = rel2abs(tmp.cov, purity = jab$purity, ploidy = jab$ploidy,
                         field = field, field.ncn = 'ncn')

    y1 = pmax(5, max(jabd$segstats$cn)*1.1)
    jabd$gtrack$y1 = y1
    jabd.simple$gtrack$y1 = y1

    td.cov = gTrack(tmp.cov, y.field = 'cn', col = alpha('black', 0.2), name = 'Cov', y1 = (y1 + jab$gamma)/jab$beta)

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
            plotted = tryCatch(plot(c(td.cov, jabd$gtrack), links = jun), error = function(e) return(1))
        } else {
            plotted = tryCatch(plot(c(jabd$agtrack, td.cov, jabd$gtrack), links = jun), error = function(e) return(1))
        }

        if (!is.null(plotted)){
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
                               error = function(e) return(1))
        } else {
            plotted = tryCatch(plot(c(jabd.simple$agtrack, td.cov, jabd.simple$gtrack), links = jun),
                               error = function(e) return(1))
        }

        if (!is.null(plotted)){
            if (verbose){
                jmessage("Something wrong with plotting JaBbA simplified results. Please try it later.")
            }
        }

        dev.off()
    }

    ## annotate loose ends
    

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
                           pp.method = "sequenza",
                           ra.file = NULL,
                           verbose = FALSE,
                           force.seqlengths = NULL,
                           purity =  NA,
                           ploidy = NA,
                           field = 'ratio',
                           mc.cores = 1,
                           max.chunk = 1e8,
                           max.na = -1,
                           ## subsample = NULL,
                           ab.exclude = NULL,
                           ab.force = NULL,
                           lp = FALSE){
    loose.ends = GRanges()

    if (!is.null(ra)){
        this.ra = ra
    } else
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
            tmp.ra = read.junctions(ra.file, seqlengths = hg_seqlengths(), get.loose = T, geno = geno)
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

    ab.exclude = intersect(seq_along(this.ra), ab.exclude)
    ab.force = setdiff(intersect(seq_along(this.ra), ab.force), ab.exclude)
    remaining = setdiff(seq_along(this.ra), ab.exclude)
    this.ra = this.ra[remaining]
    ab.force = which(remaining %in% ab.force)
    
    ## if we don't have normal segments then coverage file will be our bible for seqlengths
    if (is.character(cov.file))
    {
        if (grepl('\\.rds$', cov.file))
            this.cov = readRDS(cov.file)
        else
        {
            this.cov = rtracklayer::import(cov.file)
            names(values(this.cov)) = field
            ## field = 'score';
        }
    }
    else {
        this.cov = cov.file
    }

    ## now make sure we have the "best" seqlengths
    .fixsl = function(sl, gr) {sl[seqlevels(gr)] = pmax(GenomeInfoDb::seqlengths(gr), sl[seqlevels(gr)]); return(sl)}

    if (is.null(force.seqlengths)){
        sl = .fixsl(GenomeInfoDb::seqlengths(this.ra), this.cov)
    } else {
        sl = .fixsl(force.seqlengths, this.cov)
    }

    if (!is.null(nseg.file))
    {
        if (is.character(nseg.file))
        {
            if (file.exists(nseg.file))
            {
                stopifnot(file.size(nseg.file)>0)
                if (grepl('\\.rds$', nseg.file, ignore.case = TRUE))
                {
                    nseg = readRDS(nseg.file)
                }
                else if (grepl('(\\.txt$)|(\\.tsv$)|(\\.csv$)', nseg.file))
                {
                    nseg = dt2gr(fread(nseg.file))
                }
                else
                {
                    nseg = rtracklayer::import(nseg.file)
                }
            } else {
                jmessage('Did not find nseg file! Ignore.')
                nseg = nseg.file = NULL
            }

        }
        else {
            nseg = nseg.file
        }
        if (!is.null(nseg)){
            sl = .fixsl(sl, nseg)
        }
    }

    ## make sure all sl's are equiv
    if (is.character(seg.file)){
        if (file.exists(seg.file))
        {
            if (grepl('\\.rds$', seg.file, ignore.case = TRUE))
            {
                this.seg = readRDS(seg.file)
            }
            else if (grepl('(\\.txt$)|(\\.tsv$)|(\\.csv$)', seg.file))
            {
                this.seg = dt2gr(fread(seg.file))
            }
            else
            {
                this.seg = rtracklayer::import(seg.file)
            }
        }
        this.seg = gr.fix(this.seg, sl, drop = T)[, c()]
    } else {
        this.seg = seg.file
    }

    if (length(loose.ends>0))
    {
        if (verbose)
        {
            jmessage('Adding loose ends from vcf file to seg file')
        }
        this.seg = grbind(this.seg, gr.fix(loose.ends, sl, drop = T))
    }

    this.ra = gr.fix(this.ra, sl, drop = T)

    ## DONE: add segmentation to isolate the NA runs
    ## there were a lot of collateral damage because of bad segmentation
    ## na.runs = streduce(
    ##     this.cov[which(is.na(values(this.cov)[, field]))], 1e4
    ## ) %Q% (width>1e5)

    ## this.kag.old = karyograph(this.ra, this.seg)
    ## this.kag = karyograph(this.ra, grbind(this.seg, na.runs))
    this.kag = karyograph(this.ra, this.seg)
    if (length(this.kag$tile)>5e4){
        jmessage("WARNING: big karyograph > 50000 nodes, may take longer to finish.")
    }

    ## NA junctions thrown out
    ## if (length(this.kag$junctions)>0){
    ##     na.cov = this.cov[which(is.na(values(this.cov)[, field]))]
    ##     abe = data.table(as.data.frame(this.kag$ab.edges[,1:2,1, drop=F]))
    ##     colnames(abe) = c("from", "to")
    ##     incident.nodes = abe[, c(from, to)]
    ##     incident.gr = this.kag$tile[incident.nodes]
    ##     incident.gr$nafrac = incident.gr %O% na.cov
    ##     ov.na.runs = which(incident.gr %^% na.runs)
    ##     ov.na.bins = which(incident.gr$nafrac > 0.2)
    ##     na.ix = incident.nodes[union(ov.na.runs, ov.na.bins)]
    ##     ab.exclude = union(ab.exclude,
    ##                        abe[, which(from %in% na.ix | to %in% na.ix)])
    ##     ab.exclude = setdiff(ab.exclude, ab.force)
    ##     if (length(ab.exclude)>0){
    ##         new.jix = setdiff(setNames(seq_along(this.ra), seq_along(this.ra)), ab.exclude)
    ##         ab.force = which(new.jix %in% ab.force)
    ##         this.ra = this.ra[new.jix]
    ##         this.kag = karyograph(this.ra, grbind(this.seg, na.runs))
    ##         if(verbose){
    ##             jmessage("Filtered out ", length(ab.exclude), " junctions in NA enriched regions")
    ##         }
    ##     }
    ## }
    
    if (is.null(nseg.file)){
        warning('No normal copy number values supplied so defaulting to 2 for all segments.')
        this.kag$segstats$ncn = 2
    }

    hets.gr = NULL
    if (!is.null(het.file))
    {
        if (is.character(het.file))
        {
            if (grepl(".rds$", het.file)){
                hets = readRDS(het.file)
            } else {
                hets = fread(het.file)
            }
        }
        else
            hets = het.file

        if (!is.data.table(hets))
            hets = as.data.table(hets)

        if (verbose)
        {
            jmessage('loaded hets')
        }

        if (inherits(hets, "data.frame")){
            if (!is.null(hets$alt.count.n) & !is.null(hets$ref.count.n)){
                ## old format, apply het filter ourselves
                hets$ref.frac.n = hets$alt.count.n / (hets$alt.count.n + hets$ref.count.n)
                hets.gr = dt2gr(hets[pmin(ref.frac.n, 1-ref.frac.n) > 0.2 & (ref.count.n + alt.count.n)>=2, ])
                hets.gr$alt = hets.gr$alt.count.t
                hets.gr$ref = hets.gr$ref.count.t
            } else if (grepl("chrom", colnames(hets)[1], ignore.case = T) &
                       grepl("(pos)?|(start)?", colnames(hets)[2], ignore.case = T) &
                       any(grepl("af$", colnames(hets), ignore.case = T))){
                ## PCAWG BAF format: Chromosome Pos BAF
                colnames(hets)[1:2] = c("chr", "pos")
                baf.col = grep("af$", colnames(hets),
                               value= TRUE, ignore.case = T)[1]
                baf = hets[, baf.col, with = FALSE][[1]]
                hets[, ":="(ref = 1 - baf, alt = baf)]
                hets.gr = dt2gr(hets)
            } else {## new, standard format, with $alt and $ref field
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


        hets.gr = hets.gr[which(hets.gr %^% this.kag$tile)]
        if (!is.null(hets.gr) & length(hets.gr)>0){
            ## save hets object for later
            saveRDS(hets.gr, paste(dirname(out.file), 'hets.gr.rds', sep = '/'))
        } else {
            if (verbose){
                jmessage("None of the provided (if any) germline heterozygosity site overlap the segments, ignore hets.")
            }
        }
    }


    if (length(hets.gr)>0){
        ## pretend we don't have hets at all
        this.kag$segstats = segstats(this.kag$tile,
                                     this.cov,
                                     field = field,
                                     prior_weight = 1,
                                     max.chunk = max.chunk,
                                     ## subsample = subsample,
                                     asignal = hets.gr,
                                     afields = c('ref', 'alt'),
                                     mc.cores = mc.cores,
                                     verbose = verbose,
                                     max.na = max.na,
                                     lp = lp)
    } else {
        this.kag$segstats = segstats(this.kag$tile,
                                     this.cov,
                                     field = field,
                                     prior_weight = 1,
                                     max.chunk = max.chunk,
                                     ## subsample = subsample,
                                     mc.cores = mc.cores,
                                     verbose = verbose,
                                     max.na = max.na,
                                     lp = lp)
    }

    this.kag$segstats$ncn = 2

    if (!is.null(nseg.file)){
        if (is.null(nseg$cn)){
            warning('Normal seg file does not have "cn" met data field. USING the default 2!!')
            this.kag$segstats$ncn = 2
        } else {
            ## check if ncn has chr prefix that is inconsistent with coverage and seg
            this.kag$segstats$ncn = round(gr.val(this.kag$segstats, nseg, 'cn')$cn)
            this.kag$segstats$mean[is.na(this.kag$segstats$ncn)] = NA ## remove segments for which we have no normal copy number
        }
    }

    ## filter ra here
    ## ## 6/15 temp fix for sd on short segments, which we overestimate for now
    ## cov.thresh = pmin(1e5, median(width(this.cov)))
    ## #    jmessage('!!!!!!!!!!! cov.thresh for fix.sd is', cov.thresh, '\n')
    ## fix.sd  = width(this.kag$segstats)<(3*cov.thresh)
    ## #    this.kag$segstats$mean[make.na] = NA
    ## this.kag$segstats$sd[fix.sd] = sqrt(this.kag$segstats$mean[fix.sd])
    ## #      if (is.character(tryCatch(png(paste(out.file, '.ppgrid.png', sep = ''), height = 500, width = 500), error = function(e) 'bla')))
    ## ss.tmp = this.kag$segstats[width(this.kag$segstats)>1e4, ] ## don't use ultra short segments
    ss.tmp = this.kag$segstats %Q% (nbins>10) ## don't use ultra short segments

    purity = as.numeric(purity)
    ploidy = as.numeric(ploidy)
    if (!is.na(purity) & !is.na(ploidy) & length(purity)==1 & length(ploidy)==1) ## purity and ploidy are completely set
    {
        pp = data.table(purity = purity, ploidy = ploidy)
    } else {
        if (grepl(pp.method, "sequenza")){
            use.sequenza = TRUE
            use.ppurple = FALSE
            use.ppgrid = FALSE
        } 
        ## temporarily deprecate Ppurple
        # else if (grepl(pp.method, "ppurple")){
        #     use.ppurple = TRUE
        #     use.sequenza = FALSE
        #     use.ppgrid = FALSE
        # } 
        else if (grepl(pp.method, "ppgrid")){
            use.ppgrid = TRUE
            use.ppurple = FALSE
            use.sequenza = FALSE
        } else {
            use.sequenza = TRUE
            use.ppurple = FALSE
            use.ppgrid = FALSE
            if (verbose){
                jmessage("Cannot recognize the choice of purity-ploidy estimation method. Try the default 'Sequenza'.")
            }
        }

        ## only allow ppurple when hets.gr is absent
        if (!exists("hets.gr")){
            hets.gr = NULL
        }

        if (is.null(hets.gr)){
            use.ppgrid = TRUE
            use.ppurple = use.sequenza = FALSE
        } else if (length(hets.gr)==0){
            use.ppgrid = TRUE
            use.ppurple = use.sequenza = FALSE
        } else if (!all(c("ref.count.t",
                          "alt.count.t",
                          "ref.count.n",
                          "alt.count.n") %in% colnames(values(hets.gr)))){
            use.sequenza = FALSE
            use.ppurple = FALSE ## temp deprecate Ppurple
            use.ppgrid = !use.ppurple
        }

        # if (use.ppurple)
        # {
        #     jmessage("Using Ppurple to estimate purity ploidy")
        #     if (is.na(purity))
        #     {
        #         purity = seq(0, 1, 0.1)
        #     }

        #     if (is.na(ploidy))
        #     {
        #         ploidy = seq(1, 6, 0.2)
        #     }
        #     this.cov$y = values(this.cov)[, field]

        #     if (verbose)
        #     {
        #         jmessage('Computing purity and ploidy with Ppurple')
        #     }

        #     max.chunk = 1e3
        #     numchunks = ceiling(length(ss.tmp)/max.chunk)
        #     if (numchunks>length(purity)*length(ploidy)){
        #         pp = Ppurple::ppurple(cov = this.cov, hets = hets.gr, seg = ss.tmp,
        #                               purities = purity, ploidies = ploidy,
        #                               verbose = verbose,
        #                               mc.cores = mc.cores,
        #                               ## numchunks = numchunks,
        #                               ignore.sex = TRUE)
        #     } else {
        #         pp = Ppurple::ppurple(cov = this.cov, hets = hets.gr, seg = ss.tmp,
        #                               purities = purity, ploidies = ploidy,
        #                               verbose = verbose,
        #                               mc.cores = mc.cores,
        #                               numchunks = numchunks,
        #                               ignore.sex = TRUE)
        #     }
        # } else 
        if (use.sequenza) {
            jmessage("Using Sequenza to estimate purity ploidy")
            if (is.na(purity))
            {
                purity = seq(0, 1, 0.01)
            }

            if (is.na(ploidy))
            {
                ploidy = seq(1, 6, 0.1)
            }
            ## read in the segmentation and heterozygosity site read counts
            sqz.seg = gr2dt(ss.tmp)[strand=="+"]
            setnames(sqz.seg,
                     old = c("seqnames", "start", "end"),
                     new = c("chrom", "start.pos", "end.pos"))

            sites = gr2dt(hets.gr)           
            
            ## prepare input file to run w/ segment.breaks
            sites[, adjusted.ratio := ((ref.count.t + alt.count.t) / (ref.count.n + alt.count.n))]
            sites[, depth.normal := (ref.count.n + alt.count.n)]
            sites[, depth.tumor := (ref.count.t + alt.count.t)]
            sites[, good.reads := 2]
            sites[, zygosity.normal := "het"]
            sites[, alt.frac.t := alt.count.t / (ref.count.t + alt.count.t)]
            setnames(sites,
                     old = c("seqnames", "start", "alt.frac.t"),
                     new = c("chromosome", "position", "Bf"))
            sites = sites[which(Bf <= 0.5)] # They only include BAF w/ values 0-0.5
            if (exists("nseg") && !is.null(nseg)){
                ## only running w/ diploid autosomes
                ## to avoid situations like HCC1143BL
                good.chr = union(as.character(seqnames(nseg %Q% (cn==2))), "X")
                sites = sites[which(chromosome %in% good.chr)]
            } else {
                ## only running w/ chr1-22 and X
                sites = sites[which(chromosome %in% seqlevels(ss.tmp))]
            }
            sites[, Af := 1-Bf]
            ## zchoo Wednesday, Apr 21, 2021 02:02:54 PM
            ## re-factor sites chromosome column to exclude Y as empty factor level
            common.chr = intersect(as.character(sites$chromosome), as.character(sqz.seg$chrom))

            sites = sites[as.character(chromosome) %in% common.chr,]
            sqz.seg = sqz.seg[as.character(chrom) %in% common.chr,]
            
            sites[, chromosome := factor(as.character(chromosome), levels = common.chr)]
            sqz.seg[, chrom := factor(as.character(chrom), levels = common.chr)]
            sqz.seg[, chromosome := chrom]

            ## xtYao Tuesday, Nov 26, 2019 04:43:57 PM: new sequenza expectation, we should freeze their code at this version
            seg.s1 = sequenza::segment.breaks(sites, breaks = sqz.seg, weighted.mean = FALSE)
            ## twalradt Thursday, Apr 26, 2018 02:58:23 PM They wrote '10e6' in their documentation
            seg.filtered  = seg.s1[(seg.s1$end.pos - seg.s1$start.pos) > 1e6, ]
            ## get the genome wide mean of the normalized depth ratio:
            weights.seg <- 150 + round((seg.filtered$end.pos -
                                        seg.filtered$start.pos) / 1e6, 0)
            avg.depth.ratio <- mean(sites$adjusted.ratio) # mean(gc.stats$adj[,2])


            if (verbose){
                jmessage("Starting BAF model fit")
            }

            ## run the BAF model fit
            CP = sequenza::baf.model.fit(Bf = seg.filtered$Bf,
                                         depth.ratio = seg.filtered$depth.ratio,
                                         weight.ratio = weights.seg,
                                         weight.Bf = weights.seg,
                                         avg.depth.ratio = avg.depth.ratio,
                                         cellularity = purity,
                                         ploidy = ploidy)

            confint = sequenza::get.ci(CP)

            pp = data.table(ploidy = confint$max.ploidy,
                            purity = confint$max.cellularity)
        } else if (use.ppgrid){
            jmessage("Using ppgrid to estimate purity ploidy")
            pdf(paste(out.file, '.ppgrid.pdf', sep = ''), height = 10, width = 10)
            if (!is.null(het.file))
            {
                pp = ppgrid(ss.tmp,
                            verbose = verbose,
                            plot = F,
                            mc.cores = mc.cores,
                            purity.min = ifelse(is.na(purity[1]), 0, purity[1]),
                            purity.max = ifelse(is.na(purity[length(purity)]), 1, purity[length(purity)]),
                            ploidy.min = ifelse(is.na(ploidy[1]), 1.2, ploidy[1]),
                            ploidy.max = ifelse(is.na(ploidy[length(ploidy)]), 6, ploidy[length(ploidy)]),
                            allelic = TRUE)
            } else {
                pp = ppgrid(ss.tmp,
                            verbose = verbose,
                            plot = F,
                            mc.cores = mc.cores,
                            purity.min = ifelse(is.na(purity[1]), 0, purity[1]),
                            purity.max = ifelse(is.na(purity[length(purity)]), 1, purity[length(purity)]),
                            ploidy.min = ifelse(is.na(ploidy[1]), 1.2, ploidy[1]),
                            ploidy.max = ifelse(is.na(ploidy[length(ploidy)]), 6, ploidy[length(ploidy)]),
                            allelic = FALSE)
            }
        } else {
            jerror("Need purity ploidy estimates to start JaBbA.")
        }
    }

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
    ## cn is the copy number b4 rounding
    this.kag$segstats$cnmle = rel2abs(this.kag$segstats,
                                      purity = this.kag$purity,
                                      ploidy = this.kag$ploidy,
                                      field = 'mean')
    this.kag$segstats$cn = pmax(round(this.kag$segstats$cnmle), 0)
    ## this.kag$ab.exclude = ab.exclude
    this.kag$ab.force = ab.force
    saveRDS(this.kag, out.file) ## DONE

    ## TODO make these plots more helpful to the users
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
    }
    plot(c(gTrack(gr.fix(sample(this.cov, pmin(length(this.cov), 5e4)), this.kag$segstats), y.field = field, col = alpha('black', 0.3)),
           gTrack(this.kag$segstats, y.field = 'mean', angle = 0, col = 'gray10', border = alpha('black', 0.2))), links = this.kag$junctions, y1 = y1)
    dev.off()
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

    ## sampling random loci to plot not segments
    segsamp = pmin(sample(tmp$mean, 1e6, replace = T, prob = width(tmp)), xlim[2])
    hist(
        pmax(xlim[1], pmin(xlim[2], segsamp)),
        1000, xlab = 'Segment intensity',
        main = sprintf('Purity: %s Ploidy: %s Beta: %s Gamma: %s', kag$purity, kag$ploidy, round(kag$beta,2), round(kag$gamma,2)),
        xlim = c(pmax(0, xlim[1]), pmin(xlim[2], max(segsamp, na.rm = T))))
    abline(v = 1/kag$beta*(0:1000) + kag$gamma/kag$beta, col = alpha('red', 0.3), lty = c(4, rep(2, 1000)))
}

#' @name ramip_stub
#' @rdname internal
#' @noRd
ramip_stub = function(kag.file,
                      out.file,
                      mc.cores = 1,
                      max.threads = Inf,
                      mem = 16,
                      tilim = 1200,
                      slack.prior = 0.001,
                      gamma = NA,
                      beta = NA,
                      customparams = F,
                      purity.min = NA, purity.max = NA,
                      ploidy.min = NA, ploidy.max = NA,
                      init = NULL,
                      mipstart = NULL,
                      use.gurobi = FALSE,
                      epgap = 1e-4,
                      verbose = FALSE,
                      edge.nudge = 0,
                      ab.force = NULL,
                      ab.exclude = NULL,
                      loose.penalty.mode = "boolean",
                      dyn.tuning = TRUE,
                      debug.ix = NULL,
                      lp = FALSE,
                      ism = FALSE,
                      tfield = NULL,
                      fix.thres = -1,
                      min.bins = 5)
{
    outdir = normalizePath(dirname(kag.file))
    this.kag = readRDS(kag.file)

    ## if (ios.null(this.kag$gamma) | is.null(this.kag$beta))
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
    nnaix = which(rowSums(is.na(this.kag$ab.edges[,1:2, 1, drop = FALSE]))==0)
    nna.abe = this.kag$ab.edges[nnaix, , 1, drop = FALSE]
    ## adj.nudge[] = 1*edge.nudge[nnaix] ## if edge.nudge is length ab.edges, then corresponding edges will be nudged
    adj.nudge[nna.abe[,1:2,1]] = 1*edge.nudge

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
        nnaix = which(rowSums(is.na(this.kag$ab.edges[,1:2,1, drop = FALSE]))==0)
        ab.exclude = intersect(ab.exclude, nnaix)
        adj.ub = this.kag$adj*0
        ## adj.ub[rbind(this.kag$ab.edges[ab.exclude, ,1])[, 1:2, drop = FALSE]] = -adj
        ## 1.ub[rbind(this.kag$ab.edges[ab.exclude, ,2])[, 1:2, drop = FALSE]] = -1
        adj.ub[rbind(this.kag$ab.edges[ab.exclude, ,1])[, 1:2, drop = FALSE]] = 0.1
        adj.ub[rbind(this.kag$ab.edges[ab.exclude, ,2])[, 1:2, drop = FALSE]] = 0.1
        saveRDS(adj.ub, paste0(outdir, "/adj.ub.rds"))
    } else {
        adj.ub = NULL
    }

    ## if mipstart is not given, construct the naive solution
    ## if mipstart is given (a gGnome or JaBBA) object
    ## here we create an mipstart "adj" matrix for the new graph
    ## by looking up the junctions in the current graph in the old object
    if (is.null(mipstart)){
        if (file.exists(paste0(outdir, "/mipstart.rds"))){
            jmessage("Using existing mipstart in the current directory")
            mipstart = readRDS(paste0(outdir, "/mipstart.rds"))
        } else {
            jmessage("Adjusting the kag (naive solution) as mipstart (initial solution).")
            ndt = gr2dt(this.kag$segstats)[, ":="(cnmle = cn)][, id := seq_along(this.kag$segstats)]
            setkey(ndt, "id")
            adj = this.kag$ab.adj ## logical mat
            es = data.table(Matrix::which(adj!=0, arr.ind=T))
            if (nrow(es)>0){
                es[, ":="(elb = adj.lb[cbind(row, col)],
                          eub = adj.ub[cbind(row, col)])]
                es[, eub := ifelse(eub!=0, 0, Inf)]
                cn.out.lb = es[, .(out.lb = sum(elb)), keyby=row]
                cn.in.lb = es[, .(in.lb = sum(elb)), keyby=col]
                cn.out.ub = es[, .(out.ub = sum(eub)), keyby=row]
                cn.in.ub = es[, .(in.ub = sum(eub)), keyby=col]
                ## edge lb and ub
                ndt = merge(ndt, cn.out.lb, all.x = TRUE, by.x = "id", by.y = "row")
                ndt = merge(ndt, cn.in.lb, all.x = TRUE, by.x = "id", by.y = "col")
                ndt = merge(ndt, cn.out.ub, all.x = TRUE, by.x = "id", by.y = "row")
                ndt = merge(ndt, cn.in.ub, all.x = TRUE, by.x = "id", by.y = "col")
                pl = this.kag$ploidy
                ## make the best segment CN solution
                ndt[, cn := cnmle]
                ndt[is.na(cn), cn := ceiling(pl)]
                ndt[, ":="(cn.lb = pmax(in.lb, out.lb),
                           cn.ub = pmin(in.ub, out.ub))]
                ## if any violation
                ## ndt[cn > cn.ub, cn := cn.ub]
                ndt[cn < cn.lb, cn := cn.lb]
                if (ndt[, any(is.na(cn) | cn<cn.lb, na.rm=T)]){
                    jerror("Infeasible bounds!!")
                }
                
                ## construct the adj
                es[, ":="(so.cn = ndt[.(row), cn],
                          si.cn = ndt[.(col), cn])]
                es[, ":="(cn = pmax(elb, pmin(eub, pmin(so.cn, si.cn, na.rm=T), na.rm=T), na.rm=T))]
                es[is.na(cn), cn := 0] ## shouldn't be any tho
                mipstart = list(segstats = this.kag$segstats,
                                adj = sparseMatrix(es$row, es$col, x = es$cn,
                                                   dims = dim(this.kag$adj)))
                mipstart$segstats$cn = ndt[, cn] ## use my new cn
                saveRDS(mipstart, paste0(outdir, "/mipstart.rds"))
            } else {
                mipstart = NULL
            }
        }
    }

    if (!is.null(mipstart))
    {
        if (verbose)
            jmessage('Applying mipstarts from previous jabba solution')

        ## for mipstart graph
        if ("loose" %in% colnames(values(mipstart$segstats))){
            not.loose = which(mipstart$segstats$loose==FALSE)
            mgre = suppressWarnings(
                gr.end(mipstart$segstats[not.loose], 1, ignore.strand = FALSE)
            )
            mgrs = suppressWarnings(
                gr.start(mipstart$segstats[not.loose], 1, ignore.strand = FALSE)
            )
            mij = Matrix::which(mipstart$adj[not.loose, not.loose, drop=FALSE] != 0,
                                arr.ind = TRUE)
        } else {
            mgre = suppressWarnings(gr.end(mipstart$segstats,1, ignore.strand = FALSE))
            mgrs = suppressWarnings(gr.start(mipstart$segstats,1, ignore.strand = FALSE))
            mij = Matrix::which(mipstart$adj!=0, arr.ind = TRUE)
        }

        mgend = gr.string(mgre)
        mgstart = gr.string(mgrs)
        mijs = data.table(mstr = paste(mgend[mij[,1]], mgstart[mij[,2]]),
                          cn = mipstart$adj[mij])
        setkey(mijs, mstr)       

        ## for new (i.e. this) graph
        gre = suppressWarnings(gr.end(this.kag$segstats, 1, ignore.strand = FALSE))
        grs = suppressWarnings(gr.start(this.kag$segstats, 1, ignore.strand = FALSE))
        gend = gr.string(gre)
        gstart = gr.string(grs)
        ij = Matrix::which(this.kag$adj!=0, arr.ind = TRUE)
        ijs = data.table(i = ij[,1],
                         j = ij[,2], gstr = paste(gend[ij[,1]], gstart[ij[,2]]),
                         mipstart = as.numeric(NA))
        ijs$is.ref = GenomicRanges::shift(gre[ij[,1]], ifelse(as.logical(strand(gre)[ij[,1]]=="+"), 1, -1)) == grs[ij[,2]]
        setkey(ijs, gstr)

        mcn = mipstart$segstats$cn[gr.match(this.kag$segstats, mipstart$segstats)]
        ## for all remaining ref edges just pick the cn as the cn of the segment that it matches in mipstart
        ijs[gstr %in% mijs$mstr, mipstart := mijs[list(gstr), cn]]
        ijs[is.na(mipstart) & is.ref==TRUE, mipstart := pmin(mcn[i], mcn[j])]
        ijs[is.na(mipstart), mipstart := 0] ## all remaining are 0

        mipstart = sparseMatrix(ijs$i, ijs$j, x = ijs$mipstart,
                                dims = dim(this.kag$adj))
    }

    if (lp) {
        ra.sol = jbaLP(kag.file = kag.file,
                       verbose = verbose,
                       tilim = tilim,
                       epgap = epgap,
                       lambda = 1/slack.prior,
                       ism = ism,
                       tfield = tfield,
                       max.mem = mem,
                       min.bins = min.bins, ## should this always be 1?
                       fix.thres = fix.thres)

    } else {
        ra.sol = jbaMIP(this.kag$adj,
                        this.kag$segstats,
                        beta = this.kag$beta,
                        gamma = this.kag$gamma,
                        tilim = tilim,
                        slack.prior = slack.prior,
                        mipemphasis = 0,
                        mipstart = mipstart, ## make mipstart if not provided
                        adj.lb = adj.lb,
                        epgap = epgap,
                        adj.ub = adj.ub,
                        use.gurobi = use.gurobi,
                        mc.cores = mc.cores,
                        adj.nudge = adj.nudge,
                        outdir = outdir,
                        cn.ub = rep(500, length(this.kag$segstats)),
                        use.L0 = loose.penalty.mode == 'boolean',
                        verbose = verbose,
                        dyn.tuning = dyn.tuning,
                        debug.ix = debug.ix)
    }
    saveRDS(ra.sol, out.file)

    ## add optimization status logging for LP
    if (lp) {
        opt.report = data.table(status = ra.sol$status,
                                obj = ra.sol$obj,
                                epgap = ra.sol$epgap)
    } else {
        opt.report =
            do.call(`rbind`,
                    lapply(seq_along(ra.sol$sols),
                           function(cl){
                               x = ra.sol$sols[[cl]]
                               if (inherits(x$nll.cn, "Matrix") |
                                   inherits(x$nll.cn, "matrix")){
                                   nll.cn = x$nll.cn[1, 1]
                               } else {
                                   nll.cn = NA
                               }
                               width.tot = sum(width(x$segstats %Q% (strand=="+"))/1e6)
                               data.table(cl = cl,
                                          obj = ifelse(is.null(x$obj), NA, x$obj),
                                          width.tot = width.tot,
                                          status = ifelse(is.null(x$status), NA, x$status),
                                          nll.cn = nll.cn,
                                          nll.opt = x$nll.opt,
                                          gap.cn = x$gap.cn,
                                          epgap = ifelse(is.null(x$epgap), NA, x$epgap),
                                          converge = ifelse(is.null(x$converge), NA, x$converge))
                           }))
    }
    saveRDS(opt.report, paste0(outdir, "/opt.report.rds"))
    if (verbose){
        jmessage("Recording convergence status of subgraphs")
    }

    if (customparams)
    {
        system(paste('rm', param.file))
        Sys.setenv(ILOG_CPLEX_PARAMETER_FILE='')
    }
}

##############################
#' @name segstats
#' @title segstats
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
#' ## @param subsample number between 0 and 1 with which to subsample per segment for coverage (useful for superdense coverage eg 50 bases to avoid correlations between samples due to read overlap)
#' @param mc.cores number of cores to run on (default 1)
#' @param lp (logical) if running LP use binstats-style loess smoothing
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
                    max.na = -1,
                    na.thresh = 0.2,
                    verbose = FALSE,
                    ## subsample = NULL, ## number between 0 and 1 to subsample per segment for coverage (useful for dense coverage)
                    mc.cores = 1,
                    nsamp_prior = 1e3, ## number of data samples to estimate alpha / beta prior value
                    ksamp_prior = 100, ## size of data samples to estimate alpha / beta prior values
                    lp = FALSE
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

            ## asignal.df = as.data.frame(asignal)
            ## alpha_high = vaggregate(high.count ~ ix, asignal.df, .postalpha)[as.character(seq_along(target))]
            ## beta_high = vaggregate(high.count ~ ix, asignal.df, .postbeta)[as.character(seq_along(target))]
            ## alpha_low = vaggregate(low.count ~ ix, asignal.df, .postalpha)[as.character(seq_along(target))]
            ## beta_low = vaggregate(low.count ~ ix, asignal.df, .postbeta)[as.character(seq_along(target))]

            asignal.dt = gr2dt(asignal)
            asignal.dt[,
                       ":="(alpha_high = .postalpha(high.count),
                            beta_high = .postbeta(high.count),
                            alpha_low = .postalpha(low.count),
                            beta_low = .postbeta(low.count)),
                       by=ix]
            asignal.dt = asignal.dt[!duplicated(ix), ]
            setkey(asignal.dt, ix)

            target$mean_high = asignal.dt[list(seq_along(target)), alpha_high / beta_high]
            target$sd_high = asignal.dt[list(seq_along(target)), sqrt(alpha_high / (beta_high)^2)]
            target$mean_low = asignal.dt[list(seq_along(target)), alpha_low / beta_low]
            target$sd_low = asignal.dt[list(seq_along(target)), sqrt(alpha_low / (beta_low)^2)]
        }
        else
            jerror('One or more of the afields ', paste(afields, collapse = ', '), ' not found as meta data columns of asignal')
    }

    if (!is.null(signal))
    {
        if (!(field %in% names(values(signal))))
            jerror('Field not found in signal GRanges')

        binwidth = as.numeric(names(sort(
            table(width(sample(signal, 1000, replace=TRUE))), decreasing = TRUE
        )[1]))
        
        utarget = unique(gr.stripstrand(target)) ## strand-agnostic
        if (is.null(names(utarget))){
            names(utarget) = as.character(seq_along(utarget))
        }

        ## start mapping signal to segments
        map = gr2dt(gr.findoverlaps(utarget, signal))
        map[, target.name := names(utarget)[query.id]]
        map[, target.width := width(utarget)[query.id]]
        setkey(map, "target.name")
        mapped = unique(map[, target.name])
        ## these are the segments without a overlapping coverage point
        unmapped = setdiff(names(utarget), mapped) 
        map = map[names(utarget)]
        ## map the value of field
        map[, val := values(signal)[, field][subject.id]]

        ## xtYao ## Monday, Feb 15, 2021 11:10:32 PM
        ## Here explicitly set the infinite coverage values to NA
        ## Otherwise they will make the raw.var NAN
        map[is.infinite(val), val := NA_real_]
        
        ## target$raw.sd = target$sd
        ## map = gr.tile.map(utarget, signal, verbose = T, mc.cores = mc.cores)
        ## sample mean and sample var
        ## sample.mean = sapply(vall, mean, na.rm = TRUE)
        .geom.mean = function(x, na.rm = TRUE){
            exp(mean(log(x[!is.infinite(log(x))]), na.rm=na.rm))
        }

        ## four type of means, all recorded but only use arithmetic
        ## sample.art.mean = sapply(vall, mean, na.rm = TRUE)
        ## sample.median = sapply(vall, median, na.rm = TRUE)
        map[, raw.mean := .geom.mean(val), by = target.name]
        map[, raw.var := var(val, na.rm = TRUE), by = target.name] ## na.rm  = TRUE!!

        ## summarize valid bins per node
        map[, good.bin := !is.na(val)]
        ## valid.signal.ix = which(!is.na(values(signal)[, field]))
        ## utarget$nbins = utarget %N% signal[valid.signal.ix]
        ## target$nbins = sapply(vall, function(x) sum(!is.na(x)))[
        ##     as.character(abs(as.numeric(names(target))))
        ## ]
        ## utarget$nbins.tot = pmax(ceiling(width(utarget)/binwidth), utarget$nbins)
        ## target$nbins.tot = sapply(map, length)[as.character(abs(as.numeric(names(target))))]
        ## utarget$nbins.nafrac = 1 - utarget$nbins/utarget$nbins.tot
        ## utarget$wbins.nafrac = 1 - do.call(gUtils::`%o%`, list(utarget, signal[valid.signal.ix]))/width(utarget)

        ## map[is.na(map)] = NA_real_
        ## sample.geom.mean = sapply(vall, .geom.mean, na.rm = TRUE) ## change to geometric mean???
        ## sample.trim.mean = sapply(vall,
        ##                           function(x){
        ##                               if (sum(!is.na(x))>20){
        ##                                   return(mean(x, trim=0.05, na.rm = TRUE))
        ##                               } else {
        ##                                   return(mean(x, na.rm = TRUE))
        ##                               }
        ##                           })

        ## DEBUGGING: replace arithmetic mean with geometric mean        
        ## sample.mean = sample.art.mean 
        ## sample.mean = sample.geom.mean

        ## sample.var = sapply(vall, var, na.rm = TRUE) ## computing sample variance for each segment
        ix = map[!is.na(raw.mean) & !is.na(raw.var), unique(target.name)]

        if (length(ix)>0){
            target.mdat = map[
               ,.(raw.mean = raw.mean[1],
                  raw.var = raw.var[1],
                  target.name = target.name[1],
                  nbins = sum(good.bin),
                  nbins.tot = .N,
                  nbins.nafrac = 1 - sum(good.bin)/.N,
                  wbins.nafrac = 1 - sum(width[which(good.bin==TRUE)], na.rm = TRUE)/target.width[1]),
                  ## wbins.nafrac = 1 - sum(width[which(good.bin==TRUE)])/sum(width)),
                keyby = target.name]
            values(utarget) = cbind(
                values(utarget),
                target.mdat[names(utarget), .(raw.mean, raw.var, nbins, nbins.tot, nbins.nafrac, wbins.nafrac)]
            )
            ## target$mean[ix] = sample.mean[ix]
            ## target$var[ix] = sample.var[ix]
        } else {
            ## jmessage("Abort: No valid coverage present anywhere!")
            jerror("No valid coverage present anywhere!")
        }
        utarget$mean = utarget$raw.mean

        ## val = values(signal)[, field]
        ## val[is.infinite(val)] = NA
        ## ## val[which(signal$good.prop<0.9)] = NA
        ## vall = lapply(map, function(x) val[x])
        ## vall = vall[match(gr.stripstrand(target), utarget)]

        ## final clean up
        ## target$art.mean = sample.art.mean
        ## target$median = sample.median
        ## target$geom.mean = sample.geom.mean
        ## target$trim.mean = sample.trim.mean
        ## target$raw.mean = target$mean
        ## target$raw.var = target$var
        
        ## map = gr.tile.map(utarget, signal, verbose = T, mc.cores = mc.cores)
        ## val = values(signal)[, field]
        ## val[is.infinite(val)] = NA
        ## vall = lapply(map, function(x) val[x])
        ## vall = vall[match(gr.stripstrand(target), utarget)]

        ## ## sample mean and sample var
        ## sample.mean = sapply(vall, mean, na.rm = TRUE)
        ## sample.var = sapply(vall, var, na.rm = TRUE) ## computing sample variance for each segment
        ## ix = !is.na(sample.mean) & !is.na(sample.var)

        ## target$mean = NA;
        ## if (any(ix)){
        ##     target$mean[ix] = sample.mean[ix]
        ##     target$var[ix] = sample.var[ix]
        ## } else {
        ##     jmessage("Abort: No valid coverage present anywhere!")
        ##     jerror("No valid coverage present anywhere!")
        ## }

        ## target$nbins = sapply(vall, function(x) sum(!is.na(x)))[
        ##     as.character(abs(as.numeric(names(target))))
        ## ]
        ## target$nbins.tot = sapply(map, length)[as.character(abs(as.numeric(names(target))))]
        ## target$nbins.nafrac = 1-target$nbins/target$nbins.tot

        ## ## final clean up
        ## target$raw.mean = target$mean

        ## target$good.prop = (target+1e5) %O% good.bin
        utarget$bad = FALSE
        ## if the user didn't give the max.na, we infer it
        if (!is.numeric(max.na) || !between(max.na, 0, 1)){
            if (verbose){
                jmessage("No `max.na` argument found, inferring it for you now...")
            }
            ## gather the values of nafrac
            ## colnames(values(target))
            ## nafrac = gr2dt(utarget[which(!duplicated(gr.stripstrand(target[, c()])))])[
            ##   , .(seqnames, start, end, tile.id = 1:.N, nbins.nafrac)]
            
            nafrac = utarget[mapped]$wbins.nafrac
            if (var(nafrac)>0){
                ## dat = nafrac[!is.na(nbins.nafrac), cbind(nbins.nafrac)]
                dat = cbind(nafrac)
                ## rownames(dat) = nafrac[!is.na(nbins.nafrac), tile.id]
                ## rownames(dat) = names(utarget)
                rownames(dat) = mapped
                km2 = stats::kmeans(dat, center=2)
                ## telll which part is good/bad
                good = which.min(km2$centers)
                max.na = mean(max(dat[which(km2$cluster==good)]),
                              min(dat[which(km2$cluster!=good)]))
                ## TODO: max.na cannot end up equal to 0!!
                if (verbose){
                    jmessage("The suggested `max.na` is at ", max.na)
                }
            } else {
                max.na = 0
                if (verbose){
                    jmessage("WARNING: your coverage input has no NAs, allowing all of the data...")
                }
            }
        }
        
        ## FIXME: sometimes we'd throw away 1-bin not bad nodes because its variance is NA
        if (length(bad.nodes <- which((utarget$wbins.nafrac >= max.na) | (is.na(utarget$wbins.nafrac))))>0)
        {
            utarget$max.na = max.na ## what about really small segs in a good "environment"
            utarget$bad[bad.nodes] = TRUE
            utarget$mean[bad.nodes] = NA_real_
            
            ## target$sd[bad.nodes] = NA
            if (verbose)
            {
                na.wid = sum(width(utarget[mapped] %Q% which(bad==TRUE)))/1e6
                jmessage("Definining coverage good quality nodes as >=", (1 - max.na)*100, "% bases covered by non-NA and non-Inf values in +/-100KB region")
                jmessage("Hard setting ", na.wid,
                         " Mb of the genome to NA that didn't pass our quality threshold")
                if (na.wid > (sum(as.double(seqlengths(target[mapped])/1e6))/2)){
                    jmessage("WARNING: more than half of the mapped reference genome is set to NA, and ignored by JaBbA!!")
                }
            }
        }

        ## ## loess var estimation
        ## ## i.e. we fit loess function to map segment mean to variance across the sample
        ## ## the assumption is that such a function exists
        ## loess var estimation
        ## i.e. we fit loess function to map segment mean to variance across the sample
        ## the assumption is that such a function exists
        ##        target$nbins = sapply(map, length)[as.character(abs(as.numeric(names(target))))]        
        MINBIN = 1 ## enough data so variance~mean function can be estimated
        tmp = data.table(var = utarget$raw.var,
                         mean = utarget$mean,
                         nbins = utarget$nbins,
                         bad = utarget$bad)[var>0 & nbins>MINBIN & !bad & !is.na(mean) & !is.na(var), ]

        ## if (lp) {
        ##     browser()
        ##     bins.gr = gr.tile(signal, 5e5)
        ##     bins.gr$id = 1:length(bins.gr)
        ##     bins.dt = as.data.table(signal %$% bins.gr) %>% setnames(field, "cn")
        ##     tmp = bins.dt[, .(mean = mean(cn, na.rm = TRUE),
        ##                           var = var(cn, na.rm = TRUE),
        ##                           nbins = sum(!is.na(cn), na.rm = TRUE))][nbins > MINBIN & var > 0,]
        ## }
            
                                  
            

        ## xtYao ## Thursday, Feb 18, 2021 02:07:22 PM
        ## To prevent extreme outlying variances, limit traing data to variance between 0.05 and 0.95 quantile
        middle.mean = tmp[, which(between(mean, quantile(mean, 0.05, na.rm = TRUE), quantile(mean, 0.95, na.rm = TRUE)))]
        middle.var = tmp[, which(between(var, quantile(var, 0.05, na.rm = TRUE), quantile(var, 0.95, na.rm = TRUE)))]

        if (verbose)
        {
            jmessage('Using loess to fit mean to variance relationship in segments with greater than ', MINBIN, ' bins')
        }

        if (nrow(tmp)<10)
        {
            warning(sprintf('Could not find enough (>=10) segments with more than %s bins for modeling mean to variance relationship in data.  Data might be hypersegmented.', MINBIN))
        }
        
        ## overdispersion correction
        ##lmd = tmp[, lm(var ~ mean)]
        ## loe = tmp[, loess(var ~ mean, weights = nbins, span = 2)]
        ## loe = tmp[, loess(var ~ mean, weights = nbins, span = 5)] 
        ## No don't, this is stupid
        ## xtYao ## Wednesday, Feb 17, 2021 02:56:29 PM
        ## Switch to "surface='direct'" for LOESS as it extrapolates
        ## Also, tune up the span parameter to reduce overfitting
        ## loe = tmp[, loess(var ~ mean, weights = nbins, span = 2, control = loess.control(surface = "direct"))]

        ## tmp[, predict.var := predict.lm(lmd, newdata = data.table(mean))] 
        ## tmp[, predict.var := predict(loe, newdata = mean)]
        ## tmp[, predict.var := predict(loe2, newdata = mean)]

        ## tmp = tmp[order(mean)]
        ## ppdf(print(
        ##     ggplot(tmp) +
        ##     geom_point(aes(x = mean, y = var)) +
        ##     geom_line(aes(x = mean, y = predict.var), col = "red") +
        ##     geom_line(aes(x = mean, y = predict.var2), col = "purple") +
        ##     theme_pub()
        ## ))
        ## inferring segment specific variance using loess fit of mean to variance per node
        ## using loess fit as the prior sample var
        ## get the bayesian point estimator (expectation of posterior distribution)
        ## with conjugate prior of scaled inverse chi-sq
        ## min allowable var

        loe.middle.i = tmp[intersect(middle.var, middle.mean), loess(var ~ mean, weights = nbins, span = 5)]
        loe = loe.middle.i
        utarget$loess.var = predict(loe, utarget$mean)
        ## hyperparameter neu, same unit as sample size,
        ## the larger the more weight is put on prior
        ## neu = median(target$nbins, na.rm=T)## neu = 5 ## let's start with this
        neu = utarget$nbins
        ## neu = 434 ## let's start with this
        neu.post = utarget$nbins + neu
        utarget$tau.sq.post = (neu * utarget$loess.var + (utarget$nbins - 1) * utarget$raw.var)/ (neu + utarget$nbins)
        utarget$post.var = try(neu.post * utarget$tau.sq.post / (neu.post - 2))
        utarget$var = utarget$post.var

        min.var = min(tmp$var, na.rm = TRUE)
        ## max.var = min(tmp$var, na.rm = TRUE)

        ## xtYao ## Tuesday, Feb 16, 2021 03:31:04 PM
        ## There could be good small segments with a valid mean without a var
        ## fill them in with just the LOESS prediction        
        miss.var = which(is.na(utarget$var) & !is.na(utarget$mean))
        utarget$var[miss.var] = utarget$loess.var[miss.var]

        ## plot the prediction
        ## tdt = gr2dt(target)
        ## tdt[, more_than_20bins := ifelse(nbins>20, ">20bins", "<=20bins")]
        ## tdt[, more_than_5bins := ifelse(nbins>5, ">5bins", "<=5bins")]

        ## ppdf(print(
        ##     tdt[order(mean)] %>%
        ##     ggplot() +
        ##     geom_point(aes(x = mean, y = raw.var, size = nbins), color = "grey75", alpha = 0.5) +
        ##     geom_line(aes(x = mean, y = loess.var), color = "salmon", lwd = 2) +
        ##     geom_hline(yintercept = 0, lty = "dashed") +
        ##     theme_pub(20) +
        ##     facet_wrap(~ more_than_5bins, nrow = 1, scales = "free")
        ## ), width = 12)

        ## wtf = target %Q% which(strand=="+" & raw.var>8000)
        
        ## clean up NA values which are below or above the domain of the loess function which maps mean -> variance
        ## basically assign all means below the left domain bnound of the function the variance of the left domain bound
        ## and analogously for all means above the right domain bound
        na.var = is.na(utarget$var)
        rrm = range(utarget$mean[!na.var])
        rrv = pmax(predict(loe, rrm), min.var)
        utarget$var[utarget$mean<=rrm[1]] = rrv[1]
        utarget$var[utarget$mean>=rrm[2]] = rrv[2]
        ## no negative var!
        utarget$var[utarget$var<=0] = rrv[1]

        ## if running LP, update loess.var as well
        ## if (lp) {
        ##     jmessage("Running LP mode, filling in loess.var")
        ##     utarget$loess.var[utarget$mean<=rrm[1]] = rrv[1]
        ##     utarget$loess.var[utarget$mean>=rrm[2]] = rrv[2]
        ##     utarget$loess.var[utarget$var<=0] = rrv[1]
        ## }

        ## fill in loess.var generally
        utarget$loess.var[utarget$mean<=rrm[1]] = rrv[1]
        utarget$loess.var[utarget$mean>=rrm[2]] = rrv[2]

        #' zchoo Wednesday, Sep 22, 2021 10:05:59 AM
        ## check specifically for loess.var < 0
        utarget$loess.var[utarget$loess.var<=0] = rrv[1] 

        pdf("var.mean.loess.pdf")
        ## all points training
        ## loe = tmp[, loess(var ~ mean, weights = nbins, span = 5)]
        ## plot(x = tmp$mean, y = tmp$var, pch = 19, cex = 0.5)
        ## lines(x = sort(tmp$mean), y = predict(loe, sort(tmp$mean)), col = "red")
        ## middle means only
        ## loe.middle.mean = tmp[middle.mean, loess(var ~ mean, weights = nbins, span = 5)]
        ## plot(x = tmp$mean, y = tmp$var, pch = 19, cex = 0.5, xlim = tmp[, quantile(mean, c(0.05, 0.95))])
        ## lines(x = sort(tmp[, mean]), y = predict(loe.middle.mean, sort(tmp[, mean])), col = "orange")
        ## middle vars only
        ## loe.middle.var = tmp[middle.var, loess(var ~ mean, weights = nbins, span = 5)]
        ## plot(x = tmp$mean, y = tmp$var, pch = 19, cex = 0.5, xlim = tmp[, quantile(mean, c(0.05, 0.95))])
        ## lines(x = sort(tmp[, var]), y = predict(loe.middle.var, sort(tmp[, mean])), col = "salmon")
        ## union
        ## loe.middle.u = tmp[union(middle.var, middle.mean), loess(var ~ mean, weights = nbins, span = 5)]
        ## plot(x = tmp$mean, y = tmp$var, pch = 19, cex = 0.5, xlim = tmp[, quantile(mean, c(0.05, 0.95))])
        ## lines(x = sort(tmp[, var]), y = predict(loe.middle.u, sort(tmp[, mean])), col = "salmon")
        ## intersect        
        ## pcols = .get_density(tmp$mean, tmp$var)
        ## cr = colorRamp(c("#f7f7f7", "#2166ac"))
        plot(x = tmp$mean, y = tmp$var, pch = 19, cex = 0.5, xlim = tmp[, quantile(mean, c(0.05, 0.95))], ylim = tmp[, quantile(var, c(0.05, 0.95))]## ,
             ## col = col2hex(t(cr(pcols)))
             )
        lines(x = sort(tmp[, var]), y = predict(loe.middle.i, sort(tmp[, mean])), col = "salmon")
        plot(x = tmp[nbins>100]$mean, y = tmp[nbins>100]$var, pch = 19, cex = 0.5, xlim = tmp[, quantile(mean, c(0.05, 0.95))], ylim = tmp[, quantile(var, c(0.05, 0.95))])
        lines(x = sort(tmp[nbins>100, var]), y = predict(loe.middle.i, sort(tmp[nbins>100, mean])), col = "salmon")
        plot(utarget$raw.var, utarget$var)
        dev.off()


        ## xtYao ## Tuesday, Feb 16, 2021 09:53:50 AM
        ## CHECK: no good node should have NA var
        if (any(is.na(utarget$var) & !is.na(utarget$mean))){
            jerror("Some segments with valid mean do not have a variance.")
        }

        if (any(utarget$var<=0, na.rm = TRUE)){
            jerror("Some segments have non-positive variance.")
        }

        ## computing sd / sem for each utarget
        utarget$sd = sqrt((2*utarget$var)/utarget$nbins)

        var.ratio = max(utarget$var,na.rm = TRUE)/min(utarget$var, na.rm = TRUE)

        if ((var.ratio)>1e7)
        {
            warning('Ratio of highest and lowest segment variances exceed 1e7. This could result from very noisy bin data and/or extreme hypersegmentation.  Downstream optimization results may be unstable.')
        }

        ## browser()
        ## ## debug
        ## library(ggplot2)
        ## tdt = gr2dt(utarget)
        ## tdt[, in.tmp := raw.var>0 & nbins>MINBIN & !bad]
        ## tdt = tdt[order(raw.mean)]

        ## ppdf({
        ##     print(
        ##         tdt %>%
        ##         ggplot(aes(x = raw.mean, y = raw.var)) +
        ##         geom_point(aes(size = nbins, color = in.tmp, shape = bad)) +
        ##         geom_line(aes(x = raw.mean, y = loess.var), color = "red") + 
        ##         scale_color_viridis(alpha = 0.5, discrete = TRUE, begin = 0.2, end = 0.8, option = "magma") +
        ##         geom_vline(xintercept = rrm[1], lty = "dashed") +
        ##         geom_vline(xintercept = rrm[2], lty = "dashed") +
        ##         geom_hline(yintercept = rrv[1], lty = "dashed") +
        ##         geom_hline(yintercept = rrv[2], lty = "dashed") +
        ##         geom_vline(xintercept = tmp[, quantile(mean, 0.05)], lty = "dotted") +
        ##         geom_vline(xintercept = tmp[, quantile(mean, 0.95)], lty = "dotted") +
        ##         geom_hline(yintercept = tmp[, quantile(var, 0.05)], lty = "dotted") +
        ##         geom_hline(yintercept = tmp[, quantile(var, 0.95)], lty = "dotted") +
        ##         scale_x_continuous(trans = "log10") +
        ##         scale_y_continuous(trans = "log10") +
        ##         theme_pub()
        ##     )
        ##     print(
        ##         tdt[(in.tmp)] %>%
        ##         ggplot(aes(x = raw.mean, y = raw.var)) +
        ##         geom_point(aes(size = nbins)) +
        ##         geom_line(aes(x = raw.mean, y = loess.var), color = "red") +
        ##         geom_vline(xintercept = rrm[1], lty = "dashed") +
        ##         geom_vline(xintercept = rrm[2], lty = "dashed") +
        ##         geom_hline(yintercept = rrv[1], lty = "dashed") +
        ##         geom_hline(yintercept = rrv[2], lty = "dashed") +
        ##         geom_vline(xintercept = tmp[, quantile(mean, 0.05)], lty = "dotted") +
        ##         geom_vline(xintercept = tmp[, quantile(mean, 0.95)], lty = "dotted") +
        ##         geom_hline(yintercept = tmp[, quantile(var, 0.05)], lty = "dotted") +
        ##         geom_hline(yintercept = tmp[, quantile(var, 0.95)], lty = "dotted") +
        ##         ## scale_color_viridis(alpha = 0.5, discrete = TRUE, begin = 0.2, end = 0.8, option = "magma") + 
        ##         scale_x_continuous(trans = "log10") +
        ##         scale_y_continuous(trans = "log10") +
        ##         theme_pub()
        ##     )
        ##     print(
        ##         tdt %>%
        ##         ggplot(aes(x = raw.mean, y = raw.var)) +
        ##         geom_point(aes(size = nbins, color = in.tmp, shape = bad)) +
        ##         geom_line(aes(x = raw.mean, y = loess.var), color = "red") + 
        ##         scale_color_viridis(alpha = 0.5, discrete = TRUE, begin = 0.2, end = 0.8, option = "magma") +
        ##         geom_vline(xintercept = rrm[1], lty = "dashed") +
        ##         geom_vline(xintercept = rrm[2], lty = "dashed") +
        ##         geom_hline(yintercept = rrv[1], lty = "dashed") +
        ##         geom_hline(yintercept = rrv[2], lty = "dashed") +
        ##         geom_vline(xintercept = tmp[, quantile(mean, 0.05)], lty = "dotted") +
        ##         geom_vline(xintercept = tmp[, quantile(mean, 0.95)], lty = "dotted") +
        ##         geom_hline(yintercept = tmp[, quantile(var, 0.05)], lty = "dotted") +
        ##         geom_hline(yintercept = tmp[, quantile(var, 0.95)], lty = "dotted") +
        ##         ## scale_x_continuous(trans = "log10") +
        ##         ## scale_y_continuous(trans = "log10") +
        ##         theme_pub()
        ##     )
        ##     print(
        ##         tdt[(in.tmp)] %>%
        ##         ggplot(aes(x = raw.mean, y = raw.var)) +
        ##         geom_point(aes(size = nbins)) +
        ##         geom_line(aes(x = raw.mean, y = loess.var), color = "red") +
        ##         geom_vline(xintercept = rrm[1], lty = "dashed") +
        ##         geom_vline(xintercept = rrm[2], lty = "dashed") +
        ##         geom_hline(yintercept = rrv[1], lty = "dashed") +
        ##         geom_hline(yintercept = rrv[2], lty = "dashed") +
        ##         geom_vline(xintercept = tmp[, quantile(mean, 0.05)], lty = "dotted") +
        ##         geom_vline(xintercept = tmp[, quantile(mean, 0.95)], lty = "dotted") +
        ##         geom_hline(yintercept = tmp[, quantile(var, 0.05)], lty = "dotted") +
        ##         geom_hline(yintercept = tmp[, quantile(var, 0.95)], lty = "dotted") +
        ##         ## scale_color_viridis(alpha = 0.5, discrete = TRUE, begin = 0.2, end = 0.8, option = "magma") + 
        ##         ## scale_x_continuous(trans = "log10") +
        ##         ## scale_y_continuous(trans = "log10") +
        ##         theme_pub()
        ##     )
        ##     print(
        ##         tdt[raw.mean>25 & nbins>2] %>%
        ##         ggplot(aes(x = raw.mean, y = raw.var)) +
        ##         geom_point(aes(size = nbins)) +
        ##         geom_line(aes(x = raw.mean, y = loess.var), color = "red") +
        ##         geom_vline(xintercept = rrm[1], lty = "dashed") +
        ##         geom_vline(xintercept = rrm[2], lty = "dashed") +
        ##         geom_hline(yintercept = rrv[1], lty = "dashed") +
        ##         geom_hline(yintercept = rrv[2], lty = "dashed") +
        ##         geom_vline(xintercept = tmp[, quantile(mean, 0.05)], lty = "dotted") +
        ##         geom_vline(xintercept = tmp[, quantile(mean, 0.95)], lty = "dotted") +
        ##         geom_hline(yintercept = tmp[, quantile(var, 0.05)], lty = "dotted") +
        ##         geom_hline(yintercept = tmp[, quantile(var, 0.95)], lty = "dotted") +
        ##         ## scale_color_viridis(alpha = 0.5, discrete = TRUE, begin = 0.2, end = 0.8, option = "magma") + 
        ##         ## scale_x_continuous(trans = "log10") +
        ##         ## scale_y_continuous(trans = "log10") +
        ##         theme_pub()
        ##     )
        ##     print(
        ##         tdt %>%
        ##         ggplot(aes(x = raw.var, y = var)) +
        ##         geom_point(aes(size = nbins, color = in.tmp)) +
        ##         scale_x_continuous(trans = "log10") +
        ##         scale_y_continuous(trans = "log10") +
        ##         theme_pub()
        ##     )
        ## })

        ## finally copy all metadata from utarget to target
        values(target) = values(utarget)[gr.match(target, utarget), ]
    }
    return(target)
}

#' @name jmessage
#' @rdname internal
jmessage = function(..., pre = 'JaBbA')
    message(pre, ' ', paste0(as.character(Sys.time()), ': '), ...)

#' @name jwarning
#' @rdname internal
jwarning = function(..., pre = 'JaBbA', call. = FALSE)
    warning(paste0(pre, ' ', paste0(as.character(Sys.time()), ': '), ...), call. = call.)

#' @name jerror
#' @rdname internal
jerror = function(..., pre = 'JaBbA', call. = TRUE)
    stop(paste0(pre, ' ', paste0(as.character(Sys.time()), ': '), ...), call. = call.)

#' @name jbaLP
#' @title jbaLP
#' 
#' @details 
#'
#' LP analog of jbaMIP
#'
#' @param kag.file (character) path to karyograph
#' @param gg.file (character) path to gGraph 
#' @param kag (karyograph object) karyograph (list)
#' @param gg (gGraph object) gGraph
#' @param cn.field (character) column in karyograph with CN guess, default cnmle
#' @param var.field (character) column in karyograph with node variance estimate, default loess.var
#' @param bins.field (character) column in karyograph containing number of bins, default nbins
#' @param tfield (character) column in junction metadata containing junction tier (default tier)
#' @param min.var (numeric) min allowable variance default 1e-5
#' @param min.bins (numeric) min allowable bins default 1
#' @param lambda (numeric) slack penalty, default 100
#' @param L0 (logical) default TRUE
#' @param M (numeric) max copy number, default 1e3
#' @param verbose (numeric) 0 (nothing) 1 (everything  MIP) 2 (print MIP), default 2 print MIP
#' @param tilim (numeric) default 1e3
#' @param ism (logical) whether to add infinite site assumption constraints. default TRUE
#' @param epgap (numeric) default 1e-3
#' @param max.mem (numeric) maximum memory in GB
#' @param fix.thres (numeric) multiple of lambda above which to fix nodes
#' @param return.type (character) either "gGraph" or "karyograph"
#' @param require.convergence (logical) warn if not converged? default TRUE
#' @param max.epgap.thresh (numeric) above this value, all node and edge CNs are NA (default 0.5)
#' @param nodefileind (numeric) one of 0, 1, 2, 3 (for storing CPLEX tree node files), default 3
#'
#' @returnn
#' karyograph with modified segstats/adj. Adds fields epgap, cl, ecn.in, ecn.out, eslack.in, eslack.out to $segstats and edge CNs to $adj
#' 
#' @author Marcin Imielinski, Zi-Ning Choo
jbaLP = function(kag.file = NULL,
                 gg.file = NULL,
                 kag = NULL,
                 gg = NULL,
                 cn.field = "cnmle",
                 var.field = "loess.var",
                 bins.field = "nbins",
                 tfield = "tier",
                 min.var = 1,
                 min.bins = 1,
                 lambda = 100,
                 L0 = TRUE,
                 M = 1e3,
                 verbose = 2,
                 tilim = 1e3,
                 ism = TRUE,
                 epgap = 1e-3,
                 max.mem = 16,
                 fix.thres = -1,
                 round.thresh = 0.25,
                 return.type = "karyograph",
                 require.convergence = FALSE,
                 max.epgap.thresh = 0.5,
                 nodefileind = 3)
{
    if (is.null(kag.file) & is.null(kag) & is.null(gg.file) & is.null(gg)) {
        stop("one of kag, kag.file, gg.file, gg must be supplied")
    }
    if (!is.null(kag.file) & !is.null(kag)) {
        warning("both kag.file and kag supplied. using kag.")
    }
    if (!is.null(kag)) {
        if (verbose) {
            message("using supplied karyograph")
            kag.gg = gG(jabba = kag)
            beta = kag$beta
        }
    } else {
        if (!is.null(kag.file) && file.exists(kag.file)) {
            if (verbose) {
                message("reading karyograph from file")
            }
            kag = readRDS(kag.file)
            kag.gg = gG(jabba = kag)
            beta = kag$beta
        } else if (!is.null(gg)) {
            if (verbose) {
                message("using supplied gGraph")
            }
            kag.gg = gg$copy
            beta = gg$meta$beta
        } else if (!is.null(gg.file) && file.exists(gg.file)) {
            if (verbose) {
                message("reading gGraph from provided file")
            }
            kag.gg = readRDS(file = gg.file)
            beta = gg$meta$beta
        } else {
            stop("kag.file does not exist and kag not supplied")
        }
    }

    if (verbose) {
        message("Marking nodes with cn contained in column: ", cn.field)
    }
    
    if (is.null(values(kag.gg$nodes$gr)[[cn.field]])) {
        stop("karyograph must have field specified in cn.field")
    }
    kag.gg$nodes$mark(cn  = values(kag.gg$nodes$gr)[[cn.field]])

    if (verbose) {
        message("Computing node weights using variance contained in column: ", var.field)
    }
    
    if (is.null(values(kag.gg$nodes$gr)[[var.field]]) | is.null(values(kag.gg$nodes$gr)[[bins.field]])) {
        warning("karyograph missing var.field. setting weights to node widths")
        wts = width(kag.gg$nodes$gr)
    } else {
        
        ## process variances
        vars = values(kag.gg$nodes$gr)[[var.field]]
        
        if (is.null(beta)) {
            warning("no $beta provided in karyograph/gg metadata")
        } else if (is.na(beta)) {
            warning("NA value for $beta provided in karyograph/gg metadata")
        } else {
            vars = beta * beta * vars
        }

        ## make sure there are no negative variances and that variance is at least CN
        ## this is because there are a larger number of points with low CN and few high CN points
        ## resulting in low estimated variance for high CN nodes even if raw variance was high
        vars = pmax(pmax(vars, kag.gg$nodes$dt$cn), min.var)
        sd = sqrt(vars)

        ## process bins
        bins = values(kag.gg$nodes$gr)[[bins.field]]
        bins = ifelse(bins < min.bins, NA, bins)

        ## compute node weights
        wts = bins / (sd * sqrt(2)) ## for consistency with Laplace distribution
        wts = ifelse(is.infinite(wts) | is.na(wts) | wts < 0, NA, wts)

    }
    kag.gg$nodes$mark(weight = wts)

    ## check for edge rewards
    if (!is.null(kag.gg$edges$dt$reward)) {
        if (verbose) {
            message("Checking edge rewards...")
        }
        erewards = kag.gg$edges$dt[, reward]
        reward.ix = which(erewards != 0)
        if (verbose) {
            message("Number of nonzero rewards: ", length(reward.ix))
        }
        erewards[reward.ix] = lambda / 2 ## set to lambda? or half lambda?
        kag.gg$edges$mark(reward = erewards)
    } else {
        if (verbose) {
            message("Rewards not supplied on edges!")
        }
    }
    
    ## no edge CNs
    kag.gg$edges$mark(cn = NULL)
    kag.gg$nodes[cn > M]$mark(cn = NA, weight = NA)

    ## browser()

    ## add lower bounds depending on ALT junction tier
    if (tfield %in% colnames(kag.gg$edges$dt)) {
        lbs = ifelse(kag.gg$edges$dt[, ..tfield] == 1, 1, 0)
        kag.gg$edges$mark(lb = lbs)
    }

    if (verbose) {
        message("Starting LP balance on gGraph with...")
        message("Number of nodes: ", length(kag.gg$nodes))
        message("Number of edges: ", length(kag.gg$edges))
    }

    ## check for duplicate breakpoints in karyograph junctions
    dup.junctions = detect_duplicate_breakpoints(kag.gg$junctions, tfield = tfield, verbose = verbose)

    ## reset ISM if there are duplicate breakpoints
    if (length(dup.junctions)) {
        if (verbose) {
            message("Number of overlapping Tier 1 junctions: ", length(dup.junctions))
        }
        if (ism) {
            warning("ISM set to TRUE with duplicate Tier 1 junctions, resetting to FALSE to avoid infeasibility")
            ism = FALSE
        }
    } else {
        if (verbose) {
            message("No duplicate Tier 1 junctions detected")
        }
    }

    ## check for heavy nodes to fix
    if (fix.thres > 0) {

        if (fix.thres < 4) {
            warning("Small value for fix.thres selected. Resetting to 4, the minimum recommended value")
            fix.thres = 4
        }
        
        if (verbose) {
            message("Checking for heavy nodes to fix")
        }

        penalty.dt = kag.gg$nodes$dt[, .(node.id, cn, weight)]

        ## compute the difference between 'optimal' CN and 'next-best' CN
        penalty.dt[, penalty := weight * abs(1 - 2 * (cn - floor(cn)))]

        ## compute the 'best' CN - this is just CNMLE
        penalty.dt[, best.cn := round(cn)]

        ## get node ids and CNs of nodes with penalty > lambda * fix.thres
        penalty.dt[, fixed := ((penalty > fix.thres * lambda) & (best.cn >= 0))]
        penalty.dt[fixed == TRUE, lb := best.cn]
        penalty.dt[fixed == TRUE, ub := best.cn]
        penalty.dt[, new.cn := cn]
        penalty.dt[fixed == TRUE, new.cn := best.cn]
        penalty.dt[, new.weight := weight]
        penalty.dt[fixed == TRUE, new.weight := NA]

        if (verbose) {
            message("Number of fixed heavy nodes: ", penalty.dt[fixed == TRUE, .N])
        }

        kag.gg$nodes$mark(cn = penalty.dt$new.cn,
                          unfixed.cn = penalty.dt$cn,
                          weight = penalty.dt$new.weight,
                          unfixed.weight = penalty.dt$weight,
                          fixed = penalty.dt$fixed)

        nfix = penalty.dt[fixed == TRUE, node.id]
    } else {

        nfix = NULL

    }

    if (verbose) {
        message("Grabbing available memory...")
    }

    gc.dat = gc()
    mem.mb = sum(gc.dat[, 2])

    tm = (max.mem * 1e3 - 2 * mem.mb) - 1e3 ## 1 gb buffer - better is to call in balance

    if (verbose) {
        message("Currently used: ", mem.mb, " Mb")
        message("Allowed: ", max.mem * 1e3, " Mb")
    }

    if (tm <= 0) {
        stop("Not enough memory to continue")
    }

    if (verbose) {
        message("Treemem: ", tm, " Mb")
    }

    res = balance(kag.gg,
                  debug = TRUE,
                  lambda = lambda,
                  L0 = L0,
                  verbose = verbose,
                  tilim = tilim,
                  epgap = epgap,
                  lp = TRUE,
                  ism = ism,
                  trelim = tm, ## max.mem * 1e3,
                  nfix = nfix,
                  nodefileind = 3)
    
    bal.gg = res$gg
    sol = res$sol

    if (return.type == "gGraph") {
        return(bal.gg)
    }

    ## check for convergence
    if (require.convergence) {
        if (verbose) {
            message("Checking for convergence to epgap below: ", epgap)
        }
        if (sol$epgap > epgap) {
           warning("Optimization did not converge! Reached epgap: ", sol$epgap)
        }
        if (sol$epgap > max.epgap.thresh) {
            warning("Very high epgap, marking CN as NA")
            bal.gg$nodes$mark(cn = NA)
            bal.gg$edges$mark(cn = NA)
        }
    }
    
    ## just replace things in the outputs
    ## this can create weird errors if the order of kag and bal.gg isn't the same
    out = copy(kag)
    new.segstats = bal.gg$gr
    nnodes = length(new.segstats)

    new.segstats$cl = 1 ## everything same cluster
    new.segstats$epgap = sol$epgap ## add epgap from genome-side opt
    new.segstats$status = sol$status ## solution status to node metadata
    new.segstats$obj = bal.gg$meta$obj ## objective
    
    ## weighted adjacency
    adj = sparseMatrix(i = bal.gg$sedgesdt$from, j = bal.gg$sedgesdt$to,
                       x = bal.gg$sedgesdt$cn, dims = c(nnodes, nnodes))
    
    ## add the necessary columns
    new.segstats$ecn.in = Matrix::colSums(adj, na.rm = TRUE)
    new.segstats$ecn.out = Matrix::rowSums(adj, na.rm = TRUE)
    new.segstats$eslack.in = new.segstats$cn - new.segstats$ecn.in
    new.segstats$eslack.out = new.segstats$cn - new.segstats$ecn.out

    ## NA the telomeric segments?
    qtips = gr.end(si2gr(seqlengths(bal.gg$nodes))) ## location of q arm tips
    term.in = c(which(start(bal.gg$nodes$gr) == 1), ## beginning of chromosome
                -which(bal.gg$nodes$gr %^% qtips)) ## flip side of chromosome end
    term.out = -term.in ## out is reciprocal of in
    telo.in = which(new.segstats$snode.id %in% term.in)
    telo.out = which(new.segstats$snode.id %in% term.out)
    new.segstats$eslack.in[telo.in] = NA
    new.segstats$eslack.out[telo.out] = NA
    
    out$adj = adj

    ## add metadata
    out$segstats = new.segstats
    out$status = sol$status
    out$epgap = sol$epgap
    out$obj = bal.gg$meta$obj
    
    return(out)
}

#' @name detect_duplicate_breakpoints
#' @title detect_duplicate_breakpoints
#'
#' @details
#'
#' Identifies tier 1 junctions sharing breakpoints
#' These will cause MIP to be infeasible if ISM = TRUE
#' 
#' @param juncs (Junction) junction object
#' @param tfield (character) tier field default 'tier'
#' @param verbose (logical) default FALSE
#'
#' @return Junction containing junctions with overlapping breakpoints
#'
#' if tfield is not in junction metadata or all are tier 2 then empty junctions are returned
detect_duplicate_breakpoints = function(juncs, tfield = "tier", verbose = FALSE) {

    if (!tfield %in% colnames(juncs$dt)) {
        warning("tfield missing from metadata")
        return (jJ())
    }

    if (all(is.na(juncs$dt[, ..tfield]))) {
        warning("all entries in tfield NA")
        return (jJ())
    }

    if (!any(juncs$dt[, ..tfield] == 1, na.rm = TRUE)) {
        if (verbose) {
            message("detected no tier 1 junctions")
        }
        return (jJ())
    }

    if (all(juncs$dt$type != "ALT")) {
        if (verbose) {
            message("no ALT junctions detected!")
        }
        return (jJ())
    }

    cols = c(tfield, "seqnames", "start", "end", "strand", "type", "edge.id")
    juncs.dt = as.data.table(stack(juncs$grl))[, ..cols]

    ## select just tier 1 alt edges
    t1.dt = juncs.dt[which(juncs.dt[, ..tfield]==1),][type == "ALT"]

    if (nrow(t1.dt)==0) {
        return (jJ())
    }

    ## add breakpoint and n.unique annotations
    t1.dt[, bp := paste0(seqnames, ":", start, "-", end, strand)]
    t1.dt[, n.unique := length(unique(edge.id)), by = bp]

    ## check if no conflicts
    if (all(t1.dt$n.unique <= 1)) {
        return (jJ())
    }

    ## identify conflicts
    conflict.dt = t1.dt[n.unique > 1,]
    return(juncs[edge.id %in% conflict.dt$edge.id])
}


#' @name jbaMIP
#' @title jbaMIP
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
#' $slack.prior input slack.prior
############################################
jbaMIP = function(adj, # binary n x n adjacency matrix ($adj output of karyograph)
                  segstats, # n x 1 GRanges object with "mean" and "sd" value fields
                  mipstart = NULL, ## sparse adjacency matrix of mipstarts (0 = NA, 0+eps + 0, k>=1 = k)
########### optional args
                  beta, # beta guess
                  gamma, # gamma guess
                  field.ncn = 'ncn', # will use this field to take into account normal copy number in transformation of relative to integer copy number
                  tilim = 20, mipemphasis = 0, epgap = 1e-4, # MIP params
                  ploidy.normal = NULL, ## usually inferred from ncn field but can be entered for subgraph analysis
                  partition = T, ## whether to partition the problem into MIP subproblems depending on the relationships of the segment standard deviation and the value of the slack.prior
                  cn.fix = rep(NA, length(segstats)), ## vector of NA's and (integer) values to which to "fix" copy states, only non NA's are incorporated
                  cn.lb = cn.fix,
                  cn.ub = cn.fix,
                  loose.ends = c(), ## integer vector specifies indices of "loose ends", slack won't be penalized at these vertices
                  adj.lb = 0*adj, # lower bounds for adjacency matrix
                  adj.ub = NULL,
                  adj.nudge = 0*adj, # linear objective function coefficients for edges (only which(adj!=0) components considered)
                  na.node.nudge = TRUE,
                  use.L0 = FALSE,
                  use.gurobi = FALSE, # otherwise will use cplex
                  nsolutions = 1,
                  verbose = F,
                  debug = F,
                  outdir = NULL,
                  mc.cores = 1, ## only matters if partition = T
                  slack.prior = 1,
                  tuning = FALSE, ## whether to invoke CPLEX auto parameter tuning
                  dyn.tuning = TRUE,
                  debug.ix = c(96, 21, 621, 1122, 179, 28, 363, 180, 56, 239, 333),
                  ... # passed to optimizer
                  )
{
    if (length(segstats) != nrow(adj))
        jerror('length(segstats) !=  nrow(adj)')

    if (is.null(adj.lb))
        adj.lb = 0*adj

    ## save the naive solutions
    segstats$kag.cn = segstats$cn

    ## wrapper that calls jbaMIP recursively on subgraphs after "fixing"
    if (partition)
    {
        ## transform means from data space into copy number space
        m = rel2abs(segstats, gamma = gamma, beta = beta, field = 'mean', field.ncn = field.ncn)

        ## transform sds from data space into copy number space (only need to multiply by beta)
        segstats$sd = segstats$sd * beta

        cnmle = round(m) ## MLE estimate for CN
        residual.min = ((m-cnmle)/(segstats$sd))^2
        residual.other =
            apply(cbind(
            (m-cnmle-1)/segstats$sd,
            (m-cnmle+1)/segstats$sd
            )^2,
            1, min)
        ## penalty for moving to closest adjacent copy state
        residual.diff = residual.other - residual.min

        ## we fix nodes for which the penalty for moving to non (locally) optimal copy state
        ## is greater than k / slack.prior penalty (where k is some copy difference
        ## since each node has 4 loose ends
        ## fix = as.integer(which(residual.diff>(4/slack.prior) &
        ##                        cnmle >= 0))
        fix = as.integer(which(residual.diff>(2/slack.prior) &
                               cnmle >= 0))
        ## save the fixing threshold
        segstats$m = m
        segstats$cnmle = cnmle
        segstats$residual.min = residual.min
        segstats$residual.other = residual.other
        segstats$residual.diff = residual.diff

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
        unfix = as.numeric(setdiff(seq_along(segstats), fix))
        tmp.adj = adj
        if (!is.null(adj.ub)){
            tmp.ix = Matrix::which(adj.ub != 0, arr.ind=TRUE)
            tmp.adj[tmp.ix] = 0
        }
        G = graph(as.numeric(t(Matrix::which(tmp.adj != 0, arr.ind=T))), n = length(segstats), directed = T)
        V(G)$name = as.numeric(V(G)) ##  seq_along(V(G)) ## igraph vertex naming is a mystery

        if (length(fix)>0)
            G.unfix = induced.subgraph(G, unfix) + vertices(c(paste('from', fix), paste('to', fix)))
        else
            G.unfix = induced.subgraph(G, unfix)

        if (length(fix)>0 & length(unfix)>0)
            node.map = structure(c(unfix, fix, fix),
                                 names = c(as.character(unfix),
                                           paste('from', fix),
                                           paste('to', fix)))
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
        cl = igraph::clusters(G.unfix, 'weak')
        cll = split(V(G.unfix)$name, cl$membership) ## keep augmented graph names, use node.map later

        ## combine components with their reverse complement components
        ## (only intervals that have a (fixed node free) path from their positive to their negative strand
        ## will be part of the same component .. all other intervals will be separated from their
        ## reverse complement.  However, in the MIP we always optimize
        ## over both strands, and thus must merge components with their reverse complement
        pos.ix = which( as.logical( strand(segstats)=='+') )
        neg.ix = which( as.logical( strand(segstats)=='-') )

        ## maps segments and reverse complements
        seg.map = c(seq_along(pos.ix), suppressWarnings(pos.ix[match(segstats[neg.ix], gr.flipstrand(segstats[pos.ix]))]))

        cll.m = sapply(cll, function(x) paste(sort(seg.map[node.map[x]]), collapse = ' '))
        dup.ix = match(cll.m, unique(cll.m))
                                        #      cll = lapply(split(seq_along(dup.ix), dup.ix), function(x) sort(unique(do.call('c', cll[x]))))
        cll = lapply(split(seq_along(dup.ix), dup.ix), function(x) c(cll[[x[1]]], cll[[x[2]]]))

        ord.ix = order(-sapply(cll, length))
        cll = cll[ord.ix]

        if (verbose)
        {
            jmessage('Partitioned graph into ', length(cll), ' connected components with the size of the highest 10 components being:\n',
                     paste(sapply(cll[1:min(10, length(cll))], length), collapse = ','), '')
        }

        cn.fix = ifelse(seq_along(segstats) %in% fix, cnmle, NA)

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

        sols = parallel::mclapply(seq_along(cll), function(k, args)
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
            if (k<=6){
                saveRDS(args, paste0(outdir, "/.args.", k,".rds"))
            }

            if (!is.null(debug.ix)){
                if (k %in% debug.ix){
                    saveRDS(args, paste0(outdir, "/.args.", k,".rds"))
                }
            }
            
            if (verbose)
                jmessage('Junction balancing subgraph ', k, ' of ',
                         length(cll), ' which has ', length(uix), ' nodes comprising ',
                         round(sum(as.numeric(width(segstats[uix])))/2/1e6, 2), ' MB and ',
                         length(unique(seqnames((segstats[uix])))),
                         ' chromosomes, including chrs ',
                         paste(names(sort(-table(as.character(seqnames((segstats[uix])))))[1:min(4,length(unique(seqnames((segstats[uix])))))]), collapse = ', '))

            if (dyn.tuning){
                ##
                ## New tilim, epgap interplay
                ## the given tilim is tilim.short
                ## tilim.short = args$tilim
                tilim.long = args$tilim
                tilim.short = pmax(tilim.long/10, 10)
                ## tilim.long = 9 * tilim.short ## total time not exceeding 10 times user-defined tilim
                ## epgap.low = args$epgap
                epgap.high = args$epgap
                epgap.low = pmax(epgap.high/100, 1e-4)
                ## epgap.high = pmin(10 * epgap.low, 0.3) ## permissive bound for hard problems
                this.args = args
                this.args$tilim = tilim.short
                this.args$epgap = epgap.low
                this.args$tuning = FALSE
                if (verbose){
                    jmessage("Starting initial run for subgraph ", k,
                             ", aiming epgap at ", this.args$epgap,
                             ", within time limit of ", this.args$tilim)
                }
                out = do.call('jbaMIP', this.args)

                if (!is.na(out$status)){
                    if (out$status %in% c(101, 102)){
                        out$converge = 1
                        jmessage("Subgraph ", k, " converged quickly.")
                    } else {
                        ## prolong the tilim to tilim.long
                        jmessage("Subgraph ", k, " needs prolonged running.")
                        if (out$epgap > epgap.high){
                            ## this.args$tuning = TRUE
                            ## ## Harder prob, try if
                            ## this.args$tilim = tilim.long
                            ## this.args$epgap = epgap.high
                            ## ## this.args$segstats = out$segstats
                            ## this.args$mipstart = out$adj

                            ## if (verbose){
                            ##     jmessage("Starting prolonged run with tuning for subgraph ", k,
                            ##              ", aiming epgap at ", this.args$epgap,
                            ##              ", within time limit of ", this.args$tilim)
                            ## }

                            ## run a round of L1 mode
                            jmessage("Using half the time limit on L1 mode optimization")
                            this.args$tilim = tilim.long/2
                            this.args$epgap = epgap.high
                            this.args$mipstart = out$adj
                            this.args$use.L0 = FALSE
                            out = do.call('jbaMIP', this.args)
                            if (k==1){
                                saveRDS(out, paste0(outdir,"/.tmp.sol.1.rds"))
                            }

                            ## then run a round of L0
                            this.args$use.L0 = TRUE
                            this.args$mipstart = out$adj
                            out = do.call('jbaMIP', this.args)
                            
                            ## converge value:
                            if (out$status %in% c(101, 102)){
                                jmessage("Subgraph ", k, " roughly converged after prolonged session.")
                                out$converge = 3
                            } else {
                                jmessage("Subgraph ", k, " VERYHARD.")
                                out$converge = 4
                            }
                        } else {
                            ## this.args$tilim = tilim.long
                            ## this.args$mipstart = out$adj
                            ## out = do.call('jbaMIP', this.args)
                            out$converge = 2
                            jmessage("Subgraph ", k, " roughly converged.")
                        }
                    }
                } else {
                    jmessage("Subgraph ", k, " has no data to optimize.")
                }
            } else {
                out = do.call('jbaMIP', args)
            }

            
            if (k<=6){
                saveRDS(out, paste0(outdir, "/.sol.", k,".rds"))
            }
            if (!is.null(debug.ix)){
                if (k %in% debug.ix){
                    saveRDS(out, paste0(outdir, "/.sol.", k,".rds"))
                }
            }

            gc() ## garbage collect .. not sure why this needs to be done

            return(out)
        }, args, mc.cores = mc.cores, mc.preschedule = FALSE)

        ## saveRDS(sols, "raw.sols.rds")
        out = list()
        ## scalar fields --> length(cluster) vector
        for (f in c('residual', 'nll.cn', 'nll.opt', 'gap.cn', 'slack.prior')){
            out[[paste('component', f, sep = '')]] = sapply(sols, function(x) x[[f]])
        }

        ## length 2 fields --> length(cluster) x 2 matrix
        for (f in c('ploidy.constraints', 'beta.constraints')){
            out[[paste('component', f, sep = '')]] = do.call('rbind', lapply(sols, function(x) x[[f]]))
        }

        ## adjacency matrix
        out$adj = 0 * adj
        for (i in seq_along(sols))
        {
            ix1 = as.numeric(rownames(sols[[i]]$adj))
            out$adj[ix1, ix1] = out$adj[ix1, ix1] + sols[[i]]$adj
        }

        ## segstats
        sol.ix = lapply(sols, function(x) as.numeric(rownames(x$adj)))

        out$segstats = do.call('grbind', lapply(sols, function(x) x$segstats))[match(seq_along(segstats), unlist(sol.ix))]

        ## annotate segstats keep to keep track and "fixed nodes"
        out$segstats$fixed = seq_along(out$segstats) %in% fix
        out$segstats$cn.fix = cn.fix
        out$segstats$cl  = NA
        out$segstats$id = seq_along(out$segstats)

        ## keep track of which clusters segments originated
        sol.ixul = munlist(sol.ix)
        tmp = vaggregate(sol.ixul[,1], by = list(sol.ixul[,3]), FUN = paste, collapse = ',')
        out$segstats$cl = NA
        out$segstats$cl[as.numeric(names(tmp))] = tmp

        out$segstats$epgap = NA
        out$segstats$epgap[sol.ixul[,3]] = rep(sapply(sols, '[[', "epgap")[sol.ixul[,1]])

        out$purity = 2/(2+gamma)
        v = out$segstats$cn; w = as.numeric(width(out$segstats))
        out$ploidy = sum((v*w)[!is.na(v)]) / sum(w[!is.na(v)])
        out$beta = beta;
        out$gamma = gamma;

        target.less = Matrix::rowSums(adj, na.rm = T)==0
        source.less = Matrix::colSums(adj, na.rm = T)==0
        out$segstats$eslack.out[!target.less] = out$segstats$cn[!target.less] - Matrix::rowSums(out$adj)[!target.less]
        out$segstats$eslack.in[!source.less] =  out$segstats$cn[!source.less] - Matrix::colSums(out$adj)[!source.less]

        out$segstats$ecn.out =  Matrix::rowSums(out$adj)
        out$segstats$ecn.in =  Matrix::colSums(out$adj)

        out$segstats$edges.in = sapply(seq_along(out$segstats),
                                       function(x) {ix = Matrix::which(adj[,x]!=0); paste(ix, '(', out$adj[ix,x], ')', '->', sep = '', collapse = ',')})
        out$segstats$edges.out = sapply(seq_along(out$segstats),
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

    ## take into account (variable) normal cn
    segstats$ncn = rep(2, length(segstats))
    if (!is.null(field.ncn))
        if (field.ncn %in% names(values(segstats)))
            segstats$ncn = values(segstats)[, field.ncn]

    sid = .sid(segstats)
    names(segstats) = sid

    edges = Matrix::which(adj!=0, arr.ind = T)
    if (nrow(edges)>0)
    {
        rownames(edges) = .esid(edges, sid)
    }

    ##
    ## Setting up MIP variables (tracked in varmeta data.table)
    ##

    ## varmeta will keep track of all variables
    ## (i.e. interval, edge, source.slack, target.slack, residual, source.slack.indicator, target.slack.indicator)
    ## id = actual column index in the final matrix
    ## pid = integer specifying parent of variable, which is either row of edges matrix (for edges) and the index of the parent interval in segstats (for everything else)
    ## (all the code below assumes that pid in varmeta are in order from 1 to .N and non missing
    ## for each variable of a given type e.g. varmeta[type == 'residual', identical(pid, 1:.N)])
    ## psid = parent signed id, so that both strands of the same edges / segstats parent,
    ## have the same abs(psid), as.character(psid) also indexes names of the respective edges / segstats object
    if (length(fix.ix <- which(!is.na(cn.fix)))>0){
        cn.lb[fix.ix] = cn.fix[fix.ix]
        cn.ub[fix.ix] = cn.fix[fix.ix]
    }

    varmeta = .varmeta(segstats,
                       edges,
                       adj.lb = adj.lb,
                       adj.ub = adj.ub,
                       cn.lb = cn.lb,
                       cn.ub = cn.ub,
                       gamma = gamma,
                       beta = beta,
                       use.L0 = use.L0)

    ##
    ## Set up MIP constraints (tracked in consmeta)
    ## each constraint has a unique label and we track it's sense
    ## and right hand side, and store its formula
    constraints = .constraints(varmeta,
                               segstats,
                               edges,
                               ploidy.normal = ploidy.normal,
                               use.L0 = use.L0)

    if (is.null(constraints))
    {
        ## if constraints are NULL, then
        ## there are no segments with non NA mean, so we return NA solution
        sol = list();
        sol$residual = NA;
        sol$beta = beta;
        sol$gamma = gamma;
        sol$purity = NA;
        sol$ploidy = NA;
        sol$adj = adj*NA;
        sol$nll.cn = NA;
        sol$nll.opt = NA;
        sol$gap.cn = NA;
        sol$segstats = segstats[, c('mean', 'sd')];
        sol$segstats$cn = NA;
        sol$segstats$ecn.in = NA;
        sol$segstats$ecn.out = NA;
        segstats$ncn = NA;
        sol$segstats$edges.out = sol$segstats$edges.in = rep('', length(segstats));
        sol$segstats$eslack.in = NA;
        sol$segstats$eslack.out = NA;
        sol$slack.prior = slack.prior;
        sol$status = NA;
        sol$epgap = NA;
        sol$converge = NA
        return(sol)
    }

    ## pull constraints data.table and Amat constraints matrix
    consmeta = constraints$consmeta
    Amat = constraints$Amat

    ##
    ## set up objective function
    ##

    ## quadratic portion of objective function
    ## a.k.a. "noise penalty"
    Qobj = Zero = sparseMatrix(1, 1, x = 0, dims = c(nrow(varmeta), nrow(varmeta)))
    s.ix = varmeta[type == 'residual' & dup == FALSE, id]
    noisep = (1/segstats[varmeta[s.ix, pid]]$sd)^2
    noisep = ifelse(is.infinite(noisep), NA, noisep)
    ## remove any infinite noise penalty, eg if sd = 0
    noisep = ifelse(is.na(noisep), 0, noisep) ## set all NA noise penalty segments to 0
    Qobj[cbind(s.ix, s.ix)] = noisep

    ## linear portion of objective function
    ## a.k.a.  "slack penalty"
    cvec = Zero[,1]

    if (use.L0)
    {
        if (verbose)
        {
            jmessage('Applying L0 slack penalty')
        }

        slack.ix = varmeta[type %in% c('source.slack.indicator', 'target.slack.indicator') & !dup, id]
        cvec[slack.ix] = 1/slack.prior
    } else {
        if (verbose)
        {
            jmessage('Applying L1 slack penalty')
        }
        slack.ix = varmeta[type %in% c('source.slack', 'target.slack') & !dup, id]
        cvec[slack.ix] = 1/slack.prior
    }

    ## let any specified "loose ends" have unpenalized slack
    if (length(loose.ends)>0)
    {
        cvec[c(es.s.ix[loose.ends], es.t.ix[loose.ends])] = 0
    }

    if (verbose>1)
    {
        jmessage(sprintf('Total mass on cn portion of objective function: %s. Total mass on edge slack: %s',
                         sum(Qobj[cbind(s.ix, s.ix)]),
                         sum(cvec[slack.ix])))
    }

    if (nrow(edges)>0)
    {
        ## Future TODO: weigh the edges in objective functions
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

        e.ix = varmeta[type == 'edge', id]
        e.pix = varmeta[type == 'edge', pid]
        cvec[e.ix] = e.penalty-adj.nudge[edges[e.pix, , drop = FALSE]] ### reward each edge use in proportion to position in edge nudge
    }

    if (!is.null(mipstart))
    {
        varmeta$mipstart = .mipstart(mipstart, segstats, edges, varmeta, consmeta, Amat, beta, gamma, use.L0)
    }

    ## cap astronomical Qobj values so that CPLEX / gurobi does not freak out about large near-infinite numbers
    ## astronomical = value that is 1e8 higher than lowest value
    qr = range(setdiff(Matrix::diag(Qobj), 0))
    CPLEX.INFIN = 1e9
    Qobj[cbind(1:nrow(Qobj), 1:nrow(Qobj))] =
        pmin(CPLEX.INFIN*qr[1], Qobj[cbind(1:nrow(Qobj), 1:nrow(Qobj))])

    ## run MIP
    if (use.gurobi) # translate into gurobi
    {
        if (verbose)
        {
            jmessage('Running gurobi!')
        }

        model = list()
       model$A = Amat
        model$rhs = varmeta$b;
        model$sense = c('E'='=', 'G'='>=', 'L'='<=')[varmeta$sense]
        model$Q = Qobj;
        model$obj = cvec;
        model$lb = varmeta$lb;
        model$ub = varmeta$ub;
        model$vtype = varmeta$vtype;
        model$modelsense = 'min';

        if (!is.null(varmeta$mipstart))
            model$start = varmeta$mipstart
        else
        {
            model$start = rep(NA, nrow(varmeta));
            mu_hat = as.vector(round(((segstats$mean-gamma)/beta)));
            model$start[varmeta[type == 'interval', id]] = mu_hat;
            model$start[varmeta[type == 'residual', id]] = segstats$mean-(mu_hat*beta+gamma);
            model$start[varmeta[label = 'gamma', id]] = varmeta[label = 'gamma', lb]
            model$start[varmeta[label = 'beta', id]] = varmeta[label = 'beta', lb]
            model$start[is.infinite(model$start)] = NA;
        }
        sol = gurobi::gurobi(model, params = c(list(TimeLimit=tilim), list(...)));
        sol$xopt = sol$x;
    }
    else
    {
        control = c(list(...), list(trace = ifelse(verbose>=2, 1, 0), tilim = tilim, epgap = epgap, mipemphasis = mipemphasis))
        if (!is.null(mipstart)) ## apply mipstart if provided
            control$mipstart = varmeta$mipstart

        if (verbose)
            jmessage('Running CPLEX with relative optimality gap threshold ', epgap)

        if (tuning){
            tuning.control = control
            ## add a few more controls to the parameter tuning function
            tuning.control$tuning.display = 2L
            tuning.control$tuning.rep = 3L
            tuning.control$tuning.tilim = 30
            sol = Rcplex2(cvec = cvec,
                          Amat = Amat,
                          bvec = consmeta$b,
                          sense = consmeta$sense,
                          Qmat = Qobj,
                          lb = varmeta$lb,
                          ub = varmeta$ub,
                          n = nsolutions,
                          objsense = "min",
                          vtype = varmeta$vtype,
                          control = tuning.control,
                          tuning = TRUE)
        } else {
            sol = Rcplex2(cvec = cvec,
                          Amat = Amat,
                          bvec = consmeta$b,
                          sense = consmeta$sense,
                          Qmat = Qobj,
                          lb = varmeta$lb,
                          ub = varmeta$ub,
                          n = nsolutions,
                          objsense = "min",
                          vtype = varmeta$vtype,
                          control = control,
                          tuning = FALSE)
        }
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

    sol.l = lapply(sol.l, function(sol)
    {
        v.ix = varmeta[type == 'interval', id]
        e.ix = varmeta[type == 'edge', id]
        s.ix = varmeta[type == 'residual', id]
        es.t.ix = varmeta[type == 'target.slack', id]
        es.s.ix = varmeta[type == 'source.slack', id]
        beta.ix = varmeta[label == 'beta', id]
        gamma.ix = varmeta[label == 'gamma', id]
        vcn = round(sol$xopt[v.ix])
        ecn = round(sol$xopt[e.ix])
        sol$residual = round(sol$xopt[s.ix])
        sol$beta = sol$xopt[beta.ix]
        sol$gamma = sol$xopt[gamma.ix]
        sol$purity = 2/(2+sol$gamma)
        sol$ploidy = (vcn%*%width(segstats))/sum(as.numeric(width(segstats)))
        sol$adj = adj*0;
        
        sol$nll.cn = ((sol$xopt[s.ix]%*%Qobj[s.ix, s.ix])%*%sol$xopt[s.ix])[1,1]

        if (sum(!is.na(segstats$mean))>0)
            sol$nll.opt = pp.nll(segstats[!is.na(segstats$mean)], gamma = sol$gamma, beta = sol$beta, field = 'mean', field.ncn = field.ncn)$NLL
        else
            sol$nll.opt = NA

        ## supposed to be how far away from naive MLE is the optima
        sol$gap.cn = as.numeric(1 - sol$nll.opt / sol$nll.cn)
        sol$adj[edges] = ecn;
        sol$segstats = segstats
        sol$segstats$cn = round(vcn)
        sol$segstats$ecn.in = round(Matrix::colSums(sol$adj))
        sol$segstats$ecn.out = round(Matrix::rowSums(sol$adj))
        sol$segstats$edges.in = sapply(seq_along(sol$segstats),
                                       function(x) {ix = Matrix::which(adj[,x]!=0); paste(ix, '(', sol$adj[ix,x], ')', '->', sep = '', collapse = ',')})
        sol$segstats$edges.out = sapply(seq_along(sol$segstats),
                                        function(x) {ix = Matrix::which(adj[x, ]!=0); paste('->', ix, '(', sol$adj[x,ix], ')', sep = '', collapse = ',')})

        sol$segstats$eslack.in = round(sol$xopt[es.t.ix])
        sol$segstats$eslack.out = round(sol$xopt[es.s.ix])
        sol$eslack.in = round(sol$xopt[es.t.ix])
        sol$eslack.out = round(sol$xopt[es.s.ix])
        sol$slack.prior = slack.prior

        return(sol)
    });

    sol.l = sol.l[order(sapply(sol.l, function(x) x$obj))]

    if (length(sol.l)==1)
        sol.l = sol.l[[1]]

    return(sol.l)
}

#' @name .sid
#' @title labels segstats with sid in jbaMIP
#' @description
#'
#' Assigns stranded id "sid"
#'
#' @rdname internal
#' @noRd
#'
.sid = function(segstats)
{
    ## Book-keep vertices and their reverse complements
    ##

    ## map intervals to their reverse complement to couple their copy number (and edge variables)
    pos.ix = which( as.logical( strand(segstats)=='+') )
    neg.ix = which( as.logical( strand(segstats)=='-') )

    ## "original vertices"
    og.ix = pos.ix

    ## map flipping positive to negative vertices
    rev.ix = match(segstats, gr.flipstrand(segstats))

    ## use rev.ix to label all reverse complement pairs
    rcpairs = igraph::clusters(graph.edgelist(cbind(seq_along(rev.ix), rev.ix)), 'weak')$membership
    sid = ifelse(duplicated(rcpairs), -1, 1)*rcpairs

    if (any(is.na(rev.ix)))
        jerror('Input genome graph malformed, some nodes missing their exact reverse complement')

    ## "duplicates" of og.ix i.e. revcomp vertices
    dup.ix = suppressWarnings(neg.ix[match(segstats[og.ix], gr.flipstrand(segstats[neg.ix]))])

    if (!identical(segstats$mean[og.ix] , segstats$mean[dup.ix]) & !identical(segstats$sd[og.ix] , segstats$sd[dup.ix]))
        jerror('Segstats mean or sd not identical for all pos / neg strand interval pairs: check segstats computation')
    return(sid)
}


#' @name .esid
#' @title labels segstats with sid in jbaMIP
#' @description
#'
#' Assigns stranded id "sid"
#'
#' @rdname internal
#' @noRd
#'
.esid = function(edges, sid)
{
    edges.dt = as.data.table(edges)
    edges.dt[, es := paste(sid[row], sid[col])]
    edges.dt[, res := paste(-sid[col], -sid[row])]
    erev.ix = match(edges.dt$es, edges.dt$res)
    if (any(is.na(erev.ix)))
        jerror('Input genome graph malformed, some edges missing their exact reverse complement')

    rcepairs = igraph::clusters(graph.edgelist(cbind(seq_along(erev.ix), erev.ix)), 'weak')$membership
    edges.dt[, esid := ifelse(duplicated(rcepairs), -1, 1)*rcepairs]
    return(edges.dt$esid)
}

#' @name constraints
#' @title generates constraints matrix Amat and constraints tracker data.table consmeta
#' @description
#'
#' @rdname internal
#'
#' @param varmeta data.table output of .varmeta tracking all variables
#' @param segstats GRanges named with "sid" (ie output of sid
#' @param edges m x 2 matrix of edges
#' @param ploidy.normal numeric normal (non-tumor) ploidy estimate (default 2)
#' @param use.L0 logical flag whether to set up L0 constraints
#' @noRd
#'
.constraints = function(varmeta, segstats, edges, ploidy.normal, use.L0)
{
    consmeta = data.table() ## store meta data about constraints to keep track
    Zero = sparseMatrix(1, 1, x = 0, dims = c(nrow(varmeta), nrow(varmeta)))
    pid.nna = which(!is.na(segstats$mean) & !is.na(segstats$sd))
    mu.all = NA
    ncn = segstats$ncn

    if (length(pid.nna)==0)
        return(NULL)

    ## query varmeta for relevant indices
    v.ix = varmeta[type == 'interval', id]
    gamma.ix = varmeta[label == 'gamma', id]
    beta.ix = varmeta[label == 'beta', id]
    v.ix.c = varmeta[type == 'interval',][pid.nna, ][psid>0, id]
    s.ix = varmeta[type == "residual", ][varmeta[v.ix.c, ]$pid, id]

    if (is.null(ploidy.normal)){
        ploidy.normal =
            sum(width(segstats)[v.ix.c]*ncn[v.ix.c]) /
            sum(as.numeric(width(segstats))[v.ix.c])
    }

    ## weighted mean across vertices contributing to mean
    mu.all = (width(segstats)[v.ix.c] %*% segstats$mean[v.ix.c]) /
        sum(as.numeric(width(segstats)[v.ix.c]))

    ## set copy number constraints
    Acn = Zero[rep(1, length(v.ix.c)+1), ]
    Acn[cbind(seq_along(v.ix.c), v.ix.c)] = 1;
    Acn[cbind(seq_along(v.ix.c), s.ix)] = 1
    ## taking into account (normal) variable cn
    Acn[cbind(seq_along(v.ix.c), gamma.ix)] = ncn[v.ix.c]/2 
    Acn[cbind(seq_along(v.ix.c), beta.ix)] = -segstats$mean[v.ix.c]

    ## ## final "conservation" constraint
    ## Acn[length(v.ix.c)+1, v.ix] = width(segstats)/sum(as.numeric(width(segstats)));
    ## ##  Acn[length(v.ix.c)+1, s.ix[length(s.ix)]] = 1
    ## ##      Acn[length(v.ix.c)+1, gamma.ix] = 1; ## replacing with below
    ## Acn[length(v.ix.c)+1, gamma.ix] = ploidy.normal/2; ## taking into account (normal) variable cn
    ## Acn[length(v.ix.c)+1, beta.ix] = -mu.all;

    bcn = rep(0, nrow(Acn)) ##

    consmeta = rbind(consmeta,
                     data.table(type = 'Copy',
                                label = paste('Copy', 1:nrow(Acn)),
                                sense = 'E',
                                b = bcn,
                                stringsAsFactors = F))

    ## dup constraints on vertices
    ## constrain every vertex to get the same copy number as its reverse complement
    og.ix = varmeta[type == 'interval' & dup == FALSE, id]
    dup.ix = varmeta[type == 'interval', ][match(-varmeta$psid[og.ix], psid), id]
    Dcn = Zero[rep(1, length(dup.ix)),, drop = F];
    Dcn[cbind(1:nrow(Dcn), dup.ix)] = 1
    Dcn[cbind(1:nrow(Dcn), og.ix)] = -1
    dcn = rep(0, nrow(Dcn))
    sensedcn = rep("E", nrow(Dcn))

    consmeta = rbind(consmeta,
                     data.table(type = 'Dup',
                                label = paste('Dup', 1:nrow(Dcn)),
                                sense = 'E',
                                b = dcn,
                                stringsAsFactors = F))

    ## dup constraints on (reverse complement) edge.slack
    ## (these make sure that reverse complement edge.slacks are
    ## given the same solution as their reverse complement)
    sslack.ix = varmeta[type %in% c('source.slack') & psid>0, id]
    sslack.dup.ix = varmeta[type %in% c('target.slack') & psid<0, ][match(-varmeta$psid[sslack.ix], psid), id]
    tslack.ix = varmeta[type %in% c('target.slack') & psid>0, id]
    tslack.dup.ix = varmeta[type %in% c('source.slack') & psid<0, ][match(-varmeta$psid[tslack.ix], psid), id]
    slack.ix = c(sslack.ix, tslack.ix)
    dup.slack.ix = c(sslack.dup.ix, tslack.dup.ix)
    Ecn = Zero[rep(1, length(slack.ix)),, drop = F];
    Ecn[cbind(1:nrow(Ecn), slack.ix)] = -1
    Ecn[cbind(1:nrow(Ecn), dup.slack.ix)] = 1
    ecn = rep(0, nrow(Ecn))
    Dcn = rbind(Dcn, Ecn)
    dcn = c(dcn, ecn);

    consmeta = rbind(consmeta,
                     data.table(type = 'EdgeSlack',
                                label = paste('EdgeSlack', 1:nrow(Ecn)),
                                sense = 'E',
                                b = ecn,
                                stringsAsFactors = F))

    Acn = rbind(Acn, Dcn)

    Aineq = NULL
    bineq = NULL

    if (nrow(edges)>0)
    {
        ## add edge consistency criteria
        ## for every node that is source of an edge
        ## ensure that sum of weights on outgoing edges
        ## = node weight
        ## do the same for nodes that are targets of edges

        ## gather up ids of edges sources and sinks
        e.ix = varmeta[type == "edge", id]
        e.pix = varmeta[type == "edge", pid]
        so.ix = varmeta[type == 'interval', ][edges[e.pix,1], id]
        ta.ix = varmeta[type == 'interval', ][edges[e.pix,2], id]

        ## unique sources and sink ids of edges
        uso.ix = unique(so.ix)
        uta.ix = unique(ta.ix)

        ## encode these in Bs and Bt where every row
        ## is a junction balance constraint
        Bs = Zero[rep(1, length(uso.ix)), , drop = F]
        Bt = Zero[rep(1, length(uta.ix)), , drop = F]
        Bs[cbind(1:nrow(Bs), uso.ix)] = 1
        Bt[cbind(1:nrow(Bt), uta.ix)] = 1
        Bs[cbind(match(so.ix, uso.ix), e.ix)] = -1
        Bt[cbind(match(ta.ix, uta.ix), e.ix)] = -1

        ## add slack edges for loose ends
        uso.pid = varmeta[uso.ix, pid]
        uta.pid = varmeta[uta.ix, pid]
        uso.lix = varmeta[type == 'source.slack', ][uso.pid, id]
        uta.lix = varmeta[type == 'target.slack', ][uta.pid, id]

        Bs[cbind(1:nrow(Bs), uso.lix)] = -1
        Bt[cbind(1:nrow(Bt), uta.lix)] = -1

        B = rbind(Bs, Bt)

        ## edge duplicate constraints!
        tmp.e = varmeta[type=="edge"]
        ## exact fold-back junctions don't have duplicate
        fb.ix = tmp.e[, which(is.na(match(psid, -psid)))]
        good.e.ix = setdiff(seq_len(nrow(tmp.e)), fb.ix)
        if (length(good.e.ix)==0){
            consmeta =
                rbind(consmeta,
                      data.table(type = 'EdgeSource',
                                 label = paste('EdgeSource', 1:nrow(Bs)),
                                 sense = 'E',
                                 b = 0),
                      data.table(type = 'EdgeTarget',
                                 label = paste('EdgeTarget', 1:nrow(Bt)),
                                 sense = 'E',
                                 b = 0))
            Aed = B
        } else {
            Aedup = Zero[
                seq_len(tmp.e[good.e.ix, sum(dup)]),
              , drop=FALSE]
            reid = tmp.e[, match(psid, -psid)]
            tmp.e[, rid := tmp.e[reid, id]]
            ijs = tmp.e[good.e.ix][dup==FALSE, .(j1 = id, j2 = rid)][, i := 1:nrow(Aedup)]
            Aedup[ijs[, cbind(i, j1)]] = 1
            Aedup[ijs[, cbind(i, j2)]] = -1
            consmeta =
                rbind(consmeta,
                      data.table(type = 'EdgeSource',
                                 label = paste('EdgeSource', 1:nrow(Bs)),
                                 sense = 'E',
                                 b = 0),
                      data.table(type = 'EdgeTarget',
                                 label = paste('EdgeTarget', 1:nrow(Bt)),
                                 sense = 'E',
                                 b = 0),
                      data.table(type = 'EdgeDup',
                                 label = paste('EdgeDup', 1:nrow(Aedup)),
                                 sense = 'E',
                                 b = 0))
            ## populate linear constraints
            Aed = rbind(B, Aedup);
        }
        Amat = rbind(Acn, Aed, Aineq);
    }
    else
    {
        Amat = rbind(Acn, Aineq);
    }

    if (use.L0)
    {
        ## set up indicator constraints
        M = min(c(max(varmeta[type == 'interval', ub], na.rm = TRUE)+1, 1e10))
        if (M>1e10){
            jwarning('Using extremely high copy number upper bounds (above 10000) is not recommended for this model')
        }

        ## only make indicator constraints for non dup source and target
        ## slack variables
        es.s.nz.ix = varmeta[type == 'source.slack.indicator', id]
        es.t.nz.ix = varmeta[type == 'target.slack.indicator', id]
        es.s.nz.pid = varmeta[type == 'source.slack.indicator', pid]
        es.t.nz.pid = varmeta[type == 'target.slack.indicator', pid]
        es.s.ix = varmeta[type == 'source.slack', ][es.s.nz.pid, id]
        es.t.ix = varmeta[type == 'target.slack', ][es.t.nz.pid, id]

        nz.len = length(es.s.nz.ix) + length(es.t.nz.ix)
        vtype = varmeta[, vtype]
        ub = varmeta[, ub]
        lb = varmeta[, lb]

        ## add constraints to make them indicators of loose ends
        ## constraint 1: loose - 0.1 * loose.bool > 0
        ## constraint 2: loose - cn.ub * loose.bool < 0
        new.i = rep(seq_len(nz.len * 2), 2)
        new.j = c(rep(c(es.s.nz.ix, es.t.nz.ix), 2),
                  rep(c(es.s.ix, es.t.ix), 2))
        new.x = c(rep(-0.1, nz.len),
                  rep(-M, nz.len),
                  rep(1, (length(es.s.ix) + length(es.t.ix)) * 2))
        Amat.boolean = sparseMatrix(i = new.i,
                                    j = new.j,
                                    x = new.x,
                                    dims = c(nz.len * 2, nrow(varmeta)))
        b.boolean = rep(0, nz.len * 2)
        sense.boolean = rep(c("G", "L"), each = nz.len)
        consmeta = rbind(consmeta,
                         data.table(type = "SlackBoolean",
                                    label = paste0('SlackBoolean', seq_len(nz.len * 2)),
                                    sense = sense.boolean,
                                    b = b.boolean))
        Amat = rbind(Amat, Amat.boolean)
    }

    colnames(Amat) = varmeta$label
    ## consmeta$formula = paste(arrstring(Amat), '=', consmeta$b)


    return(list(consmeta = consmeta, Amat = Amat))
}


#' @name mipstart
#' @title generates mipstart solution from initial interval by interval mipstart matrix
#' @description
#'
#' @rdname internal
#' @param mipstart n x n matrix of edge solutions
#' @param segstats GRanges named with "sid" (ie output of sid
#' @param edges m x 2 matrix of edges
#' @param varmeta data.table output of .varmeta tracking all variables
#' @param consmeta data.table $consmeta output of .constraints tracking all constraints
#' @param Amat constraint matrix
#' @param beta scalar beta setting
#' @param gamma scalar gamma setting
#' @param use.L0 logical flag whether to set up L0 constraints
#' @noRd
#'
.mipstart = function(mipstart, segstats, edges, varmeta, consmeta, Amat, beta, gamma, use.L0)
{
    ## mips.dt = data.table(edges)
    mips.dt = as.data.table(Matrix::which(mipstart>0, arr.ind = TRUE))
    setnames(mips.dt, c("row", "col"))
    mips.dt[, cn := mipstart[cbind(row, col)]]
    setkeyv(mips.dt, c("row", "col"))

    ## convert everything to data.tables
    consmeta[, id := 1:.N]
    setkey(consmeta, "id")
    varmeta[, id := 1:.N]
    varmeta[type=="interval", mipstart := segstats$cn]
    setkey(varmeta, "id")

    ## NODES
    ## recompute
    ## first mipstart the node copy number n_hat by rounding
    cno = mips.dt[, list(cn = sum(cn, na.rm = TRUE)),  keyby = 'row']
    cni = mips.dt[, list(cn = sum(cn, na.rm = TRUE)),  keyby = 'col']
    varmeta[type == 'interval', in.e.lb := cni[list(pid), cn]]
    varmeta[type == 'interval', out.e.lb := cno[list(pid), cn]]

    ## dial down the edge copy number when there can't be enough node
    dial.down = varmeta[out.e.lb>ub | in.e.lb>ub, pid]
    dial.down.rc = varmeta[type=="interval"][psid %in% varmeta[type=="interval"][pid %in% dial.down, -psid], pid]
    ## dial.down.to = setkey(varmeta[dial.down, .(pid, ub)], "pid")
    ## mips.dt[row %in% dial.down, cn := dial.down.to[.(row), ub]]
    ## mips.dt[col %in% dial.down, cn := dial.down.to[.(col), ub]]
    mips.dt[row %in% dial.down, cn := 0]
    mips.dt[col %in% dial.down.rc, cn := 0]
    mips.dt[col %in% dial.down, cn := 0]
    mips.dt[row %in% dial.down.rc, cn := 0]
    
    ## recompute
    cno = mips.dt[, list(cn = sum(cn, na.rm = TRUE)),  keyby = 'row']
    cni = mips.dt[, list(cn = sum(cn, na.rm = TRUE)),  keyby = 'col']
    varmeta[type == 'interval', in.e.lb := cni[list(pid), cn]]
    varmeta[type == 'interval', out.e.lb := cno[list(pid), cn]]

    varmeta[type == "interval", mipstart := pmin(ub, pmax(mipstart, lb, out.e.lb, in.e.lb, na.rm = TRUE), na.rm=TRUE)]

    varmeta[type == 'interval' & is.na(mipstart), mipstart := pmax(ceiling(segstats$cn[pid]), 0)]
    nna.ix = which(!is.na(segstats$cn))
    pl = sum(segstats$cn[nna.ix] * width(segstats[nna.ix])/1e6)/sum(width(segstats)[nna.ix]/1e6)
    varmeta[type == 'interval' & is.na(mipstart), mipstart := ceiling(pl)]

    ## final cleanup
    varmeta[type == 'interval', mipstart := pmin(ub, pmax(mipstart, lb, out.e.lb, in.e.lb, na.rm = TRUE), na.rm=TRUE)]
    ## check
    node.meta = setkey(varmeta[type=="interval"], "pid")
    mips.dt[, ":="(row.mips = node.meta[.(row), mipstart],
                   col.mips = node.meta[.(col), mipstart])]

    ## EDGES
    varmeta[
        type == 'edge',
        mipstart := mips.dt[
            list(as.data.table(edges[pid,,drop = FALSE])), cn
        ]]
    varmeta[type == 'edge' & is.na(mipstart), mipstart := 0]
    varmeta[, mipstart := pmax(pmin(mipstart, ub), lb)]


    ## SLACKS
    ## which is just the difference between the copy number at
    ## c_i and the incoming / outgoing edges
    e_hat = varmeta[type == 'edge', mipstart]
    n_hat = varmeta[type == 'interval', mipstart]
    ## Bs stores the constraints c_i - - slack_s_i - sum_j \in Es(i) e_j for all i in six
    Bs = Amat[consmeta[type == 'EdgeSource', id],]
    Bs.interval = Bs[, varmeta[type == "interval", id]]
    Bs.interval.ij = data.table(Matrix::which(Bs.interval != 0, arr.ind=T))
    six = Bs.interval.ij[, col, by=row][order(row), col]

    s_slack_hat = rep(0, length(n_hat))
    if (length(varmeta[type == "edge", id])>0)
        s_slack_hat[six] = Bs[, varmeta[type == "edge", id]] %*% e_hat + n_hat[six]
    s_slack_hat[is.na(s_slack_hat)] = 0
    varmeta[type == 'source.slack', mipstart := s_slack_hat]

    ## Bt stores the constraints c_i - - slack_t_i - sum_j \in Et(i) e_j for all i in tix
    Bt = Amat[consmeta[type == 'EdgeTarget', id],]
    Bt.interval = Bt[, varmeta[type == "interval", id]]
    Bt.interval.ij = data.table(Matrix::which(Bt.interval != 0, arr.ind=T))
    tix = Bt.interval.ij[, col, by=row][order(row), col]

    t_slack_hat = rep(0, length(n_hat))
    if (length(varmeta[type == "edge", id])>0)
        t_slack_hat[tix] = Bt[, varmeta[type == "edge", id]] %*% e_hat + n_hat[tix]
    t_slack_hat[is.na(t_slack_hat)] = 0
    varmeta[type == 'target.slack', mipstart := t_slack_hat]

    ## negative slacks can happen when the number of incoming and outgoing
    ## edges exceed the number of nodes ..
    ## this will only happen when lb is provided to certain edges
    ## in this case, we adjust the mipstart by adding cn to nodes and slacks
    ## so that the copy number of the parent nodes associated
    ## with the negative slacks to make everything non-negative
    bad.slacks = varmeta[type %in% c('source.slack', 'target.slack') & mipstart<0, pid]
    slacks.to.adjust = varmeta[type %in% c('source.slack', 'target.slack'), ][pid %in% bad.slacks, ]
    if (nrow(slacks.to.adjust)>0)
    {
        cn.to.adjust = slacks.to.adjust[, min(mipstart), keyby = pid]
        nodes.to.adjust = varmeta[type == 'interval', ][slacks.to.adjust$pid, ]
        nodes.to.adjust$adjust = -cn.to.adjust[list(nodes.to.adjust$pid), V1]
        slacks.to.adjust$adjust = -cn.to.adjust[list(slacks.to.adjust$pid), V1]
        varmeta[nodes.to.adjust$id, ]$mipstart =   nodes.to.adjust$mipstart + nodes.to.adjust$adjust
        varmeta[slacks.to.adjust$id, ]$mipstart = slacks.to.adjust$mipstart + slacks.to.adjust$adjust
    }

    varmeta[label == 'beta', mipstart := beta]
    varmeta[label == 'gamma', mipstart := gamma]

    ## compute residual as difference between rounded and "mean" value
    n_hat = varmeta[type == 'interval', mipstart]
    mu_hat = ifelse(!is.na(segstats$mean), segstats$mean*beta-segstats$ncn/2*gamma, 0)
    eps_hat = mu_hat - n_hat
    varmeta[type == 'residual', mipstart := eps_hat]

    if (use.L0)
    {
        ssid = varmeta[type == 'source.slack.indicator', pid]
        varmeta[type == 'source.slack.indicator', ]$mipstart = sign(varmeta[type == 'source.slack', ][ssid, mipstart]>0)

        tsid = varmeta[type == 'target.slack.indicator', pid]
        varmeta[type == 'target.slack.indicator', ]$mipstart = sign(varmeta[type == 'target.slack', ][tsid, mipstart]>0)
    }

    ## final check on variable upper bound
    varmeta[, mipstart := pmin(mipstart, ub)]
    varmeta[, mipstart := pmax(mipstart, lb)]

    ## sanity check
    b.mipstart = Amat %*% varmeta[, mipstart]
    consmeta[, mipstart := b.mipstart[, 1, drop=TRUE]]

    return(varmeta$mipstart)
}

#' @name varmeta
#' @title makes varmeta in jbaMIP
#' @description
#'
#'
#' assumes that segstats and rows of edges are named with sid
#'
#'   varmeta will keep track of all variables
#' (i.e. interval, edge, source.slack, target.slack, residual, source.slack.indicator, target.slack.indicator)
#' id = actual column index in the final matrix
#'  pid = integer for edges the row in the edges matrix, for everything else the index of the parent interval in segstats
#'  psid = parent signed id, so that both strands of the same edges / segstats parent, have the same abs(psid), as.character(psid) also indexes names of the respective edges / segstats object
#' @rdname internal
#' @param segstats GRanges named with "sid" (ie output of sid
#' @param edges Matrix with columns named row and col, rows labeled with esid
#' @noRd
#'
.varmeta = function(segstats,
                    edges,
                    adj.lb,
                    adj.ub,
                    cn.lb,
                    cn.ub,
                    gamma,
                    beta,
                    use.L0)
{
    ##
    ## Setting up MIP variables (tracked in varmeta data.table)
    ##
    ## TODO: make gGraph conform to this
    sid = as.numeric(names(segstats)) ## assume that segstats has signed integer sid

    varmeta = data.table()
    v.ix = nrow(varmeta) +  seq_along(segstats)

    ## adding intervals
    varmeta = data.table(
        id = v.ix,
        pid = seq_along(sid),
        psid = sid,
        label = paste0('interval', ifelse(sid>0, sid, paste0(-sid, 'r'))),
        type = 'interval',
        dup = sid<0,
        vtype = 'I',
        lb = pmax(0, cn.lb, na.rm = TRUE),
        ub = pmin(cn.ub, Inf, na.rm = TRUE),
        stringsAsFactors = F)

    ## adding residuals: interval cn - cnmle
    s.ix = nrow(varmeta) + v.ix
    varmeta = rbind(varmeta,
                    data.table(
                        id = s.ix,
                        pid = seq_along(sid),
                        psid = sid,
                        label = paste0('residual', ifelse(sid>0, sid, paste0(-sid, 'r'))),
                        dup = sid<0,
                        vtype = 'C',
                        lb = -Inf,
                        ub = Inf,
                        type = 'residual',
                        stringsAsFactors = F))

    if (nrow(edges)>0)
    {
        if (any(tmpix <- adj.ub<0)){ ## tmp fix for negative adj.ub if any exist
            adj.ub[tmpix] = 0.1 ## ie make them >0 <1 effectively constraining to 0
        }

        MAX.EUB = ifelse(max(cn.ub[!is.na(cn.ub)])>0, max(cn.ub[!is.na(cn.ub)]), Inf)
        if (is.null(adj.ub)){
            eub = rep(MAX.EUB, nrow(edges))
        } else {
            eub = ifelse(adj.ub[edges]==0, MAX.EUB, round(pmax(adj.ub[edges], 0)))
        }        
        edges.dt = as.data.table(edges)[, esid := as.numeric(rownames(edges))]
        edges.dt[, lb := pmax(adj.lb[edges], 0)]
        edges.dt[, ub := eub]

        ## adding edges
        e.ix = nrow(varmeta) + (1:nrow(edges.dt))
        varmeta = rbind(varmeta,
                        data.table(
                            id = e.ix,
                            pid = 1:nrow(edges.dt),
                            psid = edges.dt$esid,
                            label = paste0('edge',
                                           ifelse(edges.dt$esid>0,
                                                  edges.dt$esid,
                                                  paste0(-edges.dt$esid, 'r'))),
                            type = 'edge',
                            dup = edges.dt$esid<0,
                            vtype = 'I',
                            lb = edges.dt$lb,
                            ub = edges.dt$ub,
                            stringsAsFactors = F))
    }

    ## adding gamma and beta
    gb.ix = nrow(varmeta) + 1:2
    varmeta = rbind(varmeta,
                    data.table(
                        id = gb.ix,
                        pid = rep(1,2),
                        psid = rep(1, 2),
                        label = c('gamma', 'beta'),
                        type = 'global',
                        dup = FALSE,
                        vtype = 'C',
                        lb = c(gamma, beta),
                        ub = c(gamma, beta),
                        stringsAsFactors = F))

    ## adding slacks
    es.s.ix = nrow(varmeta)+(seq_along(v.ix)) ## adding "source slack" variable
    varmeta = rbind(varmeta,
                    data.table(
                        id = es.s.ix,
                        pid = seq_along(sid),
                        psid = sid,
                        label = paste0('source.slack', ifelse(sid>0, sid, paste0(-sid, 'r'))),
                        type = 'source.slack',
                        dup = sid<0,
                        vtype = 'I',
                        lb = 0,
                        ub = Inf,
                        stringsAsFactors = F))

    es.t.ix = nrow(varmeta)+(seq_along(v.ix)) ## adding "target slack"
    varmeta = rbind(varmeta,
                    data.table(
                        id = es.t.ix,
                        pid = seq_along(sid),
                        psid = sid,
                        label = paste0('target.slack', ifelse(sid>0, sid, paste0(-sid, 'r'))),
                        type = 'target.slack',
                        dup = sid<0,
                        vtype = 'I',
                        lb = 0,
                        ub = Inf, stringsAsFactors = F))

    if (use.L0)
    {
        ## add loose end binary indicator variables for all (non dup) loose end
        ## source and target slack variables
        lix = varmeta[type %in% c('source.slack') & !dup, id]
        ilix = nrow(varmeta) + seq_along(lix)
        varmeta = rbind(varmeta,
                        data.table(
                            id = ilix,
                            pid = varmeta$pid[lix],
                            psid = varmeta$psid[lix],
                            label = paste0('SourceSlackIndicator',
                                           ifelse(varmeta$psid[lix]>0,
                                                  varmeta$psid[lix],
                                                  paste0(varmeta$psid[lix], 'r'))),
                            type = 'source.slack.indicator',
                            dup = varmeta$psid[lix]<0,
                            vtype = 'B',
                            lb = 0,
                            ub = 1,
                            stringsAsFactors = FALSE))

        lix = varmeta[type %in% c('target.slack') & !dup, id]
        ilix = nrow(varmeta) + seq_along(lix)
        varmeta = rbind(varmeta,
                        data.table(
                            id = ilix,
                            pid = varmeta$pid[lix],
                            psid = varmeta$psid[lix],
                            label = paste0('TargetSlackIndicator', ifelse(varmeta$psid[lix]>0, varmeta$psid[lix], paste0(varmeta$psid[lix], 'r'))),
                            type = 'target.slack.indicator',
                            dup = varmeta$psid[lix]<0,
                            vtype = 'B',
                            lb = 0,
                            ub = 1,
                            stringsAsFactors = FALSE))
    }

    return(varmeta)
}


####################################################################
#' @name JaBbA.digest
#' @title JaBbA.digest
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
JaBbA.digest = function(jab, kag, verbose = T, keep.all = T)
{
    if (any(dim(jab$adj) != dim(kag$adj)))
        jerror('JaBbA and karyograph object mismatch')

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
        sinks$id = nrow(adj) + seq_along(sinks)
        sinks$loose = TRUE
        sinks$right = as.logical(strand(sinks)=='+')

        source.ix = which(jab$segstats$eslack.in>0)
        sources = gr.start(jab$segstats[source.ix],ignore.strand = FALSE)
        sources$cn = jab$segstats$eslack.in[source.ix]
        sources$partner.id = source.ix
        sources$id = nrow(adj) + length(sinks) + seq_along(sources)
        sources$loose = TRUE
        sources$right = as.logical(strand(sources)=='+')

        nlends = length(sources) + length(sinks)

        ## pad original matrix with new nodes
        adj.plus = rbind(cbind(adj, sparseMatrix(1,1,x = 0, dims = c(nrow(adj), nlends))),
                         cbind(sparseMatrix(1,1,x = 0, dims = c(nlends, ncol(adj))), sparseMatrix(1,1,x = 0, dims = c(nlends, nlends))))

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

    out = list()

    ## now we have augmented adjacency matrix with loose ends, let's simplify the structure
    ## by collapsing all simple paths
    collapsed = collapse.paths(adj, verbose = verbose)

    ## new segstats formed by reducing "collapsed' sets
    segstats$set.id = collapsed$map
    out$segstats = gr.fix(
        gr.simplify(
            segstats[unlist(lapply(collapsed$sets, sort))],
            val = rep(seq_along(collapsed$sets), sapply(collapsed$sets, length))),
        segstats)

    tmp.ss = gr.string(gr.stripstrand(out$segstats), other.cols = 'loose')
    check1 = all(table(match(tmp.ss, tmp.ss)))
    ## check2 = identical(1:length(collapsed$sets), sort(out$segstats$set.id))

    ##   if (!check1 | !check2) ## quick sanity check to make sure we didn't screw up collapsing
    if (!check1) ## quick sanity check to make sure we didn't screw up collapsing
        jerror('collapse yielded funny / missing segments')
    else
        out$segstats = out$segstats[match(seq_along(collapsed$sets), out$segstats$set.id), c('cn')]

    ## out$segstats$og.ix = sapply(collapsed$sets, function(x) paste(sort(x), collapse = ','))

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
    out$edges = data.table(from = edge.ix[,1], to = edge.ix[,2], cn = out$adj[edge.ix])

    estr = paste(edge.ix[,1], edge.ix[,2])
    abestr = paste(out$ab.edges[,1,1:2], out$ab.edges[,2,1:2])

    if (nrow(out$edges)>0)
    {
        out$edges$type = 'reference'
        if (any(ix <- estr %in% abestr))
            out$edges$type[ix] = 'aberrant'
        if (any(ix <- out$segstats$loose[out$edges[, from]] | out$segstats$loose[out$edges[, to]]))
            out$edges$type[ix] = 'loose'

        out$edges$col = ifelse(out$edges$type == 'aberrant', ifelse(out$edges$cn>0, alpha('red', 0.4), alpha('purple', 0.3)), alpha('gray', 0.2))

        loose.ix = which(out$edges$type == 'loose')
        out$edges$h = 1
    }

    if (length(loose.ix)>0)
    {
        seg.map = match(out$segstats, gr.flipstrand(out$segstats)) ## maps segs to rev comp
        ## maps edges to rev comp
        ed.map = match(paste(out$edges[, from], out$edges[, to]),
                       paste(seg.map[out$edges[, to]], seg.map[out$edges[, from]]))
        temp.ix = which(ed.map>(seq_along(ed.map)));
        ed.id = ed.map
        ed.id[temp.ix] = temp.ix
        ## edges whose rev comp has higher id we name with their index,
        ## and we name their rev comp with their index
        out$edges$col[loose.ix] = alpha('blue', 0.6)
        rh = 0.5 + runif(length(loose.ix)/2)
        out$edges$h[loose.ix] = rh[match(ed.id[loose.ix], unique(ed.id[loose.ix]))]
        ## out$edges$h = ifelse(out$edges$type == 'loose', rand(nrow(out$edges)), 1)
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

    out$segstats$edges.in = sapply(seq_along(out$segstats),
                                   function(x) {ix = which(out$adj[,x]!=0); paste(ix, '(', out$adj[ix,x], ')', '->', sep = '', collapse = ',')})
    out$segstats$edges.out = sapply(seq_along(out$segstats),
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
jabba.alleles = function(jab,
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
        jwarning('count fields not found in meta data of het.sites input, trying BAF...')
        if (!(baf.field %in% names(values(het.sites))))
            jerror('BAF field not found in meta data of het.sites input either!')
        else{
            ## outputs are re.seg$low and re.seg$high
            ## test deviations of observed BAF from expected by beta distribution
            if (verbose)
                jmessage('Processing', length(het.sites),
                         'het sites using fields', baf.field, '\n')

        }
    } else {
        ## jerror('count fields not found in meta data of het.sites input')

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


        highs = split(het.sites$high.count, het.sites$ix)[as.character(seq_along(re.seg))]
        lows = split(het.sites$low.count, het.sites$ix)[as.character(seq_along(re.seg))]

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
        re.seg$low = sapply(seq_along(re.seg), function(i)
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
    high$parent = low$parent = seq_along(jab$segstats)
    high$type = 'high'
    low$type = 'low'
    high$id = seq_along(jab$segstats)
    low$id = length(jab$segstats) + seq_along(jab$segstats)
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

    aadj.final = rbind(
        cbind(aadj.phph, aadj.phunph),
        cbind(aadj.unphph, aadj.unphunph)
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

    tmp.string = gr.string(asegstats, mb = F, other.cols = 'type'); tmp.string2 = gr.string(gr.flipstrand(asegstats), mb = F, other.cols = 'type')
    asegstats$flip.ix = match(tmp.string, tmp.string2)
    asegstats$phased = !unphased

    asegstats.final$edges.in = sapply(seq_along(asegstats.final),
                                      function(x) {ix = Matrix::which(aadj.final[,x]!=0); paste(ix, '(', aadj.final[ix,x], ')', '->', sep = '', collapse = ',')})
    asegstats.final$edges.out = sapply(seq_along(asegstats.final),
                                       function(x) {ix = Matrix::which(aadj.final[x, ]!=0); paste('->', ix, '(', aadj.final[x,ix], ')', sep = '', collapse = ',')})

    asegstats.final$slack.in = asegstats.final$cn - Matrix::colSums(aadj.final)
    asegstats.final$slack.out = asegstats.final$cn - Matrix::rowSums(aadj.final)

    asegstats.final$new.ind = asegstats.final$phased.out = asegstats.final$phased.in = asegstats.final$id = NULL
    asegstats.final$tile.id = as.integer(factor(gr.string(gr.stripstrand(asegstats.final), mb = F, other.cols = 'type')))

    m = sparseMatrix(seq_along(asegstats.final), asegstats.final$parent, x = 1);

    hh = rep(het.sites[, c()], 2)
    hh$count = c(het.sites$high.count, het.sites$low.count)
    hh$type = rep(c('high', 'low'), each = length(het.sites))

    hh$ywid = 0.5
    atd = c(
        gTrack(hh, angle = 0, y.field = 'count', y0 = 0,
               colormaps = list(type = c('high' = alpha('red', 0.3), 'low' = alpha('blue', 0.3))), name = 'hets', y.quantile = 0.001, lwd.border = 2),
        gTrack(asegstats.final, angle = 0, y.field = 'cn', y0 = 0,
               colormaps = list(type = c('high' = alpha('red', 0.3), 'low' = alpha('blue', 0.3), 'total' = alpha('purple', 0.3))), name = 'alleles')
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
        return(cbind(ix = unlist(lapply(seq_along(x), function(y) rep(y, length(x[[y]])))),
                     iix = unlist(lapply(seq_along(x), function(y) seq_along(x[[y]]))),
                     unlist(x)))
    else if (force.rbind)
        return(cbind(ix = unlist(lapply(seq_along(x), function(y) rep(y, nrow(x[[y]])))),
                     iix = unlist(lapply(seq_along(x), function(y) 1:nrow(x[[y]]))),
                     do.call('rbind', x)))
    else if (force.cbind)
        return(t(rbind(ix = unlist(lapply(seq_along(x), function(y) rep(y, ncol(x[[y]])))),
                       iix = unlist(lapply(seq_along(x), function(y) 1:ncol(x[[y]]))),
                       do.call('cbind', x))))
}



## cplex set max threads (warning can only do once globally per machine, so be wary of multiple hosts running on same machine)
##
.cplex_customparams = function(out.file, numthreads = 0, nodefileind = NA, treememlim = NA)
{
    param_lines = "CPLEX Parameter File Version 12.6.0.0"

    param_lines = c(param_lines, paste("CPX_PARAM_THREADS", numthreads, sep = '\t')) ## 12.6.0.0 version
    ## param_lines = c(param_lines, paste("CPXPARAM_Threads", numthreads, sep = '\t'))

    if (!is.na(nodefileind))
        param_lines = c(param_lines, paste("CPX_PARAM_NODEFILEIND", nodefileind, sep = '\t'))

    if (!is.na(treememlim))
    {
        ## #      param_lines = c(param_lines, paste("CPX_PARAM_WORKDIR", getwd(), sep = '\t'))
        param_lines = c(param_lines, paste("CPX_PARAM_TRELIM", treememlim, sep = '\t'))
        ## param_lines = c(param_lines, paste("CPXPARAM_MIP_Limits_TreeMemory", treememlim, sep = '\t'))
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
#' @noRd
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
            for (i in seq_along(all.pos))
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
    out = split(m[,2], m[,1])[as.character(seq_along(query))]
    names(out) = as.character(seq_along(query))
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
#' @noRd
##################################
vaggregate = function(...)
{
    out = aggregate(...);
    return(structure(out[,ncol(out)], names = do.call(paste, lapply(names(out)[1:(ncol(out)-1)], function(x) out[,x]))))
}

#' @name write.tab
#' @rdname internal
#' @noRd
write.tab = function(x, ..., sep = "\t", quote = F, row.names = F)
{
    if (!is.data.frame(x))
        x = as.data.frame(x)

    write.table(x, ..., sep = sep, quote = quote, row.names = row.names)
}



#' @name alpha
#' @noRd
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
#' @noRd
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
#' @noRd
#####################################################
all.paths = function(A, all = F, ALL = F, sources = c(), sinks = c(), source.vertices = sources, sink.vertices = sinks, verbose = FALSE,...)
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

    ij = which(A!=0, arr.ind = T)
    B = sparseMatrix(c(ij[,1], ij[,2]), rep(1:nrow(ij), 2), x = rep(c(-1, 1), each = nrow(ij)), dims = c(nrow(A), nrow(ij)))
    I = diag(rep(1, nrow(A)))

    source.vertices = setdiff(match(source.vertices, node.ix), NA)
    sink.vertices = setdiff(match(sink.vertices, node.ix), NA)

    B2 = cbind(B, I[, source.vertices, drop = FALSE], -I[, sink.vertices, drop = FALSE])

    if (verbose)
        cat(sprintf('Computing paths for %s vertices and %s edges\n', nrow(B2), ncol(B2)))

    K = convex.basis(B2, verbose = verbose, ...)

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
        tmp.cix = cbind(unlist(lapply(seq_along(out$cycles), function(x) rep(x, length(out$cycles[[x]])))), unlist(out$cycles))
        out$cycles = out$cycles[!duplicated(as.matrix(sparseMatrix(tmp.cix[,1], tmp.cix[,2], x = 1)))]
    }

    if (length(out$paths)>0)
    {
        tmp.pix = cbind(unlist(lapply(seq_along(out$paths), function(x) rep(x, length(out$paths[[x]])))), unlist(out$paths))
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
#' @noRd
###############################################
collapse.paths = function(G, verbose = T)
{
    if (inherits(G, 'igraph'))
        G = G[,]

    out = G!=0

    ## ##  if (verbose)
    ## ##     cat('graph size:', nrow(out), 'nodes\n')

    ## first identify all nodes with exactly one parent and child to do initial collapsing of graph
    singletons = which(Matrix::rowSums(out)==1 & Matrix::colSums(out)==1)

    ## #  if (verbose)
    ## #      cat('Collapsing simple paths..\n')

    sets = split(1:nrow(G), 1:nrow(G))
    if (length(singletons)>0)
    {
        tmp = out[singletons, singletons]
        cl = igraph::clusters(graph(as.numeric(t(Matrix::which(tmp, arr.ind = TRUE))), n = nrow(tmp)), 'weak')$membership
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
    todo[Matrix::rowSums(out)==1 | Matrix::colSums(out)==1] = TRUE ## could also be 1/0 or 0/1!!!

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
    map[unlist(sets)] = match(rep(seq_along(sets), slen), ix)
    out = out[ix, ix]
    colnames(out) = rownames(out) = NULL

    return(list(adj = out, map = map, sets = split(seq_along(map), map)))
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
#' @noRd
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


#' @name convex.basis
#' @rdname internal
#' @description
#' convex.basis
#'
#' Outputs a matrix K of the convex basis of matrix A
#'
#' i.e. each column x = K[,i] is a minimal solution (with respect to sparsity) to
#' Ax = 0, x>=0
#'
#' @noRd
convex.basis = function(A, interval = 80, chunksize = 100, maxchunks = Inf,
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
                    ## coeff = c(-A_i[i, indpairs[,2]], A_i[i, indpairs[,1]])  ## dealing with Matrix ghost
                    coeff = c(-A_i[i, ][indpairs[,2]], A_i[i, ][indpairs[,1]])  ##
                    combs = sparseMatrix(pix, ix, x = coeff, dims = c(nrow(indpairs), nrow(K_last)))
                    combs[cbind(pix, ix)] = coeff;

                    H = combs %*% K_last;

                    ## remove duplicated rows in H (with respect to sparsity)
                    H = H[which(!duplicated(as.matrix(H)>ZERO)), , drop = FALSE];

                    ## remove rows in H that have subsets in H (with respect to sparsity) ..
                    if ((as.numeric(nrow(H))*as.numeric(nrow(H)))>maxchunks)
                    {
                        print('Exceeding maximum number of chunks in convex.basis computation')
                        jerror('Exceeding maximum number of chunks in convex.basis computation')
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
                                jerror('Exceeding maximum number of chunks in convex.basis computation')
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
                                jerror('Exceeding maximum number of chunks in convex.basis computation')
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

        A_i = A %*% t(K_i)
    }

    return(t(K_i))
}

############################################################
#' read.junctions: parse junction data from various common formats
#'
#' @name read.junctions
#'
#' @description Parsing various formats of structural variation data into junctions.
#'
#' @usage read.juncs(rafile,
#' keep.features = T,
#' seqlengths = hg_seqlengths(),
#' chr.convert = T,
#' geno=NULL,
#' flipstrand = FALSE,
#' swap.header = NULL,
#' breakpointer = FALSE,
#' seqlevels = NULL,
#' force.bnd = FALSE,
#' skip = NA)
#'
#' @param rafile path to the junctions file. See details for the compatible formats.
#' @param keep.features \code{logical}, if TRUE preserve meta data from the input
#' @param seqlengths a named \code{numeric} vector containing reference contig lengths
#' @param chr.convert \code{logical}, if TRUE strip "chr" prefix from contig names
#' @param geno \code{logical}, whether to parse the 'geno' fields of VCF
#' @param flipstrand \code{logical}, if TRUE will flip breakpoint strand
#' @param swap.header path to the alternative VCF header file
#' @param breakpointer \code{logical}, if TRUE will parse as breakpointer output
#' @param seqlevels vector for renaming the chromosomes
#' @param force.bnd if TRUE overwrite all junction "type" to "BND"
#' @param skip \code{numeric} lines to skip
#'
#' @details
#' A junction is a unordered pair of strand-specific genomic locations (breakpoints). Within a given
#' reference genome coordinate system, we call the direction in which coordinates increase "+". A breakpoint
#' is a width 1 (\code{start==end})genomic range with \code{strand} specified, and "+" means the side with larger
#' coordinate is fused with the other breakpoint in a junction.
#'
#' \code{rafile} must be one of the following formats:
#' 1) Some VCF (variant call format). We currently support the VCF output
#' from a number of structural variation detection methods, namely
#' SvABA (https://github.com/walaj/svaba),
#' DELLY (https://github.com/dellytools/delly),
#' LUMPY (https://github.com/arq5x/lumpy-sv),
#' novoBreak (https://sourceforge.net/projects/novobreak/). In theory,
#' VCF defined with BND style should be compatible but be cautious
#' when using the output from other methods since
#' no universal data definition is adopted by the community yet.
#' 2) BEDPE (http://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format)
#' 3) Textual output from Breakpointer
#' (http://archive.broadinstitute.org/cancer/cga/breakpointer)
#' 4) R serialized object storing junctions (.rds)
#'
#' @section Warning:
#' We assume the orientation definition in the input is consistent with ours. Check with
#' the documentation of your respective method to make sure. If the contrary, use
#' \code{flipstrand=TRUE} to reconcile.
#'
#' @return a \code{GRangesList} of the junctions
#'
#' @importFrom VariantAnnotation readVcf info geno
#' 
#' @export
###########################################################
read.junctions = function(rafile,
                          keep.features = T,
                          ## seqlengths = hg_seqlengths(),
                          seqlengths = NULL,
                          chr.convert = T,
                          geno=FALSE,
                          flipstrand = FALSE,
                          swap.header = NULL,
                          breakpointer = FALSE,
                          seqlevels = NULL,
                          force.bnd = FALSE,
                          skip = NA,
                          get.loose = FALSE){
    if (is.null(rafile)){
        return(GRangesList())
    } else if (inherits(rafile, "GRangesList")){
        return(verify.junctions(rafile))
    } else if (inherits(rafile, "Junction")){
        return(verify.junctions(rafile$grl))
    } else if (is.character(rafile)){
        if (!file.exists(rafile)){
            return(NULL)
        }
        if (grepl('.rds$', rafile)){
            ra = readRDS(rafile)
            ## validity check written for "junctions" class
            if (inherits(ra, "Junction")){
                ra = ra$grl
            }            
            return(verify.junctions(ra))
        } else if (grepl('bedpe(\\.gz)?$', rafile)){
            ra.path = rafile
            cols = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'name', 'score', 'str1', 'str2')

            ln = readLines(ra.path)
            if (is.na(skip)){
                nh = min(c(Inf, which(!grepl('^((#)|(chrom1))', ln))))-1
                if (is.infinite(nh)){
                    nh = 1
                }
            } else{
                nh = skip
            }

            if ((length(ln)-nh)==0){
                ## if (get.loose){
                ##     return(list(junctions = GRangesList(GRanges(seqlengths = seqlengths))[c()], loose.ends = GRanges(seqlengths = seqlengths)))
                ## }
                ## else{
                return(GRangesList(GRanges(seqlengths = seqlengths))[c()])
                ## }
            }

            if (nh ==0){
                rafile = fread(rafile, header = FALSE)
            } else {

                rafile = tryCatch(fread(ra.path, header = FALSE, skip = nh), error = function(e) NULL)
                if (is.null(rafile)){
                    rafile = tryCatch(fread(ra.path, header = FALSE, skip = nh, sep = '\t'), error = function(e) NULL)
                }

                if (is.null(rafile)){
                    rafile = tryCatch(fread(ra.path, header = FALSE, skip = nh, sep = ','), error = function(e) NULL)
                }

                if (is.null(rafile)){
                    jerror('Error reading bedpe')
                }
            }

            if (nrow(rafile)==0)
                return(GRangesList())
            ## this is not robust enough! there might be mismatching colnames
            setnames(rafile, seq_along(cols), cols)
            rafile[, str1 := ifelse(str1 %in% c('+', '-'), str1, '*')]
            rafile[, str2 := ifelse(str2 %in% c('+', '-'), str2, '*')]
            ## converting bedpe to 1-based coordinates for verify.junctions()
            rafile[, `:=`(
                end1 = ifelse(start1==end1-1, start1, end1),
                end2 = ifelse(start2==end2-1, start2, end2)
            )]
        } else if (grepl('(vcf$)|(vcf.gz$)', rafile)){

          if (!is.null(seqlengths) && all(!is.na(seqlengths)))
            {
              vcf = VariantAnnotation::readVcf(
                                         rafile, genome = Seqinfo(
                                                   seqnames = names(seqlengths),
                                                   seqlengths = as.vector(seqlengths)))
            }
          else ## get seqlengths from vcf
          {
            vcf = VariantAnnotation::readVcf(
                                         rafile)
          }
            ## vgr = rowData(vcf) ## parse BND format
            vgr = DelayedArray::rowRanges(vcf) ## this range is identical to using read_vcf
            ## old.vgr = read_vcf(rafile, swap.header = swap.header, geno=geno)
            mc = data.table(as.data.frame(mcols(vgr)))

            ## append the INFO
            info.dt = data.table(
                as.data.frame(VariantAnnotation::info(vcf))
            )
            mc = cbind(mc, info.dt)
            values(vgr) = mc

            if (!('SVTYPE' %in% colnames(mc))) {
                jwarning('Vcf not in proper format.  Is this a rearrangement vcf?')
                return(GRangesList());
            }

            if (any(w.0 <- (width(vgr)<1))){
                jwarning("Some breakpoint width==0.")
                ## right bound smaller coor
                ## and there's no negative width GR allowed
                vgr[which(w.0)] = gr.start(vgr[which(w.0)]) %-% 1
            }

            ## BND format doesn't have duplicated rownames
            if (any(duplicated(names(vgr)))){
                names(vgr) = NULL
            } 

            ## no events
            if (length(vgr) == 0){
                return (GRangesList())
            }

            ## local function that turns old VCF to BND
            .vcf2bnd = function(vgr){
                if (!"END" %in% colnames(values(vgr))){
                    jerror("Non BND SV should have the second breakpoint coor in END columns!")
                }

                if (!"CHR2" %in% colnames(values(vgr)) | any(is.na(vgr$CHR2))){
                    vgr$CHR2 = as.character(seqnames(vgr))
                }

                bp2 = data.table(as.data.frame(mcols(vgr)))
                bp2[, ":="(seqnames=CHR2, start=as.numeric(END), end=as.numeric(END))]
                bp2.gr = dt2gr(bp2, seqlengths = seqlengths(vgr))
                mcols(bp2.gr) = mcols(vgr)

                if (!is.null(names(vgr)) & !anyDuplicated(names(vgr))){
                    jid = names(vgr)
                } else {
                    jid = seq_along(vgr)
                }
                names(vgr) = paste(paste0("exp", jid), "1", sep=":")
                names(bp2.gr) = paste(paste0("exp", jid), "2", sep=":")

                vgr=resize(c(vgr, bp2.gr), 1)

                if (all(grepl("[_:][12]$",names(vgr)))){
                    ## row naming same with Snowman
                    nm <- vgr$MATEID <- names(vgr)
                    ix <- grepl("1$",nm)
                    vgr$MATEID[ix] = gsub("(.*?)(1)$", "\\12", nm[ix])
                    vgr$MATEID[!ix] = gsub("(.*?)(2)$", "\\11", nm[!ix])
                    vgr$SVTYPE="BND"
                }
                return(verify.junctions(vgr))
            }

            ## TODO: Delly and Novobreak
            ## fix mateids if not included
            if (!"MATEID" %in% colnames(mcols(vgr))) {
                ## TODO: don't assume every row is a different junction
                ## Novobreak, I'm looking at you.
                ## now delly...
                ## if SVTYPE is BND but no MATEID, don't pretend to be
                if (length(fake.bix <- which(values(vgr)$SVTYPE=="BND"))!=0){
                    values(vgr[fake.bix])$SVTYPE = "TRA"
                }

                ## add row names just like Snowman
                if (all(names(vgr)=="N" | ## Novobreak
                        is.null(names(vgr)) |
                        all(grepl("^DEL|DUP|INV|BND", names(vgr)))) ## Delly
                    ){
                    ## otherwise if all "N", as Novobreak
                    ## or starts with DEL|DUP|INV|BND, as Delly
                    ## expand and match MATEID
                    vgr=.vcf2bnd(vgr)
                }
            } else if (any(is.na(mid <- as.character(vgr$MATEID)))){
                ## like Lumpy, the BND rows are real BND but blended with non-BND rows
                ## treat them separately
                if (is.null(vgr$CHR2)){
                    vgr$CHR2 = as.character(NA)
                }

                names(vgr) = gsub("_", ":", names(vgr))
                vgr$MATEID = sapply(vgr$MATEID, function(x) gsub("_", ":", x))

                values(vgr) = data.table(as.data.frame(values(vgr)))

                ## break up the two junctions in one INV line!
                if ("STRANDS" %in% colnames(mc) & any(ns <- sapply(vgr$STRANDS, length)>1)){
                    ## first fix format errors, two strand given, but not comma separeted
                    ## so you'd have taken them as single
                    if (any(fuix <- sapply(vgr[which(!ns)]$STRANDS, stringr::str_count, ":")>1)){
                        which(!ns)[fuix] -> tofix
                        vgr$STRANDS[tofix] = lapply(vgr$STRANDS[tofix],
                                                    function(x){
                                                        strsplit(gsub("(\\d)([\\+\\-])", "\\1,\\2", x), ",")[[1]]
                                                    })
                        ns[tofix] = TRUE
                    }

                    ## for the one line two junction cases
                    ## split into two lines
                    vgr.double = vgr[which(ns)]
                    j1 = j2 = vgr.double
                    st1 = lapply(vgr.double$STRANDS, function(x)x[1])
                    st2 = lapply(vgr.double$STRANDS, function(x)x[2])
                    j1$STRANDS = st1
                    j2$STRANDS = st2
                    vgr.double = c(j1, j2)
                    names(vgr.double) = dedup(names(vgr.double))
                    vgr = c(vgr[which(!ns)], vgr.double)
                }

                mid <- as.logical(sapply(vgr$MATEID, length))
                vgr.bnd = vgr[which(mid)]
                mc.bnd = data.table(as.data.frame(values(vgr.bnd)))
                mc.bnd$MATEID = as.character(mc.bnd$MATEID)

                vgr.nonbnd = vgr[which(!mid)]
                if (length(loose.ix <- which(vgr.nonbnd$FILTER=="LOOSEEND"))>0){
                    ## Non-BND rows contains loose ends
                    vgr.loose = vgr.nonbnd[loose.ix]
                    vgr.nonbnd = vgr.nonbnd[setdiff(seq_along(vgr.nonbnd), loose.ix)]
                }
                
                if (length(vgr.nonbnd)>0){
                    vgr.nonbnd = .vcf2bnd(vgr.nonbnd)
                    mc.nonbnd = data.table(as.data.frame(values(vgr.nonbnd)))
                    vgr = c(vgr.bnd[,c()], vgr.nonbnd[,c()])
                    values(vgr) = rbind(mc.bnd, mc.nonbnd)
                }
            }

            ## sanity check
            if (!any(c("MATEID", "SVTYPE") %in% colnames(mcols(vgr)))){
                jerror("MATEID or SVTYPE not included. Required")
            }

            vgr$mateid = vgr$MATEID
            ## what's this???
            vgr$svtype = vgr$SVTYPE

            if (!is.null(info(vcf)$SCTG)){
                vgr$SCTG = info(vcf)$SCTG
            }

            if (force.bnd){
                vgr$svtype = "BND"
            }

            if (sum(vgr$svtype == 'BND')==0){
                jwarning('Vcf not in proper format.  Will treat rearrangements as if in BND format')
            }

            if (!all(vgr$svtype == 'BND')){
                jwarning(sprintf('%s rows of vcf do not have svtype BND, treat them as non-BND!',
                                sum(vgr$svtype != 'BND')))
            }

            bix = which(vgr$svtype == "BND")
            vgr = vgr[bix]
            alt <- sapply(vgr$ALT, function(x) x[1])

            ## Determine each junction's orientation
            if ("CT" %in% colnames(mcols(vgr))){
                jmessage("CT INFO field found.")
                if ("SVLEN" %in% colnames(values(vgr))){
                    ## proceed as Novobreak
                    ## ALERT: overwrite its orientation!!!!
                    del.ix = which(vgr$SVTYPE=="DEL")
                    dup.ix = which(vgr$SVTYPE=="DUP")
                    vgr$CT[del.ix] = "3to5"
                    vgr$CT[dup.ix] = "5to3"
                }

                ## also, Delly is like this
                ori = strsplit(vgr$CT, "to")
                iid = sapply(strsplit(names(vgr), ":"), function(x)as.numeric(x[2]))
                orimap = setNames(c("+", "-"), c("5", "3"))
                strd = orimap[sapply(seq_along(ori), function(i) ori[[i]][iid[i]])]
                strand(vgr) = strd
                vgr.pair1 = vgr[which(iid==1)]
                vgr.pair2 = vgr[which(iid==2)]
            } else if ("STRANDS" %in% colnames(mcols(vgr))){
                ## TODO!!!!!!!!!!!!!!!
                ## sort by name, record bp1 or bp2
                jmessage("STRANDS INFO field found.")
                iid = sapply(strsplit(names(vgr), ":"), function(x)as.numeric(x[2]))
                vgr$iid = iid
                vgr = vgr[order(names(vgr))]
                iid = vgr$iid

                ## get orientations
                ori = strsplit(substr(unlist(vgr$STRANDS), 1, 2), character(0))
                orimap = setNames(c("+", "-"), c("-", "+"))

                ## map strands
                strd = orimap[sapply(seq_along(ori), function(i) ori[[i]][iid[i]])]
                strand(vgr) = strd

                vgr.pair1 = vgr[which(iid==1)]
                vgr.pair2 = vgr[which(iid==2)]
            }
            else if (any(grepl("\\[|\\]", alt))){
                jmessage("ALT field format like BND")
                ## proceed as Snowman
                vgr$first = !grepl('^(\\]|\\[)', alt) ## ? is this row the "first breakend" in the ALT string (i.e. does the ALT string not begin with a bracket)
                vgr$right = grepl('\\[', alt) ## ? are the (sharp ends) of the brackets facing right or left
                vgr$coord = as.character(paste(seqnames(vgr), ':', start(vgr), sep = ''))
                vgr$mcoord = as.character(gsub('.*(\\[|\\])(.*\\:.*)(\\[|\\]).*', '\\2', alt))
                vgr$mcoord = gsub('chr', '', vgr$mcoord)

                ## add extra genotype fields to vgr
                if (all(is.na(vgr$mateid)))
                    if (!is.null(names(vgr)) & !any(duplicated(names(vgr)))){
                        jwarning('MATEID tag missing, guessing BND partner by parsing names of vgr')
                        vgr$mateid = paste(gsub('::\\d$', '', names(vgr)),
                        (sapply(strsplit(names(vgr), '\\:\\:'), function(x) as.numeric(x[length(x)])))%%2 + 1, sep = '::')
                    }
                    else if (!is.null(vgr$SCTG))
                    {
                        jwarning('MATEID tag missing, guessing BND partner from coordinates and SCTG')
                        ## require(igraph)
                        ucoord = unique(c(vgr$coord, vgr$mcoord))
                        vgr$mateid = paste(vgr$SCTG, vgr$mcoord, sep = '_')

                        if (any(duplicated(vgr$mateid)))
                        {
                            warning('DOUBLE WARNING! inferred mateids not unique, check VCF')
                            bix = bix[!duplicated(vgr$mateid)]
                            vgr = vgr[!duplicated(vgr$mateid)]
                        }
                    }
                    else{
                        jerror('Error: MATEID tag missing')
                    }

                vgr$mix = as.numeric(match(vgr$mateid, names(vgr)))

                pix = which(!is.na(vgr$mix))

                vgr.pair = vgr[pix]

                if (length(vgr.pair)==0){
                    jerror('Error: No mates found despite nonzero number of BND rows in VCF')
                }

                vgr.pair$mix = match(vgr.pair$mix, pix)

                vix = which(seq_along(vgr.pair)<vgr.pair$mix)
                vgr.pair1 = vgr.pair[vix]
                vgr.pair2 = vgr.pair[vgr.pair1$mix]

                ## now need to reorient pairs so that the breakend strands are pointing away from the breakpoint

                ## if "first" and "right" then we set this entry "-" and the second entry "+"
                tmpix = vgr.pair1$first & vgr.pair1$right
                if (any(tmpix)){
                    strand(vgr.pair1)[tmpix] = '-'
                    strand(vgr.pair2)[tmpix] = '+'
                }

                ## if "first" and "left" then "-", "-"
                tmpix = vgr.pair1$first & !vgr.pair1$right
                if (any(tmpix)){
                    strand(vgr.pair1)[tmpix] = '-'
                    strand(vgr.pair2)[tmpix] = '-'
                }

                ## if "second" and "left" then "+", "-"
                tmpix = !vgr.pair1$first & !vgr.pair1$right
                if (any(tmpix)){
                    strand(vgr.pair1)[tmpix] = '+'
                    strand(vgr.pair2)[tmpix] = '-'
                }

                ## if "second" and "right" then "+", "+"
                tmpix = !vgr.pair1$first & vgr.pair1$right
                if (any(tmpix)){
                    strand(vgr.pair1)[tmpix] = '+'
                    strand(vgr.pair2)[tmpix] = '+'
                }

                pos1 = as.logical(strand(vgr.pair1)=='+') ## positive strand junctions shift left by one (i.e. so that they refer to the base preceding the break for these junctions
                if (any(pos1)){
                    start(vgr.pair1)[pos1] = start(vgr.pair1)[pos1]-1
                    end(vgr.pair1)[pos1] = end(vgr.pair1)[pos1]-1
                }

                pos2 = as.logical(strand(vgr.pair2)=='+') ## positive strand junctions shift left by one (i.e. so that they refer to the base preceding the break for these junctions
                if (any(pos2)){
                    start(vgr.pair2)[pos2] = start(vgr.pair2)[pos2]-1
                    end(vgr.pair2)[pos2] = end(vgr.pair2)[pos2]-1
                }
            }

            ra = grl.pivot(GRangesList(vgr.pair1[, c()], vgr.pair2[, c()]))

            ## ALERT: vgr has already been subsetted to only include BND rows
            ## bix is the original indices, so NOT compatible!
            ## this.inf = values(vgr)[bix[pix[vix]], ]
            if (exists("pix") & exists("vix")){
                this.inf = values(vgr)[pix[vix], ]
            }
            if (exists("iid")){
                this.inf = values(vgr[which(iid==1)])
            }

            if (is.null(this.inf$POS)){
                this.inf = cbind(data.frame(POS = ''), this.inf)
            }
            if (is.null(this.inf$CHROM)){
                this.inf = cbind(data.frame(CHROM = ''), this.inf)
            }

            if (is.null(this.inf$MATL)){
                this.inf = cbind(data.frame(MALT = ''), this.inf)
            }

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

            if (is.null(values(ra)$TIER)){
                ## baseline tiering of PASS vs non PASS variants
                ## ALERT: mind the naming convention by diff programs
                ## TODO: make sure it is compatible with Delly, Novobreak, Meerkat
                ## Snowman/SvABA uses "PASS"
                ## Lumpy/Speedseq uses "."
                values(ra)$tier = ifelse(values(ra)$FILTER %in% c(".", "PASS"), 2, 3)
            } else {
                values(ra)$tier = values(ra)$TIER
            }

            if (geno==TRUE){
                ## expand into a list of GRLs, named by the sample name in the VCF
                geno.dt = data.table(
                    data.table(as.data.frame(VariantAnnotation::geno(vcf)$GT))
                )
                if (ncol(geno.dt)>1) {
                    cnms = colnames(geno.dt)
                    single.ra = ra
                    ra = lapply(setNames(cnms, cnms),
                                function(cnm){
                                    this.ra = copy(single.ra)
                                    this.dt = data.table(as.data.frame(values(this.ra)))
                                    this.geno = geno.dt[[cnm]]
                                    this.dt[
                                      , tier := ifelse(
                                            tier==2, ifelse(grepl("1", this.geno), 2, 3), 3)]
                                    values(this.ra) = this.dt
                                    return(this.ra)
                                })
                    loose=FALSE ## TODO: temporary until we figure out how
                }
            }

            if (!get.loose | is.null(vgr$mix)){
                return(verify.junctions(ra))
            } else {
                npix = is.na(vgr$mix)
                ## these are possible "loose ends" that we will add to the segmentation
                vgr.loose = vgr[npix, c()] 

                ## NOT SURE WHY BROKEN
                tmp =  tryCatch( values(vgr)[bix[npix], ],
                                error = function(e) NULL)
                if (!is.null(tmp)){
                    values(vgr.loose) = tmp
                } else{
                    values(vgr.loose) = cbind(vcf@fixed[bix[npix], ], info(vcf)[bix[npix], ])
                }

                return(list(junctions = verify.junctions(ra), loose.ends = vgr.loose))
            }
        } else{
            rafile = data.table::fread(rafile)
        }
    } else if (is.na(rafile)){
        return(GRangesList())
    }

    if (is.data.table(rafile)){
        rafile = as.data.frame(rafile)
    }

    if (nrow(rafile)==0){
        out = GRangesList()
        values(out) = rafile
        return(verify.junctions(out))
    }
    
    ## flip breaks so that they are pointing away from junction
    if (flipstrand) {
        rafile$str1 = ifelse(rafile$strand1 == '+', '-', '+')
        rafile$str2 = ifelse(rafile$strand2 == '+', '-', '+')
    }

    if (!is.null(seqlevels)) ## convert seqlevels from notation in tab delim file to actual
    {
        rafile$chr1 = seqlevels[rafile$chr1]
        rafile$chr2 = seqlevels[rafile$chr2]
    }


    if (is.null(rafile$str1)){
        rafile$str1 = rafile$strand1
    }

    if (is.null(rafile$str2)){
        rafile$str2 = rafile$strand2
    }

    if (!is.null(rafile$pos1) & !is.null(rafile$pos2)){
        if (breakpointer){
            rafile$pos1 = rafile$T_BPpos1
            rafile$pos2 = rafile$T_BPpos2
        }

        if (!is.numeric(rafile$pos1)){
            rafile$pos1 = as.numeric(rafile$pos1)
        }

        if (!is.numeric(rafile$pos2)){
            rafile$pos2 = as.numeric(rafile$pos2)
        }

        ## clean the parenthesis from the string

        rafile$str1 <- gsub('[()]', '', rafile$str1)
        rafile$str2 <- gsub('[()]', '', rafile$str2)

        ## goal is to make the ends point <away> from the junction where - is left and + is right
        if (is.character(rafile$str1) | is.factor(rafile$str1)){
            rafile$str1 = gsub('0', '-', gsub('1', '+', gsub('\\-', '1', gsub('\\+', '0', rafile$str1))))
        }

        if (is.character(rafile$str2) | is.factor(rafile$str2)){
            rafile$str2 = gsub('0', '-', gsub('1', '+', gsub('\\-', '1', gsub('\\+', '0', rafile$str2))))
        }


        if (is.numeric(rafile$str1)){
            rafile$str1 = ifelse(rafile$str1>0, '+', '-')
        }

        if (is.numeric(rafile$str2)){
            rafile$str2 = ifelse(rafile$str2>0, '+', '-')
        }

        rafile$rowid = 1:nrow(rafile)

        bad.ix = is.na(rafile$chr1) | is.na(rafile$chr2) | is.na(rafile$pos1) | is.na(rafile$pos2) | is.na(rafile$str1) | is.na(rafile$str2) | rafile$str1 == '*'| rafile$str2 == '*' | rafile$pos1<0 | rafile$pos2<0

        rafile = rafile[which(!bad.ix), ]

        if (nrow(rafile)==0){
            return(GRangesList())
        }

        seg = rbind(data.frame(chr = rafile$chr1, pos1 = rafile$pos1, pos2 = rafile$pos1, strand = rafile$str1, ra.index = rafile$rowid, ra.which = 1, stringsAsFactors = F),
                    data.frame(chr = rafile$chr2, pos1 = rafile$pos2, pos2 = rafile$pos2, strand = rafile$str2, ra.index = rafile$rowid, ra.which = 2, stringsAsFactors = F))

        if (chr.convert){
            seg$chr = gsub('chr', '', gsub('25', 'M', gsub('24', 'Y', gsub('23', 'X', seg$chr))))
        }

        out = seg2gr(seg, seqlengths = seqlengths)[, c('ra.index', 'ra.which')];
        out = split(out, out$ra.index)
    } else if (!is.null(rafile$start1) & !is.null(rafile$start2) & !is.null(rafile$end1) & !is.null(rafile$end2)){
        ra1 = gr.flipstrand(GRanges(rafile$chr1, IRanges(rafile$start1, rafile$end1), strand = rafile$str1))
        ra2 = gr.flipstrand(GRanges(rafile$chr2, IRanges(rafile$start2, rafile$end2), strand = rafile$str2))
        out = grl.pivot(GRangesList(ra1, ra2))
    }

    if (keep.features){
        values(out) = rafile[, ]
    }

    ## if (!is.null(pad)){
    ##     out = ra.dedup(out, pad = pad)
    ## }
    out = verify.junctions(out)
    if (!get.loose){
        return(out)
    } else{
        return(list(junctions = out, loose.ends = GRanges()))
    }

    return(out)
    ## return(new("junctions", out))
}

#' @name filter_oob_junctions
#' @rdname internal
#' @details
#'
#' Remove any out-of-range junctions
#'
#' A junction will be removed if:
#' - an endpoint exceeds seqlength of chromosome
#' - the start point is less than 1
#'
#' @param ra GRangesList object to be verified
#' @return GRangesList
filter_oob_junctions = function(ra) {
    pivoted.ra = grl.pivot(ra)
    sl = seqlengths(ra)
    bp1 = pivoted.ra[[1]]
    bp2 = pivoted.ra[[2]]
    start.oob = start(bp1) < 1 | start(bp2) < 1
    end.oob = start(bp1) > sl[as.character(seqnames(bp1))] | start(bp2) > sl[as.character(seqnames(bp2))]
    names(end.oob) = NULL
    names(start.oob) = NULL
    start.oob[which(is.na(start.oob))] = FALSE
    end.oob[which(is.na(end.oob))] = FALSE
    if (any(start.oob, na.rm = TRUE)) {
        jwarning("Some junction breakpoints are < 1 and will be removed")
        jwarning("Number of affected junctions: ", sum(start.oob, na.rm = T))
    }
    if (any(end.oob, na.rm = TRUE)) {
        jwarning("Some junction breakpoints are > seqlength and will be removed")
        jwarning("Number of affected junctions: ", sum(end.oob, na.rm = T))
    }
    return(ra[which(start.oob == FALSE & end.oob == FALSE, useNames = FALSE)])
}


#' @name verify.junctions
#' @rdname internal
#' @details
#' Produce error if the input does not meet any of the four following criteria
#' 1) is a GRangesList
#' 2) every element is length 2
#' 3) all elements are strand-specific
#' 4) all elements are width 1
#' @param ra rearrangement junctions object to be verified
#' @return the input if meets all four criteria
verify.junctions = function(ra){
    .ra.stop = function(ra){
        stopifnot(inherits(ra, "GRangesList"))
        if (length(ra)>0){
            stopifnot(all(elementNROWS(ra)==2))
            stopifnot(all(as.character(strand(unlist(ra))) %in% c("+", "-")))
            ## stopifnot(all(as.numeric(width(unlist(ra)))==1)) ## cancel this requirement but need to handle it later
        }
    }
    if (inherits(ra, "GRangesList")){
        .ra.stop(ra)
    } else if (inherits(ra, "list")) {
        for (i in seq_along(ra)){
            .ra.stop(ra[[i]])
        }
    }
    ## stopifnot(inherits(ra, "GRangesList"))
    ra = filter_oob_junctions(ra)
    return(ra)
}

#' @name karyograph
#' @title karyograph
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
                      label.edges = FALSE)
{
    if (length(junctions)>0)
    {
        bp.p = grl.pivot(junctions)
        bp1 = suppressWarnings(gr.end(gr.fix(bp.p[[1]]), 1, ignore.strand = F))
        bp2 = suppressWarnings(gr.start(gr.fix(bp.p[[2]]), 1, ignore.strand = F))


        ## #' mimielinski Sunday, Aug 06, 2017 06:46:15 PM
        ## #' fix added to handle strange [0 0] junctions outputted
        ## #' by Snowman ... which failed to match to any tile
        ## #' todo: may want to also handle junctions that
        ## #' fall off the other side of the chromosome
        #' this should now be unnecessary due to junction filtering
        ## end(bp1) = pmax(1, end(bp1))
        ## start(bp1) = pmax(1, start(bp1))
        ## end(bp1) = pmax(1, end(bp1))
        ## start(bp1) = pmax(1, start(bp1))

        if (any(as.logical(strand(bp1) == '*') | as.logical(strand(bp2) == '*')))
            jerror('bp1 and bp2 must be signed intervals (i.e. either + or -)')

        if (length(bp1) != length(bp2))
            jerror('bp1 and bp2 inputs must have identical lengths')

        ## #    if (sum(width(reduce(bp1))) != sum(width(bp1)) | sum(width(reduce(bp2))) != sum(width(bp2)))
        ## #      jerror('bp1 or bp2 cannot have duplicates / overlaps (with respect to location AND strand)')

        values(bp1)$bp.id = seq_along(bp1);
        values(bp2)$bp.id = seq_along(bp1)+length(bp1);

        pgrid = sgn1 = c('-'=-1, '+'=1)[as.character(strand(bp1))]
        sgn2 = c('-'=-1, '+'=1)[as.character(strand(bp2))]
### HACK HACK to force seqlengths to play with each other if malformedo
        tmp.sl = GenomeInfoDb::seqlengths(grbind(bp1, bp2))
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
                jwarning('Empty input given, producing empty output')
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
    else {
        tbp = NULL
    }

    if (length(junctions)>0){
        if (length(tbp)>0)
            g = gaps(gr.stripstrand(sort(c(bp1[, c()], bp2[, c()], tbp[, c()]))))
        else
            g = gaps(gr.stripstrand(sort(c(bp1[, c()], bp2[, c()]))))
    } else {
        g = gaps(gr.stripstrand(sort(tbp)))
    }

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
    tile$is.tel = start(tile)==1 | end(tile) == GenomeInfoDb::seqlengths(tile)[as.character(seqnames(tile))]
    values(tile)$tile.id = seq_along(tile);

    ## find "breakpoint" i.e. bp associated intervals, i.e. width 1 intervals that end with a bp1 or bp2 location
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

        ## # collect all pairwise adjacencies implied by breakpoints
        ## # eg imagine a|bp1|b
        ## #            c|bp2|d
        ## # "+" bp point to the right (eg b or d), "-" bp point to the left (a or c)

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
        ab.pairs[pp,1] = -ab.pairs[pp,1] ## # ++ breakpoints --> (-b)d adjacency
        ab.pairs[mm,2] = -ab.pairs[mm,2] ## -- breakpoints --> a(-c) adjacency
        ab.pairs[mp, ] = -ab.pairs[mp, ] ## +- breakpoints --> (-b)(-c) adjacency

        ## # clean up adj pairs
        ## # remove any that have crossed a chromosome boundary from their breakpoint
        ## # this will occur in cases of badly formed breakpoint input (eg breakpoints that point outward
        ## # from their telomeres)
        edge.id = rep(1:nrow(ab.pairs), 2)
        ab.pairs = rbind(ab.pairs, cbind(-ab.pairs[,2], -ab.pairs[,1]));
        ab.pairs.bpid = c(ab.pairs.bpid, ab.pairs.bpid)

        ## # build "aberrant" adjacency matrix representing directed graph of edges connecting
        ## # <signed> nodes.
        ## # note: indices of matrix represent edge labels
        adj.ab = Matrix( 0,
                        nrow = 2*length(tile),
                        ncol = 2*length(tile),
                        dimnames = rep( list( as.character(c(seq_along(tile), -(seq_along(tile))))), 2))
        tmp.ix = cbind(match(as.character(ab.pairs[,1]), rownames(adj.ab)),
                       match(as.character(ab.pairs[,2]), colnames(adj.ab)))

        ## added !is.na filters 8/4/2020 MI to take care of literal edge cases where edges point off the chromosome eg negative
        ## coordinates
        tix = !duplicated(tmp.ix) & !is.na(tmp.ix[,1]) & !is.na(tmp.ix[,2])
        adj.ab[tmp.ix[tix, , drop = F]] = ab.pairs.bpid[tix]
    }
    else
    {
        ab.pairs.bpid = edge.id = c()
        ab.pairs = matrix(nrow = 0, ncol = 2);
        adj.ab = Matrix(FALSE, nrow = 2*length(tile), ncol = 2*length(tile),
                        dimnames = rep(list(as.character(c(seq_along(tile), -(seq_along(tile))))), 2))
    }

    ## # build reference adjacency matrix (representing consecutive segments on the reference genome)
    ## # note: indices of matrix represent edge labels
    seg.ix = seq_along(tile)
    ref.pairs = cbind(seg.ix[1:(length(seg.ix)-1)],
                      seg.ix[2:(length(seg.ix))])
    ## # ref.pairs = ref.pairs[ref.pairs[,1]>0 & ref.pairs[,2]!=length(tile), ]
    ref.pairs = ref.pairs[which(as.character(seqnames(tile[ref.pairs[,1]])) ==
                                as.character(seqnames(tile[ref.pairs[,2]]))), ]

    ## XT fix 08/12: edge.id could be length 0 when no aberrant junction is used
    ## we should still make the ref edges in that case
    ## if (nrow(ref.pairs)>0 & length(edge.id)>0)
    if (nrow(ref.pairs)>0)
    {
        edge.id = c(edge.id, max(edge.id) + rep(1:nrow(ref.pairs), 2))
        ref.pairs = rbind(ref.pairs, cbind(-ref.pairs[,2], -ref.pairs[,1])) # reverse ref pairs
        adj.ref = Matrix(0, nrow = 2*length(tile), ncol = 2*length(tile),
                         dimnames = rep(list(as.character(c(seq_along(tile), -(seq_along(tile))))), 2))
        adj.ref[cbind(match(as.character(ref.pairs[,1]), rownames(adj.ref)),
                      match(as.character(ref.pairs[,2]), colnames(adj.ref)))] = nrow(ab.pairs)+1:nrow(ref.pairs)
    }
    else
    {
        adj.ref = Matrix(FALSE, nrow = 2*length(tile), ncol = 2*length(tile),
                         dimnames = rep(list(as.character(c(seq_along(tile), -(seq_along(tile))))), 2))
    }

    ## current tile is partition of genome only in positive orientation + dummy intervals for breakpoints
    ## output tile is forward partition and followed by reverse partition
    tmp.nm = as.character(c(seq_along(tile), -(seq_along(tile))))
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
    if (length(ab.ix)>0)
      {
        E(G)$eid[ab.ix] = edge.id[adj.ab[cbind(E(G)$from[ab.ix], E(G)$to[ab.ix])]]
        E(G)$eid[!ab.ix] = edge.id[adj.ref[cbind(E(G)$from[!ab.ix], E(G)$to[!ab.ix])]]
      }
    values(tile) = values(tile)[, c('tile.id', 'is.tel')]
    tile$ab.source = seq_along(tile) %in% E(G)$from[ab.ix]
    tile$ab.target = seq_along(tile) %in% E(G)$to[ab.ix]

    ## # important: map input ra to aberrant graph edges, i.e. ab.edges matrix with $from $to and $edge.ix columns
    ## # and one row for each aberrant edge
    ab.edges = array(NA, dim = c(length(junctions), 3, 2), dimnames = list(NULL, c('from', 'to', 'edge.ix'), c('+', '-')))
    dupped = duplicated(ab.pairs.bpid) ## but there are still duplicated ones
    ab.edges[,1:2,1] = cbind(match(ab.pairs[!dupped,1], names(tile)),
                             match(ab.pairs[!dupped,2], names(tile)))
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


## FIXME: problem too large

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
#' @noRd
#########################################
jabba2vcf = function(jab, fn = NULL, sampleid = 'sample', hg = NULL, include.loose = TRUE, include.cn0 = TRUE, cnv = FALSE)
{
    if (is.null(hg))
        hg = tryCatch(skidb::read_hg(), error = function(e) NULL)

    vcffields = c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'GENO')

    ## convert all aberrant connections into pairs of VCF rows
    if (!cnv)
    {
      jix = which(!is.na(jab$ab.edges[,3,1])) ## these are the only junctions with breaks in the reconstruction
      if (!include.cn0) ## remove from jix
      {
        jcn = jab$adj[jab$ab.edges[jix, 1:2, 1]]
        jix = jix[jcn>0]
        message('Removing cn=0')
      }
      abs = rbind(jab$ab.edges[jix,1:2,1])
      rabs = rbind(jab$ab.edges[jix,1:2,2])
      rcix = match(jab$segstats, gr.flipstrand(jab$segstats)) ## map of seg to its reverse complement
      
      adj.ref = jab$adj ## reference graph has reference copy numbers, we obtain by zeroing out all ab.edges and loose end edges
      adj.ref[rbind(jab$ab.edges[jix,1:2,1])] = 0
      adj.ref[rbind(jab$ab.edges[jix,1:2,2])] = 0
      
      ## #' xtYao #' Wednesday, Mar 20, 2019 11:09:46 AM
      ## Fix missing $
      if (any(jab$segstats$loose))
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
            gr1$ID = paste(sampleid, '_seg', jab$segstats$tile.id[abs[,1]], ifelse(as.logical(strand(gr1)=='+'), '_R', '_L'), sep = '')

            gr2 = gr.start(jab$segstats[abs[,2]], ignore.strand = F)[, 'cn']
            gr2$jid = jix
            gr2$nid = abs[,2]
            gr2$acn = jcn
            gr2$rcn = Matrix::colSums(adj.ref[,gr2$nid, drop = FALSE])
            gr2$ID = paste(sampleid, '_seg', jab$segstats$tile.id[abs[,2]], ifelse(as.logical(strand(gr2)=='+'), '_L', '_R'), sep = '')

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
                             ";JABID=", abs[,1], ";RJABID=", rcix[abs[,1]], ";JUNCID=", seq_along(gr1), sep = '')
            gr2$INFO = paste("SVTYPE=BND", ";MATEID=", gr2$mid, ";CNADJ=", gr2$acn, ";CNRADJ=", gr2$rcn, ";CN=", gr2$cn,
                             ";JABID=", abs[,2], ";RJABID=", rcix[abs[,2]], ";JUNCID=", seq_along(gr2), sep = '')

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
        if (length(lix)>0 & include.loose)
        {
            ## loose ends should be width 1, but just in case
            if (is.element("passed", colnames(values(jab$segstats)))){
                ## with le quality
                gr.loose = gr.start(jab$segstats[lix, c('passed', 'cn')]) 
                gr.loose$QUAL = ifelse(gr.loose$passed, "PASSED", "FAILED")
            } else {
                ## without
                gr.loose = gr.start(jab$segstats[lix, c('cn')]) 
                gr.loose$QUAL = '.'
            }

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
            gr.loose$ID = paste(sampleid, '_looseend', jab$segstats$tile.id[pcid], ifelse(isp, '_L', '_R'), sep = '')
            gr.loose$mid = NA
            gr.loose$REF = tryCatch(as.character(ffTrack::get_seq(hg, gr.stripstrand(gr.loose))), error = function(e) 'N')

            ## again, same rationale as above, parent = left loose end, child = right loose end
            gr.loose$ALT = ifelse(isp, paste('.', gr.loose$REF, sep = ''), paste(gr.loose$REF, '.', sep = ''))
            gr.loose$FORMAT = "GT:CN:RCN:SCN"
            gr.loose$GENO = paste(ifelse(gr.loose$rcn>0, '0/1', '1'), gr.loose$acn, gr.loose$rcn, gr.loose$cn, sep = ":")
            gr.loose$CHROM = as.character(seqnames(gr.loose))
            gr.loose$POS = start(gr.loose)
            gr.loose$FILTER = "LOOSEEND"
            
            gr.loose$INFO = paste("SVTYPE=BND", ";CNADJ=", gr.loose$acn, ";CNRADJ=", gr.loose$rcn, ";CN=", gr.loose$cn, ";JABID=", lix, ";RJABID=", rcix[lix], sep = '')
            gr.loose = gr.loose[, vcffields]
        }
        else
            gr.loose = GRanges()

        ## make header
        sl = GenomeInfoDb::seqlengths(jab$segstats)
        header = '##fileformat=VCFv4.2'
        header = c(header, sprintf('##fileDate=%s', format(Sys.Date(), '%Y%m%d')))
        header = c(header, '##source=JaBbAV0.1')
        if (inherits(hg, "BSgenome"))
        {
            header = c(header, sprintf("##reference=%s", BSgenome::sourceUrl(hg)))
            header = c(header, unlist(mapply(function(x, y) sprintf('##contig=<ID=%s,length=%s,assembly=%s,species="%s">', x, y, BSgenome::providerVersion(hg), BSgenome::organism(hg)), names(sl), sl)))
        }
        else
      {
        if (!is.null(hg))
        {
          header = c(header, sprintf("##reference=%s", BSgenome::filename(hg)['rds']))
        }
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
                out = sapply(seq_along(unames), function(i) {
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
        sl = GenomeInfoDb::seqlengths(jab$segstats)
        header = '##fileformat=VCFv4.2'
        header = c(header, sprintf('##fileDate=%s', format(Sys.Date(), '%Y%m%d')))
        header = c(header, '##source=JaBbAV0.1 CNV')

        if (inherits(hg, "BSgenome"))
        {
            header = c(header, sprintf("##reference=%s", BSgenome::sourceUrl(hg)))
            header = c(header, unlist(mapply(function(x, y) sprintf('##contig=<ID=%s,length=%s,assembly=%s,species="%s">', x, y, BSgenome::providerVersion(hg), BSgenome::organism(hg)), names(sl), sl)))
        }
        else if (inherits(hg, "ffTrack"))
        {
            header = c(header, sprintf("##reference=%s", BSgenome::filename(hg)['rds']))
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

#' @name read_vcf
#' @rdname internal
#' @title read_vcf
#'
#' @description
#'
#' wrapper around variantAnnotation reads VCF into granges or data.table format
#'
#' @author Marcin Imielinski
#' @noRd
read_vcf = function(fn,
                    hg = "hg19",
                    swap.header = NULL,
                    verbose = FALSE,
                    add.path = FALSE,
                    tmp.dir = '~/temp/.tmpvcf',
                    ...)
{
    in.fn = fn

    if (verbose){
        cat('Loading', fn, '\n')
    }

    ##   if (!is.null(swap.header))
    ##   {
    ##     if (!file.exists(swap.header))
    ##       jerror(sprintf('Swap header file %s does not exist\n', swap.header))

    ##     system(paste('mkdir -p', tmp.dir))
    ##     tmp.name = paste(tmp.dir, '/vcf_tmp', gsub('0\\.', '', as.character(runif(1))), '.vcf', sep = '')
    ##     if (grepl('gz$', fn))
    ##       system(sprintf("zcat %s | grep '^[^#]' > %s.body", fn, tmp.name))
    ##     else
    ##       system(sprintf("grep '^[^#]' %s > %s.body", fn, tmp.name))

    ##     system(sprintf("cat %s.header %s.body > %s", tmp.name, tmp.name, tmp.name))
    ##     vcf = VariantAnnotation::readVcf(tmp.name, hg, ...)
    ##     system(sprintf("rm %s %s.body %s.header", tmp.name, tmp.name, tmp.name))
    ##   }
    ## else

    ## QUESTION: why isn't genome recognized when ... is provided?
    args = list(...)
    ## vcf = VariantAnnotation::readVcf(file = fn, hg, ...)
    vcf = VariantAnnotation::readVcf(file = fn, genome = hg)

    out = granges(vcf)

  if (!is.null(values(out)))
    values(out) = cbind(values(out), VariantAnnotation::info(vcf))
  else
    values(out) = VariantAnnotation::info(vcf)

    if (add.path)
        values(out)$path = in.fn

    return(out)
}

#' @name write_vcf
#' @title write_vcf
#'
#' @description
#'
#' writes any GRanges vars into vcf using columns of vars to guide choice of common fields like
#' $FILTER
#' $GT
#' $REF
#' $ALT
#'
#' and adding all other fields to INFO
#'
#' @noRd
#' @author Marcin Imielinski
write_vcf = function(vars, filename, sname = "mysample", info.fields = setdiff(names(values(vars)), c("FILTER", "GT", "REF", "ALT")))
{
    
    genoh = DataFrame(row.names = 'GT', Number = 1, Type = 'Float', Description = 'Genotypes')

    for (field in names(values(vars))) ## clean up vars of weird S4 data structures that are not compatible with before
    {
        tmp = tryCatch(as.character(values(vars)[, field]), error = function(e) NULL)
        if (is.null(tmp))
        {
            values(vars)[, field] = NULL
            warning(paste('Could not process field', field, "due to S4 conversion issues, discarding"))
        }
        is.num = !all(is.na(as.numeric(tmp)))
        if (!is.num)
            values(vars)[, field] = tmp
    }

    info.fields = intersect(info.fields, names(values(vars)))

    if (length(info.fields)==0) # dummy field to keep asVCF happy
    {
        info.fields = "DM"
        vars$DM = '.'
    }

    is.num = sapply(info.fields, function(x) !suppressWarnings(all(is.na(as.numeric(as.character(values(vars)[, x]))))))
    infoh = DataFrame(
        row.names = info.fields, Number = 1,
        Type = ifelse(is.num, 'Float', 'String'),
        Description = paste('Field', info.fields))

    if (is.null(vars$REF))
        vars$REF = vars$refbase

    if (is.null(vars$ALT))
        vars$ALT = vars$altbase

    if (is.null(vars$REF))
        vars$REF =  "N"

    if (is.null(vars$ALT))
        vars$ALT =  "X"

    if (is.null(vars$FILTER))
        vars$FILTER =  "PASS"

    ## vcf = asVCF(vars)
    vr = VRanges(seqnames(vars), ranges(vars), ref = vars$REF, alt = vars$ALT, sampleNames = rep(sname, length(vars)))
    names(vr) = names(vars)
    vcf = asVCF(vr)


    for (field in info.fields)
        info(vcf)[[field]] = values(vars)[[field]]

##    geno(vcf)$DP = vars$DP; geno(vcf)$AD = vars$AD; geno(vcf)$FT = vars$FT

    info(header(vcf)) = infoh

    if (is.null(vars$FILTER))
        filt(vcf) = rep('PASS', length(vars))
    else
        filt(vcf) = vars$FILTER

    ## xtYao ## Thursday, Feb 25, 2021 02:55:16 PM
    ## use the names if they are there!!
    if (!is.null(names(vars))){
        rownames(vcf) = names(vars)
    } else {
        rownames(vcf) = vars$assembly.coord ## WHY, WHY, WHY??????
    }
    
    geno(header(vcf)) = genoh

    geno(vcf)$GT = vcf$GT

    writeVcf(vcf, filename)
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
#' @noRd
levapply = function(x, by, FUN = 'order')
{
    if (!is.list(by))
        by = list(by)

    f = factor(do.call('paste', c(list(sep = '|'), by)))
    ixl = split(seq_along(x), f);
    ixv = lapply(ixl, function(y) x[y])
    res = structure(unlist(lapply(ixv, FUN)), names = unlist(ixl))
    out = rep(NA, length(x))
    out[as.numeric(names(res))] = res;
    return(out)
}

#' @name chr2num
#' @description
#' Strips chromosome numbers
#'
#' @param x character vector to strip
#' @param xy logical flag specifying whether to keep X and Y or convert to 23 and 24
#' @noRd
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
    bps = unname(grl.unlist(juncs))
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
#' @noRd
sv.size = function(juncs,
                   ...){
    bps = gUtils::grl.pivot(juncs)
    return(IRanges::distance(bps[[1]], bps[[2]], ...))
}

#' @name reciprocal.cycles
#' @rdname internal
#' @description
#' Returns indices (subset of seq_along(junc) corresponding to cycles of (quasi) reciprocal cycles
#' @param juncs GRangesList of junctions
#' @param mc.cores parallel
#' @param ignore.strand usually TRUE
#' @return numerical vector of the same length, Inf means they r not facing each other
#' @noRd
reciprocal.cycles = function(juncs, paths = FALSE, thresh = 1e3, mc.cores = 1, verbose = FALSE, chunksize = 1e3)
{
    bp = grl.unlist(juncs)[, c("grl.ix", "grl.iix")]

    ix = split(seq_along(bp), ceiling(runif(length(bp))*ceiling(length(bp)/chunksize)))
    ixu = unlist(ix)
    eps = 1e-9
    ij = do.call(rbind, split(seq_along(bp), bp$grl.ix))
    adj = sparseMatrix(1, 1, x = FALSE, dims = rep(length(bp), 2))

    ## matrix of (strand aware) reference distances between breakpoint pairs
    adj[ixu, ] = do.call(rbind, parallel::mclapply(ix,
                                         function(iix)
                                         {
                                             if (verbose)
                                                 cat('.')
                                             tmpm = gr.dist(bp[iix], gr.flipstrand(bp), ignore.strand = FALSE)+eps
                                             tmpm[is.na(tmpm)] = 0
                                             tmpm[tmpm>thresh] = 0
                                             tmpm = as(tmpm>0, 'Matrix')
                                         },
                                         mc.cores = mc.cores))
    if (verbose)
        cat('\n')

    adj = adj | Matrix::t(adj) ## symmetrize

    ## bidirected graph --> skew symmetric directed graph conversion
    ## split each junction (bp pair) into two nodes, one + and -
    ## arbitrarily call each bp1-->bp2 junction is "+" orientation
    ## then all odd proximities adjacent to bp1 will enter the "+"
    ## version of that junction and exit the "-" version

    ## new matrix will be same dimension as adj
    ## however the nodes will represents + and -
    ## orientation of junctions
    ## using the foollowing conversion

    ## i.e.
    ## bp2 --> bp1 + +
    ## bp2 --> bp2 + -
    ## bp1 --> bp1 - +
    ## bp1 --> bp2 - -

    ## we'll use the same indices just to keep things confusing
    junpos = bp1 = bp$grl.iix == 1
    junneg = bp2 = bp$grl.iix == 2

    adj2 = adj & FALSE ## clear out adj for new skew symmetric version
    adj2[junpos, junpos] = adj[bp2, bp1]
    adj2[junpos, junneg] = adj[bp2, bp2]
    adj2[junneg, junpos] = adj[bp1, bp1]
    adj2[junneg, junneg] = adj[bp1, bp2]

    ## strongly connected components consists of (possibly nested) cycles
    cl = split(seq_along(bp), igraph::clusters(graph.adjacency(adj2), 'strong')$membership)

    ## choose only clusters with length > 1
    cl = cl[S4Vectors::elementNROWS(cl)>1]
    cl = cl[order(S4Vectors::elementNROWS(cl))]


    jcl = lapply(cl, function(x) unique(sort(bp$grl.ix[x])))
    jcls = sapply(jcl, paste, collapse = ' ')
    jcl = jcl[!duplicated(jcls)]

    if (paths)
    {
        adj3 = adj2

        ## remove all cycles and enumerate remaining paths > 1
        adj3[unlist(jcl), unlist(jcl)] = FALSE
        sinks = which(rowSums(adj3)==0)
        sources = which(colSums(adj3)==0)

        cl2 = split(seq_along(bp), igraph::clusters(graph.adjacency(adj3), 'weak')$membership)
        cl2 = cl2[S4Vectors::elementNROWS(cl2)>1]

        if (any(ix <- S4Vectors::elementNROWS(cl2)>2))
        { ## only need to do this for connected components that have 3 or more junctions
            cl3 = do.call(c, parallel::mclapply(cl2[ix], function(x)
            {
                tmp.adj = adj3[x, x]
                lapply(all.paths(tmp.adj, sources = sources, sinks = sinks)$paths, function(i) x[i])
            }, mc.cores = mc.cores))

            cl2 = c(cl2[!ix], cl3)
        }
        jcl2 = lapply(cl2, function(x) unique(sort(bp$grl.ix[x])))
        jcls2 = sapply(jcl2, paste, collapse = ' ')
        jcl2 = jcl2[!duplicated(jcls2)]

        return(list(cycles = jcl, paths = jcl2))
    }

    return(jcl)
}

#' @name ra.merge
#' Merges rearrangements represented by \code{GRangesList} objects
#'
#' Determines overlaps between two or more piles of rearrangement junctions (as named or numbered arguments) +/- padding
#' and will merge those that overlap into single junctions in the output, and then keep track for each output junction which
#' of the input junctions it was "seen in" using logical flag  meta data fields prefixed by "seen.by." and then the argument name
#' (or "seen.by.ra" and the argument number)
#'
#' @param ... GRangesList representing rearrangements to be merged
#' @param pad non-negative integer specifying padding
#' @param ind  logical flag (default FALSE) specifying whether the "seen.by" fields should contain indices of inputs (rather than logical flags) and NA if the given junction is missing
#' @param ignore.strand whether to ignore strand (implies all strand information will be ignored, use at your own risk)
#' @return \code{GRangesList} of merged junctions with meta data fields specifying which of the inputs each outputted junction was "seen.by"
#' @examples
#'
#' # generate some junctions
#' gr1 <- GRanges(1, IRanges(1:10, width = 1), strand = rep(c('+', '-'), 5))
#' gr2 <- GRanges(1, IRanges(4 + 1:10, width = 1), strand = rep(c('+', '-'), 5))
#' ra1 = split(gr1, rep(1:5, each = 2))
#' ra2 = split(gr2, rep(1:5, each = 2))
#'
#' ram = ra.merge(ra1, ra2)
#' values(ram) # shows the metadata with TRUE / FALSE flags
#'
#' ram2 = ra.merge(ra1, ra2, pad = 5) # more inexact matching results in more merging
#' values(ram2)
#'
#' ram3 = ra.merge(ra1, ra2, ind = TRUE) #indices instead of flags
#' values(ram3)
#' @noRd
ra.merge = function(..., pad = 0, ind = FALSE, ignore.strand = FALSE){
    ra = list(...)
    ra = ra[which(!sapply(ra, is.null))]
    ra = ra[which(!sapply(ra, function(x) {length(x)==0}))] ## filter zero length junctions
    nm = names(ra)
    if (is.null(nm)){
        nm = paste('ra', seq_along(ra), sep = '')
    }
    nm = paste('seen.by', nm, sep = '.')
    if (length(nm)==0){
        return(NULL)
    }
    out = ra[[1]]
    values(out) = cbind(as.data.frame(matrix(FALSE, nrow = length(out), ncol = length(nm), dimnames = list(NULL, nm))), values(out))

    if (!ind){
        values(out)[, nm[1]] = TRUE
    } else{
        values(out)[, nm[1]] = seq_along(out)
    }

    if (length(ra)>1){
        for (i in 2:length(ra)){
            this.ra = ra[[i]]
            if (length(this.ra)>0){
                values(this.ra) = cbind(as.data.frame(matrix(FALSE, nrow = length(this.ra), ncol = length(nm), dimnames = list(NULL, nm))), values(this.ra))
                ovix = ra.overlaps(out, this.ra, pad = pad, ignore.strand = ignore.strand)

                if (!ind){
                    values(this.ra)[[nm[i]]] = TRUE
                } else{
                    values(this.ra)[[nm[i]]] = seq_along(this.ra)
                }

                if (!ind){
                    if (!all(is.na(ovix))){
                        values(out)[, nm[i]][ovix[,1]] = TRUE
                    }
                } else{
                    values(out)[, nm[i]] = NA
                    if (!all(is.na(ovix))){
                        values(out)[, nm[i]][ovix[,1]] = ovix[,1]
                    }
                }
                ## which are new ranges not already present in out, we will add these
                if (!all(is.na(ovix))){
                    nix = setdiff(seq_along(this.ra), ovix[,2])
                } else{
                    nix = seq_along(this.ra)
                }

                if (length(nix)>0){
                    val1 = values(out)
                    val2 = values(this.ra)
                    if (ind){
                        val2[, nm[1:(i-1)]] = NA
                    }
                    else{
                        val2[, nm[1:(i-1)]] = FALSE
                    }
                    values(out) = NULL
                    values(this.ra) = NULL
                    out = grl.bind(out, this.ra[nix])
                    d1 = as.data.table(val1)
                    d2 = as.data.table(val2[nix, ])
                    ## prevent column class mismatches
                    c1 = data.table(class = sapply(d1, class), cnm = colnames(d1))
                    c2 = data.table(class = sapply(d2, class), cnm = colnames(d2))
                    conflict = merge(c1, c2, by = "cnm")[class.x != class.y]
                    for (r in seq_len(nrow(conflict))){
                        d1[[conflict[r, cnm]]] = as(d1[[conflict[r, cnm]]], conflict[r, c(class.y)])
                    }
                    values(out) = rbind(d1, d2, fill = TRUE)
                }
            }
        }
    }
    return(out)
}


####################################################################
#' ppgrid
#'
#' least squares grid search for purity and ploidy modes
#'
#' @param segstats GRanges object of intervals with meta data fields "mean" and "sd" (i.e. output of segstats function)
#' @param allelic logical flag, if TRUE will also look for mean_high, sd_high, mean_low, sd_low variables and choose among top solutions from top copy number according to the best allelic fit
#' @param purity.min min purity value allowed
#' @param purity.max max purity value allowed
#' @param ploidy.min min ploidy value allowed
#' @param ploidy.max max ploidy value allowed
#' @param ploidy.step grid length of ploidy values
#' @param purity.step grid length of purity values
#' @param plot whether to plot the results to file
#' @param verbose print intermediate outputs
#' @param mc.cores integer number of cores to use (default 1)
#' @return data.frame with top purity and ploidy solutions and associated gamma and beta values, for use in downstream jbaMI
############################################
ppgrid = function(segstats,
                  allelic = FALSE,
                  purity.min = 0.01,
                  purity.max = 1.0,
                  ploidy.step = 0.01,
                  purity.step = 0.01,
                  ploidy.min = 1.2, # ploidy bounds (can be generous)
                  ploidy.max = 6,
                  plot = F,
                  verbose = F,
                  mc.cores = 1){
    if (verbose)
        jmessage('setting up ppgrid matrices .. \n')

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
            jwarning('If allelic = TRUE then must have meta data fields mean_high, mean_low, sd_high, sd_low in input segstats')
            allelic = FALSE
        }

    if (is.null(segstats$mean))
        jerror('segstats must have field $mean')

    segstats = segstats[!is.na(segstats$mean) & !is.na(segstats$sd)]

    if (!is.null(segstats$ncn))
        segstats = segstats[segstats$ncn==2, ]

    ## if (is.null(segstats$ncn))
    ##     ncn = rep(2, length(mu))
    ## else
    ##     ncn = segstats$ncn

    if (any(tmpix <-is.infinite(segstats$mean) | is.infinite(segstats$sd)))
    {
      segstats$sd[tmpix] = segstats$mean[tmpix] = NA
    }
    segstats = segstats[!is.na(segstats$mean) & !is.na(segstats$sd), ]
    if (length(segstats)==0)
      jerror('No non NA segments provided')

    mu = segstats$mean
    w = as.numeric(width(segstats))
    Sw = sum(as.numeric(width(segstats)))
    sd = segstats$sd
    m0 = sum(as.numeric(mu*w))/Sw

    if (verbose)
        cat(paste(c(rep('.', length(purity.guesses)), '\n'), collapse = ''))

    NLL = matrix(unlist(parallel::mclapply(seq_along(purity.guesses), function(i)
    {
        if (verbose)
            cat('.')
        nll = rep(NA, length(ploidy.guesses))
        for (j in seq_along(ploidy.guesses))
        {
            alpha = purity.guesses[i]
            tau = ploidy.guesses[j]
            gamma = 2/alpha - 2
            beta = (tau + gamma)/m0 
            v = pmax(0, round(beta*mu-gamma))
            nll[j] = sum((v-beta*mu+gamma)^2/((sd)^2))
        }
        return(nll)
    }, mc.cores = mc.cores)), nrow = length(purity.guesses), byrow = T)

    dimnames(NLL) = list(as.character(purity.guesses), as.character(ploidy.guesses))

  if (verbose) {
      cat('\n')
  }
    ## rix = as.numeric(rownames(NLL))>=purity.min & as.numeric(rownames(NLL))<=purity.max
    ## cix = as.numeric(colnames(NLL))>=ploidy.min & as.numeric(colnames(NLL))<=ploidy.max
    ## NLL = NLL[rix, cix, drop = FALSE]

    a = rep(NA, nrow(NLL));
    b = rep(NA, ncol(NLL)+2)
    b.inf = rep(Inf, ncol(NLL)+2)
    #  a = rep(Inf, nrow(NLL));
    #  b = rep(Inf, ncol(NLL)+2)
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
        M = (NLLc < NLLul &
             NLLc < NLLuc &
             NLLc < NLLur &
             NLLc < NLLcl &
             NLLc < NLLcr &
             NLLc < NLLll &
             NLLc < NLLlc &
             NLLc < NLLlr)[-c(1, nrow(NLLc)),
                           -c(1, ncol(NLLc)),
                           drop = FALSE]
    else if (ncol(NLL)==1) ## one column, only go up and down
        M = (NLLc < NLLuc & NLLc < NLLlc)[-c(1, nrow(NLLc)), -c(1, ncol(NLLc)), drop = FALSE]
    else ## only row, only go left right
        M = (NLLc < NLLcl & NLLc < NLLcr)[-c(1, nrow(NLLc)), -c(1, ncol(NLLc)), drop = FALSE]

    if (length(M)>1)
    {
        ix = Matrix::which(M, arr.ind=T);
        if (nrow(ix)>1)
        {
            C = hclust(d = dist(ix), method = 'single')
            cl = cutree(C, h = min(c(nrow(NLL), ncol(NLL), 2)))
            minima = ix[vaggregate(1:nrow(ix), by = list(cl), function(x) x[which.min(NLL[ix[x, drop = FALSE]])]), , drop = FALSE]
        }
        else if (nrow(ix) == 0) {
            ## if NLL is monotonically increaing or dereasing, minima will not be found
            ## in this case, return the local minimum of NLL over the tested grid
            minima = Matrix::which(NLL == min(NLL, na.rm = T), arr.ind = T)
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
          if (verbose)
            {
              jmessage(sprintf('Evaluating alleles for solution %s of %s\n', i, nrow(out)))
            }
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

            for (j in seq_along(vlow.mle))
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
    out = out.all;
    out = out[order(out$group, !out$keep, out$NLL), ]
    out$rank = NA
    out$rank[out$keep] = 1:sum(out$keep)
    out$keep = out$i = out$j = NULL
    rownames(out) = NULL
    return(out)
}



####################
#' @name arrstring
#' @title arrstring
#'
#' @description
#' string representation of row array as linear combination of nonzero entries
#' of that row using column names as variables
#'
#' @param A array
#' @param sep separator to use between table elements
#' @return character representation of table
#' @author Marcin Imielinski
####################
arrstring = function(A, x = NULL, sep = ', ', sep2 = '_', signif = 3, dt = FALSE)
{
  if (is.null(dim(A)))
  {
    A = rbind(A)
  }

  if (is.null(colnames(A)))
  {
    colnames(A) = paste0('V', 1:ncol(A))
  }

  if (is.null(x))
  {
    x = colnames(A)
  }
  else
  {
    x = signif(x, signif)
  }

  tmp = as.data.table(Matrix::which(A!=0, arr.ind = TRUE))
  tmp[,  y := A[cbind(row, col)]][  , x:= x[col]]

  str = tmp[, paste(signif(y,signif), x[y!=0], sep = '*', collapse = ' + '), keyby = row][list(1:nrow(A)), V1]

  return(str)
}




#' filter.loose
#'
#' analyze coverage surrounding given loose ends to evaluate quality
#' @param gg gGraph of JaBbA model
#' @param cov.rds character path to binned coverage data
#' @param l data.table of loose ends to evaluate
#' @param purity optional, fractional purity of sample, default assumes 1
#' @param ploidy optional, ploidy of sample, default infers from gg
#' @param field optional, column name in cov.rds, default="ratio"
#' @param PTHRESH optional, threshold for GLM p-value for calling true positive loose ends, default=3.4e-7 provides consanguinity with large dataset bonferroni correction
#' @param max.epgap optional, threshold for JaBbA MIQP convergence epgap. Values over this is considered incomplete optimization and disregarded.
#' @param verbose optional, default=FALSE
#' @return data.table containing a row for every input loose end and logical column `true.pos` indicating whether each loose end has passed all filters (TRUE) or not (FALSE)
#' @author Julie Behr
filter.loose = function(gg, cov, l, purity=NULL, ploidy=NULL, field="ratio", PTHRESH=3.4e-7, max.epgap = 1e-3, verbose=F){
    ## ## load coverage and beta (coverage CN fit)
    ## if(verbose) message("Loading coverage bins")
    ## cov = readRDS(cov.rds)
    ## cov = gr.sub(cov, "chr", "")
    if(!(field %in% colnames(values(cov)))) stop("must provide field in cov.rds")
    if(field != "ratio") cov$ratio = values(cov)[, field]
    if(!("tum.counts" %in% colnames(values(cov)))){
        yf = ifelse("reads.corrected" %in% colnames(values(cov)), "reads.corrected", field)
        cov$tum.counts = values(cov)[, yf]
    }
    if(!("norm.counts" %in% colnames(values(cov)))){
        cov$norm.counts = 1 ## dummy to make it flat
    }
    if(is.null(purity)){
        if (is.null(purity <- gg$meta$purity)){
            jerror("Purity must be given.")
        }
    }  ## purity = 1
    if(is.null(ploidy)) {
        ploidy = weighted.mean(gg$nodes$gr$cn, gg$nodes$gr %>% width, na.rm=T)
    }
    ## remove bins with infinite values up front
    cov = cov[which(is.finite(cov$ratio) & !is.na(cov$ratio))]
    ratios = cov$ratio
    beta = mean(ratios[is.finite(ratios)], na.rm=T) * purity/(2*(1-purity) + purity * ploidy)
    segs = gg$nodes$gr
    ## segs = gr.sub(segs, "chr", "")
    ## l = gr.sub(l, "chr", "")

    if (!is(l, "GRanges")){
        try({l = dt2gr(l)})
        if (inherits(l, "try-error")){
            jerror("l must be a GRanges or a data.table that can be converted to a GRanges.")
        }
    }

    ## identify nodes flanking each loose end, extending up to 100kb away
    o = gr.findoverlaps(segs, l+1)
    segs = segs[o$query.id]; segs$leix = l[o$subject.id]$leix
    sides = gr.findoverlaps(segs, l+1e5, by="leix")
    values(sides) = cbind(values(sides), values(l[sides$subject.id]))
    sides$fused = !is.na(gr.match(sides, l, by="leix"))
    sides$wid = width(sides)

    ## gather coverage bins corresponding to fused & unfused sides of loose ends
    if(verbose) jmessage("Overlapping coverage with loose end fused and unfused sides")

    rel = gr.findoverlaps(cov, sides)
    values(rel) = cbind(values(cov[rel$query.id]), values(sides[rel$subject.id]))
    qq = 0.05
    rel = gr2dt(rel)[, ":="(
                   in.quant.r = ratio >= quantile(ratio, qq, na.rm=T) & ratio <= quantile(ratio, 1-qq, na.rm=T),
                   good.cov=sum(is.na(tum.counts))/.N < 0.1 & sum(is.na(norm.counts))/.N < 0.1 & sum(is.na(ratio))/.N < 0.1 & wid > 5e4
               ), by=.(subject.id, fused)]

    # dealing with tiny flanking nodes
    # if a flanking node has 2 or less bins then in.quant.r will be TRUE (unless ratio is NA)
    rel[, in.quant.r := ifelse(.N > 2, in.quant.r, !is.na(ratio)), by = .(subject.id, fused)]

    rel[, lxxx := leix]
    variances = rel[(in.quant.r), var(ratio), keyby=.(fused, lxxx)]

    ## only run if variances has > 0 rows:
    if (!nrow(variances)) {
        return(data.table())
    }
    variances[, side := ifelse(fused, "f_std", "u_std")]
    variances[, std := sqrt(V1)]
    ## if unfused loose ends are present, a column should still be made
    variances[, side := factor(side, levels = c("f_std", "u_std"))]
    vars = dcast.data.table(variances, lxxx ~ side, value.var="std", fill = NA, drop = FALSE)
    rel[is.na(in.quant.r), in.quant.r := FALSE]
    rel[, tum.median := median(tum.counts[in.quant.r]), by=.(lxxx)]
    rel[, norm.median := median(norm.counts[in.quant.r]), by=.(lxxx)]
    rel[, tum.res := tum.counts - tum.median]
    rel[, norm.res := norm.counts - norm.median]
    tum.ks = rel[(in.quant.r), tryCatch(dflm(ks.test(tum.res[fused], tum.res[!fused])), error = function(e) dflm(ks.test(tum.res, tum.res))), by=lxxx][, p, by=lxxx]
    norm.ks = rel[(in.quant.r), tryCatch(dflm(ks.test(norm.res[fused], norm.res[!fused])), error = function(e) dflm(ks.test(norm.res, norm.res))), by=lxxx][, p, by=lxxx]
    pt1 = merge(vars, merge(tum.ks, norm.ks, by="lxxx", suffixes=c("_tum", "_norm"),all=T), by="lxxx", all=T)
    pt1$lxxx = as.character(pt1$lxxx); setkey(pt1, lxxx)
    pt1[, n_fdr := p.adjust(p_norm, "bonferroni")]
    pt1[, t_fdr := p.adjust(p_tum, "bonferroni")]

    rel[, ":="(
        tumor.mean.fused = mean(tum.counts[fused], na.rm=T),
        tumor.mean.unfused = mean(tum.counts[!fused], na.rm=T),
        normal.mean.fused = mean(norm.counts[fused], na.rm=T),
        normal.mean.unfused = mean(norm.counts[!fused], na.rm=T)
    ), by=leix]

    ## evaluate waviness across bins per loose end
    rel[, good.cov := all(good.cov), by=subject.id]
    if(verbose) message("Calculating waviness around loose end")
    rel[, waviness := max(.waviness(start[fused], ratio[fused]), .waviness(start[!fused], ratio[!fused]), na.rm=T), by=subject.id]

    ## prep glm input matrix
    if(verbose) message("Prepping GLM input matrix")
    glm.in = melt.data.table(rel[(in.quant.r),], id.vars=c("leix", "fused"), measure.vars=c("tum.counts", "norm.counts"), value.name="counts")[, tumor := variable=="tum.counts"]
    glm.in[, ix := 1:.N, by=leix]
    rel2 = copy(glm.in)
    setnames(glm.in, "leix", "leix2")

    ## calculate residuals from glm 
    rel2[, residual := .mod(glm.in[leix2==leix[1],]), by=leix]

    ## evaluate KS-test on residuals and calculate effect size
    ## effect will be from KS test on residuals
    ## estimate is replaced with difference of median coverage
    if(verbose) message("Running KS-test on fused vs unfused sides of loose ends")
    res = rel2[(tumor), tryCatch(dflm(ks.test(residual[fused], residual[!fused])), error=function(e) dflm(ks.test(residual, residual))), by=leix]
    est = rel2[, median(counts), keyby=.(fused, tumor, leix)][, V1[tumor] / V1[!tumor], keyby=.(leix, fused)][, V1[fused]-V1[!fused], keyby=leix]
    res$leix = as.character(res$leix); setkey(res, leix)
    res[as.character(est$leix), estimate := est$V1]

    ## combine relevant fields from each test
    test = rel2[(tumor), mean(counts), keyby=.(fused, leix)][, V1[fused]-V1[!fused], keyby=leix]; setnames(test, "V1", "testimate")
    test$leix = as.character(test$leix); setkey(test, leix)
    nest = rel2[!(tumor), mean(counts), keyby=.(fused, leix)][, V1[fused]-V1[!fused], keyby=leix]; setnames(nest, "V1", "nestimate")
    nest$leix = as.character(nest$leix); setkey(nest, leix)
    cnl = rev(rev(colnames(rel))[1:6])
    cns = c("name", "method", "estimate", "effect", "p")
    le.class = cbind(gr2dt(l[rel[!duplicated(subject.id), subject.id]]), rel[!duplicated(subject.id), cnl, with=F], res[rel[!duplicated(subject.id), as.character(leix)], cns, with=F], nest[rel[!duplicated(subject.id), as.character(leix)], "nestimate"], test[rel[!duplicated(subject.id), as.character(leix)], "testimate"])[, effect.thresh := beta]

    le.class[, f.std := pt1[.(as.character(leix)), f_std]]
    le.class[, u.std := pt1[.(as.character(leix)), u_std]]
    le.class[, ":="(n_fdr = pt1[.(as.character(leix)), n_fdr], t_fdr = pt1[.(as.character(leix)), t_fdr])]
    le.class[, bon := p.adjust(p, "bonferroni")]

    ## correct p values
    if(verbose) message("Identifying true positives")
    if(!("epgap" %in% colnames(le.class))){
        le.class[, passed := !is.na(p) & p < PTHRESH & estimate > (0.6*effect.thresh) & testimate > (0.6*effect.thresh) & waviness < 2 & abs(nestimate) < (0.6*effect.thresh)]
    } else {
        le.class[, passed := !is.na(p) & p < PTHRESH & estimate > (0.6*effect.thresh) & testimate > (0.6*effect.thresh) & waviness < 2 & abs(nestimate) < (0.6*effect.thresh) &
                       epgap < max.epgap]
    }
    le.class[, true.pos := passed & estimate > f.std & estimate > u.std & n_fdr > 0.05 & t_fdr < 0.01]
    ## le.class$passed = NULL
    
    gc()
    return(le.class)
}

#' @name dflm
#' @title dflm
#' @description
#' @noRd
#'
#' Formats lm, glm, or fisher.test outputs into readable data.table
#' @author Marcin Imielinski
dflm = function(x, last = FALSE, nm = '')
{
    if (is.null(x))
        out = data.frame(name = nm, method = as.character(NA), p = as.numeric(NA), estimate = as.numeric(NA), ci.lower = as.numeric(NA),  ci.upper = as.numeric(NA), effect = as.character(NA))
    else if (any(c('lm', 'betareg') %in% class(x)))
    {

        coef = as.data.frame(summary(x)$coefficients)
        colnames(coef) = c('estimate', 'se', 'stat', 'p')
        if (last)
            coef = coef[nrow(coef), ]
        coef$ci.lower = coef$estimate - 1.96*coef$se
        coef$ci.upper = coef$estimate + 1.96*coef$se
        if (!is.null(summary(x)$family))
        {
            fam = summary(x)$family$family
            if (summary(x)$family$link %in% c('log', 'logit'))
            {
                coef$estimate = exp(coef$estimate)
                coef$ci.upper= exp(coef$ci.upper)
                coef$ci.lower= exp(coef$ci.lower)
            }
        }
        else
            fam = 'Unknown'

        if (!last)
            nm = paste(nm, rownames(coef))
        
        out = data.frame(name = nm,
                         method = fam,
                         stat = coef$stat,
                         p = signif(coef$p, 3),
                         estimate = coef$estimate,
                         ci.lower = coef$ci.lower,
                         ci.upper = coef$ci.upper,
                         effect = paste(signif(coef$estimate, 3), ' [',
                                        signif(coef$ci.lower,3),'-',
                                        signif(coef$ci.upper, 3), ']',
                                        sep = ''))
    }
    else if (class(x) == 'htest')
    {
        if (is.null(x$estimate))
            x$estimate = x$statistic
        if (is.null(x$conf.int))
            x$conf.int = c(NA, NA)
        out = data.table(name = nm, method = x$method, estimate = x$estimate, ci.lower = x$conf.int[1], ci.upper = x$conf.int[2], effect = paste(signif(x$estimate, 3), ' [',  signif(x$conf.int[1],3),'-', signif(x$conf.int[2], 3), ']', sep = ''), p = x$p.value)
    }
    else if (class(x) == 'polr')
    {
        coef = coef(summary(x)) %>% as.data.frame
        nm = paste(nm, rownames(coef))
        coef = as.data.table(coef)
        setnames(coef, c('estimate', 'se', 't'))
        out = data.table(name = nm) %>% cbind(coef)
        out$p =  pnorm(abs(out$t), lower.tail = FALSE) * 2
        out[, ci.lower := estimate-1.96*se]
        out[, ci.upper := estimate+1.96*se]
        out[, effect := paste(signif(estimate, 3), ' [',  signif(ci.lower,3),'-', signif(ci.upper, 3), ']', sep = '')]
    }
    else
    {
        out = data.frame(name = nm, method = x$method, p = signif(x$p.value, 3), estimate = x$estimate, ci.lower = x$conf.int[1], ci.upper = x$conf.int[2], effect = paste(signif(x$estimate, 3), ' [',  signif(x$conf.int[1],3),'-', signif(x$conf.int[2], 3), ']', sep = ''))
    }

    out$effect = as.character(out$effect)
    out$name = as.character(out$name)
    out$method = as.character(out$method)
    rownames(out) = NULL
    return(as.data.table(out))
}


#' @name .waviness
#'
#' quantifies autocorrelation in coverage
.waviness = function(x, y, min.thresh = 5e3, max.thresh = 10e4, spar = 0.5, smooth = TRUE, filter = rep(FALSE, length(x)), trim = 10) {
    if(length(x)==0) return(NA)
    dat = data.table(x, y)[order(x), ]
    dat[, lag := x-min(x)]
    fdat = dat[!is.na(y), ][!is.infinite(y), ]
    ## autocorrelation
    if(nrow(fdat)==0) return(NA)
    fdat[, ac := as.numeric(acf(c(y, y), plot = FALSE, lag.max = length(y))$acf[-1])]
    if (smooth) { ## smoothing the autocorrelation gets rid of some more noise
        fdat = fdat[!is.na(lag) & !is.na(ac),]
        if (nrow(fdat[!is.na(lag) & !is.na(ac),])<4)
            return(NA)
        fdat$ac = predict(smooth.spline(fdat$lag, fdat$ac, spar = spar), fdat$lag)$y
    }
    return(fdat[lag>min.thresh & lag<max.thresh, sum(ac^2)])
}

#' @name .mod
#'
#' fit data table with glm and return residuals
.mod = function(dt){
    mod = dt[, glm(counts ~ tumor + fused + ix, family='gaussian')]
    res = dt$counts - predict(mod, dt, type='response')
    return(res)
}

#' @name .mod2
#'
#' fit data table with glm and return residuals
.mod2 = function(dt){
    mod = dt[, glm(counts ~ tumor + ix, family='gaussian')]
    res = dt$counts - predict(mod, dt, type='response')
    return(res)
}
