#' JaBbA: Junction Balance Analysis
#' 
#' @importFrom gUtils streduce si2gr seg2gr rrbind ra.overlaps ra.duplicated parse.gr hg_seqlengths grl.unlist grl.pivot grl.in grl.eval grl.bind grbind gr2dt gr.val gr.tile.map gr.tile
#' @importFrom gUtils gr.stripstrand gr.sum gr.string gr.start gr.end gr.simplify gr.setdiff gr.sample gr.reduce gr.rand gr.quantile gr.nochr
#' @importFrom gUtils gr.match gr.in gr.flipstrand gr.fix gr.findoverlaps gr.duplicated gr.dist gr.disjoin gr.breaks dt2gr "%^%" "%Q%" "%&%" "%$%"
#' @importFrom S4Vectors mcols mcols<- values values<- elementMetadata elementMetadata<-
#' @importFrom BiocGenerics width
#' @importFrom gGnome gG balance
#' @importFrom GenomicRanges GRanges GRangesList values split match setdiff
#' @importFrom gTrack gTrack
#' @importFrom igraph graph induced.subgraph V E graph.adjacency clusters
#' @importFrom optparse make_option OptionParser parse_args print_help
#' @importFrom data.table data.table as.data.table setnames setkeyv fread setkey between
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
"_PACKAGE"