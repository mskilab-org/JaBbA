
library(jabbadevtest)   ### package
library(gUtils)
library(testthat)


junctions = system.file("extdata", "junctions.vcf", package = 'jabbadevtest')
coverage = system.file("extdata", "coverage.txt", package = 'jabbadevtest')
hets = system.file("extdata", "hets.txt", package = 'jabbadevtest')



library(Rcplex)

cvec <- c(1,2,3)
Amat <- matrix(c(-1,1,1,-1,3,-1),byrow=TRUE,nc=3)
bvec <- c(20,-30)
ub <- c(40,Inf,Inf)
res <- Rcplex(cvec,Amat,bvec,ub=ub,objsense="max",sense=c('L','G'))
print(res)


test_that('testing Rcplex works', {
    
    expect_equal(res$xopt[1], 40)
    expect_equal(res$xopt[3], 42.5)
    expect_equal(res$obj, 202.5)
    expect_equal(res$extra$slack[1], 0)

})


##set.seed(42);

##jab = JaBbA(junctions = junctions, coverage = coverage, tilim = 10, verbose = 1, overwrite = TRUE)

##print(jab)

## JaBbA

## jabba_stub

## karyograph_stub

## .plot_ppfit

## ramip_stub

## segstats 

## jmessage

## jbaMIP

## JaBbA.digest

## jbaMIP.process

## jabba.alleles

## pp.nll

## munlist
test_that('testing munlist() works', {
    
    expect_equal(dim(munlist(gr2dt(example_genes)$start))[1], 18812)
    expect_equal(dim(munlist(gr2dt(example_genes)$start))[2], 3)

})

## .correct.slack

## .cplex_customparams
test_that('testing alpha() works', {
    
    ## not sure how to test this properly
    expect_error(.cplex_customparams('out.txt'), NA) ## check no error

})

## gr.tile.map 
## test_that('testing gr.tile.map() works', {
## 
##     > gr.tile.map(example_genes, example_dnase)
##     Error in m[, 2] : subscript out of bounds
## 
## })


## vaggregate  ### this should be replaced! data.table


## write.tab
test_that('testing write.tab() works', {

    expect_error(write.tab(example_genes), NA) ## check works

})



## alpha()
test_that('testing alpha() works', {

    expect_match(alpha('blue', 1), '#0000FFFF')

})




## rel2abs()


## abs2rel()

## adj2inc()



## mmatch()


## all.paths()



## collapse.paths()


## sparse_subset()
test_that('testing sparse_subset() works', {

    ## function(A, B, strict = FALSE, chunksize = 100, quiet = FALSE)
    Amat = matrix(c(20, 40, 30, 10), nrow=2, ncol=2, byrow = TRUE) 
    Bmat = matrix(c(250, 450, 350, 150), nrow=2, ncol=2, byrow = TRUE)
    expect_equal(dim(sparse_subset(Amat, Bmat, quiet=TRUE))[1], 2)
    expect_equal(dim(sparse_subset(Amat, Bmat, quiet=TRUE))[2], 2)
    expect_equal(dim(sparse_subset(Amat, Bmat, chunksize=5, strict=TRUE, quiet=TRUE))[1], 2)
    expect_equal(dim(sparse_subset(Amat, Bmat, chunksize=5, strict=TRUE, quiet=TRUE))[2], 2)


})




## convex.basis()
test_that('testing convex.basis() works', {

    example_matrix = matrix(c(2, 4, 3, 1, 5, 7), nrow=2, ncol=3, byrow = TRUE)  
    expect_equal(as.logical(convex.basis(example_matrix)), NA)

})



## read.junctions()

## karyograph()



## jabba2vcf()

## chromoplexy()

## read_vcf()

## levapply()

test_that('testing levapply() works', {

    gr = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
    foo = levapply(width(gr), as.character(seqnames(gr)), function(x) if (length(x)>1) cumsum(c(0, x[1:(length(x)-1)])) else return(0))
    expect_equal(foo, c(0, 3, 6))

})


## chr2num

test_that('testing chr2num() works', {

    expect_equal(chr2num('chr1'), 1)
    expect_equal(chr2num('chrX'), 23)
    expect_equal(chr2num('chrY'), 24)
    expect_match(chr2num('chr1', xy=TRUE), '1')
    expect_match(chr2num('chrX', xy=TRUE), 'X')
    expect_match(chr2num('chrY', xy=TRUE), 'Y')
    expect_equal(as.logical(chr2num('chrZ')), NA)
    ### hmmmmm
    expect_match(chr2num('chrZ', xy=TRUE), 'Z')

})
