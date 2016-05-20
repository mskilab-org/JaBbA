General Vignette (Detailing some functions)
===========================================

segstats
~~~~~~~
Calculate the **mean** and **standard deviation of GRanges object when supplied the tumor depth and normal depth** 


.. sourcecode:: r
    

    ## upload coverage data for 10X HCC1143 cell line
    ocovh = readRDS("~/Desktop/Projects/10X/files/HCC1143/cov/cov.rds")
    ## tile the segments based on size of chromosomes
    cov = gr.tile(seqlengths(ocovh, 5e3)
    
    ## calculate mean and standard deviation
    cov = segstats(cov, ocovh, field  = 'ratio')


::

    ## Error: <text>:7:1: unexpected symbol
    ## 6: ## calculate mean and standard deviation
    ## 7: cov
    ##    ^




