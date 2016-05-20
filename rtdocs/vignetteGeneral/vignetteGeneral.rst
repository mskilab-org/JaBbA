General Vignette (Detailing some functions)
===========================================

segstats
~~~~~~~
Calculate the **mean** and **standard deviation of GRanges object when supplied the tumor depth and normal depth** 


.. sourcecode:: r
    

    ## upload coverage data for 10X HCC1143 cell line
    ocovh = readRDS("~/Desktop/Projects/10X/files/HCC1143/cov/cov.rds")
    ## tile the segments based on size of chromosomes
    cov = gr.tile(seqlengths(ocovh), 5e3)
    
    ## calculate mean and standard deviation
    cov = segstats(cov, ocovh, field  = 'ratio')


::

    ## Iteration 100000 of 1285699 
    ## Iteration 200000 of 1285699 
    ## Iteration 300000 of 1285699 
    ## Iteration 400000 of 1285699 
    ## Iteration 500000 of 1285699 
    ## Iteration 600000 of 1285699 
    ## Iteration 700000 of 1285699 
    ## Iteration 800000 of 1285699 
    ## Iteration 900000 of 1285699 
    ## Iteration 1000000 of 1285699 
    ## Iteration 1100000 of 1285699 
    ## Iteration 1200000 of 1285699 
    ## Iteration 100000 of 703949 
    ## Iteration 200000 of 703949 
    ## Iteration 300000 of 703949 
    ## Iteration 400000 of 703949 
    ## Iteration 500000 of 703949 
    ## Iteration 600000 of 703949 
    ## Iteration 700000 of 703949 
    ## Iteration 100000 of 699153 
    ## Iteration 200000 of 699153 
    ## Iteration 300000 of 699153 
    ## Iteration 400000 of 699153 
    ## Iteration 500000 of 699153 
    ## Iteration 600000 of 699153 
    ## Iteration 100000 of 688218 
    ## Iteration 200000 of 688218 
    ## Iteration 300000 of 688218 
    ## Iteration 400000 of 688218 
    ## Iteration 500000 of 688218 
    ## Iteration 600000 of 688218 
    ## Iteration 100000 of 593544 
    ## Iteration 200000 of 593544 
    ## Iteration 300000 of 593544 
    ## Iteration 400000 of 593544 
    ## Iteration 500000 of 593544 
    ## Iteration 100000 of 553117 
    ## Iteration 200000 of 553117 
    ## Iteration 300000 of 553117 
    ## Iteration 400000 of 553117 
    ## Iteration 500000 of 553117 
    ## Iteration 100000 of 521763 
    ## Iteration 200000 of 521763 
    ## Iteration 300000 of 521763 
    ## Iteration 400000 of 521763 
    ## Iteration 500000 of 521763 
    ## Iteration 100000 of 461903 
    ## Iteration 200000 of 461903 
    ## Iteration 300000 of 461903 
    ## Iteration 400000 of 461903 
    ## Iteration 100000 of 409629 
    ## Iteration 200000 of 409629 
    ## Iteration 300000 of 409629 
    ## Iteration 400000 of 409629 
    ## Iteration 100000 of 395810 
    ## Iteration 200000 of 395810 
    ## Iteration 300000 of 395810 
    ## Iteration 100000 of 307470 
    ## Iteration 200000 of 307470 
    ## Iteration 300000 of 307470 
    ## Iteration 100000 of 1263347 
    ## Iteration 200000 of 1263347 
    ## Iteration 300000 of 1263347 
    ## Iteration 400000 of 1263347 
    ## Iteration 500000 of 1263347 
    ## Iteration 600000 of 1263347 
    ## Iteration 700000 of 1263347 
    ## Iteration 800000 of 1263347 
    ## Iteration 900000 of 1263347 
    ## Iteration 1000000 of 1263347 
    ## Iteration 1100000 of 1263347 
    ## Iteration 1200000 of 1263347 
    ## Iteration 100000 of 324668 
    ## Iteration 200000 of 324668 
    ## Iteration 300000 of 324668 
    ## Iteration 100000 of 244111 
    ## Iteration 200000 of 244111 
    ## Iteration 100000 of 258397 
    ## Iteration 200000 of 258397 
    ## Iteration 100000 of 1029717 
    ## Iteration 200000 of 1029717 
    ## Iteration 300000 of 1029717 
    ## Iteration 400000 of 1029717 
    ## Iteration 500000 of 1029717 
    ## Iteration 600000 of 1029717 
    ## Iteration 700000 of 1029717 
    ## Iteration 800000 of 1029717 
    ## Iteration 900000 of 1029717 
    ## Iteration 1000000 of 1029717 
    ## Iteration 100000 of 994002 
    ## Iteration 200000 of 994002 
    ## Iteration 300000 of 994002 
    ## Iteration 400000 of 994002 
    ## Iteration 500000 of 994002 
    ## Iteration 600000 of 994002 
    ## Iteration 700000 of 994002 
    ## Iteration 800000 of 994002 
    ## Iteration 900000 of 994002 
    ## Iteration 100000 of 940462 
    ## Iteration 200000 of 940462 
    ## Iteration 300000 of 940462 
    ## Iteration 400000 of 940462 
    ## Iteration 500000 of 940462 
    ## Iteration 600000 of 940462 
    ## Iteration 700000 of 940462 
    ## Iteration 800000 of 940462 
    ## Iteration 900000 of 940462 
    ## Iteration 100000 of 888680 
    ## Iteration 200000 of 888680 
    ## Iteration 300000 of 888680 
    ## Iteration 400000 of 888680 
    ## Iteration 500000 of 888680 
    ## Iteration 600000 of 888680 
    ## Iteration 700000 of 888680 
    ## Iteration 800000 of 888680 
    ## Iteration 100000 of 825873 
    ## Iteration 200000 of 825873 
    ## Iteration 300000 of 825873 
    ## Iteration 400000 of 825873 
    ## Iteration 500000 of 825873 
    ## Iteration 600000 of 825873 
    ## Iteration 700000 of 825873 
    ## Iteration 800000 of 825873 
    ## Iteration 100000 of 760630 
    ## Iteration 200000 of 760630 
    ## Iteration 300000 of 760630 
    ## Iteration 400000 of 760630 
    ## Iteration 500000 of 760630 
    ## Iteration 600000 of 760630 
    ## Iteration 700000 of 760630 
    ## Iteration 100000 of 729422 
    ## Iteration 200000 of 729422 
    ## Iteration 300000 of 729422 
    ## Iteration 400000 of 729422 
    ## Iteration 500000 of 729422 
    ## Iteration 600000 of 729422 
    ## Iteration 700000 of 729422 
    ## Iteration 100000 of 805552 
    ## Iteration 200000 of 805552 
    ## Iteration 300000 of 805552 
    ## Iteration 400000 of 805552 
    ## Iteration 500000 of 805552 
    ## Iteration 600000 of 805552 
    ## Iteration 700000 of 805552 
    ## Iteration 800000 of 805552



::

    ## Error in gr.tile.map(utarget, signal, verbose = T, mc.cores = mc.cores): could not find function "munlist"




