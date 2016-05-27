#!/bin/bash

Rscript -e 'source("https://bioconductor.org/biocLite.R"); 
biocLite("BiocInstaller"); 
biocLite("BSgenome.Hsapiens.UCSC.hg19"); 
biocLite("GenomicRanges"); 
install.packages("devtools"); 
devtools:: install_github("mskilab/gUtils"); 
devtools::install_github("jimhester/covr"); 
devtools::install_github("mskilab/gTrack"); 

library(gTrack);'
