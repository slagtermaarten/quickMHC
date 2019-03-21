#!/bin/zsh

olddir=$(pwd)
cd ~/libs/quickMHC
Rscript -e 'Rcpp::compileAttributes()'
R CMD build .
R CMD INSTALL ./quickMHC_0.1.0.0.tar.gz
cd $olddir
