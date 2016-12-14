#!/bin/bash
set -euo pipefail
set +x

# what are we?
echo "$(uname -a)"

. /mnt/software/Modules/current/init/bash
# the order of R/3.2.3 module with gcc/4.8.4 matters
module load R/3.2.3-bauhaus
module load gcc/4.8.2
module load git/2.8.3

# deployment

export R_LIBS=/pbi/dept/consensus/bauhaus/rlibs
rm -rf $R_LIBS
mkdir -p $R_LIBS

echo "$(env | grep R_LIBS)"
echo "$(R -q -e 'print(.libPaths())')"

R -q -e 'install.packages(c("Rcpp", "RcppEigen"), repos="http://cran.us.r-project.org")'
R -q -e 'install.packages(c("argparse", "argparser", "data.table", "devtools", "digest", "doMC", "doParallel", "dplyr", "feather", "ggplot2", "", "gridExtra", "gtools", "logging", "plyr", "tidyr", "wesanderson", "xml2"), repos="http://cran.us.r-project.org")'
R -q -e 'library(devtools) ; install_github("PacificBiosciences/pbbamr")'
R -q -e 'library(devtools) ; install_git("http://bitbucket.nanofluidics.com:7990/scm/sat/unitem-sat.git")'
