#!/bin/bash
source /mnt/software/Modules/current/init/bash

# The latter line here is a directive to the R module load: don't let
# user's R environment interfere with our scripts.  Developers may
# want to comment out these lines.
unset R_LIBS
export R_IGNORE_USER_LIBS=1

module purge
module use /pbi/dept/primary/modulefiles
module load smrtanalysis/mainline
module load gfftools/dalexander
module load R/3.2.3-internal
module load ninja

THISDIR=$(cd "$(dirname "$0")" && pwd)
cd $THISDIR
export PATH=$THISDIR/scripts:$PATH
export NINJA_STATUS="[%f/%t] [Elapsed: %e] "

ninja -j 999 -v -k 1 | tee ninja.log
