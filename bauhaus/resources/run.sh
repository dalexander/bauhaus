#!/bin/bash
source /mnt/software/Modules/current/init/bash

module purge
module use /pbi/dept/primary/modulefiles
module load smrtanalysis/mainline
module load gfftools/dalexander
module load R/3.2.3-bauhaus
module load ninja

THISDIR=$(cd "$(dirname "$0")" && pwd)
cd $THISDIR
export PATH=$THISDIR/scripts:$PATH
export NINJA_STATUS="[%f/%t] [Elapsed: %e] "

ninja -j 999 -v -k 1 | tee ninja.log
