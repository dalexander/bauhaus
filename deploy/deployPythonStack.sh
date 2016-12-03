#!/bin/bash

# Module system
. /mnt/software/Modules/current/init/bash


module purge
module load gcc/4.8.4 python/2.7.9 zlib/1.2.5 openblas/0.2.14

set -euo pipefail
set -x

BAUHAUS_VE=${BAUHAUS_VE:-/mnt/secondary/Share/MLG/Bauhaus/VE-Bauhaus}
rm -rf $BAUHAUS_VE

set +u
/mnt/software/v/virtualenv/13.0.1/virtualenv.py $BAUHAUS_VE
source $BAUHAUS_VE/bin/activate
set -u

pip install -r requirements.txt
python setup.py install
