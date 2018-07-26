#!/bin/tcsh

source /afs/cern.ch/na61/Releases/SHINE/pro64/scripts/env/lxplus_64bit_slc6.csh
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/afs/cern.ch/user/m/munger/Prop/lib
set PROG=/afs/cern.ch/user/m/munger/Prop/bin/createPropMatrixFile

rsync --partial --verbose --progress --stats --recursive ${BASEDIR}/ .

set minDists = "0 1 5 10 50 100"

foreach minDist ($minDists)

  $PROG $EVOLUTION $minDist $PIONFILE1 $PIONFILE2 $PHOTONFIELD/*.root
  mv propMatrix_${EVOLUTION}_nu.root $OUTDIR/${PHOTONFIELD}_${EVOLUTION}_${minDist}_nu.root
  
end
