#!/bin/tcsh

source /afs/cern.ch/na61/Releases/SHINE/pro/scripts/env/lxplus_64bit_slc6.csh
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/afs/cern.ch/user/m/munger/Prop/lib

/afs/cern.ch/user/m/munger/Prop/bin/createPropMatrixFile $EVOLUTION $PIONFILE $BASEDIR/$PHOTONFIELD/*.root
mv propMatrix_${EVOLUTION}.root $OUTDIR/${PHOTONFIELD}_${EVOLUTION}.root
mv propMatrix_${EVOLUTION}_nu.root $OUTDIR/${PHOTONFIELD}_${EVOLUTION}_nu.root
