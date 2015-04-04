#!/bin/tcsh

cd /afs/cern.ch/user/m/munger/Prop
source bin/setvars.csh
cd -
/afs/cern.ch/user/m/munger/Prop/bin/createPropMatrixFile $EVOLUTION $PIONFILE $BASEDIR/$PHOTONFIELD/*.root
mv propMatrix_${EVOLUTION}.root $OUTDIR/${PHOTONFIELD}_${EVOLUTION}.root
mv propMatrix_${EVOLUTION}_nu.root $OUTDIR/${PHOTONFIELD}_${EVOLUTION}_nu.root
