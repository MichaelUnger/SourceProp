#!/bin/tcsh 

set options = "eMm40z10 eMm40z20  eMm40z30  eMm40z40 eMm40z50 eMm20z10 eMm20z20 eMm20z30 eMm20z40 eMm20z50 eM00z10 eM00z20 eM00z30 eM00z40 eM00z50 eMp20z10 eMp20z20 eMp20z30 eMp20z40 eMp20z50 eMp40z10 eMp40z20 eMp40z30 eMp40z40  eMp40z50"


setenv BASEDIR  /afs/cern.ch/work/u/unger/crp
setenv PIONFILE $BASEDIR/pionDecay.root
setenv PHOTONFIELD CRPropaG12

setenv OUTDIR /afs/cern.ch/work/m/munger/crp

set queue=8nh

set counter = 0
foreach option ($options)
  setenv EVOLUTION $option
  set jName = mc$counter
  set logName = $OUTDIR/logs/${PHOTONFIELD}_${option}
  bsub -q $queue -oo $logName.lsf.out -eo $logName.lsf.err -J $jName $PWD/bin/runProp.csh
  @ counter ++
end
