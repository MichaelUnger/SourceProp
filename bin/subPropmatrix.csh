#!/bin/tcsh

set options = "eUniform eUniformCutAt3 eAGN eSFR2 eAAGHRW05 eM10 eM15 eM20 eM25 eM30 eM35 eM40 eM45 eM50"

set options = "eUniform"

setenv BASEDIR  /afs/cern.ch/work/u/unger/crp
setenv PIONFILE $baseDir/pionDecay.root
setenv PHOTONFIELD Test

setenv OUTDIR /afs/cern.ch/work/m/munger/crp

set queue=8nh

set counter = 0
foreach option ($options)
  setenv EVOLUTION $option
  set jName mc$counter
  set logName $OUTDIR/logs/${PHOTONFIELD}_${option}
  bsub -q $queue -oo $logName.lsf.out -eo $logName.lsf.err -J $jName runProp.csh
  @ counter ++
end
