#!/bin/tcsh 

set options = "uniform uniformCutAt3 AGN SFR2 AAGHRW05 M10 M15 M20 M25 M30 M35 M40 M45 M50"

#set options = "uniform"

setenv BASEDIR  /afs/cern.ch/work/u/unger/crp
setenv PIONFILE $BASEDIR/pionDecay.root
setenv PHOTONFIELD Test

setenv OUTDIR /afs/cern.ch/work/m/munger/crp

set queue=1nh

set counter = 0
foreach option ($options)
  setenv EVOLUTION $option
  set jName = mc$counter
  set logName = $OUTDIR/logs/${PHOTONFIELD}_${option}
  bsub -q $queue -oo $logName.lsf.out -eo $logName.lsf.err -J $jName $PWD/bin/runProp.csh
  @ counter ++
end
