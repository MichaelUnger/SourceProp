#!/bin/tcsh 

set options = "SFR2 uniform uniformCutAt3 AGN AAGHRW05 M10 M15 M20 M25 M30 M35 M40 M45 M50"

set options = "uniform"

setenv EXEDIR /home/mu495/Software/Prop
setenv PIONFILE $EXEDIR/bin/pionDecay.root
setenv PHOTONFIELD Test #CRPropaG12
setenv OUTDIR /scratch/mu495/Matrices
setenv INDIR /scratch/mu495/

set queue=s48

set counter = 0
foreach option ($options)
  setenv EVOLUTION $option
  set jName = mc$counter
  set logName = $OUTDIR/logs/${PHOTONFIELD}_${option}
  qsub -q $queue $EXEDIR/bin/runProp.csh
#-oo $logName.lsf.out -eo $logName.lsf.err -J $jName 
  @ counter ++
end
