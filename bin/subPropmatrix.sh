#!/bin/tcsh

#options = "SFR2 uniform uniformCutAt3 AGN AAGHRW05 M10 M15 M20 M25 M30 M35 M40 M45 M50"

option="uniform"
export EXEDIR=/home/mu495/Software/Prop
export PIONFILE=$EXEDIR/bin/pionDecay.root
export PHOTONFIELD=Test #CRPropaG12
export OUTDIR=/scratch/mu495/Matrices
export INDIR=/scratch/mu495/

set queue=s48

set counter = 0
#foreach option ($options)
  export EVOLUTION=$option
  jName=mc$counter
  logName=$OUTDIR/logs/${PHOTONFIELD}_${option}
  qsub -q $queue -o $logName.out -e $logName.err $EXEDIR/bin/runProp.csh
# -J $jName
  @ counter ++
#end
