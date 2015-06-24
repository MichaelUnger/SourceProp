#!/bin/bash

declare -a options=("SFR2" "uniform" "uniformCutAt3" "AGN" "AAGHRW05" "M10" "M15" "M20" "M25" "M30" "M35" "M40" "M45" "M50")

option="uniform"
export EXEDIR=/home/mu495/Software/Prop
export PIONFILE=$EXEDIR/bin/pionDecay.root
export PHOTONFIELD=Test #CRPropaG12
export OUTDIR=/scratch/mu495/Matrices
export INDIR=/scratch/mu495/

queue=s48

for option in "${options[@]}"
do
  export EVOLUTION=$option
  logName=$OUTDIR/logs/${PHOTONFIELD}_${option}
  echo qsub -q $queue -o $logName.out -e $logName.err $EXEDIR/bin/runProp.sh
#  qsub -V -q $queue -o $logName.out -e $logName.err $EXEDIR/bin/runProp.sh
done
