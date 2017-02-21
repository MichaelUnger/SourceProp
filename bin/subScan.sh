#!/bin/bash

#function calc() {
#    awk "BEGIN { print "$*" }"
#}
#eps0_1=`calc $T1/1000`

temperatures=("10" "15" "20" "25" "30" "35" "40" "45" "50" "60" "70" "80" \
    "90" "100" "110" "120" "130" "140" "150" "175" "200" "225" "250" "275" \
    "300" "350" "400" "450" "500" "750")

temperatures=("10"  "150" "500")
temperatures=("10")

nTemperature=${#temperatures[@]}
echo number of temperatures $nTemperature
((nTemperature--))

export EXEDIR=/home/mu495/Software/Prop
export OUTDIR=/scratch/mu495/Matrices
export DATADIR=/scratch/mu495/

queue=s48

i=0
for temperature in "${temperatures[@]}"
do
  for j in $(seq $i $nTemperature)
  do
    export T1=$temperature
    export T2=${temperatures[$j]}
    echo $T1 $T2
    logName=$OUTDIR/logs/${PHOTONFIELD}_${option}
    echo qsub -q $queue -o $logName.out -e $logName.err $EXEDIR/bin/runProp.sh
    qsub -V -q $queue -o $logName.out -e $logName.err $EXEDIR/bin/runProp.sh
  done
  ((i++))
done
