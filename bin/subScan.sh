#!/bin/bash

temperatures=("10" "15" "20" "25" "30" "35" "40" "45" "50" "60" "70" "80" \
    "90" "100" "110" "120" "130" "140" "150" "175" "200" "225" "250" "275" \
    "300" "350" "400" "450" "500" "750")

temperatures=("10"  "150" "500")
temperatures=("10")

export EXEDIR=/home/mu495/Software/Prop
export OUTDIR=/archive/mu495/Fit
export DATADIR=$EXEDIR/Data
export EVO=SFR2
export IRB=G12
export PRODNAME=Test_${EVO}_${IRB}

queue=s48

nTemperature=${#temperatures[@]}
echo number of temperatures $nTemperature
((nTemperature--))

i=0
for temperature in "${temperatures[@]}"
do
  for j in $(seq $i $nTemperature)
  do
    export T1=$temperature
    export T2=${temperatures[$j]}
    echo $T1 $T2
    export FILEBASE=${PRODNAME}_${T1}_${T2}
    logName=$OUTDIR/$FILEBASE
    echo qsub -q $queue -o $logName.out -e $logName.err $EXEDIR/bin/runFit.sh
    qsub -V -q $queue -o $logName.out -e $logName.err $EXEDIR/bin/runFit.sh
  done
  ((i++))
done
