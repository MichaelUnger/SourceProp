#!/bin/bash

temperatures=("10" "15" "20" "25" "30" "35" "40" "45" "50" "60" "70" "80" \
    "90" "100" "110" "120" "130" "140" "150" "175" "200" "225" "250" "275" \
    "300" "350" "400" "450" "500" "750" "1000" "2000" "3000" "4000"  "5000" \
    "6000" "7000" "8000" "9000")

#temperatures=("250")

export EXEDIR=/home/mu495/Software/Prop
export OUTDIR=/scratch/mu495/Fit
export DATADIR=$EXEDIR/Data
export EVO=SFR2
export IRB=G12
export PRODNAME=FixMass_${EVO}_${IRB}

singleTemperature=0
queue=s48

nTemperature=${#temperatures[@]}
echo number of temperatures $nTemperature
((nTemperature--))

i=0
for temperature in "${temperatures[@]}"
do
  for j in $(seq $i $nTemperature)
  do
    export T1=${temperatures[$j]}
    export T2=$temperature
    echo $T1 $T2
    export FILEBASE=${PRODNAME}
    logName=/home/mu495/Logfiles/${FILEBASE}_${T1}_${T2}
    echo qsub -q $queue -o $logName.out -e $logName.err $EXEDIR/bin/runFit.sh
    qsub -V -q $queue -o $logName.out -e $logName.err $EXEDIR/bin/runFit.sh
    if [ $singleTemperature -ne 0 ]
    then
        break
    fi
  done
  ((i++))
done
