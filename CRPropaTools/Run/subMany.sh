#!/bin/sh

if [ $# != 2 ]; then
   echo usage: $0 firstJob lastJob
   exit 1
fi

for i in $(seq $1 $2); do
   export JOBNUM=$1
   qsub -q s48 run.sh
done


