#!/bin/bash
#PBS -l nodes=1:ppn=1,walltime=1:00:00

function calc() {
    awk "BEGIN { print "$*" }"
}
#eps0_1=`calc $T1/1000`


if [ "$HOSTNAME" != "bowman" ]
then
    module load NYUAD/2.0
    module load devel
    module load gcc/4.9.1
    module load python/2.7.9
    module load openmpi/1.8.3
    module load gsl/1.16
    . $EXEDIR/bin/setvars_butinah.sh
    . /home/mu495/Software/ROOT/root_v5.34.24_install/bin/thisroot.sh
else
    export EXEDIR=/home/munger/Mag/Prop
    export DATADIR=/home/munger/Mag/Prop/data
    echo running on laptop
    cd /tmp/run
fi

rootCmd="root.exe -b -l -q -x $EXEDIR/macros/fitWrapper.C"

jobdir=/scratch/mu495/tmp/$PBS_JOBID
mkdir $jobdir
cd $jobdir

echo executing $PBS_JOBID at $PWD on $HOSTNAME
echo "# $FILEBASE $PBS_JOBID at $PWD on $HOSTNAME" > common.txt
echo "OutDir $PWD" >> common.txt
echo "DataDir $DATADIR" >> common.txt
echo "evolution $EVO" >> common.txt
echo "IRB $IRB" >> common.txt
if [ "$T1" -ne "T2" ]
then
    echo "par lgfPhoton    0 0.1 0 0 0" >> common.txt
fi


export FITFILE=$PWD/Fit.txt

cat common.txt > $FITFILE

$rootCmd
ls
mv Fit.root $OUTDIR/${FILEBASE}.root
mv Fit.pdf $OUTDIR/${FILEBASE}.pdf
mv Fit_nu.pdf $OUTDIR/${FILEBASE}_nu.pdf
mv $FITFILE $OUTDIR/${FILEBASE}.txt
rm common.txt
rmdir $jobdir
