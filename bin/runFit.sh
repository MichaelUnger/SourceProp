#!/bin/bash
#PBS -l nodes=1:ppn=1,walltime=24:00:00

function calc() {
    awk "BEGIN { print "$*" }"
}


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

#jobdir=$PBS_JOBTMP
jobdir=/state/partition1/$PBS_JOBID
mkdir $jobdir
if [ $? -ne 0 ] ; then
   echo cannot create $jobdir
   exit 1
fi
cd $jobdir

masses=("56"  "40" "28" "14" "4")
#masses=("56")
for mass in "${masses[@]}"
do
    echo "########################## Astart = $mass #################"
    echo executing $PBS_JOBID at $PWD on $HOSTNAME
    echo "# $FILEBASE $PBS_JOBID at $PWD on $HOSTNAME" > common.txt
    echo "OutDir $PWD" >> common.txt
    echo "DataDir $DATADIR" >> common.txt
    echo "evolution $EVO" >> common.txt
    echo "IRB $IRB" >> common.txt
#    echo "par gammaInj   -1.1 0.1 -3 -0.8 0" >> common.txt
#    echo "par deltaEsc   -0.9 0.1 -1.01 -0.2 0" >> common.txt
    echo "mass $mass 0.1 1 56 0 1" >> common.txt

    if [ "$T1" -ne "$T2" ]
    then
        echo "par lgfPhoton    0 0.1 0 0 0" >> common.txt
    fi


    export FITFILE=$PWD/Fit.txt

    sigmas=("0"  "1" "2")
    for sigma in "${sigmas[@]}"
    do
        echo "########################### MBB$sigma #######################"
        cat common.txt > $FITFILE
        echo "PhotonMBB $T1 $sigma" >> $FITFILE
        suffix="MBB1_${T1}_${sigma}_A${mass}"
        if [ "$T1" -ne "$T2" ]
        then
            echo "PhotonMBB $T2 $sigma" >> $FITFILE
            suffix="MBB2_${T1}_${T2}_${sigma}_A${mass}"
        fi
        $rootCmd > ${FILEBASE}_${suffix}.log 2>&1
        mv Fit.root ${FILEBASE}_${suffix}.root
        mv Fit.pdf ${FILEBASE}_${suffix}.pdf
        mv Fit_nu.pdf ${FILEBASE}_${suffix}_nu.pdf
        mv $FITFILE ${FILEBASE}_${suffix}.txt
    done


    echo "########################### BPL #######################"
    cat common.txt > $FITFILE
    eps0_1=`calc $T1/1000`
    echo "PhotonBPL $eps0_1 52 2.0" >> $FITFILE
    suffix="BPL1_${T1}_A${mass}"
    if [ "$T1" -ne "$T2" ]
    then
        eps0_2=`calc $T2/1000`
        echo "PhotonBPL $eps0_2 52 2.0" >> $FITFILE
        suffix="BPL2_${T1}_${T2}_A${mass}"
    fi
    $rootCmd > ${FILEBASE}_${suffix}.log 2>&1
    mv Fit.root ${FILEBASE}_${suffix}.root
    mv Fit.pdf ${FILEBASE}_${suffix}.pdf
    mv Fit_nu.pdf ${FILEBASE}_${suffix}_nu.pdf
    mv $FITFILE ${FILEBASE}_${suffix}.txt
    rm common.txt
done
gzip *.log
gzip *.txt
tarballName=${FILEBASE}_${T1}_${T2}.tar
tar -cvf $tarballName *.root *.pdf *.gz
mv $tarballName ${OUTDIR}/
rm *.root *.pdf *.gz
rm -f core
ls -la
rmdir $jobdir
