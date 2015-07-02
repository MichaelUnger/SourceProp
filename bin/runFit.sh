#!/bin/bash
#PBS -l nodes=1:ppn=1,walltime=1:00:00

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

echo executing at $PWD on $HOSTNAME

rootCmd="root.exe -b -l -q -x $EXEDIR/macros/fitWrapper.C"

export FITFILE=$PWD/Fit.txt
echo "# executing at $PWD on $HOSTNAME" > $FITFILE
echo "OutDir $PWD" >> $FITFILE
echo "DataDir $DATADIR" >> $FITFILE
echo "evolution $EVO" >> $FITFILE
echo "IRB $IRB" >> $FITFILE

#cd $EXEDIR
$rootCmd
ls
mv Fit.root $OUTDIR/${FILEBASE}.root
mv Fit.pdf $OUTDIR/${FILEBASE}.pdf
mv Fit_nu.pdf $OUTDIR/${FILEBASE}_nu.pdf
mv $FITFILE $OUTDIR/${FILEBASE}.txt

#source /afs/cern.ch/na61/Releases/SHINE/pro/scripts/env/lxplus_64bit_slc6.csh
#setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/afs/cern.ch/user/m/munger/Prop/lib

#$EXEDIR/bin/createPropMatrixFile $EVOLUTION $PIONFILE $INDIR/$PHOTONFIELD/*.root
#mv propMatrix_${EVOLUTION}.root $OUTDIR/${PHOTONFIELD}_${EVOLUTION}.root
#mv propMatrix_${EVOLUTION}_nu.root $OUTDIR/${PHOTONFIELD}_${EVOLUTION}_nu.root
