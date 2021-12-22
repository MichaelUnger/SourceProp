#!/bin/bash
#PBS -l nodes=1:ppn=1,walltime=15:00:00
#PBS -N crp
#PBS -e localhost:/scratch/mu495/logs/${PBS_JOBNAME}_${PBS_JOBID}.e
#PBS -o localhost:/scratch/mu495/logs/${PBS_JOBNAME}_${PBS_JOBID}.o

JOBNUMBER=`echo ${PBS_JOBID} | awk -F. '{print $1}'`
WORKDIR=/state/partition1/mu495_${PBS_JOBID}
mkdir $WORKDIR
cd $WORKDIR
echo running on $HOSTNAME
module  load cmake/gnu/2.8.12
module load python/gnu/2.7.5
export PATH=$PATH:/home/mu495/Software/swig/swig-3.0.2_install/bin
. /home/mu495/Software/ROOT/root_v5.34.24_install/bin/thisroot.sh
. /home/mu495/Software/CRPropaNew/src/CRPropaTools/Run/setpaths_butinah.sh
python /home/mu495/Software/CRPropaNew/src/CRPropaTools/Run/run1D.py
mv crpropa.root /scratch/mu495/CRPropa/crp_${JOBNUMBER}.root
cd
rm -rf $WORKDIR
