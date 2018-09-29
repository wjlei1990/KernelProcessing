#!/bin/bash

#PBS -A GEO111
#PBS -N psf
#PBS -j oe
#PBS -o job_psf.$PBS_JOBID.log
#PBS -l walltime=0:30:00
#PBS -l nodes=24

NPROC=384
runcmd=mpirun

cd $PBS_O_WORKDIR
echo "pwd: `pwd`"

solver_file="/lustre/atlas2/geo111/proj-shared/Wenjie/PSF/model.M23/M23/solver_data.bp"
outputfn="/lustre/atlas2/geo111/proj-shared/Wenjie/PSF/model.M23/psf_direction.bp"

if [ ! -f $solver_file ]; then
  echo "Missing solver file: $solver_file"
  exit
fi

echo Compute Gaussian perturbation for PSF
echo `date`
echo "$runcmd -n $NPROC ../../bin/xgauss_psf $solver_file $outputfn"
echo
$runcmd -n $NPROC ../../bin/xgauss_psf $solver_file $outputfn

echo
echo done successfully: `date`
