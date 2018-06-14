#!/bin/sh

#PBS -A GEO111
#PBS -N sd_direction
#PBS -j oe
#PBS -o sd_direction.$PBS_JOBID.log
#PBS -l walltime=0:10:00
#PBS -l nodes=24

NPROC=384
runcmd=mpirun

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

# dir for this iteration, kernels_precond.bp
gradient="/lustre/atlas2/geo111/proj-shared/Wenjie/kernels_M24/new/10/kernels_precond.bp"

# dir for this iteration, cg_direction.bp
direction="/lustre/atlas2/geo111/proj-shared/Wenjie/kernels_M24/new/10/sd_direction.bp"

if [ ! -f $gradient ]; then
  echo Missing grad file: $gradient
  exit
fi

echo submit compute sd direction
echo `date`
echo "$runcmd -n $NPROC ./bin/xsteepDescent $gradient $direction"
echo "----------------"
$runcmd -n $NPROC ./bin/xsteepDescent $gradient $direction
echo "----------------"
echo done successfully: `date`
