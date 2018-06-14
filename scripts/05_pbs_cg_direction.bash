#!/bin/bash

#PBS -A GEO111
#PBS -N CG_direction
#PBS -j oe
#PBS -o cg_direction.$PBS_JOBID.log
#PBS -l walltime=0:30:00
#PBS -l nodes=24

NPROC=384
runcmd=mpirun

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

solver_file="/lustre/atlas2/geo111/proj-shared/Wenjie/kernels_M24/new/10/M23/solver_data.bp"

# dir for last iteration, kernels_precond.bp
gradient_0="/lustre/atlas/proj-shared/geo111/rawdata/M23/kernel/10/kernels_precond.bp"
# dir for this iteration, kernels_precond.bp
gradient_1="/lustre/atlas2/geo111/proj-shared/Wenjie/kernels_M24/new/10/kernels_precond.bp"

# dir for last iteration, cg_direction.bp
direction_0="/lustre/atlas/proj-shared/geo111/rawdata/M23/kernel/10/direction_M23/cg_direction.bp"
# dir for this iteration, cg_direction.bp
direction_1="/lustre/atlas2/geo111/proj-shared/Wenjie/kernels_M24/new/10/cg_direction.bp"

if [ ! -f $solver_file]; then
  echo "Missing solver file: $solver_file"
  exit
fi

if [ ! -f $direction_0 ]; then
  echo Error direction 0: $direction_0
  exit
fi
if [ ! -f $gradient_0 ]; then
  echo Missing grad 0: $gradient_0
  exit
fi
if [ ! -f $gradient_1 ]; then
  echo Missing grad 1: $gradient_1
  exit
fi

echo submit compute cg direction
echo `date`
echo "aprun -n $NPROC ./bin/xcg_direction $gradient_0 $gradient_1 $direction_0 $solver_file $direction_1"
echo "----------------"
$runcmd -n $NPROC ./bin/xcg_direction $gradient_0 $gradient_1 $direction_0 $solver_file $direction_1
echo "----------------"
echo done successfully: `date`
