#!/bin/bash

#PBS -A GEO111
#PBS -N PRECOND_KERNELS
#PBS -j oe
#PBS -o precond_kernels.$PBS_JOBID.log
#PBS -l walltime=01:00:00
#PBS -l nodes=24

NPROC=384
runcmd=mpirun

input_kernel="/lustre/atlas2/geo111/proj-shared/Wenjie/kernels_M24/new/10/kernels_smooth.bp"
output_kernel="/lustre/atlas2/geo111/proj-shared/Wenjie/kernels_M24/new/10/kernels_precond.bp"

cd $PBS_O_WORKDIR
echo "pwd: `pwd`"
echo "Precondition kernels"
echo `date`
echo "============================================"

echo "mpirun -n $NPROC ./bin/xprecond_kernels $input_kernel $output_kernel"
$runcmd -n $NPROC ./bin/xprecond_kernels $input_kernel $output_kernel

echo done successfully
echo `date`
