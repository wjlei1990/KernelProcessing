#!/bin/bash

#PBS -A GEO111
#PBS -N MERGE_KERNELS
#PBS -j oe
#PBS -o merge_kernels.$PBS_JOBID.log
#PBS -l walltime=01:00:00
#PBS -l nodes=24

NPROC=384
runcmd=mpirun

cd $PBS_O_WORKDIR
echo "pwd: `pwd`"
echo "time: `date`"

echo "========================================="
echo "submit merging kernels: `date`"

kernel_dir="/lustre/atlas2/geo111/proj-shared/Wenjie/kernels_M24/new/10"
echo "kernel dir: $kernel_dir"
if [ ! -d $kernel_dir ]; then
  echo "No kernel dir!!!"
  exit
fi

outputfn="/lustre/atlas2/geo111/proj-shared/Wenjie/kernels_M24/new/10/kernels_smooth.bp"

echo "mpirun -n $NPROC ./bin/xmerge_kernels $kernel_dir $outputfn"
$runcmd -n $NPROC ./bin/xmerge_kernels $kernel_dir $outputfn
echo "done successfully: `date`"
