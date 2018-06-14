#!/bin/bash

#PBS -A GEO111
#PBS -N model_update
#PBS -j oe
#PBS -o job_model_update.$PBS_JOBID.log
#PBS -l walltime=1:00:00
#PBS -l nodes=24

NPROC=384
runcmd=mpirun

iter_new=M22

cd $PBS_O_WORKDIR
echo "pwd: `pwd`"

step_lens=(0.0150)

input_model="/lustre/atlas/proj-shared/geo111/Wenjie/kernels_M22/post_proc/10/M21/model_gll.bp"
input_solver="/lustre/atlas/proj-shared/geo111/Wenjie/kernels_M22/post_proc/10/M21/solver_data.bp"
input_kernel="/lustre/atlas/proj-shared/geo111/Wenjie/kernels_M22/post_proc/10/kernels_precond.bp"

if [ ! -f $input_model ]; then
	echo "Missing input_model: $input_model"
	exit
fi

if [ ! -f $input_solver ]; then
  echo "Missing input_solver: $input_solver"
  exit
fi

if [ ! -f $input_kernel ]; then
	echo "Missing input_kernel: $input_kernel"
	exit
fi

for sl in ${step_lens[@]}
do
  echo "=============================="
  echo "Model update at step length: $sl"
  outputdir="/lustre/atlas/proj-shared/geo111/Wenjie/kernels_M22/post_proc/10/model_$iter_new-$sl"
  if [ ! -d $outputdir ]; then
	  mkdir -p $outputdir
  fi
  rm $outputdir/*
  echo "Output dir: $outputdir"
  log_file=$outputdir/model_update.log
  echo "$runcmd -n $NPROC ./bin/xupdate_model $sl $input_model $input_solver $input_kernel $outputdir &> $log_file"
  $runcmd -n $NPROC ./bin/xupdate_model $sl $input_model $input_solver $input_kernel $outputdir
  echo "Done at: `date`"
done

echo
echo "done successfully: `date`"
