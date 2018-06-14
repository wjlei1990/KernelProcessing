#!/bin/bash

#PBS -A GEO111
#PBS -N model_update
#PBS -j oe
#PBS -o model_update.$PBS_JOBID.log
#PBS -l walltime=0:30:00
#PBS -l nodes=24

NPROC=384
runcmd=mpirun

cd $PBS_O_WORKDIR
echo "pwd: `pwd`"

iter_old=M23
iter_new=M24
#step_lens=(0.0025 0.0050 0.0075 0.0100 0.0125)
step_lens=(0.0100)

#basedir="/lustre/atlas/proj-shared/geo111/Wenjie/kernels_M23/post_proc.bb/10"
#input_model=$basedir"/"$iter_old"_model/model_gll.bp"
#input_solver="$basedir/"$iter_old"_model/solver_data.bp"
#input_kernel="$basedir/cg_direction.bp"

input_model="/lustre/atlas2/geo111/proj-shared/Wenjie/kernels_M24/new/10/M23/model_gll.bp"
input_solver="/lustre/atlas2/geo111/proj-shared/Wenjie/kernels_M24/new/10/M23/solver_data.bp"
input_kernel="/lustre/atlas2/geo111/proj-shared/Wenjie/kernels_M24/new/10/cg_direction.bp"
outputbase="/lustre/atlas2/geo111/proj-shared/Wenjie/kernels_M24/new/10/"

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
  outputdir="$outputbase/model_$iter_new-$sl"
  if [ ! -d $outputdir ]; then
	  mkdir -p $outputdir
  fi
  rm $outputdir/*
  echo "Output dir: $outputdir"
  echo "aprun -n $NPROC ./bin/xupdate_model $sl $input_model $input_solver $input_kernel $outputdir"
  $runcmd -n $NPROC ./bin/xupdate_model $sl $input_model $input_solver $input_kernel $outputdir
  echo "Done at: `date`"
done

echo
echo "done successfully: `date`"
