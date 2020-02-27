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

# ########################################################################
iter_old=M23
step_lens=(0.02)
#step_lens=(0.0064)

basedir="/lustre/atlas2/geo111/proj-shared/Wenjie/PSF/model.M23"
#input_kernel="$basedir/afar_psf/psf_direction.bp"
#outputbase="$basedir/afar_psf"

input_kernel="$basedir/pacific/psf_direction.bp"
outputbase="$basedir/pacific"

input_model="$basedir/$iter_old/model_gll.bp"
input_solver="$basedir/$iter_old/solver_data.bp"

# ########################################################################

if [ ! -f $input_model ]; then
	echo "Missing input_model: $input_model"
	exit
fi
echo "Input model file: $input_model"

if [ ! -f $input_solver ]; then
  echo "Missing input_solver: $input_solver"
  exit
fi
echo "input solver file: $input_solver"

if [ ! -f $input_kernel ]; then
	echo "Missing input_kernel: $input_kernel"
	exit
fi
echo "input kernel file: $input_kernel"

for sl in ${step_lens[@]}
do
  echo "=============================="
  echo "Model update at step length: $sl"
  outputdir="$outputbase/model_$iter_old-$sl-psf"
  if [ ! -d $outputdir ]; then
	  mkdir -p $outputdir
  fi
  rm $outputdir/*
  echo "Output dir: $outputdir"
  echo "$runcmd -n $NPROC ../../bin/xupdate_model $sl $input_model $input_solver $input_kernel $outputdir"
  $runcmd -n $NPROC ../../bin/xupdate_model $sl $input_model $input_solver $input_kernel $outputdir
  echo "Done at: `date`"
done

echo
echo "done successfully: `date`"
