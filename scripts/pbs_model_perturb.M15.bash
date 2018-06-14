#!/bin/bash

#PBS -A GEO111
#PBS -N perturb_model
#PBS -j oe
#PBS -o job_perturb.$PBS_JOBID.log
#PBS -l walltime=1:00:00
#PBS -l nodes=24

NPROC=384

cd $PBS_O_WORKDIR
echo "pwd: `pwd`"

ref_model_file="/lustre/atlas/proj-shared/geo111/Wenjie/M15_NEX256/model_gll.bp"
new_model_file="/lustre/atlas/proj-shared/geo111/Wenjie/specfem3d_globe_990cd4.M24/specfem3d_globe_990cd4/DATA/GLL/model_gll.bp"
outputfn="/lustre/atlas/proj-shared/geo111/Wenjie/kernels_M23/model_perturb/model_M23-M15.bp"

if [ ! -f $ref_model_file ]; then
  echo "Missing ref_model_file: $ref_model_file"
  exit
fi

if [ ! -f $new_model_file ]; then
  echo "Missing new_model_file: $new_model_file"
  exit
fi

echo "Time:`date`"
echo

echo "mpirun -n 384 ./bin/xmodel_perturb_ref $ref_model_file $new_model_file $outputfn"
mpirun -n 384 ./bin/xmodel_perturb_ref $ref_model_file $new_model_file $outputfn

echo
echo "All done at: `date`"
