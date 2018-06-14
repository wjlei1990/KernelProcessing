#!/bin/sh

#PBS -A GEO111
#PBS -N CG_direction
#PBS -j oe
#PBS -o job_cg_direction.$PBS_JOBID.log
#PBS -l walltime=0:30:00
#PBS -l nodes=24

NPROC=384
runcmd=mpirun

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

# the file that contains the gradient and model change for each previous iterations
input_path_file="/ccs/proj/geo111/wenjie/AdjointTomography/M23/post_proc/LBFGS.rhea/examples/lbfgs/path.txt"

solver_file="/lustre/atlas/proj-shared/geo111/Wenjie/kernels_M23/post_proc.bb/10/M22_model/model_gll.bp"

outputfn="/lustre/atlas/proj-shared/geo111/Wenjie/kernels_M23/benchmark/new/10/lbfgs_direction.bp"

if [ ! -f $input_path_file ]; then
  echo Missing input_path_file: $input_path_file
  exit
fi

if [ ! -f $solver_file]; then
  echo "Missing solver file: $solver_file"
  exit
fi

echo Compute lbfgs direction
echo "<----------------"
echo `date`
echo "$runcmd -n $NPROC ./bin/xlbfgs $input_path_file $solver_file $outputfn"

$runcmd -n $NPROC ./bin/xlbfgs $input_path_file $solver_file $outputfn
echo "---------------->"
echo done successfully: `date`
