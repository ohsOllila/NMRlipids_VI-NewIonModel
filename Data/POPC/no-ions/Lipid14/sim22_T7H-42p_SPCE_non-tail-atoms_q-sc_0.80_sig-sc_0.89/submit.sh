#!/bin/bash

#PBS -N L14_headgr80_0mM
#PBS -q gpu_long
#PBS -l walltime=167h
#PBS -l mem=20gb
#PBS -l nodes=1:ppn=1:gpu=1:nfs4:cl_doom
#PBS -l scratch=10gb
#PBS -j oe



export stob3='/storage/brno3-cerit/home/melcrj'
export storage=$stob3
export workfolder=$PBS_O_WORKDIR   # $storage/wilCleGFP_3lipids_POPS/md12_dAFED/val-met_dist_dAFED


# suggestion from https://wiki.metacentrum.cz/wiki/Pr%C3%A1ce_s_daty_v_%C3%BAloze
trap 'cp -r $SCRATCH/* $workfolder/ && clean_scratch' TERM




echo 'load modules ... '
#  also loads gcc-4.8.1 automatically
module load cmake-3.6.1 swig-3.0.8 cuda-8.0  doxygen-1.8.6 fftw-3.3.4omp gcc-4.8.4 
module load gcc-4.9.2 
module load gcc-5.3.0 

# get rid of annoying warning from openMPI
export OMPI_MCA_mpi_warn_on_fork=0



# Just a basic initial printout
echo
echo 'Basic job -- info'
echo '-----------------'
echo 'job name: ' $PBS_JOBNAME
echo 'job id: ' $PBS_JOBID
echo 'name: ' $name
echo 'running on: ' `hostname`
echo 'proc per node: ' $PBS_NUM_PPN
echo 'num nodes: ' $PBS_NUM_NODES
echo 'Req-d walltime: ' $TORQUE_RESC_TOTAL_WALLTIME
# echo 'nodelist: (PBS_NODEFILE) '
# cat $PBS_NODEFILE
echo 

cd $workfolder          || exit 1
echo 'cwd: ' `pwd`
echo


echo `env` > $infofile.env



# check hwether previous run ended/copied OK
if [ -e ___COPY_FAILED___ ] 
then
   echo previous copy had failed, ___COPY_FAILED___ file exist.
   echo Solve the problem and remove this file after that. 
   exit 3
fi


infofile=info.on.job.${PBS_JOBID%.arien*}
echo ' Info on job: ' $PBS_JOBNAME > $infofile
echo ' job ID: ' $PBS_JOBID >> $infofile
echo ' running on: ' `hostname` >> $infofile
echo ' started at: ' `date` >> $infofile


# Copy things to Scratch - not necessary, working in the "results" folder
#cp -a $workfolder/topol*.tpr  $SCRATCH/  || exit 2
#cp -a $workfolder/state*.cpt  $SCRATCH/  || exit 2


echo 'current time:  ' `date`
echo







time_scale=0.97
maxh=`echo "scale=2 ; " $TORQUE_RESC_TOTAL_WALLTIME "/ 3600 * " $time_scale | bc`


########################################
#CALCULATION
########################################

# this var must be adopted to the specific simulation
export OMP_NUM_THREADS=$PBS_NUM_PPN

# add openMM python wrapper to PYTHONPATH
export PYTHONPATH=/storage/brno3-cerit/home/melcrj/prog/BIN/openmm7.1.0/python:$PYTHONPATH

python simulate_oMM_psf_test_1x4fs_langevin_CUDA.py





########################################
# Copy things back to result folder
########################################

mv $SCRATCH/* $workfolder/
if [ $? -ne 0 ]; then 
  echo
  echo Copy output data failed. Copy them manualy from `hostname` 1>&2 
  touch ___COPY_FAILED___
  exit 1
else
  rm -rf $SCRATCH/*
fi


echo
echo 'current time:  ' `date`
echo

echo   >> $infofile
echo ' ended at:   ' `date` >> $infofile

exit 0

