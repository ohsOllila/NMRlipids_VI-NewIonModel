#!/bin/bash

#PBS -N sim22a_sim
#PBS -q gpu
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=8:ngpus=1:cl_doom=True:scratch_local=10gb:mem=20gb
#PBS -j oe



export stob3='/storage/brno3-cerit/home/melcrj'
export storage=$stob3
export workfolder=$PBS_O_WORKDIR   # $storage/wilCleGFP_3lipids_POPS/md12_dAFED/val-met_dist_dAFED

export SCRATCH=$SCRATCHDIR


# suggestion from https://wiki.metacentrum.cz/wiki/Pr%C3%A1ce_s_daty_v_%C3%BAloze
trap 'cp -r $SCRATCH/* $workfolder/ && clean_scratch' TERM




echo 'load modules ... '
#  also loads gcc-4.8.1 automatically
module load intelcdk-15 openmpi-1.8.2-gcc cuda-7.5 cmake-3.2.3 boost-1.57-gcc 



########################################
# load my build of Gromacs 5.1.2 with Plumed2
########################################

mpirun_command=""

if [ $PBS_NUM_NODES -eq 1 ]
then
	#
	# don't use MPI (intel-15 compiler)
	#
	# if "avx" is not contained in /proc/cpuinfo, load sse4.1 version of Gromacs
	if [ `grep -e 'avx' /proc/cpuinfo -c` -eq 0 ]
	then
	   # sse4.1
	   source $stob3/prog/BIN/gromacs-5.1.2+plumed2_ownFFTW_noMPI_intel15_sse4/bin/GMXRC
	   export gmx_mdrun="gmx mdrun"
	else
	   # avx_256 mdrun only
	   export PATH=$stob3/prog/BIN/gromacs-5.1.2+plumed2_ownFFTW_noMPI_intel15_avx_256/bin:$PATH
	   export gmx_mdrun="mdrun"
	fi
else
	mpirun_command="mpirun -n" $PBS_NUM_NODES
	#
	# use MPI (gcc 4.8.1 compiler)
	#
	# if "avx" is not contained in /proc/cpuinfo, load sse4.1 version of Gromacs
	if [ `grep -e 'avx' /proc/cpuinfo -c` -eq 0 ]
	then
	   # sse4.1
	   source $stob3/prog/BIN/gromacs-5.1.2+plumed2_ownFFTW_openMPI_gcc4.8.1_sse4/bin/GMXRC
	   export gmx_mdrun="gmx_mpi mdrun"
	else
	   # avx_256 mdrun only
	   export PATH=$stob3/prog/BIN/gromacs-5.1.2+plumed2_ownFFTW_openMPI_gcc4.8.1_avx_256/bin:$PATH
	   export gmx_mdrun="mdrun_mpi"
	fi
fi


# PLUMED 2.2.2 env settings
#export PLUMED_ROOT=/storage/brno3-cerit/home/melcrj/prog/plumed-2.2.2
#
export PATH=/storage/brno3-cerit/home/melcrj/prog/BIN/plumed-2.2.2/bin:$PATH
export INCLUDE=/storage/brno3-cerit/home/melcrj/prog/BIN/plumed-2.2.2/include:$INCLUDE
export LD_LIBRARY_PATH=/storage/brno3-cerit/home/melcrj/prog/BIN/plumed-2.2.2/lib:$LD_LIBRARY_PATH
# for runtime binding:
export PLUMED_KERNEL=/storage/brno3-cerit/home/melcrj/prog/BIN/plumed-2.2.2/lib/libplumedKernel.so
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
OMP_NUM_THREADS=$PBS_NUM_PPN

$mpirun_command $gmx_mdrun \
   -maxh $maxh \
   -ntomp $OMP_NUM_THREADS \
   -noappend   \
   -cpi state.cpt       \
   -s   topol.tpr       \
   -o   $SCRATCH/traj.trr  \
   -x   $SCRATCH/traj.xtc  \
   -e   $SCRATCH/ener.edr  \
   -cpo $SCRATCH/state.cpt \
   -g   $SCRATCH/md.log    \

   #-plumed plumed.dat   \

   # causes problems with multiple simulations at one node
   #-pin on






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

