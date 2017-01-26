#!/bin/bash 
# analysis of tajectories with double membranes for the use with mempot project
# uses default filenaming comming from mdrun -deffnm ....
#  assuming $indexfile file 4 levels above in the dir tree
#       and $electronsfile 5 levels above in the slipids_top dir
#       and some other stuff (like existing groups in indexfiles and so)
# Josef Joe Melcr 15/04/2015

 die() {
       echo $@
       exit 1
 }


 [ -z $1 ] && die "path not specified.   Usage: this_script.sh PATH"

 simpath=$1  # is path-to-topol in this script
 dfnm="ana"  # just a name for analysis-results files

 cwd=`pwd`
 cd $simpath

 # add mempot/scripts dir to PATH
 PATH=${PATH}":/mnt/scratch/scaledip/scripts/"
 # path to general index file -- same for all mempot methods (ion-imbal/const-field/none)
 # relative path should be conserved among trajectories
 indexfile=${cwd}"/index.ndx"
 electronsfile="../../../../../slipids_top/electrons_SlipidsVS.dat"
 trajfilename="traj_comp.xtc"



 # trajectory with pbc removed for VMD and for g_potential (where it is crucial!) and other analysis with POPC-centering
 # echo System | trjconv -s $dfnm -f $dfnm -skip 20 -pbc mol -o ${dfnm}_traj_sk20_pbc.xtc
 #echo POPC System | trjconv -s $dfnm -f $dfnm -n $indexfile -center -pbc mol -o $trajfilename


 # 17-19 should be box- x,y,z dimensions
 echo 17 18 19 | g_energy > ${dfnm}.boxdim
 calc_apl.py -f ${dfnm}.boxdim > ${dfnm}.apl
 

 # groups for desity profile plots (9)
 groups="NA   CL   Water" 
 for type in number
 do
   echo $groups | gmx density -n $indexfile -f $trajfilename -dens $type -ng 3 -sl 3500 -o ${dfnm}_densities_${type}
 done

 # correct acts on chosen groups (not charge groups) and corrects the total charge and total field intensity in all slices so that they are zero as (charge should be, field not necessary) (checked in source code, Joe)
 # super-large value for number of slices required for precision (~2pm thick slices are approx optimal; more slices=more noise and fluctuation, less slices=lower precision of integral {larger step})
 echo "System" | gmx potential -n $indexfile -f $trajfilename  -sl  3500 -o ${dfnm}_potential -of ${dfnm}_field -oc ${dfnm}_charge


 exit

