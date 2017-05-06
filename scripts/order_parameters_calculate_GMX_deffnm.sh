#!/bin/bash

# Gromacs 5.x version of ...
# bash wrapper for calculating Order Parameters of lipid bilayer 
# python script/library calcOrderParameters.py
# meant for use with NMRlipids projects
#------------------------------------------------------------
# Made by J.Melcr,  Last edit 2017/03/21
#------------------------------------------------------------

scriptdir=`dirname $0`


 die() {
       echo $@
       exit 1
 }


[ -z $1 ] && die "Default file name not specified.   Usage: this_script.sh DFNM"
dfnm=$1

traj_file_name=$dfnm".xtc" #"../traj.trr" 
traj_pbc_nonwat_file_name=$dfnm"_nonwat_pbc.xtc" #"../traj.trr" 
top_file_name=$dfnm"_last_frame_nonwat.gro"
tpr_file_name=$dfnm".tpr"
#op_def_file="../../Headgroup_Glycerol_OPs.def"
op_def_file="../../order_parameter_definitions_POPC_all.def"
op_out_file="OrdPars.dat"
top="topol.top"
f_conc=55430  # in mM/L

if ! [ -s $tpr_file_name ] 
then
    echo "We really need " $tpr_file_name " , but we can't find it!"
    exit 1
fi

# remove PBC:
! [ -s $traj_pbc_nonwat_file_name ] && echo non-water | gmx trjconv -f $traj_file_name -s $tpr_file_name -o $traj_pbc_nonwat_file_name -pbc mol

# get a non-water gro-file (topology)
if ! [ -s $top_file_name ] 
then
    if [ -s ${dfnm}.cpt ]
    then
        echo non-water | gmx trjconv -f ${dfnm}.cpt -s $tpr_file_name -o $top_file_name -pbc mol
    else
        echo "Couldn't find ${dfnm}.cpt"
        exit 1
    fi
fi

#CALCULATE ORDER PARAMETERS
python $scriptdir/calcOrderParameters.py -i $op_def_file -t $top_file_name -x $traj_pbc_nonwat_file_name -o $op_out_file && rm $traj_pbc_nonwat_file_name


#getting concentration from topol.top file (if exists)
if [ -f $top ]
then
    nwat=`grep -e "molecules" -A 10 $top | grep -e "^SOL" -e "^TIP" -e "^OPC3" | cut -d " " -f1 --complement `
    nion=`grep -e "molecules" -A 10 $top | grep -e "^NA"  -e "^CA"  | cut -d " " -f1 --complement `
    [ -z $nion ] && nion=0

    conc=`echo $f_conc "*" $nion / $nwat  | bc`
    echo conc: $conc
    sed "s/conc/$conc/" -i ${op_out_file}.line
else
    echo "Topology probably not present, can't calculate concentration."
fi
