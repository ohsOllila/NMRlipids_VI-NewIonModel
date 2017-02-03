#!/bin/bash

# bash wrapper for calculating Order Parameters of lipid bilayer
# meant for use with NMRlipids projects

scriptdir=`dirname $0`

traj_file_name="traj_comp.xtc" #"../traj.trr" 
traj_pbc_nonwat_file_name="traj_nonwat_pbc.xtc" #"../traj.trr" 
top_file_name="last_frame_nonwat.gro"
tpr_file_name="topol.tpr"
#op_def_file="../../Headgroup_Glycerol_OPs.def"
op_def_file="../../order_parameter_definitions_POPC_all.def"
op_out_file="Headgroup_Glycerol_OPs.dat"

if ! [ -s $tpr_file_name ] 
then
    echo "We really need " $tpr_file_name " , but we can't find it!"
    exit 1
fi

# remove PBC:
! [ -s $traj_pbc_nonwat_file_name ] && echo non-water | gmx trjconv -f $traj_file_name -s topol.tpr -o $traj_pbc_nonwat_file_name -pbc mol

# get a non-water gro-file (topology)
if ! [ -s $top_file_name ] 
then
    if [ -s state.cpt ]
    then
        echo non-water | gmx trjconv -f state.cpt -s topol.tpr -o $top_file_name -pbc mol
    else
        echo "Couldn't find state.cpt"
    fi
fi

#CALCULATE ORDER PARAMETERS
python $scriptdir/calcOrderParameters.py -i $op_def_file -t $top_file_name -x $traj_pbc_nonwat_file_name -o $op_out_file && rm $traj_pbc_nonwat_file_name


#getting concentration from topol.top file (if exists)
top="topol.top"
f_conc=55430  # in mM/L

if [ -f $top ]
then
    nwat=`grep -e "^SOL" -e "^TIP" $top | cut  -f1 --complement `
    nion=`grep -e "^NA" -e "^CA" $top | cut -d " " -f1 --complement `
    [ -z $nion ] && nion=0

    conc=`echo $f_conc "*" $nion / $nwat  | bc`
    echo conc: $conc
    sed "s/conc/$conc/" -i ${op_out_file}.line
else
    echo "Topology probably not present, can't calculate concentration."
fi
