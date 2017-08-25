#!/bin/bash

# openMM version of ...
# bash wrapper for calculating Order Parameters of lipid bilayer 
# python script/library calcOrderParameters.py
# meant for use with NMRlipids projects
#------------------------------------------------------------
# Made by J.Melcr,  Last edit 2017/03/21
#------------------------------------------------------------

scriptdir=`dirname $0`

traj_file_name="traj_openMM" #"../traj.trr" 
top_file_name="conf_init.gro"
top_file_name2="conf_init_PopcPalmOle.gro"
#op_def_file="../../Headgroup_Glycerol_OPs.def"
#op_def_file="../../order_parameter_definitions_POPC_all.def"
op_def_file=${scriptdir}"/order_parameter_definitions_Lipid14_POPC_all.def"
op_out_file="OrdPars.dat"
top="topol.top"
f_conc=55430  # in mM/L


# rename resnames of palmitoyl and oleoyl segments 
python $scriptdir/rename_residue_lipid14_to_PALM-POPC-OLE.py -i $top_file_name -o $top_file_name2 

#CALCULATE ORDER PARAMETERS
# can operate directly on DCD trajectory, cause it has PBC removed already
python $scriptdir/calcOrderParameters.py -i $op_def_file -t $top_file_name2 -x $traj_file_name -o $op_out_file 


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

