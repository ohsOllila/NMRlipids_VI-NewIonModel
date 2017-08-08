#!/bin/bash

# Gromacs 5.x version of ...
# bash wrapper for calculating Order Parameters of lipid bilayer 
# python script/library calcOrderParameters.py
# meant for use with NMRlipids projects
#------------------------------------------------------------
# Made by J.Melcr,  Last edit 2017/03/21
#------------------------------------------------------------

scriptdir=`dirname $0`

traj_file_name="traj_comp.xtc" #"../traj.trr" 
traj_pbc_nonwat_file_name="traj_nonwat_pbc.xtc" #"../traj.trr" 
top_file_name="last_frame_nonwat.gro"
tpr_file_name="topol.tpr"
dens_file_name="density_ca_cl_water.xvg"
#op_def_file="../../Headgroup_Glycerol_OPs.def"
op_def_file=${scriptdir}"/order_parameter_definitions_Lipid14_POPC_all.def"
op_out_file="OrdPars.dat"
top="topol.top"
f_conc=55430  # in mM/L

if ! [ -s $tpr_file_name ] 
then
    echo "We really need " $tpr_file_name " , but we can't find it!"
    exit 1
fi

# remove PBC:
# center the trajectory around POPC (should be no. 2)
# and start at 100ns
! [ -s $traj_pbc_nonwat_file_name ] && echo 2 non-water | gmx trjconv -f $traj_file_name -s topol.tpr -o $traj_pbc_nonwat_file_name -pbc mol -center -b 100000 #-n index

# get a non-water gro-file (topology)
if ! [ -s $top_file_name ] 
then
    if [ -s state.cpt ]
    then
        echo non-water | gmx trjconv -f state.cpt -s topol.tpr -o $top_file_name -pbc mol
    else
        echo "Couldn't find state.cpt"
        exit 1
    fi
fi

# get densities of cation - anion - water centered around POPC
# and start at 100ns
seq 2 5 | gmx density -sl 900 -dens number -ng 3 -f $traj_file_name -center -symm -relative -o $dens_file_name -xvg none -b 100000 

# rename resnames of palmitoyl and oleoyl segments 
python $scriptdir/rename_residue_lipid14_to_PALM-POPC-OLE.py -i $top_file_name -o $top_file_name 

#CALCULATE ORDER PARAMETERS
python $scriptdir/calcOrderParameters.py -i $op_def_file -t $top_file_name -x $traj_pbc_nonwat_file_name -o $op_out_file && rm $traj_pbc_nonwat_file_name

# getting actual bulk concentration from number density profile
python $scriptdir/get_conc_ion_bulk.py 
conc=`cut -d " " -f1 conc_ion_bulk_mmolL.dat`
# write the actual bulk concentration into the OP-output file
sed "s/conc/$conc/" -i ${op_out_file}.line

#getting nominal concentration from topol.top file (if exists)
if [ -f $top ]
then
    nwat=`grep -e "molecules" -A 10 $top | grep -e "^SOL" -e "^TIP" -e "^SPCE" -e "^OPC3" | cut -d " " -f1 --complement `
    nion=`grep -e "molecules" -A 10 $top | grep -e "^NA"  -e "^CA"  | cut -d " " -f1 --complement `
    [ -z $nion ] && nion=0

    concnom=`echo $f_conc "*" $nion / $nwat  | bc`
    echo nominal conc: $concnom
    #sed "s/conc/$conc/" -i ${op_out_file}.line
    echo $concnom > conc_ion_nominal_mmolL.dat
else
    echo "Topology probably not present, can't calculate concentration."
fi

