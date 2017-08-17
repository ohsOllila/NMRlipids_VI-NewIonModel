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

# get densities of cation - anion - water centered around POPC
# and start at 100ns
if ! [ -s $dens_file_name ] 
then
    seq 2 5 | gmx density -sl 900 -dens number -ng 3 -f $traj_file_name -center -symm -relative -o $dens_file_name -xvg none -b 100000 
fi

# calculate area and area-per-lipid
area=`bash $scriptdir/area-per-lipid_calculate.sh | grep -e "area" | cut -d " " -f3`

# getting actual bulk concentration from number density profile
python $scriptdir/get_conc_ion_bulk.py 
conc=`cut -d " " -f1 conc_ion_bulk_mmolL.dat`
conc_anion=`cut -d " " -f1 conc_anion_bulk_mmolL.dat`
conc_wat=`cut -d " " -f1 conc_wat_bulk_mmolL.dat`

#getting nominal concentration from topol.top file (if exists)
if [ -f $top ]
then
    nwat=`grep -e "molecules" -A 10 $top | grep -e "^SOL" -e "^TIP" -e "^SPCE" -e "^OPC3" | cut -d " " -f1 --complement `
    nion=`grep -e "molecules" -A 10 $top | grep -e "^NA"  -e "^CA"  | cut -d " " -f1 --complement `
    [ -z $nion ] && nion=0

    nanion=`grep -e "molecules" -A 10 $top | grep -e "^CL"  | cut -d " " -f1 --complement `
    [ -z $nanion ] && nanion=0

    concnom=`echo $f_conc "*" $nion / $nwat  | bc`
    echo nominal conc: $concnom
    #sed "s/conc/$conc/" -i ${op_out_file}.line
    echo $concnom > conc_ion_nominal_mmolL.dat
else
    echo "Topology probably not present, can't calculate concentration."
fi

surfexc=`echo scale=4 ; ($nwat "-" $nion '*' $conc_wat"/"$conc)/("2*"$area)  | bc`
surfexcani=`echo scale=4 ; ($nwat "-" $nanion '*' $conc_wat"/"$conc_anion)/("2*"$area)  | bc`

echo Relative surface excess of water--ions = $surfexc > surf_excess_wat_ions.dat
echo Relative surface excess of water--anions = $surfexcani > surf_excess_wat_anions.dat
echo Relative surface excess of water--ions = $surfexc
