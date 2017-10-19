#!/bin/bash

# Gromacs 5.x version of ...
# bash wrapper for calculating relative surface excess of a lipid bilayer 
# using
# python script/library get_conc_ion_bulk.py
# meant for use with NMRlipids projects
#------------------------------------------------------------
# Made by J.Melcr,  Last edit 2017/08/21
#------------------------------------------------------------

scriptdir=`dirname $0`

traj_file_name="traj_openMM.dcd" #"../traj.trr" 
#traj_pbc_nonwat_file_name="traj_nonwat_pbc.xtc" #"../traj.trr" 
#top_file_name="last_frame_nonwat.gro"
tpr_file_name="topol.tpr"
dens_file_name="density_ca_cl_water.xvg"
#op_def_file="../../Headgroup_Glycerol_OPs.def"
#op_def_file=${scriptdir}"/order_parameter_definitions_Lipid14_POPC_all.def"
#op_out_file="OrdPars.dat"
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
    seq 2 5 | gmx density -sl 900 -dens number -ng 3 -f $traj_file_name -center -symm -relative -o $dens_file_name -xvg none #-b 100000 
fi

# calculate area and area-per-lipid
area=`python $scriptdir/get_apl_fromDCD.py | grep -e "APL" | cut -d " " -f3`
echo area = $area nm2

# getting actual bulk concentration from number density profile
python $scriptdir/get_conc_ion_bulk.py 
conc=`cut -d " " -f1 conc_ion_bulk_mmolL.dat`
conc_anion=`cut -d " " -f1 conc_anion_bulk_mmolL.dat`
conc_wat=`cut -d " " -f1 conc_wat_bulk_mmolL.dat`


#getting nominal concentration from topol.top file (if exists)
if [ -f $top ]
then
    nwat=`grep -e "molecules" -A 10 $top | grep -e "^SOL" -e "^TIP" -e "^SPCE" -e "^OPC3" | cut -d " " -f1 --complement `
    # correct the number density of water with the number of its atoms
    # unfortunately, this does not really work, all 4-site models I employ use name SOL, however, some 3-site use it too
    nwat4=`grep -e "molecules" -A 10 $top | grep -e "^TIP4" -e "^OPC4" | cut -d " " -f1 --complement `
    if ! [ -z $nwat4 ]
    then
       watdenom=4.0
    else
       watdenom=3.0
    fi

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

# correct the number density of water with the number of its atoms
conc_wat=`echo "scale=4 ; "$conc_wat "/" $watdenom | bc `
echo water concentration = $conc_wat mmol/L

surfexc=`echo "scale=4 ; ("$nion "-" $nwat '*' $conc"/"$conc_wat")/(2*"$area")"  | bc`
surfexcani=`echo "scale=4 ; ("$nanion "-" $nwat '*' $conc_anion"/"$conc_wat")/(2*"$area")"  | bc`

echo Relative surface excess of ions  --water = $surfexc nm-2 > surf_excess_wat_ions.dat
echo Relative surface excess of ions  --water = $surfexc nm-2 
echo Relative surface excess of anions--water = $surfexcani nm-2 > surf_excess_wat_anions.dat
echo Relative surface excess of anions--water = $surfexcani nm-2 

