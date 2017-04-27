#!/bin/bash

# wrapping BASH script for calc_apl.py
#  uses GMX tools to calculate Area per Lipid (APL)
# meant for use with NMRlipids projects

scriptdir=`dirname $0`


 die() {
       echo $@
       exit 1
 }


[ -z $1 ] && die "Default file name not specified.   Usage: this_script.sh DFNM"
dfnm=$1

traj_file_name=$dfnm".xtc" #"../traj.trr" 
traj_pbc_nonwat_file_name=$dfnm"nonwat_pbc.xtc" #"../traj.trr" 
top_file_name=$dfnm"_last_frame_nonwat.gro"
tpr_file_name=$dfnm".tpr"
#op_def_file="../../Headgroup_Glycerol_OPs.def"
op_def_file="../../order_parameter_definitions_POPC_all.def"
op_out_file="OrdPars.dat"
top="topol.top"
f_conc=55430  # in mM/L


edr_file_name=$dfnm".edr"
xvg_file_name=$dfnm"_energies.xvg"
boxdim_file_name=$dfnm"_box_dim.txt"
top="lipid.top"

if ! [ -s $edr_file_name ] 
then
    echo "We really need " $edr_file_name " , but we can't find it!"
    exit 1
fi

# should give box-x,y,z
echo 17 18 19 | gmx energy -f $edr_file_name -o $xvg_file_name > $boxdim_file_name

#getting no. POPC lipids from topol.top file (if exists)
if [ -f $top ]
then
    nlip=`grep -e "molecules" -A 10 $top | grep -e "^POPC" | cut -d " "  -f1 --complement `
    nlip_per_leaflet=`echo $nlip/2.00 | bc`
else
    echo "Topology probably not present, can't get no. POPC lipids, assuming 64."
    nlip_per_leaflet=64
fi

#CALCULATE APL (per POPC only!)
python $scriptdir/calc_apl.py -f $boxdim_file_name -n $nlip_per_leaflet > apl_${nlip_per_leaflet}lip-leaflet.dat

echo "Lipids per leaflet: " $nlip_per_leaflet
cat apl_${nlip_per_leaflet}lip-leaflet.dat
