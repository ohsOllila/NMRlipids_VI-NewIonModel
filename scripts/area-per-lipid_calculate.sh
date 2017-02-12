#!/bin/bash

# wrapping BASH script for calc_apl.py
#  uses GMX tools to calculate Area per Lipid (APL)
# meant for use with NMRlipids projects

scriptdir=`dirname $0`

edr_file_name="ener.edr"
xvg_file_name="energies.xvg"
boxdim_file_name="box_dim.txt"
top="topol.top"

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
    nlip=`grep -e "molecules" -A 10 $top | grep -e "^POPC" | cut  -f1 --complement `
    nlip_per_leaflet=`echo $nlip/2.00 | bc`
else
    echo "Topology probably not present, can't get no. POPC lipids, assuming 64."
    nlip_per_leaflet=64
fi

#CALCULATE APL (per POPC only!)
python $scriptdir/calc_apl.py -f $boxdim_file_name -n $nlip_per_leaflet > apl_${nlip_per_leaflet}lip-leaflet.dat

echo "Lipids per leaflet: " $nlip_per_leaflet
cat apl_${nlip_per_leaflet}lip-leaflet.dat
