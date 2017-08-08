#!/usr/bin/env python
"""
simple python script for obtaining the actual bulk concentration of 
salts from simulation density profiles

assumes:
  existence of the file density_ca_cl...
  the first density profile is the one we're interested in
  first 50 frames are about right for averageing
"""

import numpy as np
densities = np.loadtxt("density_ca_cl_water.xvg")
ion_dens = densities[:,1]
i=5
denom_conc = 0.602
mean  = ion_dens[:10*i].mean() /denom_conc *1000.0
stdev = ion_dens[:10*i].std()  /denom_conc *1000.0
print "bulk conc. = ({: 4.0f} +/- {: 3.0f}) mmol/L".format(mean, stdev)

with open("conc_ion_bulk_mmolL.dat", "w") as f:
    f.write(" ".join([str(mean), str(stdev)]))

