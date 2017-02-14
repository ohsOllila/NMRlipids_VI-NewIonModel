import MDAnalysis as mda
import numpy as np
u = mda.Universe("conf_init.gro", "traj_openMM.dcd")
u.trajectory.ts.dimensions
x = []
y = []
for frame in u.trajectory:
    x.append(frame.dimensions[0])
    y.append(frame.dimensions[1])
nlip = np.compress(u.residues.resnames == 'POPC', u.residues.resnames).shape[0]
apl = np.mean(x)*np.mean(y)/(0.5*nlip)
print "APL = {apl}".format(apl=apl)
