import MDAnalysis as mda
import numpy as np
u = mda.Universe("conf_init.gro", "traj_openMM.dcd")
u.trajectory.ts.dimensions
x = []
y = []
for frame in u.trajectory:
    x.append(frame.dimensions[0])
    y.append(frame.dimensions[1])
nlip = np.compress(u.residues.resnames == 'POPE', u.residues.resnames).shape[0]
apl = np.mean(x)*np.mean(y)/(0.5*nlip) *3  # because the conf.gro counts every POPE molecule 3x
print "APL = {apl}".format(apl=apl)
