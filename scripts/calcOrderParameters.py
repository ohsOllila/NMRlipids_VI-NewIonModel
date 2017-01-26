#!/usr/bin/env python
"""
 calculation of order parameters from a MD trajectory

 ------------------------------------------------------------
 Made by Joe,  Last edit 2017/01/26
------------------------------------------------------------
 input:
 output: order parameters (pickled, textfile and plot)

--------------------------------------------------------
"""


# coding: utf-8

import MDAnalysis as mda
import numpy as np
#import scipy, math
import matplotlib.pyplot as plt
import cPickle
from optparse import OptionParser

try:
    import read_xvg_calc_mean as rxvg
except:
    print "Couldn't read XVG-reading library, can't process Gromacs pullx.xvg files!"

k_b = 0.0083144621  #kJ/Mol*K


#%%
class OrderParameter:
    """
    Class for storing&manipulating
    order parameter (OP) related metadata (definition, name, ...)
    and OP trajectories
    and methods to evaluate OPs.
    """
    def __init__(self, name, resname, atom_A_name, atom_B_name):
        """
        it doesn't matter which comes first,
        atom A or B, for OP calculation.
        """
        self.name = name             # name of the order parameter, a label
        self.resname = resname       # name of residue atoms are in
        self.atAname = atom_A_name
        self.atBname = atom_B_name
        for field in self.__dict__:
            if not isinstance(field, str):
                raise UserWarning, "provided name >> {} << is not a string! \n \
                Unexpected behaviour might occur.".format(field)
            else:
                if not field.strip():
                    raise RuntimeError, "provided name >> {} << is empty! \n \
                    Cannot use empty names for atoms and OP definitions.".format(field)
        self.traj = []  # for storing OPs



    def get_OPs(self, mol):
        """
        Goes through every residue in a frame
        (which is not necessary, but more robust and circumvents errors)
        and
        finds atom-name match with atoms A and B.
        Calculates OP for the pair and
        adds the OP to the OP's self.traj (which is a list)

        Works with MDAnalysis library

        Assuming 1 particular OP per 1 residue.
        """
        for res in mol.residues:
            sel = res.select_atoms("name {atA} {atB}".format(atA=self.atAname, atB=self.atBname))
            # check if we have only 2 atoms (A & B) selected
            if len(sel.atoms) != 2:
                raise UserWarning, "Selection >> name {atA} {atB} << \
                contains {nat} atoms, but should contain exactly 2!".format(
                atA=self.atAname, atB=self.atBname, nat=len(sel.atoms))
            op = self.calc_OP(sel)
            self.traj.append(op)
            print "--", mol.atoms[0].position


    def calc_OP(self, atoms):
        """
        calculates Order Parameter according to equation
        S = 1/2 * (3*cos(theta)^2 -1)
        """
        vec = atoms[1].position - atoms[0].position
        cos2 = vec[2]**2/sum(vec**2)
        S = 0.5*(3*cos2-1)
        return S


    @property
    def get_avg_std_OP(self):
        """
        Provides average and stddev of all OPs in self.traj
        """
        # convert to numpy array
        data = np.array(self.traj)
        return (data.mean(), data.std())


def read_trajs(ordPars, top, trajs):
    """
    procedure that
    creates MDAnalysis Universe with top,
    reads in trajectories trajs and then
    goes through every frame and
    evaluates each Order Parameter "S" from the list of OPs ordPars.

    top : str
        filename of a top file (e.g. conf.gro)
    trajs : list of strings
        filenames of trajectories
    """
    mol = mda.Universe(top, trajs)
    for frame in mol.trajectory:
        for op in ordPars.values():
            for res in mol.select_atoms("resname {rnm}".format(rnm=op.resname)).residues:
                sel = res.select_atoms("name {atA} {atB}".format(atA=op.atAname, atB=op.atBname))
                # check if we have only 2 atoms (A & B) selected
                if len(sel.atoms) != 2:
                    raise UserWarning, "Selection >> name {atA} {atB} << \
                    contains {nat} atoms, but should contain exactly 2!".format(
                    atA=op.atAname, atB=op.atBname, nat=len(sel.atoms))
                S = op.calc_OP(sel)
                op.traj.append(S)
        print "--", mol.atoms[0].position




#%%

if __name__ == "__main__":
    # help message is automatically provided
    # type=string, action=store is default
    parser = OptionParser()
    parser.add_option('-t', '--top',  dest='top_fname',  help='topology (gro, pdb) file name', default="confout.gro")
    parser.add_option('-x', '--traj', dest='traj_fname', help='trajectory (xtc) file name', default="traj_comp.xtc")
    opts, args = parser.parse_args()

    ordPars = {}
    beta1OP  = OrderParameter("beta1",  "POPC", "C12", "H12A")
    beta2OP  = OrderParameter("beta2",  "POPC", "C12", "H12B")
    alpha1OP = OrderParameter("alpha1", "POPC", "C11", "H11A")
    alpha2OP = OrderParameter("alpha2", "POPC", "C11", "H11B")
    op_list = (beta1OP, beta2OP, alpha1OP, alpha2OP)
    for op in op_list:
        ordPars[op.name] = op

#%%

    read_trajs(ordPars, opts.top_fname, opts.traj_fname)

#%%

    for op in ordPars.values():
        print op.name, op.get_avg_std_OP

#%%


#############################################
## ENDED HERE CONTINUE NEXT TIME...
## note: code above worked on test files
#############################################


    print "Saving the profile into a pickle-object and text file..."

    # "aliases" for plotting and saving
    #dTRAM data
    fep_dtram = dtram_obj.f_i
    fep_dtram_shift = np.copy(fep_dtram)
    fep_dtram_shift -= fep_dtram_shift.min()

    with open(opts.out_file_name+"_dTRAM-FEP_windows.pickle","w") as f : cPickle.dump(dtram_obj.f_K,f)
    with open(opts.out_file_name+"_dTRAM-FEP.pickle","w") as f : cPickle.dump(fep_dtram_shift,f)
    with open(opts.out_file_name+"_dTRAM-FEP.dat","w") as f:
        for i,c in enumerate(gridpoints):
            line = str(c)+"   "+str(fep_dtram_shift[i])+"\n"
            f.write( line )



    # In[46]:

    print "Plotting dTRAM, and WHAM data if present ... "

    # Load the result from WHAM: (generated with the same set and with the same converg. criterion 1.0E-6, 90 bins ..)
    try:
        wham_data = np.loadtxt(opts.wham_fep_fname)
        cnz_wham = wham_data[:,0]
        kT = k_b*sim_data.sims[0]['temperature']
        fep_wham = wham_data[:,1]/kT # scale by kT
        plt.plot( cnz_wham, fep_wham, '-', color='green', lw=2.0, label="WHAM" )
	try:
		fep_wham_err1 =  wham_data[:,2]/kT + fep_wham
		fep_wham_err2 = -wham_data[:,2]/kT + fep_wham
		plt.fill_between( cnz_wham, fep_wham_err1, fep_wham_err2, color='green', alpha=0.2)
	except:
		print "Couldn't read-in the errorbars for WHAM-data."
    except:
        print "Wham data missing or just something wrong happened in the process. Look into the code, lad. "

    plt.plot( gridpoints, fep_dtram_shift, '-', color='black', lw=2.0, label="dTRAM" )

    plt.legend( loc='upper left', fontsize=10 )
    plt.xlabel( r"$pitch$ [deg]", fontsize=12 )
    plt.ylabel( r"$G$ [kT]", fontsize=12 )

    plt.savefig(opts.out_file_name+"_dTRAM_WHAM.png", papertype="letter", dpi=300)

    # plotting dTRAM - Jacobi corrected (assuming Pitch-polar coordinates)
    kT = k_b*sim_data.sims[0]['temperature']
    plt.plot( gridpoints, fep_dtram_shift + kT*np.log(np.sin(gridpoints/180.0*math.pi)), '-', color='blue', lw=2.0, label="dTRAM Jacobi corrected" )

    plt.legend( loc='upper left', fontsize=10 )
    plt.xlabel( r"$pitch$ [deg]", fontsize=12 )
    plt.ylabel( r"$G$ [kT]", fontsize=12 )

    plt.savefig(opts.out_file_name+"_dTRAM_WHAM_JacobiCorr.png", papertype="letter", dpi=300)
