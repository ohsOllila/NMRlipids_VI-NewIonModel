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
import math
#import matplotlib.pyplot as plt
#import cPickle
from optparse import OptionParser
#import scipy, math

#try:
    #import read_xvg_calc_mean as rxvg
#except:
    #print "Couldn't read XVG-reading library, can't process Gromacs pullx.xvg files!"

k_b = 0.0083144621  #kJ/Mol*K
f_conc=55430  # factor for calculating concentrations of salts from numbers of ions/waters; in mM/L
bond_len_max=1.5  # in A, max distance between atoms for reasonable OP calculation
bond_len_max_sq=bond_len_max**2

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


    def calc_OP(self, atoms):
        """
        calculates Order Parameter according to equation
        S = 1/2 * (3*cos(theta)^2 -1)
        """
        vec = atoms[1].position - atoms[0].position
        d2 = np.square(vec).sum()
        if d2>bond_len_max_sq:
            raise UserWarning, "Atomic distance for atoms \
 {at1} and {at2} in residue no. {resnr} is suspiciously \
 long: {d}!\nPBC removed???".format(at1=atoms[0].name, at2=atoms[1].name, resnr=atoms[0].resid, d=math.sqrt(d2))
        cos2 = vec[2]**2/d2
        S = 0.5*(3.0*cos2-1.0)
        return S


    @property
    def get_avg_std_OP(self):
        """
        Provides average and stddev of all OPs in self.traj
        """
        # convert to numpy array
        return (np.mean(self.traj), np.std(self.traj))


def read_trajs(ordPars, top, trajs):
    """
    procedure that
    creates MDAnalysis Universe with top,
    reads in trajectories trajs and then
    goes through every frame and
    evaluates each Order Parameter "S" from the list of OPs ordPars.

    ordPars : list of OrderParameter class
       each item in this list describes an Order parameter to be calculated in the trajectory
    top : str
        filename of a top file (e.g. conf.gro)
    trajs : list of strings
        filenames of trajectories
    """
    # read-in topology and trajectory
    mol = mda.Universe(top, trajs)

    # make atom selections for each OP and store it as its attribute for later use in trajectory
    for op in ordPars.values():
        # selection = pairs of atoms, split-by residues
        selection = mol.select_atoms("resname {rnm} and name {atA} or name {atB}".format(
                                    rnm=op.resname, atA=op.atAname, atB=op.atBname)
                                    ).atoms.split("residue")
        for res in selection:
            # check if we have only 2 atoms (A & B) selected
            if res.n_atoms != 2:
                raise UserWarning, "Selection >> name {atA} {atB} << \
                contains {nat} atoms, but should contain exactly 2!".format(
                atA=op.atAname, atB=op.atBname, nat=len(sel.atoms))
        op.selection = selection

    # go through trajectory frame-by-frame
    for frame in mol.trajectory:
        for op in ordPars.values():
            for residue in op.selection:
                S = op.calc_OP(residue)
                op.traj.append(S)
        #print "--", mol.atoms[0].position




#%%

if __name__ == "__main__":
    # help message is automatically provided
    # type=string, action=store is default
    parser = OptionParser()
    parser.add_option('-i', '--inp',  dest='inp_fname',  help='input (OP definitions) file name', default="Headgroup_Glycerol_OPs.def")
    parser.add_option('-t', '--top',  dest='top_fname',  help='topology (gro, pdb) file name', default="last_frame_nonwat.gro")
    parser.add_option('-x', '--traj', dest='traj_fname', help='trajectory (xtc) file name', default="traj_nonwat_pbc.xtc")
    parser.add_option('-o', '--out',  dest='out_fname',  help='output (OPs mean&std) file name', default="Headgroup_Glycerol_OPs.dat")
    opts, args = parser.parse_args()

    # dictionary for storing of OrderParameter class instances (name-wise, of course)
    ordPars = {}
    # read-in OP definitions from the input file
    try:
        with open(opts.inp_fname,"r") as f:
            for line in f.readlines():
                if not line.startswith("#"):
                    items = line.split()
                    op_name = items[0]
                    resname = items[1]
                    atAname = items[2]
                    atBname = items[3]
                    ordPars[op_name] = OrderParameter(op_name, resname, atAname, atBname)
    except:
        raise RuntimeError, "Couldn't read input file >> {inpf} <<".format(inpf=opts.inp_fname)


#    ordPars = {}
#    beta1OP  = OrderParameter("beta1",  "POPC", "C12", "H12A")
#    beta2OP  = OrderParameter("beta2",  "POPC", "C12", "H12B")
#    alpha1OP = OrderParameter("alpha1", "POPC", "C11", "H11A")
#    alpha2OP = OrderParameter("alpha2", "POPC", "C11", "H11B")
#    op_list = (beta1OP, beta2OP, alpha1OP, alpha2OP)
#    for op in op_list:
#        ordPars[op.name] = op

#%%
    read_trajs(ordPars, opts.top_fname, opts.traj_fname)

#%%

    print "OP Name     mean     stddev "
    print "----------------------------"
    for op in ordPars.values():
        (op.avg, op.std) = op.get_avg_std_OP
        print op.name, op.avg, op.std


    try:
        with open(opts.out_fname,"w") as f:
            f.write("OP_name    resname    atom1    atom2    OP_mean   OP_stddev\n\
------------------------------------------------------------\n")
            for op in ordPars.values():
                f.write( "   ".join([op.name, op.resname, op.atAname, op.atBname, str(op.avg), str(op.std), "\n"]) )
        print "OrderParameters written to >> {fname} <<".format(fname=opts.out_fname)
    except:
        print "Problems writing main output file."


    #conc = f_conc*nion/nwat
    try:
        conc_formatted_line = "conc  {b1} 0  {b2} 0    {a1} 0  {a2} 0".format(
                              b1=ordPars['beta1'].avg, b2=ordPars['beta2'].avg,
                              a1=ordPars['alpha1'].avg, a2=ordPars['alpha2'].avg)
        print
        print "Single line format:\nconc  beta1 0  beta2 0  alpha1 0  alpha2 0"
        print conc_formatted_line
        with open(opts.out_fname+".line","w") as f:
            f.write(conc_formatted_line)
    except:
        print "Problems writing the beta-alpha single line format file."
