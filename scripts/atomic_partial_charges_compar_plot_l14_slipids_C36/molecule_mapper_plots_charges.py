# -*- coding: utf-8 -*-
"""
reorders atoms in molecules
according to the names in a reference molecule

and prints out their partial charges in plots.

IN: itp files (hard-coded in)

Originally developed to be used with NMRlipids project 6
by
J. Melcr
"""

import pmx
import pmx.forcefield
import numpy as np
import matplotlib.pyplot as plt
#import math

#%%

def fetch_atom(mol, at_name):
    """
    return the first occurence of an atom with a given name.
    return the 1st atom if searched atom is not found.
    """
    try:
        found = False
        for otherat in mol.atoms:
            if otherat.name == at_name:
                found = True
                #print "Found atom {atname}.".format(atname=at_name)
                return otherat
    except:
        print "Something went wrong during atom searching (function fetch_atom). \
        \nIs mol a pmx.forcefield.ITPfile object?"
    finally:
        if not found:
            print "Atom {atname} not found! -- will substitute it with atom no.1.".format(atname=at_name)
            # so that there's no missing space in the sequence.
            return mol.atoms[0]


def autolabel(rects, i_max=None, font_size=2.5, q_label_offset=0.03):
    """
    Attach a text label above/below each bar displaying its height
    This implementation works with a GLOBAL variable ax -- QUICK&DIRTY HACK!!
    """
    if i_max == None:
        i_max = len(rects)
    for i,rect in enumerate(rects):
        if i<i_max :
            height = rect.get_height()
            if rect.get_y() >= 0:
                ax.text(rect.get_x() + rect.get_width()/2., rect.get_y() + height + q_label_offset,
                    '%2.2f' % height,
                    ha='center', va='bottom',
                    fontsize=font_size,
                    rotation=90)
            else:
                ax.text(rect.get_x() + rect.get_width()/2., rect.get_y() - q_label_offset,
                    '-%2.2f' % height,
                    ha='center', va='top',
                    fontsize=font_size,
                    color='red',
                    rotation=90)
        else:
            break


def create_mapping_file_for_macrog(mol_d, mol_c, mapping_bx):
    """
    Returns mapping dictionary noted as mapping_dx
    (mol_d is MacRox, 'x' referres to the universal M_*_M language started in NMRlipids projects)
    assuming that molecules d and c have the same order of atoms.

    In addition: Writes the mapping_dx into a file.

    Parameters
    ------------
    mol_d, mol_c : pmx.ITPfile instances
        molecule topologies that have the same order of atoms
    mapping_bx : dictionary
        dictionary that translates molecule c to the atom-naming-language 'x'
    """

    # create mapping file forMacROg POPC from http://www.sciencedirect.com/science/article/pii/S2352340915002073
    # using atom names in the ITP files that are in the same order as for CHARMM POPC (stated in the reference)
    mapping_dx = {}
    for i,atom in enumerate(mol_d.atoms):
        mapping_dx[atom.name] = mapping_bx[mol_c.atoms[i].name]

    # write mappingFILEmacRogPOPC.txt -- works for MacRog POPC from the above reference
    with open("mappingFILEmacRogPOPC.txt","w") as f:
        f.write("# Individual atoms\n# mapping for MacRog POPC as implemented in reference http://www.sciencedirect.com/science/article/pii/S2352340915002073\n")
        for item in mapping_dx.items():
            f.write("    {x_name}    {d_name}\n".format(x_name=item[1], d_name=item[0]))
        f.write("""# Whole molecules
                M_POPC_M         POPC
                M_NA_M           NA
                M_CA_M           CA
                """)

    return mapping_dx


#%%

if __name__ == '__main__':

    # load-in the files
    mol_a = pmx.forcefield.ITPFile("l14_POPC.itp")
    mol_b = pmx.forcefield.ITPFile("slipids_POPC.itp")
    mol_c = pmx.forcefield.ITPFile("charmm_POPC.itp")
    mol_d = pmx.forcefield.ITPFile("macRog_POPC.itp")

#%%

    # translating A->B, so
    # B is a reference being a dictionary in order  FFname - Mapping_name
    # whereas A is being translated so a dictionary in opposite order is more practical
    # A: Mapping_name -> FFname
    # mappingPOPClipid14.txt -- works only for lipid14
    mapping_xa = {}
    with open("mappingPOPClipid14.txt","r") as f:
        for line in f.readlines():
            if not line.startswith("#"):
                items = line.split()
                mapping_xa[items[0]] = items[1]

    # mappingFILE.txt -- works for slipids & charmm
    mapping_bx = {}
    with open("mappingFILE.txt","r") as f:
        for line in f.readlines():
            if not line.startswith("#"):
                items = line.split()
                mapping_bx[items[1]] = items[0]

#    # mapping file for MacRog from earlier NMRlipids publication doesn't work
#    # with itp file from http://www.sciencedirect.com/science/article/pii/S2352340915002073
#
#    # mappingPOPCmacrog.txt -- works for MacRog model
#    mapping_xd = {}
#    with open("mappingPOPCmacrog.txt","r") as f:
#        for line in f.readlines():
#            if not line.startswith("#"):
#                items = line.split()
#                mapping_xd[items[0]] = items[1]


#%%
    # create dictionaries for transaltion from b->a
    mapping_ab = {}
    for key in mapping_bx.keys():
        mapping_ab[mapping_xa[mapping_bx[key]]] = key

    mapping_ba = {}
    for key in mapping_bx.keys():
        mapping_ba[key] = mapping_xa[mapping_bx[key]]

#    mapping_bd = {}
#    for key in mapping_bx.keys():
#        mapping_bd[key] = mapping_xd[mapping_bx[key]]


#%%
    # sort molecules after molecule C and record their charges
    mol_a_sorted_atoms_q = []
    mol_b_sorted_atoms_q = []
    mol_c_sorted_atoms_q = []
    mol_d_sorted_atoms_q = []
    mol_c_sorted_atoms_names = []
    mol_d_sorted_atoms_names = []
    # sort atoms after molecule C
    # following code re-uses variable otherat -- simple but DIRTY!
    for atom in mol_c.atoms:
        mol_a_sorted_atoms_q.append(fetch_atom(mol_a, mapping_ba[atom.name]).q)
        mol_b_sorted_atoms_q.append(fetch_atom(mol_b, atom.name).q)
        #mol_d_sorted_atoms_q.append(fetch_atom(mol_d, mapping_bd[atom.name]).q)
        mol_c_sorted_atoms_q.append(atom.q)
        mol_c_sorted_atoms_names.append(atom.name)
    # assuming mol D -- MacRog POPC -- atoms are already
    # CHARMM-like sorted as promised by authors.
    for atom in mol_d.atoms:
        mol_d_sorted_atoms_q.append(atom.q)


#%%
    ax_x_label_pos = []
    ax_x_labels    = []
    for i,name in enumerate(mol_c_sorted_atoms_names):
        if name.startswith(('C', 'N', 'P', 'O') ):
            ax_x_label_pos.append(i)
            ax_x_labels.append(name)

#%%
    # BAR PLOT 1
    #############
    bar_width = 0.2
    max_atom_index = 43

    fig, ax = plt.subplots()

    # add some text for labels, title and axes ticks
    ax.set_ylabel('Q')
    ax.set_title('Atomic partial charges')
    ax.set_xticks(np.array(ax_x_label_pos)+0.375)
    ax.set_xticklabels(ax_x_labels, rotation=89)
    ax.set_xlim(right=max_atom_index)  # Includes Headgroup and glycerol with carbonyls

    # plot bars
    rects1 = ax.bar(np.array(range(len(mol_c_sorted_atoms_names))),             mol_c_sorted_atoms_q, bar_width, color='navy',   linewidth=0, label="Charmm36")
    rects4 = ax.bar(np.array(range(len(mol_c_sorted_atoms_names)))+1*bar_width, mol_d_sorted_atoms_q, bar_width, color='pink',   linewidth=0, label="MacRog")
    rects2 = ax.bar(np.array(range(len(mol_c_sorted_atoms_names)))+2*bar_width, mol_b_sorted_atoms_q, bar_width, color='orange', linewidth=0, label="Slipids")
    rects3 = ax.bar(np.array(range(len(mol_c_sorted_atoms_names)))+3*bar_width, mol_a_sorted_atoms_q, bar_width, color='green',  linewidth=0, label="Lipid14")

    # assign labels
    autolabel(rects1, i_max=max_atom_index)
    autolabel(rects2, i_max=max_atom_index)
    autolabel(rects3, i_max=max_atom_index)
    autolabel(rects4, i_max=max_atom_index)

    ax.legend(loc="top left")
    plt.savefig("slipids_lipid14_charmm_macrog_charges_aligned_bars.png", dpi=300)


#%%
    # BAR PLOT 2
    #############
    fig, ax = plt.subplots()

    # add some text for labels, title and axes ticks
    ax.set_ylabel('Q')
    ax.set_title('Atomic partial charges')
    ax.set_xticks(np.array(ax_x_label_pos)+0.375)
    ax.set_xticklabels(ax_x_labels, rotation=89)

    # plot bars
    rects1 = ax.bar(np.array(range(len(mol_c_sorted_atoms_names))),             mol_c_sorted_atoms_q, bar_width, color='navy',   linewidth=0, label="Charmm36")
    rects4 = ax.bar(np.array(range(len(mol_c_sorted_atoms_names)))+1*bar_width, mol_d_sorted_atoms_q, bar_width, color='pink',   linewidth=0, label="MacRog")
    rects2 = ax.bar(np.array(range(len(mol_c_sorted_atoms_names)))+2*bar_width, mol_b_sorted_atoms_q, bar_width, color='orange', linewidth=0, label="Slipids")
    rects3 = ax.bar(np.array(range(len(mol_c_sorted_atoms_names)))+3*bar_width, mol_a_sorted_atoms_q, bar_width, color='green',  linewidth=0, label="Lipid14")

    # assign labels
    autolabel(rects1)
#    autolabel(rects2)
    autolabel(rects3)
#    autolabel(rects4)

    ax.legend()
    plt.savefig("slipids_lipid14_charmm_macrog_charges_aligned_bars_all-atoms.png", dpi=300)



#%%

