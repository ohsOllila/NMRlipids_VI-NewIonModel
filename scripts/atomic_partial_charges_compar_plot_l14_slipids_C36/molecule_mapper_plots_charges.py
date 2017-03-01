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
import math


#%%

if __name__ == '__main__':

    # load-in the files
    mol_a = pmx.forcefield.ITPFile("l14_POPC.itp")
    mol_b = pmx.forcefield.ITPFile("slipids_POPC.itp")
    mol_c = pmx.forcefield.ITPFile("charmm_POPC.itp")

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


#%%
    # create dictionaries for transaltion from b->a
    mapping_ab = {}
    for key in mapping_bx.keys():
        mapping_ab[mapping_xa[mapping_bx[key]]] = key

    mapping_ba = {}
    for key in mapping_bx.keys():
        mapping_ba[key] = mapping_xa[mapping_bx[key]]

#%%
    # sort molecules after molecule C and record their charges
    mol_a_sorted_atoms_q = []
    mol_b_sorted_atoms_q = []
    mol_c_sorted_atoms_q = []
    mol_c_sorted_atoms_names = []
    for atom in mol_c.atoms:
        mol_a_same_atom_name = mapping_ba[atom.name]
        for otherat in mol_a.atoms:
            if otherat.name == mol_a_same_atom_name:
                break
        mol_a_sorted_atoms_q.append(otherat.q)
        for otherat in mol_b.atoms:
            if otherat.name == atom.name:
                break
        mol_b_sorted_atoms_q.append(otherat.q)
        mol_c_sorted_atoms_q.append(atom.q)
        mol_c_sorted_atoms_names.append(atom.name)

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
            
    def autolabel(rects, i_max=None, font_size=4, q_label_offset=0.03):
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


    bar_width = 0.25

    fig, ax = plt.subplots()

    max_atom_index = 43

    # add some text for labels, title and axes ticks
    ax.set_ylabel('Q')
    ax.set_title('Atomic partial charges')
    ax.set_xticks(np.array(ax_x_label_pos)+0.375)
    ax.set_xticklabels(ax_x_labels, rotation=89)
    ax.set_xlim(right=max_atom_index)  # Includes Headgroup and glycerol with carbonyls

    # plot bars
    rects1 = ax.bar(range(len(mol_c_sorted_atoms_names)), mol_c_sorted_atoms_q, bar_width, color='navy', linewidth=0, label="Charmm36")
    rects2 = ax.bar(np.array(range(len(mol_c_sorted_atoms_names)))+bar_width, mol_b_sorted_atoms_q, bar_width, color='orange', linewidth=0, label="Slipids")
    rects3 = ax.bar(np.array(range(len(mol_c_sorted_atoms_names)))+2*bar_width, mol_a_sorted_atoms_q, bar_width, color='green', linewidth=0, label="Lipid14")

    # assign labels
    autolabel(rects1, i_max=max_atom_index)
    autolabel(rects2, i_max=max_atom_index)
    autolabel(rects3, i_max=max_atom_index)

    ax.legend()
    plt.savefig("slipids_lipid14_charmm_charges_aligned_bars.png", dpi=300)


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
    rects1 = ax.bar(range(len(mol_c_sorted_atoms_names)), mol_c_sorted_atoms_q, bar_width, color='navy', linewidth=0, label="Charmm36")
    rects2 = ax.bar(np.array(range(len(mol_c_sorted_atoms_names)))+bar_width, mol_b_sorted_atoms_q, bar_width, color='orange', linewidth=0, label="Slipids")
    rects3 = ax.bar(np.array(range(len(mol_c_sorted_atoms_names)))+2*bar_width, mol_a_sorted_atoms_q, bar_width, color='green', linewidth=0, label="Lipid14")

    # assign labels
    autolabel(rects1)
    autolabel(rects2)
    autolabel(rects3)

    ax.legend()
    plt.savefig("slipids_lipid14_charmm_charges_aligned_bars_all-atoms.png", dpi=300)



#%%

