# -*- coding: utf-8 -*-
"""
script for reordering atoms in molecules
according to the names of a reference molecule

Originally meant to reorder POPC atoms from CHARMM-GUI to Slipids ordering
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

    plt.xticks(ax_x_label_pos, ax_x_labels, rotation=89)
    plt.plot(mol_a_sorted_atoms_q, label="Lipid14")
    plt.plot(mol_b_sorted_atoms_q, label="Slipids")
    plt.plot(mol_c_sorted_atoms_q, label="Charmm36")
    plt.legend()
    plt.ylabel("Q")
    plt.savefig("slipids_lipid14_charmm_charges_aligned.png", dpi=300)
    plt.show()

#%%
    plt.xticks(ax_x_label_pos, ax_x_labels, rotation=89)
    plt.xlim([0,45])
    plt.plot(mol_a_sorted_atoms_q, label="Lipid14")
    plt.plot(mol_b_sorted_atoms_q, label="Slipids")
    plt.plot(mol_c_sorted_atoms_q, label="Charmm36")
    plt.legend()
    plt.ylabel("Q")
    plt.savefig("slipids_lipid14_charmm_charges_aligned_subplot.png", dpi=300)
    plt.show()

#%%
    plt.xticks(ax_x_label_pos, ax_x_labels, rotation=89)
    plt.xlim([0,45])
    plt.plot(np.abs(mol_a_sorted_atoms_q), label="Lipid14")
    plt.plot(np.abs(mol_b_sorted_atoms_q), label="Slipids")
    plt.plot(np.abs(mol_c_sorted_atoms_q), label="Charmm36")
    plt.legend()
    plt.ylabel("|Q|")
    plt.savefig("slipids_lipid14_charmm_charges_aligned_subplot_abs.png", dpi=300)
    plt.show()

#%%
    bar_width = 0.35

    fig, ax = plt.subplots()
    rects1 = ax.bar(range(len(mol_c_sorted_atoms_names)), mol_c_sorted_atoms_q, bar_width, color='r')
    rects1 = ax.bar(np.array(range(len(mol_c_sorted_atoms_names)))+bar_width, mol_b_sorted_atoms_q, bar_width, color='blue')

    # add some text for labels, title and axes ticks
    ax.set_ylabel('Q')
    ax.set_title('Atomic partial charges')
    ax.set_xticks(range(len(mol_c_sorted_atoms_names)))
    ax.set_xticklabels(mol_c_sorted_atoms_names)

    #ax.legend((rects1[0], rects2[0]), ('Men', 'Women'))


    def autolabel(rects):
        """
        Attach a text label above each bar displaying its height
        """
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                    '%d' % int(height),
                    ha='center', va='bottom')

    autolabel(rects1)
    #autolabel(rects2)
    plt.savefig("slipids_lipid14_charmm_charges_aligned_bars.png", dpi=300)
    plt.show()
